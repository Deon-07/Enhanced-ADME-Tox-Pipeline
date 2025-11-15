#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simplified ADME/Tox Pipeline with Enhanced Visualization
Comprehensive ADME/Tox analysis with improved plots and statistics.

Features:
- Enhanced molecular standardization 
- Improved visualization with better plots
- Multiple prediction models
- Comprehensive statistical analysis
- Better error handling and progress tracking
"""

import argparse
import logging
import sys
import traceback
from dataclasses import dataclass, asdict
from functools import partial
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import warnings

# Third-party imports with better error handling
try:
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.colors import LinearSegmentedColormap
    from openpyxl import load_workbook
    from openpyxl.drawing.image import Image as XLImage
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, QED, rdMolDescriptors, rdmolops, Draw
    from rdkit.Chem.FilterCatalog import FilterCatalogParams, FilterCatalog
    
    # Suppress RDKit warnings
    RDLogger.DisableLog('rdApp.*')
    warnings.filterwarnings('ignore')
    
    try:
        from rdkit.Chem import rdMolStandardize
        _HAS_RDSTANDARDIZE = True
    except ImportError:
        _HAS_RDSTANDARDIZE = False
        
except Exception as e:
    print("Missing dependencies. Please install: rdkit, pandas, numpy, matplotlib, seaborn, openpyxl")
    print(f"Error: {e}")
    sys.exit(1)

# Default configuration
DEFAULT_SMILES_COL = "smiles"
DEFAULT_ID_COL = "CID"
DEFAULT_QED_THRESHOLD = 0.55
DEFAULT_LOGS_THRESHOLD = -4.0
DEFAULT_CACO2_THRESHOLD = 1.85
DEFAULT_FSP3_THRESHOLD = 0.42
DEFAULT_CACO2_COEFFS = np.array([20.0, -0.02, -0.05, 0.30])

# Color schemes for plots
COLOR_SCHEMES = {
    'passed': '#2E8B57',  # SeaGreen
    'failed': '#DC143C',   # Crimson
    'neutral': '#4682B4',  # SteelBlue
    'warning': '#FF8C00',  # DarkOrange
}

# Enhanced logging
logger = logging.getLogger("enhanced-admet-pipeline")
logger.setLevel(logging.INFO)

class EnhancedFormatter(logging.Formatter):
    """Custom formatter with colors"""
    FORMATS = {
        logging.INFO: "\033[94m%(levelname)s:\033[0m %(message)s",
        logging.WARNING: "\033[93m%(levelname)s:\033[0m %(message)s",
        logging.ERROR: "\033[91m%(levelname)s:\033[0m %(message)s",
    }
    
    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno, "%(levelname)s: %(message)s")
        return logging.Formatter(log_fmt).format(record)

handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(EnhancedFormatter())
logger.addHandler(handler)

@dataclass
class EnhancedCompoundResult:
    """Enhanced result class with additional properties"""
    cid: str
    smiles_in: str
    smiles_std: Optional[str]
    valid: bool
    mw: Optional[float] = None
    logP: Optional[float] = None
    h_donors: Optional[int] = None
    h_acceptors: Optional[int] = None
    rot_bonds: Optional[int] = None
    tpsa: Optional[float] = None
    qed: Optional[float] = None
    fsp3: Optional[float] = None
    logS: Optional[float] = None
    logS_confidence: Optional[float] = None
    caco2: Optional[float] = None
    caco2_confidence: Optional[float] = None
    pains_matches: Optional[List[str]] = None
    filter_pass: Optional[bool] = None
    failure_reasons: Optional[List[str]] = None
    lipinski_violations: Optional[int] = None
    veber_violations: Optional[int] = None
    gsk_violations: Optional[int] = None
    synthetic_accessibility: Optional[float] = None

    def to_dict(self) -> Dict:
        d = asdict(self)
        if self.pains_matches is not None:
            d["pains_matches"] = ";".join(self.pains_matches) if self.pains_matches else ""
        if self.failure_reasons is not None:
            d["failure_reasons"] = ";".join(self.failure_reasons) if self.failure_reasons else ""
        return d

class EnhancedADMETAnalyzer:
    """Enhanced ADMET analysis with multiple models and visualization"""
    
    def __init__(self):
        self.pains_catalog = None
        
    def initialize_pains(self):
        """Initialize PAINS catalog"""
        try:
            params = FilterCatalogParams()
            params.AddCatalog(FilterCatalogParams.FilterCatalogParams.PAINS)
            self.pains_catalog = FilterCatalog(params)
        except Exception as e:
            logger.warning(f"Failed to initialize PAINS catalog: {e}")
            self.pains_catalog = None

    def enhanced_standardize_molecule(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """Enhanced molecular standardization"""
        if mol is None:
            return None
        
        try:
            # Remove salts and keep largest fragment
            frags = rdmolops.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
            if len(frags) > 1:
                frags = sorted(frags, key=lambda m: m.GetNumAtoms(), reverse=True)
                mol = frags[0]
            
            Chem.SanitizeMol(mol)
            
            # Advanced standardization if available
            if _HAS_RDSTANDARDIZE:
                try:
                    # Normalize functional groups
                    normalizer = rdMolStandardize.Normalizer()
                    mol = normalizer.normalize(mol)
                    
                    # Remove charges
                    uncharger = rdMolStandardize.Uncharger()
                    mol = uncharger.uncharge(mol)
                    
                    # Canonicalize tautomers
                    tautomer_enumerator = rdMolStandardize.TautomerEnumerator()
                    mol = tautomer_enumerator.Canonicalize(mol)
                    
                except Exception as e:
                    logger.debug(f"Advanced standardization failed: {e}")
            
            return mol
        except Exception as e:
            logger.debug(f"Standardization failed: {e}")
            return None

    def calculate_enhanced_logS(self, mol: Chem.Mol) -> Tuple[Optional[float], Optional[float]]:
        """Enhanced LogS prediction with confidence estimation"""
        try:
            # Multiple descriptor-based estimation
            logP = Descriptors.MolLogP(mol)
            mw = Descriptors.MolWt(mol)
            tpsa = rdMolDescriptors.CalcTPSA(mol)
            
            # Ensemble of models for better prediction
            model1 = 0.5 - logP - 0.01 * (mw - 100.0)  # Original model
            model2 = -0.01 * mw + 0.1 * logP - 0.002 * tpsa + 0.5  # TPSA-enhanced
            model3 = -0.008 * mw - 0.5 * logP - 0.001 * tpsa + 0.3  # Alternative
            
            logS_pred = np.mean([model1, model2, model3])
            confidence = 1.0 - (np.std([model1, model2, model3]) / 2.0)  # Simple confidence measure
            
            return float(logS_pred), float(max(0.0, min(1.0, confidence)))
        except Exception:
            return None, None

    def predict_enhanced_caco2(self, mol: Chem.Mol, mw: float, tpsa: float, 
                             logP: float, coeffs: np.ndarray) -> Tuple[Optional[float], Optional[float]]:
        """Enhanced Caco-2 prediction with confidence"""
        try:
            # Additional descriptors for better prediction
            h_bond_donor = rdMolDescriptors.CalcNumHBD(mol)
            h_bond_acceptor = rdMolDescriptors.CalcNumHBA(mol)
            
            # Extended feature vector
            features = np.array([1.0, mw, tpsa, logP, h_bond_donor, h_bond_acceptor])
            
            # Simple linear model with extended coefficients
            if len(coeffs) < len(features):
                # Pad coefficients if needed
                extended_coeffs = np.pad(coeffs, (0, len(features) - len(coeffs)), 
                                       mode='constant', constant_values=0.0)
            else:
                extended_coeffs = coeffs[:len(features)]
            
            prediction = float(np.dot(extended_coeffs, features))
            
            # Simple confidence based on molecular complexity
            complexity = (mw / 500.0 + tpsa / 140.0) / 2.0
            confidence = max(0.0, min(1.0, 1.0 - complexity))
            
            return prediction, confidence
        except Exception:
            return None, None

    def calculate_synthetic_accessibility(self, mol: Chem.Mol) -> Optional[float]:
        """Simple synthetic accessibility score"""
        try:
            mw = Descriptors.MolWt(mol)
            logP = Descriptors.MolLogP(mol)
            rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
            rings = rdMolDescriptors.CalcNumRings(mol)
            stereo_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
            
            # Simple SA score based on various factors
            sa_score = (mw / 1000.0 + abs(logP) / 10.0 + rot_bonds / 20.0 + 
                       rings / 10.0 + stereo_centers / 5.0)
            
            return float(min(1.0, sa_score))
        except Exception:
            return None

def compute_enhanced_row(args_tuple: Tuple[int, Dict], analyzer: EnhancedADMETAnalyzer,
                        smiles_col: str, id_col: str, coeffs_caco2: np.ndarray,
                        qed_threshold: float, logs_threshold: float, 
                        caco2_threshold: float, fsp3_threshold: float,
                        enable_pains: bool) -> Tuple[int, EnhancedCompoundResult]:
    """Enhanced row computation with additional metrics"""
    idx, row = args_tuple
    cid = str(row.get(id_col, "")).strip() if row.get(id_col, None) is not None else f"row{idx}"
    smiles_in = row.get(smiles_col, "")
    
    try:
        # Initial molecule creation
        mol0 = Chem.MolFromSmiles(smiles_in)
        if mol0 is None:
            return idx, EnhancedCompoundResult(cid, smiles_in, None, False, 
                                            failure_reasons=["invalid_smiles"])

        # Enhanced standardization
        mol = analyzer.enhanced_standardize_molecule(mol0)
        if mol is None:
            return idx, EnhancedCompoundResult(cid, smiles_in, None, False, 
                                            failure_reasons=["standardization_failed"])

        # Calculate basic descriptors
        mw = Descriptors.MolWt(mol)
        logP = Descriptors.MolLogP(mol)
        h_donors = rdMolDescriptors.CalcNumHBD(mol)
        h_acceptors = rdMolDescriptors.CalcNumHBA(mol)
        rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        tpsa = rdMolDescriptors.CalcTPSA(mol)
        qedv = QED.qed(mol)
        fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)
        
        # Enhanced predictions with confidence
        logS, logS_confidence = analyzer.calculate_enhanced_logS(mol)
        caco2_val, caco2_confidence = analyzer.predict_enhanced_caco2(mol, mw, tpsa, logP, coeffs_caco2)
        sa_score = analyzer.calculate_synthetic_accessibility(mol)

        # PAINS filtering
        pains_matches = []
        if enable_pains and analyzer.pains_catalog:
            try:
                matches = analyzer.pains_catalog.GetMatches(mol)
                pains_matches = [entry.GetDescription() for match in matches 
                               for entry in [match[1]] if hasattr(entry, 'GetDescription')]
            except Exception:
                pains_matches = []

        # Rule-based filtering with detailed violations
        failures = []
        lipinski_violations = 0
        veber_violations = 0
        gsk_violations = 0

        # Lipinski violations
        lip_rules = [
            (mw > 500.0, "MW>500"),
            (logP > 5.0, "LogP>5"),
            (h_donors > 5, "HBD>5"),
            (h_acceptors > 10, "HBA>10")
        ]
        lipinski_violations = sum(1 for condition, _ in lip_rules if condition)
        if lipinski_violations >= 2:
            failures.append(f"Lipinski_{lipinski_violations}")

        # Veber violations
        if rot_bonds > 10:
            veber_violations += 1
            failures.append("Veber_RotBonds")
        if tpsa > 140:
            veber_violations += 1
            failures.append("Veber_TPSA")

        # GSK violations
        if mw > 400:
            gsk_violations += 1
            failures.append("GSK_MW")
        if logP > 4.0:
            gsk_violations += 1
            failures.append("GSK_LogP")

        # Additional filters
        if enable_pains and pains_matches:
            failures.append("PAINS")
        if qedv is not None and qedv <= qed_threshold:
            failures.append("QED")
        if logS is not None and logS < logs_threshold:
            failures.append("LogS")
        if caco2_val is not None and caco2_val < caco2_threshold:
            failures.append("Caco2")
        if fsp3 is not None and fsp3 < fsp3_threshold:
            failures.append("Fsp3")

        filter_pass = len(failures) == 0
        smiles_std = Chem.MolToSmiles(mol, isomericSmiles=True) if mol is not None else None

        result = EnhancedCompoundResult(
            cid=cid,
            smiles_in=smiles_in,
            smiles_std=smiles_std,
            valid=True,
            mw=float(mw),
            logP=float(logP),
            h_donors=int(h_donors),
            h_acceptors=int(h_acceptors),
            rot_bonds=int(rot_bonds),
            tpsa=float(tpsa),
            qed=float(qedv),
            fsp3=float(fsp3),
            logS=float(logS) if logS is not None else None,
            logS_confidence=float(logS_confidence) if logS_confidence is not None else None,
            caco2=float(caco2_val) if caco2_val is not None else None,
            caco2_confidence=float(caco2_confidence) if caco2_confidence is not None else None,
            pains_matches=pains_matches,
            filter_pass=filter_pass,
            failure_reasons=failures,
            lipinski_violations=lipinski_violations,
            veber_violations=veber_violations,
            gsk_violations=gsk_violations,
            synthetic_accessibility=float(sa_score) if sa_score is not None else None
        )
        return idx, result
        
    except Exception as e:
        logger.debug(f"compute_enhanced_row exception: {e}")
        return idx, EnhancedCompoundResult(cid, smiles_in, None, False, 
                                        failure_reasons=[f"exception:{str(e)}"])

class EnhancedVisualization:
    """Enhanced visualization with multiple plot types"""
    
    @staticmethod
    def create_enhanced_histograms(df: pd.DataFrame, out_dir: Path) -> List[Path]:
        """Create enhanced histogram plots with distribution analysis"""
        out_dir.mkdir(parents=True, exist_ok=True)
        files = []
        
        # Set style
        plt.style.use('default')
        sns.set_palette("husl")
        
        properties = [
            ('mw', 'Molecular Weight', 'Da', (0, 1000)),
            ('logP', 'LogP', '', (-5, 10)),
            ('tpsa', 'TPSA', 'Å²', (0, 300)),
            ('qed', 'QED Score', '', (0, 1)),
            ('logS', 'Predicted LogS', 'log mol/L', (-10, 2)),
            ('caco2', 'Predicted Caco-2', 'log nm/s', (0, 3)),
            ('fsp3', 'Fsp3', 'Fraction', (0, 1)),
            ('synthetic_accessibility', 'Synthetic Accessibility', 'Score', (0, 1))
        ]
        
        for prop, label, unit, xrange in properties:
            if prop not in df.columns or df[prop].dropna().empty:
                continue
                
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
            
            # Histogram with KDE
            data = df[prop].dropna()
            sns.histplot(data=data, ax=ax1, kde=True, color=COLOR_SCHEMES['neutral'], alpha=0.7)
            ax1.axvline(data.mean(), color='red', linestyle='--', label=f'Mean: {data.mean():.2f}')
            ax1.axvline(data.median(), color='orange', linestyle='--', label=f'Median: {data.median():.2f}')
            ax1.set_xlabel(f'{label} ({unit})' if unit else label)
            ax1.set_ylabel('Count')
            ax1.set_title(f'Distribution of {label}')
            ax1.legend()
            ax1.set_xlim(xrange)
            
            # Box plot
            sns.boxplot(x=data, ax=ax2, color=COLOR_SCHEMES['neutral'])
            if len(data) <= 100:
                # Add data points for small datasets
                sns.stripplot(x=data, ax=ax2, color=COLOR_SCHEMES['warning'], size=4, alpha=0.7)
            ax2.set_xlabel(f'{label} ({unit})' if unit else label)
            ax2.set_title(f'Box Plot of {label}')
            
            plt.tight_layout()
            fp = out_dir / f"{prop}_enhanced_dist.png"
            fig.savefig(fp, dpi=300, bbox_inches='tight')
            plt.close(fig)
            files.append(fp)
            
        return files

    @staticmethod
    def create_rule_violation_heatmap(df: pd.DataFrame, out_dir: Path) -> Optional[Path]:
        """Create heatmap of rule violations"""
        if 'failure_reasons' not in df.columns:
            return None
            
        # Extract all failure reasons
        all_failures = set()
        for reasons in df['failure_reasons'].dropna():
            if isinstance(reasons, str):
                all_failures.update(reason.strip() for reason in reasons.split(';'))
        
        if not all_failures:
            return None
            
        # Create violation matrix
        violation_data = []
        for failure in sorted(all_failures):
            row = [failure]
            for _, compound in df.iterrows():
                reasons = compound.get('failure_reasons', '')
                if isinstance(reasons, str):
                    violation = failure in [r.strip() for r in reasons.split(';')]
                else:
                    violation = False
                row.append(1 if violation else 0)
            violation_data.append(row)
        
        if not violation_data:
            return None
            
        # Create heatmap
        fig, ax = plt.subplots(figsize=(12, 8))
        data_array = np.array([row[1:] for row in violation_data])
        
        sns.heatmap(data_array, ax=ax, cmap=['white', 'red'], cbar_kws={'label': 'Violation'},
                   xticklabels=False, yticklabels=[row[0] for row in violation_data])
        
        ax.set_ylabel('Failure Reasons')
        ax.set_xlabel('Compounds')
        ax.set_title('Rule Violation Heatmap')
        
        plt.tight_layout()
        fp = out_dir / "rule_violation_heatmap.png"
        fig.savefig(fp, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        return fp

    @staticmethod
    def create_correlation_matrix(df: pd.DataFrame, out_dir: Path) -> Optional[Path]:
        """Create enhanced correlation matrix"""
        numeric_cols = ['mw', 'logP', 'tpsa', 'qed', 'logS', 'caco2', 'fsp3', 
                       'synthetic_accessibility', 'h_donors', 'h_acceptors', 'rot_bonds']
        numeric_cols = [col for col in numeric_cols if col in df.columns]
        
        if len(numeric_cols) < 2:
            return None
            
        corr_data = df[numeric_cols].dropna()
        if corr_data.empty:
            return None
            
        # Calculate correlation matrix
        corr_matrix = corr_data.corr()
        
        # Create masked matrix for upper triangle
        mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
        
        # Create plot
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Custom diverging colormap
        cmap = sns.diverging_palette(220, 10, as_cmap=True)
        
        # Plot heatmap
        sns.heatmap(corr_matrix, mask=mask, cmap=cmap, vmax=1, vmin=-1, center=0,
                   square=True, linewidths=0.5, cbar_kws={"shrink": .8}, ax=ax,
                   annot=True, fmt=".2f", annot_kws={"size": 8})
        
        ax.set_title('Property Correlation Matrix', fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        fp = out_dir / "enhanced_correlation_matrix.png"
        fig.savefig(fp, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        return fp

    @staticmethod
    def create_property_radar_chart(df: pd.DataFrame, out_dir: Path) -> Optional[Path]:
        """Create radar chart for property distributions"""
        if df.empty:
            return None
            
        # Select properties for radar chart
        properties = ['mw', 'logP', 'tpsa', 'qed', 'fsp3']
        properties = [p for p in properties if p in df.columns]
        
        if len(properties) < 3:
            return None
            
        # Normalize properties for radar chart
        normalized_data = {}
        for prop in properties:
            data = df[prop].dropna()
            if data.empty:
                continue
            # Min-max normalization
            min_val, max_val = data.min(), data.max()
            if max_val > min_val:
                normalized_data[prop] = (data - min_val) / (max_val - min_val)
            else:
                normalized_data[prop] = data * 0 + 0.5  # Center if all same
        
        if not normalized_data:
            return None
            
        # Create radar chart
        categories = list(normalized_data.keys())
        N = len(categories)
        
        # Compute angles for radar chart
        angles = [n / float(N) * 2 * np.pi for n in range(N)]
        angles += angles[:1]  # Complete the circle
        
        fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection='polar'))
        
        # Plot percentiles
        for percentile in [10, 25, 50, 75, 90]:
            values = [np.percentile(normalized_data[prop], percentile) for prop in categories]
            values += values[:1]
            ax.plot(angles, values, 'o-', linewidth=1, label=f'{percentile}th %ile', alpha=0.7)
        
        # Add filled area for interquartile range
        q25 = [np.percentile(normalized_data[prop], 25) for prop in categories]
        q75 = [np.percentile(normalized_data[prop], 75) for prop in categories]
        q25 += q25[:1]
        q75 += q75[:1]
        ax.fill_between(angles, q25, q75, alpha=0.3, color=COLOR_SCHEMES['neutral'])
        
        # Add property labels
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(categories)
        ax.set_ylim(0, 1)
        ax.set_title('Property Distribution Radar Chart', size=14, fontweight='bold')
        ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0))
        
        plt.tight_layout()
        fp = out_dir / "property_radar_chart.png"
        fig.savefig(fp, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        return fp

    @staticmethod
    def create_comparison_plots(before_df: pd.DataFrame, after_df: pd.DataFrame, out_dir: Path) -> List[Path]:
        """Create before/after comparison plots"""
        files = []
        
        comparison_props = ['qed', 'logS', 'caco2', 'fsp3', 'mw', 'logP', 'tpsa']
        
        for prop in comparison_props:
            if prop not in before_df.columns or prop not in after_df.columns:
                continue
                
            before_data = before_df[prop].dropna()
            after_data = after_df[prop].dropna()
            
            if before_data.empty or after_data.empty:
                continue
            
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
            
            # Overlay histograms
            ax1.hist(before_data, bins=30, alpha=0.7, label='Before Filtering', color=COLOR_SCHEMES['neutral'])
            ax1.hist(after_data, bins=30, alpha=0.7, label='After Filtering', color=COLOR_SCHEMES['passed'])
            ax1.set_xlabel(prop)
            ax1.set_ylabel('Count')
            ax1.set_title(f'{prop} Distribution: Before vs After')
            ax1.legend()
            
            # Box plot comparison
            plot_data = [before_data, after_data]
            ax2.boxplot(plot_data, labels=['Before', 'After'])
            ax2.set_ylabel(prop)
            ax2.set_title(f'{prop} Box Plot: Before vs After')
            
            plt.tight_layout()
            fp = out_dir / f"comparison_{prop}.png"
            fig.savefig(fp, dpi=300, bbox_inches='tight')
            plt.close(fig)
            files.append(fp)
            
        return files

def generate_html_report(df: pd.DataFrame, output_path: Path, plot_files: List[Path]):
    """Generate an interactive HTML report"""
    try:
        # Basic statistics
        total_compounds = len(df)
        valid_compounds = df['valid'].sum() if 'valid' in df.columns else total_compounds
        passed_compounds = df['filter_pass'].sum() if 'filter_pass' in df.columns else 0
        
        # Create HTML content
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>ADME/Tox Analysis Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 40px; }}
                .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 10px; }}
                .stats {{ display: flex; justify-content: space-around; margin: 20px 0; }}
                .stat-box {{ background: white; padding: 15px; border-radius: 5px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); text-align: center; }}
                .plots {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(500px, 1fr)); gap: 20px; }}
                img {{ max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 5px; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>ADME/Tox Analysis Report</h1>
                <p>Generated on {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            </div>
            
            <div class="stats">
                <div class="stat-box">
                    <h3>Total Compounds</h3>
                    <p style="font-size: 24px; color: #333;">{total_compounds}</p>
                </div>
                <div class="stat-box">
                    <h3>Valid Compounds</h3>
                    <p style="font-size: 24px; color: #2E8B57;">{valid_compounds}</p>
                </div>
                <div class="stat-box">
                    <h3>Passed Filters</h3>
                    <p style="font-size: 24px; color: #2E8B57;">{passed_compounds}</p>
                </div>
                <div class="stat-box">
                    <h3>Success Rate</h3>
                    <p style="font-size: 24px; color: #2E8B57;">{(passed_compounds/valid_compounds*100 if valid_compounds > 0 else 0):.1f}%</p>
                </div>
            </div>
            
            <div class="plots">
        """
        
        # Add plots to HTML
        for plot_file in plot_files:
            if plot_file.suffix == '.png':
                html_content += f'<div><h3>{plot_file.stem}</h3><img src="{plot_file.name}" alt="{plot_file.stem}"></div>'
        
        html_content += """
            </div>
        </body>
        </html>
        """
        
        # Write HTML file
        report_path = output_path.parent / f"{output_path.stem}_report.html"
        with open(report_path, 'w') as f:
            f.write(html_content)
        
        logger.info(f"HTML report generated: {report_path}")
        return report_path
        
    except Exception as e:
        logger.warning(f"Failed to generate HTML report: {e}")
        return None

def main():
    """Enhanced main function with better error handling and visualization"""
    parser = argparse.ArgumentParser(description="Enhanced ADME/Tox Pipeline with Advanced Visualization")
    parser.add_argument("input", type=str, help="Input CSV file with SMILES")
    parser.add_argument("-o", "--output", type=str, default=None, help="Output Excel file")
    parser.add_argument("--smiles-col", type=str, default=DEFAULT_SMILES_COL, help="SMILES column name")
    parser.add_argument("--id-col", type=str, default=DEFAULT_ID_COL, help="ID column name")
    parser.add_argument("--no-pains", action="store_true", help="Disable PAINS filtering")
    parser.add_argument("--n-jobs", type=int, default=max(1, cpu_count() - 1), help="Number of parallel workers")
    parser.add_argument("--calibrate-caco2", type=str, default=None, help="CSV for Caco-2 calibration")
    parser.add_argument("--qed-threshold", type=float, default=DEFAULT_QED_THRESHOLD)
    parser.add_argument("--logs-threshold", type=float, default=DEFAULT_LOGS_THRESHOLD)
    parser.add_argument("--caco2-threshold", type=float, default=DEFAULT_CACO2_THRESHOLD)
    parser.add_argument("--fsp3-threshold", type=float, default=DEFAULT_FSP3_THRESHOLD)
    parser.add_argument("--generate-html", action="store_true", help="Generate HTML report")
    
    args = parser.parse_args()
    
    # Input validation
    input_path = Path(args.input)
    if not input_path.exists():
        logger.error(f"Input file not found: {input_path}")
        return 1

    if args.output:
        output_file = Path(args.output)
    else:
        output_file = input_path.with_name(input_path.stem + "_admet_enhanced.xlsx")

    logger.info(f"Input: {input_path}")
    logger.info(f"Output: {output_file}")
    
    # Read input data
    try:
        df = pd.read_csv(input_path, dtype=str)
    except Exception as e:
        logger.error(f"Failed to read input CSV: {e}")
        return 1

    # Find SMILES and ID columns
    found_smiles = None
    for c in df.columns:
        if c.lower() == args.smiles_col.lower():
            found_smiles = c
            break
    if found_smiles is None and args.smiles_col in df.columns:
        found_smiles = args.smiles_col
    if found_smiles is None:
        logger.error(f"SMILES column '{args.smiles_col}' not found. Available columns: {list(df.columns)}")
        return 1

    found_id = None
    for c in df.columns:
        if c.lower() == args.id_col.lower():
            found_id = c
            break
    if found_id is None and args.id_col in df.columns:
        found_id = args.id_col
    if found_id is None:
        logger.warning(f"ID column '{args.id_col}' not found. Using row numbers.")
        found_id = "row_number"

    # Initialize analyzer
    analyzer = EnhancedADMETAnalyzer()
    if not args.no_pains:
        analyzer.initialize_pains()

    # Use default coefficients (calibration removed for simplicity)
    coeffs_caco2 = DEFAULT_CACO2_COEFFS.copy()

    # Process compounds
    records = list(df.to_dict(orient="records"))
    indexed = list(enumerate(records))
    
    worker = partial(compute_enhanced_row,
                     analyzer=analyzer,
                     smiles_col=found_smiles,
                     id_col=found_id,
                     coeffs_caco2=coeffs_caco2,
                     qed_threshold=args.qed_threshold,
                     logs_threshold=args.logs_threshold,
                     caco2_threshold=args.caco2_threshold,
                     fsp3_threshold=args.fsp3_threshold,
                     enable_pains=not args.no_pains)

    results = []
    if args.n_jobs > 1:
        logger.info(f"Running in parallel with {args.n_jobs} workers")
        try:
            with Pool(processes=args.n_jobs) as pool:
                for idx, res in pool.imap_unordered(worker, indexed, chunksize=10):
                    results.append(res)
        except Exception as e:
            logger.warning(f"Parallel processing failed: {e}. Falling back to serial.")
            for item in indexed:
                _, r = worker(item)
                results.append(r)
    else:
        logger.info("Running serially")
        for item in indexed:
            _, r = worker(item)
            results.append(r)

    # Create output DataFrame
    out_df = pd.DataFrame([r.to_dict() for r in results])
    
    # Remove duplicates based on standardized SMILES
    if "smiles_std" in out_df.columns:
        out_df = out_df.drop_duplicates(subset=["smiles_std"], keep="first").reset_index(drop=True)

    # Write Excel output
    try:
        with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
            out_df.to_excel(writer, index=False, sheet_name="Results")
            
            # Add explanations sheet
            expl_data = [
                ("Lipinski", "Allow up to 1 violation of (MW<=500, LogP<=5, Hdonors<=5, Hacceptors<=10)"),
                ("Veber", "Rotatable bonds <=10 and TPSA <=140"),
                ("PAINS", "PAINS flagged SMARTS via RDKit FilterCatalog"),
                ("QED", f"QED threshold > {args.qed_threshold}"),
                ("LogS", f"Predicted LogS >= {args.logs_threshold} (log mol/L)"),
                ("Caco2", f"Predicted Caco-2 >= {args.caco2_threshold} (log nm/s)"),
                ("GSK", "GSK-like rule: MW <= 400 and LogP <= 4"),
                ("Fsp3", f"Fsp3 threshold >= {args.fsp3_threshold}")
            ]
            pd.DataFrame(expl_data, columns=["Filter", "Description"]).to_excel(
                writer, index=False, sheet_name="Filter_Explanations")
            
            # Add parameters sheet
            params_data = {
                "Parameter": ["QED threshold", "LogS threshold", "Caco2 threshold", "Fsp3 threshold", "PAINS enabled"],
                "Value": [args.qed_threshold, args.logs_threshold, args.caco2_threshold, args.fsp3_threshold, not args.no_pains]
            }
            pd.DataFrame(params_data).to_excel(writer, index=False, sheet_name="Parameters")
            
    except Exception as e:
        logger.error(f"Failed to write Excel file: {e}")
        return 1

    # Create enhanced visualizations
    plot_dir = output_file.parent / f"{output_file.stem}_plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    
    viz = EnhancedVisualization()
    plot_files = []
    
    # Generate various plots
    plot_files.extend(viz.create_enhanced_histograms(out_df, plot_dir))
    
    if correlation_plot := viz.create_correlation_matrix(out_df, plot_dir):
        plot_files.append(correlation_plot)
    
    if radar_plot := viz.create_property_radar_chart(out_df, plot_dir):
        plot_files.append(radar_plot)
    
    if heatmap_plot := viz.create_rule_violation_heatmap(out_df, plot_dir):
        plot_files.append(heatmap_plot)
    
    # Before/after comparison plots
    before_df = out_df[out_df["valid"] == True] if "valid" in out_df.columns else out_df
    after_df = out_df[out_df["filter_pass"] == True] if "filter_pass" in out_df.columns else pd.DataFrame()
    
    if not after_df.empty:
        comparison_plots = viz.create_comparison_plots(before_df, after_df, plot_dir)
        plot_files.extend(comparison_plots)

    # Generate HTML report if requested
    if args.generate_html:
        generate_html_report(out_df, output_file, plot_files)

    # Summary statistics
    total = len(out_df)
    valid = len(before_df) if "valid" in out_df.columns else total
    passed = len(after_df) if "filter_pass" in out_df.columns else 0
    pass_rate = (passed / valid * 100) if valid > 0 else 0

    logger.info(f"Summary: {total} total, {valid} valid, {passed} passed ({pass_rate:.1f}%)")
    logger.info(f"Results saved to: {output_file}")
    logger.info(f"Plots saved to: {plot_dir}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
