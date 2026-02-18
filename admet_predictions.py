"""
ADMET Predictions Module
Calculates key ADMET properties for drug molecules
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, QED, Lipinski
import numpy as np


def calculate_logp(mol):
    """
    Calculate lipophilicity (LogP) using Wildman-Crippen method
    
    Optimal range: 0-3 (Lipinski's Rule of Five)
    - Too low: poor membrane permeability
    - Too high: poor solubility, metabolic issues
    """
    return Crippen.MolLogP(mol)


def estimate_logs(mol):
    """
    Estimate aqueous solubility (LogS) using ESOL method (Delaney)
    
    Optimal range: > -4
    - LogS < -6: poorly soluble
    - LogS > -4: good solubility
    """
    # ESOL formula: LogS = 0.16 - 0.63*cLogP - 0.0062*MW + 0.066*RB - 0.74*AP
    logp = Crippen.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    rb = Lipinski.NumRotatableBonds(mol)
    ap = Lipinski.NumAromaticRings(mol)  # Aromatic proportion approximation
    
    logs = 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rb - 0.74 * ap
    return logs


def estimate_herg_risk(mol):
    """
    Estimate hERG cardiotoxicity risk using heuristic rules
    
    Returns probability (0-1)
    - < 0.3: Low risk
    - 0.3-0.7: Moderate risk
    - > 0.7: High risk
    
    Based on: aromatic rings, basic nitrogen, LogP, molecular weight
    """
    logp = Crippen.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    aromatic_rings = Lipinski.NumAromaticRings(mol)
    hbd = Lipinski.NumHDonors(mol)
    
    # Heuristic scoring (simplified)
    risk_score = 0.0
    
    # High LogP increases hERG binding
    if logp > 4:
        risk_score += 0.3
    elif logp > 3:
        risk_score += 0.15
    
    # Large molecules tend to bind hERG
    if mw > 400:
        risk_score += 0.2
    
    # Multiple aromatic rings increase risk
    if aromatic_rings > 2:
        risk_score += 0.25
    
    # Basic nitrogen (common in hERG binders)
    if hbd == 0:
        risk_score += 0.15
    
    return min(risk_score, 1.0)


def calculate_qed(mol):
    """
    Calculate Quantitative Estimate of Drug-likeness (QED)
    
    Range: 0-1 (higher is better)
    Combines: MW, LogP, HBD, HBA, PSA, rotatable bonds, aromatic rings, alerts
    """
    return QED.qed(mol)


def calculate_lipinski_violations(mol):
    """
    Check Lipinski's Rule of Five violations
    
    Rules:
    - MW ≤ 500 Da
    - LogP ≤ 5
    - HBD ≤ 5
    - HBA ≤ 10
    """
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    
    violations = 0
    if mw > 500:
        violations += 1
    if logp > 5:
        violations += 1
    if hbd > 5:
        violations += 1
    if hba > 10:
        violations += 1
    
    return violations


def calculate_admet_profile(smiles):
    """
    Calculate complete ADMET profile for a molecule
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        dict: ADMET properties and scores
    """
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return None
    
    logp = calculate_logp(mol)
    logs = estimate_logs(mol)
    herg_risk = estimate_herg_risk(mol)
    qed_score = calculate_qed(mol)
    lipinski_violations = calculate_lipinski_violations(mol)
    
    # Normalized scores for radar chart (0-1 scale)
    # Higher is better for all normalized scores
    logp_norm = 1 - abs(logp - 1.5) / 3.5  # Optimal at 1.5, range 0-5
    logs_norm = max(0, min(1, (logs + 6) / 4))  # Scale from -6 to -2
    herg_norm = 1 - herg_risk  # Inverse (lower risk is better)
    qed_norm = qed_score  # Already 0-1
    
    return {
        'logp': round(logp, 2),
        'logs': round(logs, 2),
        'herg_risk': round(herg_risk, 2),
        'qed': round(qed_score, 2),
        'lipinski_violations': lipinski_violations,
        'normalized': {
            'logp': round(logp_norm, 2),
            'logs': round(logs_norm, 2),
            'herg_safety': round(herg_norm, 2),
            'qed': round(qed_norm, 2)
        }
    }


def interpret_logp(logp):
    """Interpret LogP value"""
    if logp < 0:
        return "⚠️ Very hydrophilic - poor membrane permeability"
    elif logp <= 3:
        return "✅ Optimal lipophilicity - good oral bioavailability"
    elif logp <= 5:
        return "⚠️ High lipophilicity - potential solubility issues"
    else:
        return "❌ Excessive lipophilicity - metabolic and toxicity concerns"


def interpret_logs(logs):
    """Interpret LogS value"""
    if logs > -2:
        return "✅ Highly soluble - excellent formulation potential"
    elif logs > -4:
        return "✅ Good solubility - suitable for most formulations"
    elif logs > -6:
        return "⚠️ Moderate solubility - formulation challenges"
    else:
        return "❌ Poor solubility - significant formulation issues"


def interpret_herg(risk):
    """Interpret hERG risk"""
    if risk < 0.3:
        return "✅ Low cardiotoxicity risk"
    elif risk < 0.7:
        return "⚠️ Moderate cardiotoxicity risk - requires testing"
    else:
        return "❌ High cardiotoxicity risk - likely hERG liability"


def interpret_qed(qed):
    """Interpret QED score"""
    if qed > 0.7:
        return "✅ Excellent drug-likeness"
    elif qed > 0.5:
        return "✅ Good drug-likeness"
    elif qed > 0.3:
        return "⚠️ Moderate drug-likeness - optimization needed"
    else:
        return "❌ Poor drug-likeness - significant issues"
