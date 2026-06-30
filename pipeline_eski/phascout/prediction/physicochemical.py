"""
Physicochemical Analysis Module
===============================
Extracts physicochemical properties from a protein sequence to be used as ML features.
"""

from Bio.SeqUtils.ProtParam import ProteinAnalysis

def analyze_physicochemical(sequence: str) -> dict:
    """
    Calculates various physicochemical features of a protein sequence.
    
    Args:
        sequence: Amino acid sequence string.
        
    Returns:
        dict: Features including MW, pI, GRAVY, etc.
    """
    # Remove any masking characters like 'X' or '-' before analysis
    clean_seq = sequence.replace("X", "A").replace("-", "").replace("*", "")
    
    if not clean_seq:
        return {
            "molecular_weight": 0.0,
            "isoelectric_point": 0.0,
            "gravy": 0.0,
            "instability_index": 0.0,
            "aromaticity": 0.0
        }
        
    analysis = ProteinAnalysis(clean_seq)
    
    try:
        mw = analysis.molecular_weight()
    except ValueError:
        mw = 0.0
        
    try:
        pi = analysis.isoelectric_point()
    except ValueError:
        pi = 0.0
        
    try:
        gravy = analysis.gravy()
    except Exception:
        gravy = 0.0
        
    try:
        instability = analysis.instability_index()
    except Exception:
        instability = 0.0
        
    try:
        aromaticity = analysis.aromaticity()
    except Exception:
        aromaticity = 0.0

    return {
        "molecular_weight": mw,
        "isoelectric_point": pi,
        "gravy": gravy,
        "instability_index": instability,
        "aromaticity": aromaticity
    }
