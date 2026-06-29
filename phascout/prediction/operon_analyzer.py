"""
Operon Analyzer Module
======================
Determines if PHA genes (phaC, phaA, phaB, etc.) are co-located within a single operon 
by analyzing their genomic coordinates from GFF3 data.
"""

def analyze_operon(detected_genes: dict, gff_data: dict, max_distance_bp=3000) -> dict:
    """
    Analyzes genomic distances between detected PHA genes.
    
    Args:
        detected_genes: Dictionary of detected genes from hmm_scanner.
                        Example: {"phaC": {"detected": True, "protein_id": "WP_123"}, ...}
        gff_data: GFF3 coordinate dictionary from ncbi_datasets.
        max_distance_bp: Maximum allowed distance between genes to be considered an operon.
        
    Returns:
        dict: Operon analysis results containing distances and boolean flags.
    """
    results = {
        "is_class_i_operon": False,
        "distances": {}
    }
    
    if not gff_data:
        return results
        
    def get_coords(gene_name):
        gene_info = detected_genes.get(gene_name, {})
        if gene_info.get("detected"):
            prot_id = gene_info.get("protein_id")
            if prot_id and prot_id in gff_data:
                return gff_data[prot_id]
        return None
        
    phac_coords = get_coords("phaC")
    phaa_coords = get_coords("phaA")
    phab_coords = get_coords("phaB")
    
    if not phac_coords:
        return results
        
    # Function to calculate min distance between two genomic features
    def calc_dist(c1, c2):
        if not c1 or not c2:
            return None
        if c1["contig"] != c2["contig"]:
            return float('inf') # Different contigs
            
        # Distance is distance between intervals [start1, end1] and [start2, end2]
        if c1["end"] < c2["start"]:
            return c2["start"] - c1["end"]
        elif c2["end"] < c1["start"]:
            return c1["start"] - c2["end"]
        else:
            return 0 # Overlapping
            
    dist_c_a = calc_dist(phac_coords, phaa_coords)
    dist_c_b = calc_dist(phac_coords, phab_coords)
    
    if dist_c_a is not None:
        results["distances"]["phaC-phaA"] = dist_c_a
    if dist_c_b is not None:
        results["distances"]["phaC-phaB"] = dist_c_b
        
    # Operon Logic for Class I: phaC is physically close to phaA and/or phaB
    # Cupriavidus necator H16 has phaCAB operon.
    # If phaC is within max_distance of phaA OR phaB, we consider it an operon context.
    
    operon_found = False
    if dist_c_a is not None and dist_c_a <= max_distance_bp:
        operon_found = True
    if dist_c_b is not None and dist_c_b <= max_distance_bp:
        operon_found = True
        
    results["is_class_i_operon"] = operon_found
    
    return results
