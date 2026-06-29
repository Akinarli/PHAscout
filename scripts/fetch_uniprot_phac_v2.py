import urllib.request
import json
import csv
import time
from pathlib import Path
from Bio import Entrez

# Configuration
Entrez.email = "phascout@example.com"
OUTPUT_FASTA = Path("data/phac_uniprot_all.fasta")
OUTPUT_CSV = Path("data/phac_uniprot_metadata.csv")

def search_uniprot():
    # We query for phaC gene OR specific protein names.
    # We sort by annotation score to get the highest confidence TrEMBL and reviewed first.
    url = "https://rest.uniprot.org/uniprotkb/search?query=(gene:phaC+OR+protein_name:%22hydroxyalkanoic%20acid%20synthase%22)&format=json&size=500&sort=annotation_score+desc"
    print(f"Fetching from UniProt: {url}")
    
    req = urllib.request.Request(url, headers={'User-Agent': 'PHAscout'})
    try:
        response = urllib.request.urlopen(req)
        data = json.loads(response.read().decode('utf-8'))
        return data.get('results', [])
    except Exception as e:
        print("Error fetching from UniProt:", e)
        return []

def assign_class(entry):
    """Attempt to assign PhaC class based on keywords, taxonomy, and protein features."""
    tax_lineage = [t.get("scientificName", "") for t in entry.get("lineages", [])]
    if "organism" in entry:
        tax_lineage.append(entry["organism"].get("scientificName", ""))
    tax_str = " ".join(tax_lineage).lower()
    
    comments = entry.get("comments", [])
    function_text = ""
    for c in comments:
        if c.get("commentType") == "FUNCTION":
            function_text += c.get("texts", [{}])[0].get("value", "").lower()

    # Rule-based classification
    if "class iv" in function_text or "bacillus" in tax_str:
        return "Class_IV"
    elif "class iii" in function_text or "allochromatium" in tax_str or "thiocapsa" in tax_str or "synechocystis" in tax_str:
        return "Class_III"
    elif "class ii" in function_text or "pseudomonas" in tax_str or "medium-chain" in function_text or "mcl-pha" in function_text:
        # P. aeruginosa, P. putida, etc. are Class II.
        return "Class_II"
    elif "class i" in function_text or "short-chain" in function_text or "scl-pha" in function_text or "cupriavidus" in tax_str or "ralstonia" in tax_str or "azotobacter" in tax_str or "burkholderia" in tax_str or "chromobacterium" in tax_str or "aeromonas" in tax_str:
        return "Class_I" # Aeromonas is sometimes Class III but often its phaC is grouped with Class I functionally in terms of substrate SCL/MCL hybrid. We'll label Class I.
    
    return "Unknown"

def main():
    results = search_uniprot()
    print(f"Downloaded {len(results)} records.")
    
    records_to_save = []
    
    for r in results:
        uniprot_id = r["primaryAccession"]
        reviewed = r["entryType"] == "UniProtKB reviewed (Swiss-Prot)"
        sequence = r.get("sequence", {}).get("value", "")
        organism = r.get("organism", {}).get("scientificName", "")
        
        if not sequence:
            continue
            
        phac_class = assign_class(r)
        
        records_to_save.append({
            "uniprot_id": uniprot_id,
            "reviewed": reviewed,
            "organism": organism,
            "phac_class": phac_class,
            "sequence": sequence
        })
        
    # Write FASTA
    OUTPUT_FASTA.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_FASTA, "w") as f:
        for rec in records_to_save:
            f.write(f">{rec['uniprot_id']} {rec['phac_class']} {rec['organism']}\n")
            f.write(f"{rec['sequence']}\n")
            
    # Write CSV
    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["uniprot_id", "reviewed", "organism", "phac_class"])
        writer.writeheader()
        for rec in records_to_save:
            writer.writerow({
                "uniprot_id": rec["uniprot_id"],
                "reviewed": rec["reviewed"],
                "organism": rec["organism"],
                "phac_class": rec["phac_class"]
            })
            
    # Summary
    class_counts = {"Class_I": 0, "Class_II": 0, "Class_III": 0, "Class_IV": 0, "Unknown": 0}
    for rec in records_to_save:
        class_counts[rec["phac_class"]] += 1
        
    print(f"Saved {len(records_to_save)} sequences.")
    print("Class Distribution:", class_counts)

if __name__ == "__main__":
    main()
