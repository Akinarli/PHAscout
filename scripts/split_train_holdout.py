import re
import csv
from pathlib import Path
from Bio import SeqIO

def parse_cdhit_clusters(clstr_file):
    train_ids = set()
    holdout_ids = set()
    current_cluster = []
    with open(clstr_file) as f:
        for line in f:
            if line.startswith(">Cluster"):
                if current_cluster:
                    rep = next((x for x in current_cluster if x["is_rep"]), None)
                    if rep:
                        train_ids.add(rep["id"])
                        for m in current_cluster:
                            if not m["is_rep"]:
                                holdout_ids.add(m["id"])
                current_cluster = []
            else:
                is_rep = line.strip().endswith("*")
                id_match = re.search(r">(\S+)\.\.\.", line)
                if id_match:
                    current_cluster.append({"id": id_match.group(1), "is_rep": is_rep})
        if current_cluster:
            rep = next((x for x in current_cluster if x["is_rep"]), None)
            if rep:
                train_ids.add(rep["id"])
                for m in current_cluster:
                    if not m["is_rep"]:
                        holdout_ids.add(m["id"])
    return train_ids, holdout_ids

def main():
    clstr_file = Path("data/phac_clustered_80.fasta.clstr")
    meta_file = Path("data/phac_uniprot_metadata.csv")
    
    train_ids, holdout_ids = parse_cdhit_clusters(clstr_file)
    print(f"Train (Representatives): {len(train_ids)}")
    print(f"Holdout (Members): {len(holdout_ids)}")
    
    metadata = {}
    with open(meta_file, encoding='utf-8') as f:
        for row in csv.DictReader(f):
            metadata[row["uniprot_id"]] = row
            
    classes = ["Class_I", "Class_II", "Class_III", "Class_IV"]
    train_recs = {c: [] for c in classes}
    holdout_list = []
    
    for cls in classes:
        full_fasta = Path(f"data/phac_{cls.lower()}_full.fasta")
        if not full_fasta.exists(): continue
        
        for rec in SeqIO.parse(full_fasta, "fasta"):
            uniprot_id = rec.id
            meta = metadata.get(uniprot_id, {})
            
            if uniprot_id in train_ids:
                train_recs[cls].append(rec)
            elif uniprot_id in holdout_ids:
                holdout_list.append({
                    "Strain_Name": meta.get("organism", "Unknown Organism"),
                    "Assembly_Accession": uniprot_id,
                    "Expected_Result": "Positive",
                    "Expected_Class": cls,
                    "Source": "UniProt_Holdout",
                    "In_HMM_Build": "No"
                })
                
    for cls, recs in train_recs.items():
        if recs:
            out_path = Path(f"data/reference_sequences/phac_classes_clean/phac_{cls.lower()}.fasta")
            out_path.parent.mkdir(parents=True, exist_ok=True)
            SeqIO.write(recs, out_path, "fasta")
            print(f"Saved {len(recs)} train records for {cls}")
            
    csv_out = Path("data/benchmark/independent_benchmark_set.csv")
    csv_out.parent.mkdir(parents=True, exist_ok=True)
    with open(csv_out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["Strain_Name", "Assembly_Accession", "Expected_Result", "Expected_Class", "Source", "In_HMM_Build"])
        writer.writeheader()
        writer.writerows(holdout_list)
        
    print(f"Saved {len(holdout_list)} hold-out records to benchmark CSV.")

if __name__ == "__main__":
    main()
