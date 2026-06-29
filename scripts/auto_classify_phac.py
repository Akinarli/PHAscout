import pyhmmer
from Bio import SeqIO
from pathlib import Path

def main():
    hmms_dir = Path("data/hmm_profiles/phac_classes")
    fasta_in = Path("data/phac_uniprot_all.fasta")
    
    # Load HMMs
    hmms = {}
    with pyhmmer.plan7.HMMFile(hmms_dir / "phac_class_I.hmm") as f: hmms["Class_I"] = next(f)
    with pyhmmer.plan7.HMMFile(hmms_dir / "phac_class_II.hmm") as f: hmms["Class_II"] = next(f)
    with pyhmmer.plan7.HMMFile(hmms_dir / "phac_class_III.hmm") as f: hmms["Class_III"] = next(f)
    # We don't have Class IV HMM yet, but we already parsed 103 Class IVs via 'bacillus' keyword.
    
    records = list(SeqIO.parse(fasta_in, "fasta"))
    alphabet = pyhmmer.easel.Alphabet.amino()
    
    classified_records = {"Class_I": [], "Class_II": [], "Class_III": [], "Class_IV": [], "Unknown": []}
    
    for rec in records:
        desc = rec.description
        current_class = desc.split(" ")[1] if len(desc.split(" ")) > 1 else "Unknown"
        
        if current_class != "Unknown":
            classified_records[current_class].append(rec)
            continue
            
        # Run HMM
        seq_str = str(rec.seq).replace("*", "").replace("X", "")
        if not seq_str: continue
        
        digital_seq = pyhmmer.easel.TextSequence(name=rec.id.encode(), sequence=seq_str).digitize(alphabet)
        
        best_class = "Unknown"
        best_score = 0
        
        for cls_name, hmm in hmms.items():
            hits = list(pyhmmer.hmmsearch([hmm], [digital_seq]))
            if hits and len(hits[0]) > 0:
                hit = hits[0][0]
                if hit.score > best_score and hit.score > 100: # Threshold for homology
                    best_score = hit.score
                    best_class = cls_name
                    
        classified_records[best_class].append(rec)
        
    for cls_name, recs in classified_records.items():
        print(f"{cls_name}: {len(recs)}")
        # Save to full FASTA for CD-HIT
        if cls_name != "Unknown":
            out_path = Path(f"data/phac_{cls_name.lower()}_full.fasta")
            out_path.parent.mkdir(parents=True, exist_ok=True)
            SeqIO.write(recs, out_path, "fasta")

if __name__ == "__main__":
    main()
