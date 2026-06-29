import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

def calculate_identity(seq1, seq2, aligner):
    score = aligner.score(seq1, seq2)
    length = min(len(seq1), len(seq2))
    return score / length if length > 0 else 0.0

def main():
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 1
    aligner.mismatch_score = -2
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    in_fasta = Path("data/phac_uniprot_all.fasta")
    out_rep = Path("data/phac_clustered_80.fasta")
    out_clstr = Path("data/phac_clustered_80.fasta.clstr")

    records = list(SeqIO.parse(in_fasta, "fasta"))
    
    # Sort by length descending (like cd-hit)
    records.sort(key=lambda r: len(r.seq), reverse=True)
    
    clusters = []
    
    print(f"Clustering {len(records)} sequences at 80% identity...")
    
    for i, rec in enumerate(records):
        seq = str(rec.seq).replace("*", "").replace("X", "")
        if not seq:
            continue
            
        found_cluster = False
        for cluster in clusters:
            rep_seq = cluster["rep_seq"]
            ident = calculate_identity(rep_seq, seq, aligner)
            if ident >= 0.80:
                cluster["members"].append(rec)
                found_cluster = True
                break
                
        if not found_cluster:
            clusters.append({
                "rep": rec,
                "rep_seq": seq,
                "members": []
            })
            
        if (i + 1) % 50 == 0:
            print(f"Processed {i + 1}/{len(records)}...")

    print(f"Formed {len(clusters)} clusters.")
    
    # Save representatives
    out_rep.parent.mkdir(parents=True, exist_ok=True)
    with open(out_rep, "w") as f:
        for cluster in clusters:
            f.write(f">{cluster['rep'].description}\n{cluster['rep'].seq}\n")
            
    # Save .clstr file in format similar to cd-hit
    with open(out_clstr, "w") as f:
        for i, cluster in enumerate(clusters):
            f.write(f">Cluster {i}\n")
            f.write(f"0\t{len(cluster['rep'].seq)}aa, >{cluster['rep'].id}... *\n")
            for j, member in enumerate(cluster['members']):
                f.write(f"{j+1}\t{len(member.seq)}aa, >{member.id}... \n")

if __name__ == "__main__":
    main()
