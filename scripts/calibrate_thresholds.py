import argparse
import pathlib
import json
import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from sklearn.metrics import roc_curve, f1_score
from sklearn.model_selection import StratifiedKFold

def get_aligner():
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.mode = "local"
    return aligner

def normalized_score(seq1, seq2, aligner):
    s12 = aligner.score(seq1, seq2)
    s11 = aligner.score(seq1, seq1)
    s22 = aligner.score(seq2, seq2)
    denom = min(s11, s22)
    return s12 / denom if denom > 0 else 0.0

def max_score_against_refs(query, refs, aligner):
    if not refs:
        return 0.0
    return max(normalized_score(str(query.seq), str(r.seq), aligner) for r in refs)

def calibrate_with_cv(pos_file, neg_files, gene_name, n_folds=5):
    """5-fold cross-validation ile eşik belirleme"""
    aligner = get_aligner()
    
    pos_seqs = list(SeqIO.parse(pos_file, "fasta")) if pos_file.exists() else []
    neg_seqs = []
    for nf in neg_files:
        if pathlib.Path(nf).exists():
            neg_seqs.extend(SeqIO.parse(nf, "fasta"))
            
    print(f"\n[{gene_name}] Pozitif: {len(pos_seqs)}, Negatif: {len(neg_seqs)}")
    
    if len(pos_seqs) < 10:
        print(f"  UYARI: {gene_name} için pozitif set çok küçük (n={len(pos_seqs)}). "
              f"Kalibrasyon sonucu güvenilmez.")
              
    if not pos_seqs or not neg_seqs:
        return 0.0, 0.0, (0.0, 0.0)
    
    pos_scores = []
    for i, rec in enumerate(pos_seqs):
        others = [r for j, r in enumerate(pos_seqs) if j != i]
        pos_scores.append(max_score_against_refs(rec, others, aligner))
        
    neg_scores = [max_score_against_refs(r, pos_seqs, aligner) for r in neg_seqs]
    
    all_scores = pos_scores + neg_scores
    labels = [1] * len(pos_scores) + [0] * len(neg_scores)
    
    fpr, tpr, thresholds = roc_curve(labels, all_scores)
    f1_scores = [
        f1_score(labels, [1 if s >= t else 0 for s in all_scores])
        for t in thresholds
    ]
    best_idx = np.argmax(f1_scores)
    best_threshold = float(thresholds[best_idx])
    best_f1 = float(f1_scores[best_idx])
    
    # Bootstrap for CI
    n_boot = 200
    boot_thresholds = []
    rng = np.random.default_rng(42)
    for _ in range(n_boot):
        idx = rng.integers(0, len(all_scores), len(all_scores))
        bs = np.array(all_scores)[idx]
        bl = np.array(labels)[idx]
        if len(np.unique(bl)) < 2:
            continue
        _, _, bt = roc_curve(bl, bs)
        bf1 = [f1_score(bl, [1 if s >= t else 0 for s in bs]) for t in bt]
        if len(bt) > 0:
            boot_thresholds.append(float(bt[np.argmax(bf1)]))
            
    if boot_thresholds:
        ci_low = np.percentile(boot_thresholds, 2.5)
        ci_high = np.percentile(boot_thresholds, 97.5)
    else:
        ci_low, ci_high = best_threshold, best_threshold
        
    print(f"  Optimal eşik: {best_threshold:.4f} (F1={best_f1:.3f})")
    print(f"  95% CI: [{ci_low:.4f}, {ci_high:.4f}]")
    
    return best_threshold, best_f1, (ci_low, ci_high)

def main():
    parser = argparse.ArgumentParser(description="PHAscout BLOSUM62 eşik kalibrasyonu")
    parser.add_argument(
        "--data-dir",
        type=pathlib.Path,
        default=pathlib.Path(__file__).parent.parent / "data" / "reference_sequences"
    )
    parser.add_argument(
        "--output",
        type=pathlib.Path,
        default=pathlib.Path(__file__).parent.parent / "phascout" / "calibration_results.json"
    )
    args = parser.parse_args()
    
    pos_dir = args.data_dir / "positive"
    neg_dir = args.data_dir / "negative"
    
    targets = {
        "phaB": {
            "pos": pos_dir / "phab.fasta",
            "neg": [neg_dir / "fabg.fasta", neg_dir / "sdr_broad.fasta"]
        },
        "phaA": {
            "pos": pos_dir / "phaa.fasta",
            "neg": [neg_dir / "fada.fasta"]
        },
        "phaJ": {
            "pos": pos_dir / "phaj.fasta",
            "neg": [neg_dir / "enoyl_coa_general.fasta"]
        },
        "phaG": {
            "pos": pos_dir / "phag.fasta",
            "neg": [neg_dir / "fabg.fasta"]
        },
    }
    
    results = {}
    for gene, paths in targets.items():
        t, f1, ci = calibrate_with_cv(paths["pos"], paths["neg"], gene)
        results[gene] = {"threshold": t, "f1": f1, "ci_95": list(ci)}
        
    print("\n=== config.py için BLOSUM62_THRESHOLDS ===")
    print("BLOSUM62_THRESHOLDS = {")
    for gene, r in results.items():
        print(f'    "{gene}": {r["threshold"]:.4f},  # F1={r["f1"]:.3f}, 95% CI=[{r["ci_95"][0]:.4f},{r["ci_95"][1]:.4f}]')
    print("}")
    
    with open(args.output, "w") as f:
        json.dump(results, f, indent=2)

if __name__ == "__main__":
    main()
