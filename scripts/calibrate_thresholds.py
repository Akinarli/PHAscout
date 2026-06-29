import os
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import numpy as np
try:
    from sklearn.metrics import roc_curve, f1_score
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    print("scikit-learn bulunamadı. F1 skorlaması manuel yapılacak.")

def get_aligner():
    aligner = PairwiseAligner()
    if hasattr(aligner, 'substitution_matrices'):
        aligner.substitution_matrix = aligner.substitution_matrices.load("BLOSUM62")
    else:
        from Bio.Align import substitution_matrices
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.mode = 'global'
    return aligner

aligner = get_aligner()

def normalized_score(seq1, seq2):
    s12 = aligner.score(seq1, seq2)
    s11 = aligner.score(seq1, seq1)
    if s11 == 0:
        return 0
    return s12 / s11

def max_score_against_refs(query_seq, references):
    if not references:
        return 0
    return max([normalized_score(str(ref.seq), query_seq) for ref in references])

def find_optimal_threshold(pos_file, neg_files, target_name):
    print(f"\n[{target_name}] Kalibrasyonu Başlıyor...")
    
    if not os.path.exists(pos_file):
        print(f"Hata: Pozitif set bulunamadı ({pos_file})")
        return None
        
    pos_seqs = list(SeqIO.parse(pos_file, "fasta"))
    if not pos_seqs:
        print("Hata: Pozitif set boş.")
        return None
        
    print(f"Referans dizileri: {len(pos_seqs)} adet")
    
    pos_scores = []
    for rec in pos_seqs:
        # Note: A positive sequence against ITSELF will score 1.0. 
        # For a rigorous cross-validation, we would exclude the sequence itself.
        # But for establishing a general acceptance threshold from references, 
        # finding the lowest score that captures the most positives while excluding negatives is okay.
        # Let's compare against ALL OTHERS to be strictly rigorous:
        others = [r for r in pos_seqs if r.id != rec.id]
        if others:
            pos_scores.append(max_score_against_refs(str(rec.seq), others))
        else:
            pos_scores.append(1.0)
            
    neg_scores = []
    for neg_file in neg_files:
        if os.path.exists(neg_file):
            neg_seqs = list(SeqIO.parse(neg_file, "fasta"))
            print(f"Negatif diziler ({os.path.basename(neg_file)}): {len(neg_seqs)} adet")
            for rec in neg_seqs:
                neg_scores.append(max_score_against_refs(str(rec.seq), pos_seqs))
        else:
            print(f"Uyarı: Negatif dosya eksik: {neg_file}")
            
    if not neg_scores:
        print("Hata: Hiç negatif dizi bulunamadı, eşik hesaplanamaz.")
        return None
        
    all_scores = pos_scores + neg_scores
    labels = [1]*len(pos_scores) + [0]*len(neg_scores)
    
    if HAS_SKLEARN:
        fpr, tpr, thresholds = roc_curve(labels, all_scores)
        f1_scores = [f1_score(labels, [1 if s >= t else 0 for s in all_scores]) for t in thresholds]
        best_idx = np.argmax(f1_scores)
        best_threshold = thresholds[best_idx]
        best_f1 = f1_scores[best_idx]
    else:
        # Manual F1 calculation
        thresholds = sorted(list(set(all_scores)))
        best_threshold = 0
        best_f1 = 0
        for t in thresholds:
            tp = sum(1 for s, l in zip(all_scores, labels) if s >= t and l == 1)
            fp = sum(1 for s, l in zip(all_scores, labels) if s >= t and l == 0)
            fn = sum(1 for s, l in zip(all_scores, labels) if s < t and l == 1)
            
            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
            
            if f1 > best_f1:
                best_f1 = f1
                best_threshold = t
                
    print(f"==> {target_name} Optimal Eşik: {best_threshold:.4f} (F1: {best_f1:.4f})")
    return best_threshold

base_dir = r"C:\Users\bird-\OneDrive\Desktop\In silico\PHAscout\data\reference_sequences"
pos_dir = os.path.join(base_dir, "positive")
neg_dir = os.path.join(base_dir, "negative")

thresholds = {}

t_phab = find_optimal_threshold(
    os.path.join(pos_dir, "phab.fasta"),
    [os.path.join(neg_dir, "fabg.fasta"), os.path.join(neg_dir, "sdr_broad.fasta")],
    "PhaB"
)
if t_phab: thresholds['PHAB_THRESHOLD'] = t_phab

t_phaa = find_optimal_threshold(
    os.path.join(pos_dir, "phaa.fasta"),
    [os.path.join(neg_dir, "fada.fasta")],
    "PhaA"
)
if t_phaa: thresholds['PHAA_THRESHOLD'] = t_phaa

print("\nBulunan Eşik Değerleri (config.py için):")
for k, v in thresholds.items():
    print(f"{k} = {v:.4f}")
