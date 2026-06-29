import csv
import urllib.request
import argparse
from io import StringIO
from Bio import SeqIO
from phascout.detection.phac_validator import PhaCValidator
import pyhmmer
from phascout.config import PHAC_CLASS_PROFILES

def get_uniprot_seq(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    req = urllib.request.Request(url, headers={'User-Agent': 'PHAscout'})
    try:
        resp = urllib.request.urlopen(req)
        data = resp.read().decode('utf-8')
        rec = next(SeqIO.parse(StringIO(data), "fasta"))
        return str(rec.seq)
    except Exception:
        return None

def predict_class(seq_str):
    alphabet = pyhmmer.easel.Alphabet.amino()
    digital_seq = pyhmmer.easel.TextSequence(name=b"query", sequence=seq_str).digitize(alphabet)
    
    best_score = 0
    best_class = "Unknown"
    
    for cls, path in PHAC_CLASS_PROFILES.items():
        try:
            with pyhmmer.plan7.HMMFile(path) as f:
                hmm = next(f)
            hits = list(pyhmmer.hmmsearch([hmm], [digital_seq]))
            if hits and len(hits[0]) > 0:
                score = hits[0][0].score
                if score > best_score and score > 20: # GA cutoff proxy
                    best_score = score
                    best_class = cls
        except Exception:
            pass
            
    return best_class

def run_benchmark(csv_path):
    strains = []
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            strains.append(row)
            
    results = {"TP": 0, "FN": 0}
    class_results = {}
    validator = PhaCValidator()
    
    print(f"Starting Benchmark. Total Positives: {len(strains)}")
    
    for strain in strains:
        acc = strain["Assembly_Accession"] # actually UniProt ID
        exp_class = strain["Expected_Class"]
        
        seq = get_uniprot_seq(acc)
        if not seq:
            print(f"  Failed to fetch {acc}. Skipping.")
            continue
            
        pred_class = predict_class(seq)
        if pred_class == "Unknown":
            results["FN"] += 1
            print(f"  [FN] {acc} ({exp_class}) - HMM did not match.")
            continue
            
        # Validate Triad
        val_res = validator.validate_triad_hmm(seq, pred_class)
        if val_res["is_functional"]:
            results["TP"] += 1
            if exp_class not in class_results:
                class_results[exp_class] = {"correct": 0, "total": 0}
            class_results[exp_class]["total"] += 1
            if pred_class == exp_class:
                class_results[exp_class]["correct"] += 1
            else:
                print(f"  [MISCLASSIFIED] {acc} predicted {pred_class}, expected {exp_class}")
        else:
            results["FN"] += 1
            print(f"  [FN] {acc} ({exp_class}) - Found class {pred_class} but missing functional triad/box.")
            
    print("\n=== CLASS-BAZINDA PERFORMANS (Sadece Fonksiyonel Olarak Onaylananlar) ===")
    for cls, r in class_results.items():
        n = r["total"]
        acc_pct = r["correct"] / n * 100 if n > 0 else 0
        z = 1.96
        p = r["correct"] / n if n > 0 else 0
        ci = z * (p * (1-p) / n) ** 0.5 if n > 0 else 0
        print(f"  {cls}: {r['correct']}/{n} ({acc_pct:.1f}%) +- {ci*100:.1f}% (95% CI)")
        
    print("\n=== OVERALL METRICS ===")
    print(f"TP: {results['TP']}, FN: {results['FN']}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", default="data/benchmark/independent_benchmark_set.csv")
    args = parser.parse_args()
    
    run_benchmark(args.csv)
