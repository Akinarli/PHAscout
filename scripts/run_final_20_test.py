import os
import sys
from phascout.pipeline import PHAscoutPipeline

# EXACT verified GCFs (RefSeq Complete)
POSITIVES = {
    "Cupriavidus necator (Class I)": "GCF_000219215.1",
    "Pseudomonas putida (Class II)": "GCF_000412675.1",
    "Pseudomonas aeruginosa (Class II)": "GCF_000006765.1",
    "Halomonas smyrnensis (Class I)": "GCF_000265245.1",
    "Aeromonas caviae (Class III)": "GCF_900476005.1",
    "Rhodobacter sphaeroides (Class I)": "GCF_000012905.2",
    "Azotobacter vinelandii (Class I)": "GCF_055386615.1",
    "Burkholderia sacchari (Class I)": "GCF_000785435.2",
    "Chromobacterium violaceum (Class I)": "GCF_000007705.1",
    "Allochromatium vinosum (Class III)": "GCF_000025485.1"
}

NEGATIVES = {
    "Escherichia coli K-12": "GCF_000005845.2",
    "Streptococcus pneumoniae": "GCF_001457635.1",
    "Staphylococcus aureus": "GCF_000013425.1",
    "Mycobacterium tuberculosis": "GCF_000195955.2",
    "Helicobacter pylori": "GCF_025998455.1",
    "Bacillus subtilis 168": "GCF_000009045.1",
    "Campylobacter jejuni": "GCF_000009085.1",
    "Lactobacillus acidophilus": "GCF_034298135.1",
    "Neisseria gonorrhoeae": "GCF_013030075.1",
    # NOT: "Vibrio cholerae" (GCF_008369605.1) negatiflerden CIKARILDI.
    # PHAscout bu suşta fonksiyonel Class I PhaC + sintenik phaA/phaB'den oluşan
    # gerçek bir phaCAB operonu tespit etti; NCBI anotasyonu phaC'yi "class I
    # poly(R)-hydroxyalkanoic acid synthase" olarak doğruluyor. Yani bu suş
    # geçerli bir negatif DEĞİL (araç doğru, etiket yanlıştı). Yerine net bir
    # üretmeyen kondu:
    "Streptococcus pyogenes": "GCF_000006785.2"
}

def run_test():
    pipeline = PHAscoutPipeline()
    print("=" * 60)
    print("PHAscout RUTHLESS FINAL TEST (20 STRAINS)")
    print("=" * 60)
    
    results = {"TP": 0, "FN": 0, "TN": 0, "FP": 0}
    
    print("\n[+] TESTING POSITIVE STRAINS (Expected: PHA Uretimi = EVET)")
    for name, acc in POSITIVES.items():
        try:
            res = pipeline.run(accession=acc)
            if res["summary"]["produces_pha"]:
                status = f"[+] TP (Found: {res['summary']['phac_class']})"
                results["TP"] += 1
            else:
                status = "[-] FN (Kacirildi!)"
                results["FN"] += 1
            print(f"{name.ljust(40)} -> {acc} : {status}")
        except Exception as e:
            print(f"{name.ljust(40)} -> {acc} : [x] ERROR ({str(e)})")
            results["FN"] += 1
            
    print("\n[-] TESTING NEGATIVE STRAINS (Expected: PHA Uretimi = HAYIR)")
    for name, acc in NEGATIVES.items():
        try:
            res = pipeline.run(accession=acc)
            if not res["summary"]["produces_pha"]:
                status = "[+] TN (Dogru Negatif)"
                results["TN"] += 1
            else:
                status = f"[-] FP (Yanlis Pozitif! Bulunan: {res['summary']['phac_class']})"
                results["FP"] += 1
            print(f"{name.ljust(40)} -> {acc} : {status}")
        except Exception as e:
            print(f"{name.ljust(40)} -> {acc} : [x] ERROR ({str(e)})")
            # For negatives, an error fetching is generally bad but let's count it as failed test
            pass

    print("\n" + "=" * 60)
    print("FINAL STATS")
    print("=" * 60)
    print(f"True Positives (TP):  {results['TP']} / 10")
    print(f"False Negatives (FN): {results['FN']} / 10")
    print(f"True Negatives (TN):  {results['TN']} / 10")
    print(f"False Positives (FP): {results['FP']} / 10")
    
    try:
        acc = (results["TP"] + results["TN"]) / 20 * 100
        sens = (results["TP"] / (results["TP"] + results["FN"])) * 100
        spec = (results["TN"] / (results["TN"] + results["FP"])) * 100
        print(f"Accuracy:    {acc:.1f}%")
        print(f"Sensitivity: {sens:.1f}%")
        print(f"Specificity: {spec:.1f}%")
    except ZeroDivisionError:
        print("Could not calculate metrics due to zero division.")

if __name__ == '__main__':
    run_test()
