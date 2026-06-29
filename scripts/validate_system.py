"""
PHAscout Validasyon Testi (Asama 7)
=====================================
Aşama 0'da hazirlanan known_strains.csv uzerinden
tum pipeline'in gercek dunya basarimini test eder.
Sensitivity, Specificity ve False Positive oranlarini hesaplar.
"""

import sys
import os
import pandas as pd
import logging
from phascout.pipeline import PHAscoutPipeline

logging.basicConfig(level=logging.WARNING)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV_PATH = os.path.join(BASE_DIR, "data", "known_strains.csv")


def main():
    print("=" * 60)
    print("PHAscout V11 FINAL VALIDASYON TESTI (Asama 7)")
    print("=" * 60)
    
    if not os.path.exists(CSV_PATH):
        print(f"HATA: Validasyon veri seti bulunamadi: {CSV_PATH}")
        sys.exit(1)
        
    df = pd.read_csv(CSV_PATH)
    pipeline = PHAscoutPipeline()
    
    results = []
    
    tp = 0  # True Positive
    tn = 0  # True Negative
    fp = 0  # False Positive
    fn = 0  # False Negative
    
    print(f"Toplam {len(df)} sus test ediliyor...\n")
    
    for idx, row in df.iterrows():
        strain = row['Strain_Name']
        acc = row['Assembly_Accession']
        expected_res = str(row['Expected_Result']).strip()
        expected_cls = str(row['Expected_Class']).strip() if pd.notna(row['Expected_Class']) else "None"
        
        print(f"[{idx+1}/{len(df)}] Test ediliyor: {strain} ({acc})")
        
        try:
            report = pipeline.run(accession=acc)
            is_producer = report["summary"]["produces_pha"]
            detected_cls = str(report["summary"]["phac_class"])
            
            actual_res = "Positive" if is_producer else "Negative"
            
            # Confusion Matrix
            if expected_res == "Positive" and actual_res == "Positive":
                tp += 1
                status = "[+] TP"
            elif expected_res == "Negative" and actual_res == "Negative":
                tn += 1
                status = "[+] TN"
            elif expected_res == "Negative" and actual_res == "Positive":
                fp += 1
                status = "[-] FP (YANLIS ALARM)"
            elif expected_res == "Positive" and actual_res == "Negative":
                fn += 1
                status = "[-] FN (KACIRILDI)"
                
            print(f"  -> Beklenen: {expected_res}, Bulunan: {actual_res} | {status}")
            if expected_res == "Positive":
                print(f"  -> Sinif Beklenen: {expected_cls}, Bulunan: {detected_cls}")
                
            results.append({
                "Strain": strain,
                "Accession": acc,
                "Expected": expected_res,
                "Predicted": actual_res,
                "Status": status
            })
            
        except Exception as e:
            print(f"  -> HATA OLUSTU: {e}")
            
    # Istatistikler
    total = tp + tn + fp + fn
    accuracy = (tp + tn) / total if total > 0 else 0
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    
    print("\n" + "=" * 60)
    print("FINAL TEST ISTATISTIKLERI")
    print("=" * 60)
    print(f"Dogruluk (Accuracy):    {accuracy:.2%}")
    print(f"Hassasiyet (Sens):      {sensitivity:.2%}")
    print(f"Spesifiklik (Spec):     {specificity:.2%}")
    print("-" * 60)
    print(f"True Positives (TP):    {tp}")
    print(f"True Negatives (TN):    {tn}")
    print(f"False Positives (FP):   {fp}  <- Bu sifir olmali!")
    print(f"False Negatives (FN):   {fn}")
    print("=" * 60)
    
    
if __name__ == "__main__":
    main()
