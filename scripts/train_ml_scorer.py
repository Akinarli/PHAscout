import os
import joblib
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from Bio import SeqIO
from phascout.prediction.physicochemical import analyze_physicochemical

# Egitim verileri yollari
BASE_DIR = r"C:\Users\bird-\OneDrive\Desktop\In silico\PHAscout"
FASTA_DIR = os.path.join(BASE_DIR, "data", "reference_sequences", "phac_classes_clean")
MODEL_PATH = os.path.join(BASE_DIR, "data", "rf_model.pkl")

def get_dummy_operon_distance(class_name):
    # Egitim setinde gff destegi simdilik olmadigi icin, 
    # sinifa gore biyolojik beklenen mesafeleri (mock feature) ekliyoruz.
    if class_name == "Class_I":
        return 100 # Operon
    return 10000 # Daginis

def get_dummy_hmm_score(class_name):
    return 500.0 # Yüksek HMM skoru

def extract_features_from_fasta(fasta_path, class_label):
    features = []
    if not os.path.exists(fasta_path):
        return features
        
    records = list(SeqIO.parse(fasta_path, "fasta"))
    for rec in records:
        phys = analyze_physicochemical(str(rec.seq))
        
        feat = {
            "hmm_score": get_dummy_hmm_score(class_label),
            "mw": phys["molecular_weight"],
            "pi": phys["isoelectric_point"],
            "gravy": phys["gravy"],
            "operon_dist": get_dummy_operon_distance(class_label),
            "label": 1 # Bunlar PHA ureticisi (Pozitif Ornekler)
        }
        features.append(feat)
    return features

def generate_negative_features(n_samples=100):
    # Sentezlemeyen enzimler icin (Lipase vs.) sahte metrikler (ornek sentetik veriler)
    # Egitim amaciyla
    import numpy as np
    features = []
    for _ in range(n_samples):
        feat = {
            "hmm_score": np.random.uniform(0, 100),
            "mw": np.random.uniform(30000, 50000),
            "pi": np.random.uniform(4.0, 9.0),
            "gravy": np.random.uniform(0.1, 1.5), # Daha hidrofobik
            "operon_dist": 10000,
            "label": 0
        }
        features.append(feat)
    return features

def main():
    print("Ozellikler (Features) Cikariliyor...")
    data = []
    classes = {
        "phac_class_I.fasta": "Class_I",
        "phac_class_II.fasta": "Class_II",
        "phac_class_III.fasta": "Class_III",
        "phac_class_IV.fasta": "Class_IV"
    }
    
    for f_name, cls in classes.items():
        f_path = os.path.join(FASTA_DIR, f_name)
        data.extend(extract_features_from_fasta(f_path, cls))
        
    # Negatif ornekler ekle
    data.extend(generate_negative_features(len(data)))
    
    df = pd.DataFrame(data)
    if df.empty:
        print("Veri bulunamadi!")
        return
        
    X = df[["hmm_score", "mw", "pi", "gravy", "operon_dist"]]
    y = df["label"]
    
    print(f"Toplam Veri Sayisi: {len(df)} ({sum(y==1)} Pozitif, {sum(y==0)} Negatif)")
    
    print("Random Forest Modeli Egitiliyor...")
    model = RandomForestClassifier(n_estimators=100, random_state=42)
    model.fit(X, y)
    
    joblib.dump(model, MODEL_PATH)
    print(f"Model basariyla kaydedildi: {MODEL_PATH}")

if __name__ == "__main__":
    main()
