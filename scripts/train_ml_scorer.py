"""
ML Skorlayici Egitimi (GERCEK Ozelliklerle)
============================================
PhaC adaylarini, kendisiyle karistirilabilen alpha/beta-hidrolaz / SDR
proteinlerinden ayirt eden bir Random Forest yardimci-guven modeli egitir.

KRITIK: Tum ozellikler GERCEK dizilerden, cikarim aninda hesaplananlarla
AYNI yontemle uretilir. Sentetik (np.random) negatif YOK, sabit (dummy)
ozellik YOK. Boylece egitim ve cikarim ayni dagilimda calisir.

Ozellikler:
  - hmm_score : 4 PhaC sinif HMM'ine karsi en yuksek bit skoru (gercek hmmsearch)
  - mw, pi, gravy, instability, aromaticity : fizikokimyasal (ProtParam)

Pozitifler: data/reference_sequences/phac_classes_clean/phac_class_{I..IV}.fasta
Negatifler: data/reference_sequences/negative/*.fasta (FabG, FadA, SDR,
            crotonase, lipase/esteraz, epoksit hidrolaz - PhaC OLMAYAN gercek
            proteinler; ozellikle alpha/beta-hidrolazlar gercek karistiricidir)

Calistirma:
    python scripts/train_ml_scorer.py
"""

import os
import glob
import joblib
import numpy as np
import pandas as pd
import pyhmmer
from Bio import SeqIO
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import classification_report, roc_auc_score, matthews_corrcoef

from phascout.config import PHAC_CLASS_PROFILES
from phascout.prediction.physicochemical import analyze_physicochemical

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CLEAN_DIR = os.path.join(BASE_DIR, "data", "reference_sequences", "phac_classes_clean")
NEG_DIR = os.path.join(BASE_DIR, "data", "reference_sequences", "negative")
MODEL_PATH = os.path.join(BASE_DIR, "data", "rf_model.pkl")

FEATURE_ORDER = ["hmm_score", "mw", "pi", "gravy", "instability", "aromaticity"]

_ALPHABET = pyhmmer.easel.Alphabet.amino()
_BACKGROUND = pyhmmer.plan7.Background(_ALPHABET)
_CLASS_HMMS = []
for _cls, _path in PHAC_CLASS_PROFILES.items():
    if os.path.exists(_path):
        with pyhmmer.plan7.HMMFile(_path) as _f:
            _CLASS_HMMS.append(next(_f))


def best_phac_hmm_score(seq_str):
    """4 PhaC sinif HMM'ine karsi en yuksek bit skoru (gercek hmmsearch)."""
    clean = seq_str.replace("*", "").replace("X", "A")
    if not clean:
        return 0.0
    ds = pyhmmer.easel.TextSequence(name=b"q", sequence=clean).digitize(_ALPHABET)
    best = 0.0
    for hmm in _CLASS_HMMS:
        try:
            hits = list(pyhmmer.hmmsearch([hmm], [ds], background=_BACKGROUND))
            if hits and len(hits[0]) > 0:
                best = max(best, hits[0][0].score)
        except Exception:
            pass
    return float(best)


def features_for(seq_str):
    phys = analyze_physicochemical(seq_str)
    return {
        "hmm_score": best_phac_hmm_score(seq_str),
        "mw": phys["molecular_weight"],
        "pi": phys["isoelectric_point"],
        "gravy": phys["gravy"],
        "instability": phys["instability_index"],
        "aromaticity": phys["aromaticity"],
    }


def collect(files, label):
    rows = []
    for fp in files:
        for rec in SeqIO.parse(fp, "fasta"):
            s = str(rec.seq)
            if len(s) < 50:
                continue
            feat = features_for(s)
            feat["label"] = label
            rows.append(feat)
    return rows


def main():
    pos_files = [os.path.join(CLEAN_DIR, f"phac_class_{r}.fasta") for r in ["I", "II", "III", "IV"]]
    pos_files = [f for f in pos_files if os.path.exists(f)]
    neg_files = sorted(glob.glob(os.path.join(NEG_DIR, "*.fasta")))

    print("Pozitif PhaC dosyalari:", [os.path.basename(f) for f in pos_files])
    print("Negatif dosyalar:", [os.path.basename(f) for f in neg_files])
    print("Gercek ozellikler cikariliyor (hmmsearch + ProtParam)...")

    data = collect(pos_files, 1) + collect(neg_files, 0)
    df = pd.DataFrame(data).dropna()
    if df.empty:
        print("Veri bulunamadi!")
        return

    X = df[FEATURE_ORDER]
    y = df["label"]
    print(f"Toplam: {len(df)} ({int((y == 1).sum())} pozitif PhaC, {int((y == 0).sum())} negatif)")

    model = RandomForestClassifier(n_estimators=300, random_state=42, class_weight="balanced")

    # Durust performans: 5-katli capraz dogrulama (egitim setine bakmadan)
    cv_pred = cross_val_predict(model, X, y, cv=5)
    cv_prob = cross_val_predict(model, X, y, cv=5, method="predict_proba")[:, 1]
    print("\n=== 5-Katli Capraz Dogrulama (out-of-fold) ===")
    print(classification_report(y, cv_pred, target_names=["non-PhaC", "PhaC"]))
    print(f"ROC-AUC: {roc_auc_score(y, cv_prob):.4f}")
    print(f"MCC:     {matthews_corrcoef(y, cv_pred):.4f}")

    # Son model tum veriyle egitilir ve kaydedilir (feature sirasiyla birlikte)
    model.fit(X, y)
    print("\nOzellik onemleri:")
    for f, imp in sorted(zip(FEATURE_ORDER, model.feature_importances_), key=lambda x: -x[1]):
        print(f"  {f:14} {imp:.3f}")

    joblib.dump({"model": model, "feature_order": FEATURE_ORDER}, MODEL_PATH)
    print(f"\nModel kaydedildi: {MODEL_PATH}")


if __name__ == "__main__":
    main()
