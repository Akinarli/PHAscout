"""
ML Skorlayici (Yardimci Guven)
==============================
Egitilmis bir Random Forest ile, dogrulanmis bir PhaC adayinin gercek
ozelliklerine (PhaC sinif HMM bit skoru + fizikokimyasal) dayanarak bir
yardimci guven olasiligi uretir.

Ozellikler ve uretim yontemi, scripts/train_ml_scorer.py ile BIREBIR
aynidir (egitim/cikarim tutarliligi). Model bulunamazsa seffaf, kural
tabanli bir yedek (fallback) skor kullanilir.

ONEMLI KISIT (dairesellik): 'hmm_score' ozelligi, ayni pozitif dizilerden
(phac_classes_clean) insa edilen sinif HMM'lerine karsi hesaplanir. Yani model
kismen "yuksek HMM skoru = PhaC"yi ezberler ve egitimdeki capraz-dogrulama
ROC-AUC/MCC degerleri IYIMSER'dir; genelleme performansi olarak ALINTILANMAMALI.
Bu skor yalnizca YARDIMCI bir guven sinyalidir ve deterministik (triad+box)
karari ASLA gecersiz kilmaz.
"""

import os
import joblib
import logging
import pyhmmer
from phascout.config import PHAC_CLASS_PROFILES
from phascout.prediction.physicochemical import analyze_physicochemical

logger = logging.getLogger(__name__)

FEATURE_ORDER = ["hmm_score", "mw", "pi", "gravy", "instability", "aromaticity"]


class MLScorer:
    def __init__(self, model_path="data/rf_model.pkl"):
        self.model_path = model_path
        self.model = None
        self.feature_order = FEATURE_ORDER
        if os.path.exists(model_path):
            try:
                obj = joblib.load(model_path)
                if isinstance(obj, dict):  # yeni format: {model, feature_order}
                    self.model = obj["model"]
                    self.feature_order = obj.get("feature_order", FEATURE_ORDER)
                else:  # geriye uyumluluk
                    self.model = obj
                logger.info(f"ML Modeli yuklendi: {model_path}")
            except Exception as e:
                logger.error(f"ML Modeli yuklenemedi: {e}")
        else:
            logger.warning(f"ML Modeli bulunamadi: {model_path}. Fallback modu kullanilacak.")

    def predict(self, phac_result, gene_vector, operon_result, phac_seq):
        """PHA uretim yardimci-guven olasiligi (0-100) hesaplar."""
        features = self._extract_features(phac_result, phac_seq)

        result = {
            "ml_probability": 0.0,
            "features_used": features,
            "model_loaded": self.model is not None,
        }

        if not phac_result.get("is_functional") or not phac_seq:
            return result

        if self.model is not None and features is not None:
            try:
                import pandas as pd
                df = pd.DataFrame([features])[self.feature_order]
                prob = self.model.predict_proba(df)[0][1]
                result["ml_probability"] = round(float(prob * 100), 1)
            except Exception as e:
                logger.error(f"Tahmin sirasinda hata: {e}")
        else:
            # Seffaf kural tabanli yedek (model yokken)
            prob = 50.0
            if phac_result.get("is_functional"):
                prob += 20.0
            if operon_result.get("is_class_i_operon"):
                prob += 15.0
            if features and features["gravy"] < 0:  # cozunur protein
                prob += 10.0
            result["ml_probability"] = round(min(prob, 95.0), 1)

        return result

    def _extract_features(self, phac_result, phac_seq):
        if not phac_seq:
            return None
        physico = analyze_physicochemical(phac_seq)
        return {
            "hmm_score": float(phac_result.get("best_score", 0.0)),
            "mw": physico["molecular_weight"],
            "pi": physico["isoelectric_point"],
            "gravy": physico["gravy"],
            "instability": physico["instability_index"],
            "aromaticity": physico["aromaticity"],
        }
