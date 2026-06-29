"""
Machine Learning Scorer
=======================
Replaces the old heuristic logic. Uses a trained Random Forest classifier
to predict the biological probability of PHA production based on 
physicochemical, HMM, and operon features.
"""

import os
import joblib
import logging
from phascout.prediction.physicochemical import analyze_physicochemical

logger = logging.getLogger(__name__)

class MLScorer:
    def __init__(self, model_path="data/rf_model.pkl"):
        self.model_path = model_path
        self.model = None
        if os.path.exists(model_path):
            try:
                self.model = joblib.load(model_path)
                logger.info(f"ML Modeli yuklendi: {model_path}")
            except Exception as e:
                logger.error(f"ML Modeli yuklenemedi: {e}")
        else:
            logger.warning(f"ML Modeli bulunamadi: {model_path}. Fallback modu kullanilacak.")
            
    def predict(self, phac_result, gene_vector, operon_result, phac_seq):
        """
        Calculates the probability of PHA production.
        """
        features = self._extract_features(phac_result, gene_vector, operon_result, phac_seq)
        
        result = {
            "ml_probability": 0.0,
            "features_used": features,
            "model_loaded": self.model is not None
        }
        
        if not phac_result.get("is_functional") or not phac_seq:
            return result
            
        if self.model:
            try:
                import pandas as pd
                # Ensure feature order matches training
                feature_order = ["hmm_score", "mw", "pi", "gravy", "operon_dist"]
                df = pd.DataFrame([features])[feature_order]
                prob = self.model.predict_proba(df)[0][1] # Probability of class 1 (Positive)
                result["ml_probability"] = float(prob * 100)
            except Exception as e:
                logger.error(f"Tahmin sirasinda hata: {e}")
        else:
            # Fallback simple logic if model is not trained yet
            prob = 50.0
            if phac_result.get("is_functional"):
                prob += 20.0
            if operon_result.get("is_class_i_operon"):
                prob += 20.0
            if features["gravy"] < 0: # Soluble
                prob += 10.0
            result["ml_probability"] = prob
            
        return result
        
    def _extract_features(self, phac_result, gene_vector, operon_result, phac_seq):
        # HMM Score
        hmm_score = phac_result.get("best_score", 0.0)
        
        # Physicochemical
        physico = analyze_physicochemical(phac_seq)
        
        # Operon Distance (min distance to phaA or phaB)
        distances = operon_result.get("distances", {})
        min_dist = 10000 # Default large distance if no operon
        if distances:
            min_dist = min(distances.values())
            
        return {
            "hmm_score": hmm_score,
            "mw": physico["molecular_weight"],
            "pi": physico["isoelectric_point"],
            "gravy": physico["gravy"],
            "operon_dist": min_dist
        }
