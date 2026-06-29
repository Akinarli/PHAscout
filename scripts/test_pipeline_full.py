import sys
from phascout.pipeline import PHAscoutPipeline

def test_genome(pipeline, accession, expected_class):
    print(f"\n--- Test Ediliyor: {accession} (Beklenen: {expected_class}) ---")
    try:
        report = pipeline.run(accession=accession)
        phac_info = report.get("phac_analysis", {})
        ml_info = report.get("heuristic_result", {})
        
        pred_class = phac_info.get("class")
        ml_prob = ml_info.get("ml_probability", 0)
        
        print(f"Tahmin Edilen Sinif: {pred_class}")
        print(f"ML Sentez Olasiligi: %{ml_prob:.2f}")
        
        if pred_class == expected_class:
            print("[BASARILI] Siniflandirma dogru.")
        else:
            print(f"[HATA] Beklenen {expected_class}, ancak {pred_class} bulundu.")
            
    except Exception as e:
        print(f"HATA: {e}")

if __name__ == "__main__":
    pipeline = PHAscoutPipeline()
    
    test_genomes = [
        ("GCF_000009285.1", "Class_I"),   # Cupriavidus necator
        ("GCF_000007565.2", "Class_II"),  # Pseudomonas putida
        ("GCF_000832985.1", "Class_IV")   # Bacillus megaterium
    ]
    
    for acc, cls in test_genomes:
        test_genome(pipeline, acc, cls)
