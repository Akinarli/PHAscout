import csv
import sys
import os

from phascout.pipeline import PHAscoutPipeline
import logging
logging.getLogger().setLevel(logging.ERROR) # Sadece hatalari bas, consol kalabalik etmesin

def main():
    if len(sys.argv) > 1:
        csv_file = sys.argv[1]
    else:
        csv_file = 'benchmark_dataset_100.csv'
        
    if not os.path.exists(csv_file):
        print(f"Hata: {csv_file} bulunamadi.")
        sys.exit(1)

    print("PHAscout Toplu Benchmark (Sinav) basliyor...")
    pipeline = PHAscoutPipeline()
    
    results = []
    correct_count = 0
    total_count = 0

    with open(csv_file, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            acc = row['Accession']
            species = row['Species']
            expected_type = row['Expected_Type']
            expected_class = row['Expected_Class']
            
            print(f"\nTest Ediliyor: {species} ({acc})")
            
            try:
                report = pipeline.run(accession=acc)
                predicted_type_full = report.get('pha_potential', {}).get('potential', 'ERROR')
                predicted_class = report.get('phac_result', {}).get('best_class', 'None')
                if predicted_class is None:
                    predicted_class = 'None'
                    
                # predicted_type format is usually 'SCL', 'MCL', 'SCL-co-MCL', 'YOK'
                if 'YOK' in predicted_type_full or 'none' in predicted_type_full.lower():
                    predicted_type = 'None'
                elif 'SCL-co-MCL' in predicted_type_full:
                    predicted_type = 'SCL-co-MCL'
                elif 'SCL' in predicted_type_full:
                    predicted_type = 'SCL'
                elif 'MCL' in predicted_type_full:
                    predicted_type = 'MCL'
                else:
                    predicted_type = predicted_type_full
                    
                match = (predicted_type == expected_type)
                if match:
                    print(f"  -> SONUC: DOGRU! (Beklenen: {expected_type}, Bulunan: {predicted_type})")
                    correct_count += 1
                else:
                    print(f"  -> SONUC: YANLIS (Beklenen: {expected_type}, Bulunan: {predicted_type})")
                    
                total_count += 1
                
                results.append({
                    'Accession': acc,
                    'Species': species,
                    'Expected_Type': expected_type,
                    'Predicted_Type': predicted_type,
                    'Match': match,
                    'Expected_Class': expected_class,
                    'Predicted_Class': predicted_class
                })
                
            except Exception as e:
                print(f"  -> HATA: {e}")

    print("\n" + "="*50)
    print("BENCHMARK TAMAMLANDI.")
    if total_count > 0:
        accuracy = (correct_count / total_count) * 100
        print(f"Genel Dogruluk Orani: %{accuracy:.2f} ({correct_count}/{total_count})")
        
    out_csv = 'benchmark_results.csv'
    with open(out_csv, 'w', encoding='utf-8', newline='') as f:
        fieldnames = ['Accession', 'Species', 'Expected_Type', 'Predicted_Type', 'Match', 'Expected_Class', 'Predicted_Class']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    
    print(f"Detayli sonuclar {out_csv} dosyasina kaydedildi.")

if __name__ == '__main__':
    main()
