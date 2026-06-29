import csv

# 1. Parse Table 3 to find "No data" species
no_data_species = []
with open('table3_output.txt', 'r', encoding='utf-8') as f:
    lines = [line.strip() for line in f if line.strip()]

for i, line in enumerate(lines):
    if line.startswith('Halomonas ') or line.startswith('Candidatus Halomonas '):
        if i + 1 < len(lines) and lines[i+1] == 'No data':
            no_data_species.append(line)

print("--- Tabloda 'No data' (Veri Yok) Olan Turler ---")
for sp in no_data_species:
    print(sp)

print("\n--- PHAscout'un 'No data' Turlerindeki Kesifleri ---")
# 2. Check their results in benchmark_results.csv
results = {}
with open('benchmark_results.csv', 'r', encoding='utf-8-sig') as f:
    reader = csv.DictReader(f)
    for row in reader:
        results[row['Species']] = row['Predicted_Type']

for sp in no_data_species:
    pred = results.get(sp, "Test edilmedi (GCF bulunamadi)")
    print(f"{sp} -> PHAscout Sonucu: {pred}")
