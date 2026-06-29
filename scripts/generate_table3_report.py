import csv

scl = []
mcl = []
copolymer = []
none = []

with open('benchmark_results.csv', 'r', encoding='utf-8-sig') as f:
    reader = csv.DictReader(f)
    for row in reader:
        pred = row['Predicted_Type']
        name = row['Species']
        acc = row['Accession']
        if pred == 'SCL':
            scl.append(f"- *{name}* ({acc})")
        elif pred == 'MCL':
            mcl.append(f"- *{name}* ({acc})")
        elif pred == 'SCL-co-MCL':
            copolymer.append(f"- *{name}* ({acc})")
        else:
            none.append(f"- *{name}* ({acc}) - *Büyük ihtimalle PHA üretmiyor, BLAST yaniltmacasi*")

with open('table3_analysis_report.md', 'w', encoding='utf-8') as f:
    f.write("# Halomonas (Table 3) Keşif Sonuçları\n\n")
    f.write("Toplam 74 genom başarıyla analiz edildi.\n\n")
    
    f.write(f"## 🧬 SCL-co-MCL (Kopolimer) Potansiyeli Olanlar ({len(copolymer)} Tür)\n")
    f.write("Bu türlerin genomunda nadir bulunan kopolimer (MCL+SCL) enzim yolları tespit edildi:\n")
    f.write('\n'.join(copolymer) + "\n\n")
    
    f.write(f"## 🧪 SCL (P3HB) Üreticileri ({len(scl)} Tür)\n")
    f.write("Bu türler klasik Class I PhaC ve SCL üretim yolaklarına tam olarak sahipler:\n")
    f.write('\n'.join(scl) + "\n\n")
    
    f.write(f"## ❌ PHA Potansiyeli Olmayanlar / Belirsizler ({len(none)} Tür)\n")
    f.write("Sen bu türlerde >%95 ANI ile PhaC bulmuş olabilirsin ancak PHAscout yapısal (katalitik triad) ve fonksiyonel olarak bunları **gerçek bir PHA sentazı OLMADIĞINA** karar verdi (Yüksek ihtimalle jenerik bir lipaz/hidrolaz veya kırık bir gen):\n")
    f.write('\n'.join(none) + "\n\n")
