import urllib.request
import urllib.parse
import os
import json

def fetch_uniprot(query, output_file, limit=100):
    url = f"https://rest.uniprot.org/uniprotkb/search?format=fasta&query={urllib.parse.quote(query)}&size={limit}"
    print(f"Sorgulanıyor: {query}")
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'PHAscout-Bot'})
        with urllib.request.urlopen(req) as response:
            data = response.read().decode('utf-8')
            with open(output_file, "w", encoding="utf-8") as f:
                f.write(data)
            count = data.count('>')
            print(f" -> {count} adet dizi indirildi ve {output_file} konumuna kaydedildi.")
    except Exception as e:
        print(f"Hata ({query}): {e}")

base_dir = r"C:\Users\bird-\OneDrive\Desktop\In silico\PHAscout"
pos_dir = os.path.join(base_dir, "data", "reference_sequences", "positive")
neg_dir = os.path.join(base_dir, "data", "reference_sequences", "negative")

os.makedirs(pos_dir, exist_ok=True)
os.makedirs(neg_dir, exist_ok=True)

base_filter = "AND (reviewed:true) AND (taxonomy_id:2)"

queries = {
    os.path.join(pos_dir, "phab.fasta"): f'(gene:phab OR protein_name:"acetoacetyl-CoA reductase") {base_filter}',
    os.path.join(pos_dir, "phaa.fasta"): f'(gene:phaa OR protein_name:"polyhydroxyalkanoate synthesis-related beta-ketothiolase") {base_filter}',
    os.path.join(pos_dir, "phag.fasta"): f'(gene:phag OR protein_name:"3-hydroxyacyl-ACP:CoA transacylase") {base_filter}',
    os.path.join(pos_dir, "phaj.fasta"): f'(gene:phaj OR protein_name:"enoyl-CoA hydratase") {base_filter}',
    os.path.join(neg_dir, "fabg.fasta"): f'(gene:fabg) {base_filter}',
    os.path.join(neg_dir, "fada.fasta"): f'(gene:fada) {base_filter}',
    os.path.join(neg_dir, "sdr_broad.fasta"): f'(protein_name:"short-chain dehydrogenase") {base_filter}'
}

for out, q in queries.items():
    fetch_uniprot(q, out, limit=100)

print("İşlem tamamlandı!")
