import urllib.request
import urllib.parse
import os

def fetch_uniprot_query(query, output_file):
    url = f"https://rest.uniprot.org/uniprotkb/search?format=fasta&query={urllib.parse.quote(query)}&size=500"
    print(f"Sorgulanıyor: {query}")
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'PHAscout-Bot'})
        with urllib.request.urlopen(req) as response:
            data = response.read().decode('utf-8')
            with open(output_file, "w", encoding="utf-8") as f:
                f.write(data)
            count = data.count('>')
            print(f" -> {count} adet dizi indirildi: {os.path.basename(output_file)}")
            return count
    except Exception as e:
        print(f"Hata ({query}): {e}")
        return 0

base_dir = r"C:\Users\bird-\OneDrive\Desktop\In silico\PHAscout\data\reference_sequences\phac_classes"
query = '(gene:phac) AND (reviewed:true)'
fetch_uniprot_query(query, os.path.join(base_dir, "phac_all_reviewed.fasta"))
