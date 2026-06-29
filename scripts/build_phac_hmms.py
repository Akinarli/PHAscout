import urllib.request
import urllib.parse
import os
import time
import requests
import pyhmmer

# Sınıf bazlı UniProt accession numaraları (Özenle seçilmiş, literatür doğrulamalı)
CLASSES = {
    "phac_class_I": ["P23608", "P50176", "Q9ZHI2", "P45367", "Q4Q6B1", "Q939T5", "Q52728", "O66392", "P52070", "Q45028"],
    "phac_class_II": ["Q88FD4", "Q88FD3", "Q9I5V8", "Q9I5V7", "P26494", "P26495", "O68661", "O68656", "P45368", "Q8GNC0"],
    "phac_class_III": ["P45370", "Q05206", "Q8GMK8", "O52051", "Q9F9Y6"],
    "phac_class_IV": ["Q9R9V5", "Q81C10", "Q817Y3", "Q81L67", "Q81EV6"]
}

base_dir = r"C:\Users\bird-\OneDrive\Desktop\In silico\PHAscout"
seq_dir = os.path.join(base_dir, "data", "reference_sequences", "phac_classes")
hmm_dir = os.path.join(base_dir, "data", "hmm_profiles", "phac_classes")
os.makedirs(seq_dir, exist_ok=True)
os.makedirs(hmm_dir, exist_ok=True)

def fetch_sequences(class_name, ids):
    query = " OR ".join([f"accession:{acc}" for acc in ids])
    url = f"https://rest.uniprot.org/uniprotkb/search?format=fasta&query={urllib.parse.quote(query)}&size=50"
    out_file = os.path.join(seq_dir, f"{class_name}.fasta")
    print(f"[{class_name}] Diziler UniProt'tan indiriliyor...")
    req = urllib.request.Request(url, headers={'User-Agent': 'PHAscout'})
    with urllib.request.urlopen(req) as res:
        data = res.read().decode('utf-8')
        with open(out_file, "w") as f:
            f.write(data)
    return out_file

def ebi_clustalo_align(fasta_file):
    print(f"[{os.path.basename(fasta_file)}] EBI ClustalO ile hizalanıyor...")
    with open(fasta_file, "r") as f:
        seq_data = f.read()
    
    run_url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run"
    payload = {
        "email": "phascout.project@gmail.com",
        "sequence": seq_data,
        "outfmt": "fa"
    }
    r = requests.post(run_url, data=payload)
    if r.status_code != 200:
        print(f"EBI API Hatası: {r.text}")
        return None
    job_id = r.text.strip()
    
    status_url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/{job_id}"
    while True:
        time.sleep(3)
        s = requests.get(status_url)
        status = s.text.strip()
        print(f"  Durum: {status}")
        if status == "FINISHED":
            break
        elif status in ["ERROR", "FAILURE", "NOT_FOUND"]:
            print(f"Hizalama Başarısız: {status}")
            return None
            
    res_url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{job_id}/fa"
    r = requests.get(res_url)
    aln_file = fasta_file.replace(".fasta", "_aligned.fasta")
    with open(aln_file, "w") as f:
        f.write(r.text)
    return aln_file

def build_hmm(aln_file, class_name):
    print(f"[{class_name}] HMM profili pyhmmer ile oluşturuluyor...")
    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.MSAFile(aln_file, format='afa', digital=True, alphabet=alphabet) as msa_file:
        msa = msa_file.read()
        msa.name = class_name.encode()
        
        builder = pyhmmer.plan7.Builder(alphabet)
        hmm, profile, opt = builder.build_msa(msa, pyhmmer.plan7.Background(alphabet))

        
        hmm_path = os.path.join(hmm_dir, f"{class_name}.hmm")
        with open(hmm_path, "wb") as f:
            hmm.write(f)
        print(f"  --> {hmm_path} başarıyla oluşturuldu!")

for cls, ids in CLASSES.items():
    fasta_path = fetch_sequences(cls, ids)
    aln_path = ebi_clustalo_align(fasta_path)
    if aln_path:
        build_hmm(aln_path, cls)

print("Aşama 0.3 Tamamlandı!")
