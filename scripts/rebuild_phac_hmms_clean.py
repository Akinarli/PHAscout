import os
import time
from Bio import SeqIO
import requests
import pyhmmer

base_dir = r"C:\Users\bird-\OneDrive\Desktop\In silico\PHAscout"
raw_dir = os.path.join(base_dir, "data", "reference_sequences", "phac_classes")
clean_dir = os.path.join(base_dir, "data", "reference_sequences", "phac_classes_clean")
hmm_dir = os.path.join(base_dir, "data", "hmm_profiles", "phac_classes")

os.makedirs(clean_dir, exist_ok=True)
os.makedirs(hmm_dir, exist_ok=True)

# Sınıf bazlı minimum ve maksimum uzunluklar (Filtreleme için)
LENGTH_LIMITS = {
    "phac_class_I": (500, 650),
    "phac_class_II": (500, 650),
    "phac_class_III": (300, 450),
    "phac_class_IV": (300, 450),
}

def clean_fasta(class_name):
    print(f"[{class_name}] Eski FASTA dosyasından kırpılmış ve hatalı olanlar temizleniyor...")
    in_file = os.path.join(raw_dir, f"{class_name}.fasta")
    out_file = os.path.join(clean_dir, f"{class_name}_filtered.fasta")
    
    if not os.path.exists(in_file):
        print(f"HATA: Dosya bulunamadi: {in_file}")
        return None
        
    records = list(SeqIO.parse(in_file, "fasta"))
    min_len, max_len = LENGTH_LIMITS[class_name]
    
    clean_records = [r for r in records if min_len <= len(r.seq) <= max_len]
    
    for r in clean_records:
        print(f"  -> KABUL EDİLDİ: {r.id} ({len(r.seq)} AA)")
        
    if not clean_records:
        print(f"  -> UYARI: Bu sınıf için uygun uzunlukta dizi kalmadı! Filtre gevşetilmeli.")
        return None
        
    SeqIO.write(clean_records, out_file, "fasta")
    return out_file


def ebi_clustalo_align(fasta_file):
    print(f"[{os.path.basename(fasta_file)}] EBI ClustalO ile hizalaniyor...")
    with open(fasta_file, "r") as f:
        seq_data = f.read()
    
    run_url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run"
    payload = {
        "email": "phascout.clean@example.com",
        "sequence": seq_data,
        "outfmt": "fa"
    }
    r = requests.post(run_url, data=payload)
    if r.status_code != 200:
        print(f"EBI API Hatasi: {r.text}")
        return None
    job_id = r.text.strip()
    
    status_url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/{job_id}"
    while True:
        time.sleep(5)
        s = requests.get(status_url)
        status = s.text.strip()
        print(f"  Durum: {status}")
        if status == "FINISHED":
            break
        elif status in ["ERROR", "FAILURE", "NOT_FOUND"]:
            print(f"Hizalama Basarisiz: {status}")
            return None
            
    res_url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{job_id}/fa"
    r = requests.get(res_url)
    aln_file = fasta_file.replace(".fasta", "_aligned.fasta")
    with open(aln_file, "w") as f:
        f.write(r.text)
    return aln_file


def build_hmm(aln_file, class_name):
    print(f"[{class_name}] TEMIZ HMM profili pyhmmer ile olusturuluyor...")
    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.MSAFile(aln_file, format='afa', digital=True, alphabet=alphabet) as msa_file:
        msa = msa_file.read()
        msa.name = class_name.encode()
        
        builder = pyhmmer.plan7.Builder(alphabet)
        hmm, profile, opt = builder.build_msa(msa, pyhmmer.plan7.Background(alphabet))
        
        hmm_path = os.path.join(hmm_dir, f"{class_name}.hmm")
        with open(hmm_path, "wb") as f:
            hmm.write(f)
        print(f"  --> {hmm_path} basariyla olusturuldu! (Model Uzunlugu: {hmm.M} dugum)\n")

def main():
    for cls in LENGTH_LIMITS.keys():
        fasta_path = clean_fasta(cls)
        if fasta_path:
            aln_path = ebi_clustalo_align(fasta_path)
            if aln_path:
                build_hmm(aln_path, cls)
    print("ASAMA 8: Tum HMM modelleri altin standart dizilerle basariyla yenilendi!")

if __name__ == "__main__":
    main()
