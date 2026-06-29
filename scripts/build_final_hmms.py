import os
import time
from Bio import SeqIO
import requests
import pyhmmer

base_dir = r"C:\Users\bird-\OneDrive\Desktop\In silico\PHAscout"
clean_dir = os.path.join(base_dir, "data", "reference_sequences", "phac_classes_clean")
hmm_dir = os.path.join(base_dir, "data", "hmm_profiles", "phac_classes")

os.makedirs(hmm_dir, exist_ok=True)

CLASSES = ["phac_class_i", "phac_class_ii", "phac_class_iii", "phac_class_iv"]

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
        time.sleep(10)
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
    print(f"[{class_name}] HMM profili pyhmmer ile olusturuluyor...")
    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.MSAFile(aln_file, format='afa', digital=True, alphabet=alphabet) as msa_file:
        msa = msa_file.read()
        
        # Sınıf adına göre doğru ismi ver
        proper_name = "Class_" + class_name.split("_")[-1].upper()
        msa.name = proper_name.encode()
        
        builder = pyhmmer.plan7.Builder(alphabet)
        hmm, profile, opt = builder.build_msa(msa, pyhmmer.plan7.Background(alphabet))
        
        hmm_path = os.path.join(hmm_dir, f"{proper_name}.hmm") # phac_class_I.hmm degil, Class_I.hmm olmamali
        # Config dosyamiz sunu bekliyor: phac_class_I.hmm
        # PHAC_CLASS_PROFILES icinde isim boyle.
        
        real_hmm_path = os.path.join(hmm_dir, f"phac_{class_name.split('_')[-1].upper()}.hmm") 
        if class_name == "phac_class_i":
            real_hmm_path = os.path.join(hmm_dir, "phac_class_I.hmm")
        elif class_name == "phac_class_ii":
            real_hmm_path = os.path.join(hmm_dir, "phac_class_II.hmm")
        elif class_name == "phac_class_iii":
            real_hmm_path = os.path.join(hmm_dir, "phac_class_III.hmm")
        elif class_name == "phac_class_iv":
            real_hmm_path = os.path.join(hmm_dir, "phac_class_IV.hmm")

        with open(real_hmm_path, "wb") as f:
            hmm.write(f)
        print(f"  --> {real_hmm_path} basariyla olusturuldu! (Model Uzunlugu: {hmm.M} dugum)\n")

def main():
    for cls in CLASSES:
        fasta_path = os.path.join(clean_dir, f"{cls}.fasta")
        if not os.path.exists(fasta_path):
            print(f"ATLANDI: {fasta_path} bulunamadi.")
            continue
            
        aln_path = ebi_clustalo_align(fasta_path)
        if aln_path:
            build_hmm(aln_path, cls)

if __name__ == "__main__":
    main()
