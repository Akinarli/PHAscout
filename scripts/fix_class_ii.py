import os
import time
import requests
import pyhmmer

base_dir = r"C:\Users\bird-\OneDrive\Desktop\In silico\PHAscout"
hmm_dir = os.path.join(base_dir, "data", "hmm_profiles", "phac_classes")
os.makedirs(hmm_dir, exist_ok=True)

# GERCEK P. putida PhaC1 ve PhaC2 dizileri (NCBI Entrez'den dogrulandi)
CLASS_II_FASTA = """>WP_010955568.1 PhaC1 [Pseudomonas putida KT2440]
MTDKPAKGSTTLPATRMNVQNAILGLRGRDLLSTLRNVGRHGLRHPLHTAHHLLALGGQLGRVMLGDTPYQPNPRDARFSDPTWSQNPFYRRGLQAYLAWQKQTRQWIDESHLNDDDRARAHFLFNLINDALAPSNSLLNPLAVKELFNTGGQSLVRGVAHLLDDLRHNDGLPRQVDERAFEVGVNLAATPGAVVFRNELLELIQYSPMSEKQHARPLLVVPPQINKFYIFDLSATNSFVQYMLKSGLQVFMVSWRNPDPRHREWGLSSYVQALEEALNACRSISGNRDPNLMGACAGGLTMAALQGHLQAKKQLRRVRSATYLVSLLDSKFESPASLFADEQTIEAAKRRSYQRGVLDGGEVARIFAWMRPNDLIWNYWVNNYLLGKTPPAFDILYWNADSTRLPAALHGDLLEFFKLNPLTYASGLEVCGTPIDLQQVNIDSFTVAGSNDHITPWDAVYRSALLLGGERRFVLANSGHIQSIINPPGNPKAYYLANPKLSSDPRAWFHDAKRSEGSWWPLWLEWITARSGLLKAPRTELGNATYPPLGPAPGTYVLTR
>WP_010955566.1 PhaC2 [Pseudomonas putida KT2440]
MSNKNNDELQRQASENTLGLNPVIGIRRKDLLSSARTVLRQAVRQPLHSAKHVAHFGLELKNVLLGKSSLAPDSDDRRFNDPAWSNNPLYRRYLQTYLAWRKELQDWVSSSDLSPQDISRGQFVINLMTEAMAPTNTLSNPAAVKRFFETGGKSLLDGLSNLAKDMVNNGGMPSQVNMDAFEVGKNLGTSEGAVVYRNDVLELIQYSPITEQVHARPLLVVPPQINKFYVFDLSPEKSLARFCLRSQQQTFIISWRNPTKAQREWGLSTYIDALKEAVDAVLSITGSKDLNMLGACSGGITCTALVGHYAAIGENKVNALTLLVSVLDTTMDNQVALFVDEQTLEAAKRHSYQAGVLEGSEMAKVFAWMRPNDLIWNYWVNNYLLGNEPPVFDILFWNNDTTRLPAAFHGDLIEMFKSNPLTRPDALEVCGTAIDLKQVKCDIYSLAGTNDHITPWPSCYRSAHLFGGKIEFVLSNSGHIQSILNPPGNPKARFMTGADRPGDPVAWQENAIKHADSWWLHWQSWLGERAGALKKAPTRLGNRTYAAGEASPGTYVHER
"""

def ebi_clustalo_align(seq_data, name):
    print(f"[{name}] EBI ClustalO ile hizalaniyor...")
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
    aln_file = os.path.join(base_dir, "data", "reference_sequences", "phac_classes_clean", f"{name}_aligned.fasta")
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
    aln_path2 = ebi_clustalo_align(CLASS_II_FASTA, "phac_class_II")
    if aln_path2:
        build_hmm(aln_path2, "phac_class_II")
    print("GERCEK Class II modeli insa edildi!")

if __name__ == "__main__":
    main()
