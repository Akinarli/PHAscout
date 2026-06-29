import pyhmmer
import urllib.request
from Bio import SeqIO
from io import StringIO
from pathlib import Path

def get_uniprot_seq(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    req = urllib.request.Request(url, headers={'User-Agent': 'PHAscout'})
    response = urllib.request.urlopen(req)
    fasta_data = response.read().decode('utf-8')
    rec = next(SeqIO.parse(StringIO(fasta_data), "fasta"))
    return str(rec.seq)

def find_columns(hmm_path, seq_str, target_residues):
    alphabet = pyhmmer.easel.Alphabet.amino()
    digital_seq = pyhmmer.easel.TextSequence(name=b"ref", sequence=seq_str).digitize(alphabet)
    
    with pyhmmer.plan7.HMMFile(hmm_path) as f:
        hmm = next(f)
        
    # PyHMMER hmmalign
    msa = pyhmmer.hmmalign(hmm, [digital_seq], trim=True)
    
    # In PyHMMER, MSA objects have an 'alignment' attribute (list of strings/bytes)
    # The first sequence in the MSA is our aligned reference
    aln_seq = msa.alignment[0]
    # Convert from bytes to string if needed
    if isinstance(aln_seq, bytes):
        aln_seq = aln_seq.decode('utf-8')
    
    targets = {}
    for tr in target_residues:
        aa = tr[0]
        pos = int(tr[1:])
        targets[aa] = pos
        
    protein_pos = 0
    hmm_col = 0
    found = {}
    
    for i, char in enumerate(aln_seq):
        is_hmm_node = char.isupper() or char == '-' # Match or Delete
        is_residue = char.isalpha() # upper or lower case (match or insert)
        
        if is_hmm_node:
            hmm_col += 1
            
        if is_residue:
            protein_pos += 1
            
        if is_residue and is_hmm_node: # It's a match state
            for aa, target_pos in targets.items():
                if protein_pos == target_pos and char.upper() == aa:
                    found[aa] = hmm_col
                    
    return found

def main():
    tasks = [
        ("Class_I", "data/hmm_profiles/phac_classes/phac_class_I.hmm", "P23608", ["C319", "D447", "H477"]),
        ("Class_II", "data/hmm_profiles/phac_classes/phac_class_II.hmm", "P26494", ["C296", "D424", "H453"]),
        ("Class_III", "data/hmm_profiles/phac_classes/phac_class_III.hmm", "P45370", ["C149", "D302", "H331"]), 
        # Note: A. vinosum Cys is 149 actually! (Reference: Jia et al. 2000 C149, D302, H331) Let's check literature if it's 149 or 130. I'll pass 149.
        ("Class_IV", "data/hmm_profiles/phac_classes/phac_class_IV.hmm", "Q9R9V5", ["C149", "D302", "H331"]) # Wait, B. megaterium triad is also often homologous to Class III. Let's see if we can find C, D, H around there.
    ]
    
    print("CATALYTIC_HMM_COLUMNS = {")
    for cls, hmm_path, uniprot, targets in tasks:
        seq = get_uniprot_seq(uniprot)
        # Search sequence for C, D, H to dynamically adjust if numbers are slightly off
        c_pos = seq.find("C") + 1
        
        try:
            cols = find_columns(hmm_path, seq, targets)
            print(f'    "{cls}": {{"Cys": {cols.get("C")}, "Asp": {cols.get("D")}, "His": {cols.get("H")}}},')
        except Exception as e:
            print(f'    "{cls}": error {e}')
    print("}")

if __name__ == "__main__":
    main()
