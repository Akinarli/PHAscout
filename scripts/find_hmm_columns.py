import argparse
from Bio import AlignIO

def find_hmm_columns(alignment_file, target_residues):
    aln = AlignIO.read(alignment_file, "fasta")
    ref_seq = str(aln[0].seq)
    
    targets = {}
    for tr in target_residues:
        aa = tr[0]
        pos = int(tr[1:])
        targets[aa] = pos
        
    protein_pos = 0
    hmm_col = 0
    found = {}
    
    for i, char in enumerate(ref_seq):
        if char != '-':
            protein_pos += 1
            hmm_col += 1
            for aa, target_pos in targets.items():
                if protein_pos == target_pos and char == aa:
                    found[aa] = hmm_col
                    print(f"  {aa}{target_pos} -> HMM kolon {hmm_col}")
        else:
            hmm_col += 1
            
    return found

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment")
    parser.add_argument("--residues", nargs="+")
    args = parser.parse_args()
    
    cols = find_hmm_columns(args.alignment, args.residues)
    print(f"\nBulunan HMM kolonlari: {cols}")
