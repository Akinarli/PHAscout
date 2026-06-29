"""
Katalitik HMM Kolonlarini Turetme (Yeniden Uretilebilir)
=========================================================
Deneysel olarak dogrulanmis katalitik kalintilari (UniProt ACT_SITE Cys +
literatur Asp/His), her sinif icin bir referans diziye hizalayarak ve mevcut
HMM modellerine (data/hmm_profiles/phac_classes/) hmmalign ile esleştirerek
config.CATALYTIC_HMM_COLUMNS degerlerini yeniden uretir.

Calistirma:
    python scripts/derive_catalytic_columns.py
"""

import pyhmmer
import requests
from Bio import SeqIO
from io import StringIO

# (sinif, hmm_yolu, referans_uniprot, {kalinti: (aa, protein_pozisyonu_1indexed)})
# Protein pozisyonlari dizilere karsi DOGRULANMISTIR (residue identity check).
TASKS = [
    ("Class_I", "data/hmm_profiles/phac_classes/phac_class_I.hmm", "P23608",
     {"Cys": ("C", 319), "Asp": ("D", 480), "His": ("H", 508)}),
    ("Class_II", "data/hmm_profiles/phac_classes/phac_class_II.hmm", "P26494",
     {"Cys": ("C", 296), "Asp": ("D", 451), "His": ("H", 479)}),
    ("Class_III", "data/hmm_profiles/phac_classes/phac_class_III.hmm", "P45370",
     {"Cys": ("C", 149), "Asp": ("D", 302), "His": ("H", 331)}),
    ("Class_IV", "data/hmm_profiles/phac_classes/phac_class_IV.hmm", "A0A1J9SXW8",
     {"Cys": ("C", 151), "Asp": ("D", 306), "His": ("H", 335)}),
]


def get_seq(acc):
    txt = requests.get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=30).text
    return str(next(SeqIO.parse(StringIO(txt), "fasta")).seq)


def find_columns(hmm_path, seq_str, targets):
    alphabet = pyhmmer.easel.Alphabet.amino()
    ds = pyhmmer.easel.TextSequence(name=b"ref", sequence=seq_str).digitize(alphabet)
    with pyhmmer.plan7.HMMFile(hmm_path) as f:
        hmm = next(f)
    msa = pyhmmer.hmmalign(hmm, [ds], trim=False)
    aln = msa.alignment[0]
    if isinstance(aln, bytes):
        aln = aln.decode()

    protein_pos = 0
    hmm_col = 0
    found = {}
    for ch in aln:
        is_node = ch.isupper() or ch == "-"
        is_res = ch.isalpha()
        if is_node:
            hmm_col += 1
        if is_res:
            protein_pos += 1
        if is_res and is_node:
            for name, (aa, tpos) in targets.items():
                if protein_pos == tpos and ch.upper() == aa:
                    found[name] = hmm_col
    return found, hmm.M


def main():
    print("CATALYTIC_HMM_COLUMNS = {")
    for cls, hmm_path, acc, targets in TASKS:
        seq = get_seq(acc)
        # Kalinti kimligini dogrula
        for name, (aa, tpos) in targets.items():
            actual = seq[tpos - 1] if tpos - 1 < len(seq) else "?"
            assert actual == aa, f"{cls} {name}{tpos}: beklenen {aa}, bulunan {actual}"
        cols, M = find_columns(hmm_path, seq, targets)
        print(f'    "{cls}": {{"Cys": {cols.get("Cys")}, '
              f'"Asp": {cols.get("Asp")}, "His": {cols.get("His")}}},  '
              f'# {acc}, HMM M={M}')
    print("}")


if __name__ == "__main__":
    main()
