import pyhmmer
from pathlib import Path

def print_hmm_consensus(hmm_path):
    with pyhmmer.plan7.HMMFile(hmm_path) as f:
        hmm = next(f)
    print(f"--- {hmm_path} ---")
    consensus = hmm.consensus
    # find C, D, H
    # The Cys is highly conserved, usually C
    c_pos = [i+1 for i, c in enumerate(consensus) if c == 'C']
    d_pos = [i+1 for i, c in enumerate(consensus) if c == 'D']
    h_pos = [i+1 for i, c in enumerate(consensus) if c == 'H']
    print(f"C positions: {c_pos}")
    print(f"D positions: {d_pos}")
    print(f"H positions: {h_pos}")

for p in Path("data/hmm_profiles/phac_classes").glob("*.hmm"):
    print_hmm_consensus(p)
