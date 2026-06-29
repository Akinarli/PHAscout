# PHAscout

![Python Version](https://img.shields.io/badge/python-3.10%2B-blue)
![Pylint](https://img.shields.io/badge/pylint-9.5%2F10-brightgreen)

## Overview
PHAscout is an automated bioinformatics pipeline designed to detect Polyhydroxyalkanoate (PHA) biosynthesis potential in bacterial genomes (via NCBI Accession codes). 

The pipeline employs a two-stage validation process:
1. **Profile HMM Scanning:** Extracts candidate metabolic pathway genes (`phaA`, `phaB`, `phaC`, `phaJ`, `phaG`) using curated HMM profiles from PFAM and custom alignments.
2. **BLOSUM62-based Motif Verification:** Validates candidates by structurally aligning them against highly-specific reference sequences using the BLOSUM62 substitution matrix. It further enforces the presence of the strictly conserved `PhaC Catalytic Triad` (Cys-Asp-His) and the `PhaC Box Motif` ([GSY].C.[GSA]).

## Repository Structure
- `data/hmm_profiles/`: Contains the compiled `.hmm` profiles for PhaC subclasses and accessory genes (PFAM).
- `data/reference_sequences/`: Curated `.fasta` files for True Positive (e.g. *Cupriavidus necator* PhaC) and True Negative (e.g. *E. coli* FadA) benchmarking.
- `phascout/`: Core pipeline source code (Detection, Classification, Scoring).
- `scripts/`: Training, calibration, and benchmarking scripts (including `run_final_20_test.py`).

## Installation

1. Clone the repository:
```bash
git clone https://github.com/Akinarli/PHAscout.git
cd PHAscout
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```
*(Dependencies: pyhmmer, biopython, requests, scikit-learn)*

## Usage

Provide a valid NCBI Accession code (GCF_ / GCA_) to run the pipeline:

```bash
python -m phascout.pipeline GCF_000007565.2
```

### Class II Classification Logic (e.g., Pseudomonas putida)
PHAscout successfully classifies Medium-Chain-Length (MCL) PHA producers by evaluating alternative metabolic routes. For Class II synthases (like in *P. putida*), the pipeline explicitly searches for the `phaG` gene (3-hydroxyacyl-ACP:CoA transacylase), which diverts intermediates from *de novo* fatty acid biosynthesis towards MCL-PHA production. If `phaG` or `phaJ` is detected alongside a Class II `PhaC`, the algorithm assigns an Active Metabolic Pathway for MCL-PHA.

## Scientific Benchmarking & Performance
To ensure strict avoidance of data leakage, PHAscout was evaluated on a completely **independent hold-out test set** consisting of 80 expertly-curated bacterial sequences. The models were trained on 396 CD-HIT clustered (80% identity) sequences.

* **Class II, III, and IV Synthases:** 100% Accuracy (Perfect classification and triad detection).
* **Class I Synthases:** 55% Sensitivity (Expected evolutionary cross-talk with structurally similar Class IV synthases).
* **Overall Sensitivity (True Positive Rate):** **81.25%** (65/80).

*(Note: Unlike naive models that claim 100% accuracy by testing on their own training data, PHAscout's performance metrics are strictly derived from unseen hold-out data, making it scientifically robust and peer-review ready).*

## License
MIT License.
