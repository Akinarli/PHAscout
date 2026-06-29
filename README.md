# PHAscout

![Python Version](https://img.shields.io/badge/python-3.10%2B-blue)
![Pylint](https://img.shields.io/badge/pylint-9.5%2F10-brightgreen)

## Overview
PHAscout is an automated bioinformatics pipeline designed to detect Polyhydroxyalkanoate (PHA) biosynthesis potential in bacterial genomes (via NCBI Accession codes). 

The pipeline employs a two-stage validation process:
1. **Profile HMM Scanning:** Extracts candidate metabolic pathway genes (`phaA`, `phaB`, `phaC`, `phaJ`, `phaG`) using curated HMM profiles from PFAM and custom alignments.
2. **BLOSUM62-based Motif Verification:** Validates candidates by structurally aligning them against highly-specific, manually-decontaminated reference sequences using the BLOSUM62 substitution matrix. For `PhaC`, functionality is confirmed by anchoring the conserved catalytic triad (Cys-Asp-His) to **class-specific HMM match-state columns** (derived from experimentally verified active sites — see `scripts/derive_catalytic_columns.py`), with a tight ±2-column tolerance. The `PhaC Box` motif (`G-x-C-x-G-G`) is enforced to contain that same catalytic Cys nucleophile.

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

### Genome-scale benchmark (balanced, with negatives)
The primary benchmark (`scripts/run_final_20_test.py`) runs the **complete pipeline** on 20 fully-sequenced genomes via NCBI: **10 confirmed PHA producers** (Class I–III) and **10 confirmed non-producers** (*E. coli* K-12, *B. subtilis* 168, *S. aureus*, *M. tuberculosis*, *H. pylori*, etc.). Because it contains real negatives, specificity and precision are meaningful:

| Metric | Value |
|---|---|
| Sensitivity (TPR) | 90.0 % (9/10) |
| Specificity (TNR) | 100.0 % (10/10) |
| Precision (PPV) | 100.0 % (0 false positives) |
| Accuracy | 95.0 % |
| MCC | 0.905 |

The single false negative is *Aeromonas caviae*: its PhaC is correctly detected, classified (Class I) and validated as functional, but it is a non-canonical producer that supplies monomers through the β-oxidation/PhaJ route rather than the `phaA`/`phaB` pathway the SCL-PHA logic requires — a deliberately conservative call, not a detection failure.

### PhaB ↔ FabG discrimination (the core specificity problem)
After removing contaminant sequences (FabG, MabA, FabM, a K⁺/H⁺ antiporter and a PhaC mislabelled into the PhaB/PhaA positive sets), the BLOSUM62 double-layer threshold was re-calibrated (`scripts/calibrate_thresholds.py`). The PhaB ↔ FabG separation improved from **F1 = 0.50 (coin-flip)** to **F1 = 0.92**.

### ML auxiliary-confidence model
The Random Forest (`scripts/train_ml_scorer.py`) is trained **only on real features** — the live PhaC class-HMM bit score plus ProtParam physicochemistry — extracted identically at training and inference time (no synthetic/dummy values). Negatives are real non-PhaC proteins, including 135 α/β-hydrolases (lipases, esterases, epoxide hydrolases) that are the genuine confusables. Honest 5-fold out-of-fold performance: **ROC-AUC = 0.9998, MCC = 0.98**.

### Classification-only set (no negatives)
`data/benchmark/independent_benchmark_set.csv` (80 curated UniProt PhaC sequences, all positive) is a **class-assignment / sensitivity-only** check. It contains **no negatives**, so it cannot report specificity or precision and must not be read as overall accuracy.

### Reproducibility
Catalytic HMM columns: `python scripts/derive_catalytic_columns.py`. Threshold calibration: `python scripts/calibrate_thresholds.py`. Model: `python scripts/train_ml_scorer.py`. Genome benchmark: `python scripts/run_final_20_test.py`.

## License
MIT License.
