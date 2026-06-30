# PHAscout Methodology

The PHAscout pipeline is a de novo, open-source computational tool designed for the genome-wide identification and classification of polyhydroxyalkanoate (PHA) biosynthesis genes. The pipeline estimates the **genomic potential** for PHA biosynthesis (currently reported as SCL, MCL, "uncertain", or none) from an annotated proteome. It does **not** assert PHA production; genomic potential requires experimental validation. The methodology consists of five sequential analytical layers.

## 1. Data Acquisition and Preprocessing
The pipeline accepts three input types: (i) a GenBank Assembly Accession (e.g., GCF_030505695.1) — the NCBI Datasets API retrieves the proteome and GFF3 annotation; (ii) a protein FASTA; or (iii) a raw, unannotated nucleotide genome/contig FASTA. Nucleotide input is auto-detected (by alphabet composition) and gene-called with **Prodigal**; Prodigal's per-gene coordinates are reused for the operon/synteny layer, so raw genomes get the same downstream analysis as annotated proteomes. Sequences are parsed with `Biopython`.

## 2. Hidden Markov Model (HMM) Scanning
To identify putative PHA biosynthesis enzymes without relying on sequence identity tools like BLAST (which are prone to false positives due to generic hydrolases), PHAscout employs probabilistic Hidden Markov Models (HMMs). 
Profile HMMs for PHA genes (`phaC`, `phaA`, `phaB`, `phaG`, `phaJ`, `phaP`, `phaR`, `phaE`) were constructed or retrieved from the Pfam database. (`phaZ`/depolymerase is not currently scanned — no usable Pfam HMM was available — and is therefore not part of the gene vector.) The `pyhmmer` library scans the entire proteome against these profiles. Proteins exhibiting statistically significant structural homology (E-value < 1e-5) are selected as preliminary candidates.

## 3. The Double-Layer BLOSUM62 Filter
The BLOSUM62 specificity filter is applied to the **monomer-supplying enzymes** (`phaA`, `phaB`, `phaG`, `phaJ`), where Pfam membership alone cannot separate true PHA enzymes from close paralogs (e.g. PhaB vs. FabG, PhaA vs. other thiolases). Each candidate is pairwise-aligned against a Swiss-Prot reference set using a length-normalized BLOSUM62 score; candidates below the per-gene calibrated threshold (e.g. ~0.72 for PhaA, ~0.49 for PhaB) are discarded. PhaC is **not** filtered by BLOSUM here — it is validated structurally instead (Layer 4: class HMM + catalytic triad + lipase box), which is a stronger discriminator than pairwise identity.

## 4. Active Site and Catalytic Triad Verification
To prevent the false-positive identification of pseudogenes or non-functional homologs, the most critical enzyme, PHA Synthase (PhaC), is subjected to atomic-level structural verification. 
The algorithm scans the candidate's sequence for the highly conserved "Lipase Box" motif (typically `[G/S] X C X G G`). Furthermore, it verifies the spatial presence of the strictly conserved catalytic triad: Cysteine (C), Aspartate (D), and Histidine (H). Candidates lacking an intact triad are permanently rejected, ensuring that only functional synthases proceed to the next layer.

## 5. Pathway Integration and Machine Learning Verification
A functional PhaC indicates PHA biosynthesis potential, but the polymer type (SCL vs. MCL) is constrained by the PhaC class and the monomer-supplying enzymes present. The potential layer integrates the surviving enzymes with class biochemistry:
- **SCL potential:** functional Class I/III/IV PhaC **and** a monomer route — PhaA + PhaB (sugar → 3HB) or PhaJ (β-oxidation feeding an SCL synthase).
- **MCL potential:** functional Class II PhaC **and** PhaG (de novo fatty-acid route) or PhaJ (β-oxidation route).
- **Uncertain:** a functional PhaC with **no** detectable monomer-supplying gene — the pipeline abstains rather than forcing a positive call.

Note on co-polymers: PhaJ on an SCL-class synthase does **not** by itself imply an SCL-*co*-MCL (e.g. 3HHx) co-polymer — that requires an MCL-capable or broad-substrate synthase — so the current potential layer reports SCL/MCL/uncertain/none and does not emit an SCL-co-MCL call from PhaJ alone. Likewise, PHBV (3HV) potential requires a 3HV-precursor signal (odd-chain VFA feeding plus a C5-accepting thiolase), not merely PhaA/PhaB.

Finally, biophysical properties (Molecular Weight, Isoelectric Point, GRAVY, Instability Index) plus the PhaC HMM bit-score are fed to a Random Forest that provides an **auxiliary** plausibility score. This score is reported alongside the result and may lower confidence or raise a warning, but it **does not override** the deterministic, structure-based decision.
