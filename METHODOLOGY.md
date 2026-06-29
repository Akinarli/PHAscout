# PHAscout Methodology

The PHAscout pipeline is a de novo, open-source computational tool designed for the genome-wide identification and classification of polyhydroxyalkanoate (PHA) biosynthesis pathways. Operating without prior literature-based biases, the pipeline determines the specific PHA production potential (SCL, MCL, or SCL-co-MCL) directly from raw genomic or proteomic sequences. The methodology consists of five sequential analytical layers.

## 1. Data Acquisition and Preprocessing
The pipeline accepts a GenBank Assembly Accession number (e.g., GCF_030505695.1). Utilizing the NCBI Datasets API, the system automatically retrieves the complete proteome (amino acid sequences) and corresponding GFF3 annotation files. The proteome is parsed using `Biopython`, rendering thousands of uncharacterized protein sequences ready for high-throughput screening.

## 2. Hidden Markov Model (HMM) Scanning
To identify putative PHA biosynthesis enzymes without relying on sequence identity tools like BLAST (which are prone to false positives due to generic hydrolases), PHAscout employs probabilistic Hidden Markov Models (HMMs). 
Profile HMMs for essential PHA genes (`phaC`, `phaA`, `phaB`, `phaG`, `phaJ`, `phaP`, `phaR`, `phaE`, `phaZ`) were constructed or retrieved from the Pfam database. The `pyhmmer` library is utilized to scan the entire proteome against these profiles. Proteins exhibiting statistically significant structural homology (E-value < 1e-3) are selected as preliminary candidates.

## 3. The Double-Layer BLOSUM62 Filter
Preliminary HMM candidates undergo a stringent evolutionary scoring filter. A high-confidence reference database of experimentally validated PHA enzymes was compiled from the Swiss-Prot (UniProt) database. Each candidate is pairwise-aligned against this reference dataset using the BLOSUM62 substitution matrix.
This matrix evaluates evolutionary divergence by penalizing structurally disruptive amino acid substitutions. Candidates failing to meet strict baseline similarity scores (e.g., 0.71 for PhaC) are discarded as false positives (e.g., paralogous lipases or generic enoyl-CoA hydratases).

## 4. Active Site and Catalytic Triad Verification
To prevent the false-positive identification of pseudogenes or non-functional homologs, the most critical enzyme, PHA Synthase (PhaC), is subjected to atomic-level structural verification. 
The algorithm scans the candidate's sequence for the highly conserved "Lipase Box" motif (typically `[G/S] X C X G G`). Furthermore, it verifies the spatial presence of the strictly conserved catalytic triad: Cysteine (C), Aspartate (D), and Histidine (H). Candidates lacking an intact triad are permanently rejected, ensuring that only functional synthases proceed to the next layer.

## 5. Pathway Integration and Machine Learning Verification
The presence of a functional PhaC indicates PHA biosynthesis potential, but the polymer type (SCL vs. MCL) is dictated by monomer-supplying enzymes. The "Pathway Engine" integrates the surviving enzymes:
- **SCL Potential (Alpha Pathway):** Triggered by the co-occurrence of PhaA (β-ketothiolase) and PhaB (acetoacetyl-CoA reductase).
- **MCL Potential (Beta/Gamma Pathways):** Triggered by the presence of PhaG (3-hydroxyacyl-ACP thioesterase) or PhaJ (enoyl-CoA hydratase).
- **SCL-co-MCL Potential:** Triggered when both SCL and MCL monomer-supplying routes are simultaneously functional alongside a compatible PhaC class.

Finally, the biophysical properties (Molecular Weight, Isoelectric Point, GRAVY, Instability Index) of the enzymes are calculated and fed into a pre-trained **Random Forest Classifier**. This Machine Learning model acts as a secondary biological plausibility check, ensuring the biochemical parameters of the discovered enzymes align with known PHA producers, yielding a final confidence score for the predicted PHA type.
