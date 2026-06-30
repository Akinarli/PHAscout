# PHAscout

PHAscout is a de novo bioinformatics pipeline that screens bacterial/archaeal
**proteomes** for the **genomic potential** to synthesize polyhydroxyalkanoates
(PHA), classifies the PHA synthase (PhaC) class, and infers candidate
biosynthetic routes.

> ⚠️ **Scope and honesty statement**
> PHAscout reports **genomic potential**, not production. The presence of
> biosynthesis genes does **not** prove that an organism accumulates PHA — real
> production depends on gene expression, carbon source, and growth conditions,
> and **requires experimental validation** (GC-MS / NMR / FTIR / Nile-red /
> gravimetric / TEM). No claim in this repository should be read as "organism X
> produces PHA."

## What it does

Given an NCBI assembly accession (downloads the annotated proteome) or a local
protein FASTA, the pipeline:

1. Scans proteins against Pfam HMM profiles (`pyhmmer`) for phaC/A/B/J/G/P/R/E.
2. Applies a BLOSUM62 "double-layer" specificity filter to phaA/phaB/phaG/phaJ
   (separating, e.g., PhaB from FabG, PhaA from other thiolases).
3. Validates PhaC **function** via the catalytic triad (Cys/Asp/His, checked at
   class-specific HMM alignment columns) and the lipase-box motif.
4. Uses GFF-based **operon/synteny** evidence to support phaA/phaB calls.
5. Infers a PHA-type **potential** (SCL / MCL / none / "uncertain" when a
   functional synthase has no detectable monomer-supplying route).
6. Produces a JSON/text report (CLI) or an interactive dashboard (Streamlit).

A Random Forest provides an **auxiliary** confidence score only; it never
overrides the deterministic, structure-based decision.

## Current status (honest)

This is a research prototype. The trustworthy evaluation is the held-out
harness in [`benchmark/`](benchmark/), currently **14 wet-lab-labeled genomes**
(detection + PHA-type, with abstention and Wilson 95% CIs). See
[`benchmark/README.md`](benchmark/README.md) for the labeling "iron rules"
(wet-lab evidence only; annotation labels are forbidden as circular) and
[`benchmark/RESULTS_run02.md`](benchmark/RESULTS_run02.md) for the latest run.

Latest run (N=14, **not** a generalization claim — these genomes are now "burned"):
detection recall 1.00 (95% CI 0.70–1.00, n=9 positives), 0 false positives,
**3 abstentions** on hard cases, PHA-type accuracy 8/9. Specificity is reported
but **meaningless** (only n=2 negatives). These numbers must be re-validated on a
fresh held-out set before being quoted anywhere.

Known limitations are tracked in [`PROJECT_AUDIT.md`](PROJECT_AUDIT.md),
including:

- The negative set is small — specificity/precision figures are **not** yet
  trustworthy and should not be quoted as headline accuracy.
- Input must be a **pre-annotated proteome**; raw-assembly gene calling
  (Prodigal/Bakta) is not yet integrated.
- PHBV (3HV) potential requires a 3HV-precursor signal, which is being made
  explicit (see the audit).

> Earlier versions of this README reported large-scale accuracy figures
> (e.g. a "90-genome, 100% precision" benchmark) and described specific
> *Halomonas* strains as "novel producers." Those claims relied on
> annotation-based labels (circular by the project's own rules) and overstated
> genome-only evidence as proof of production; they have been removed. Any
> *Halomonas* findings are **genomic-potential candidates that require
> experimental validation**, not demonstrated producers.

## Installation

```bash
pip install -r requirements.txt
```

External tools used at runtime come from the `bioinfo` conda environment
(HMMER via `pyhmmer`, Biopython). BLAST/Bakta/Prodigal/DIAMOND/MMseqs2/MAFFT/
Snakemake are present in the environment but are **not** part of the runtime
pipeline (MAFFT/hmmbuild are used offline in `scripts/` to build HMMs).

## Usage

```bash
# By NCBI accession (downloads proteome)
python cli.py --accession GCF_000009285.1

# By local protein FASTA
python cli.py --fasta proteins.faa --out report.json

# Web dashboard
streamlit run app.py
```

## Methodology

See [METHODOLOGY.md](METHODOLOGY.md) for the algorithms (HMM scan, BLOSUM62
filter, catalytic-triad verification, pathway/potential logic).

## Provenance

The architecture and core framework were drafted with AI assistance and then
re-engineered (HMM profile fixes, threshold calibration, catalytic-triad
locking, honest benchmark harness). See `git log` for the full history.
