# PHAscout — Technical Audit

**Date:** 2026-06-30
**Auditor:** automated code+data review (read-only; no code changed)
**Scope:** whole repo on branch `fix/profile-dedup-robustness`
**Verdict in one line:** the *engine* is real, careful, and runs; the *claims wrapped around it* (README, METHODOLOGY, committed benchmark numbers) are overstated and partly describe code that no longer exists. Fix the claims and the leakage before showing this to anyone.

---

## CHANGELOG — fixes applied in this session (2026-06-30)

Small, safe edits (steps 1–3, 6, 7 of §7). Tests: **12 passed** (was 9). No
behavior change to PhaC detection; the only logic change is the PHBV caveat.

| # | Change | Files | Risk |
|---|---|---|---|
| 1 | `requirements.txt` completed (added `click, streamlit, pandas, numpy, joblib, reportlab`; moved `pylint`/`pytest` to a commented dev section) — project now installs and runs from a clean env | `requirements.txt` | none |
| 2 | Streamlit UI: replaced producer verdict **"PHA Üreticisi (Pozitif/Değil)"** with genomic-**potential** tiers (SCL/MCL / belirsiz / yok) read from the honest `pha_potential` field; softened subtitle; added "potential ≠ production, requires validation" caveat | `app.py` (labels only) | none (no logic change) |
| 3 | README rewritten: removed the unsupported "90-genome / 100% precision / flawless" numbers and the "proving novel producers" *Halomonas* claims; reframed as genomic potential + n=15 honest benchmark + explicit caveats. METHODOLOGY corrected to match code (phaZ not scanned, no SCL-co-MCL output, no phaC BLOSUM filter, E-value 1e-5, ML is auxiliary/non-overriding, PHBV needs 3HV precursor) | `README.md`, `METHODOLOGY.md` | none |
| 6 | **PHBV over-claim fixed.** `delta` pathway no longer bare-asserts "P(3HB-co-3HV)" from phaA+phaB. Product string + note now state PHBV is conditional on odd-chain VFA co-substrate + a C5-accepting thiolase that PHAscout does not detect; confidence tagged `CONDITIONAL`; intrinsic tendency reported as P(3HB). Now consistent with the `pha_type` SCL layer. Regression test added | `config.py`, `pathway_engine.py`, `tests/test_pathway_phbv.py` | low (verified end-to-end + tests) |
| 7 | Hygiene check: `.gitignore` already covers `__pycache__`/`*.pyc` and nothing is tracked (the `.pyc` in the tree are untracked cpython-3.14 artifacts) — no action needed | — | none |

### Steps 4, 5, 7 — completed (2026-06-30, second pass, with approval)

- **§3B Reference decontamination (step 4) — DONE (safe variant).** Removed *Synechocystis*
  6803 (`P73825`/`P73826`) and *S. meliloti* 1021 (`P50174`/`P50205`) from
  `positive/phaa.fasta` (9→7) and `positive/phab.fasta` (7→5). **A. vinosum was kept in the
  references** — it is the *defining* Class III sequence and `CATALYTIC_HMM_COLUMNS["Class_III"]`
  was derived from it, so removing it would require rebuilding the Class III HMM and re-deriving
  every triad column (silent breakage risk). Instead, per the benchmark's own iron rule #5
  ("burned genomes aren't held-out"), **A. vinosum was dropped from the benchmark**
  (`dataset.csv` 15→14). Net: all leakage removed for genomes that remain in the benchmark,
  **with zero change to the phaC HMMs or triad columns.**
- **§3D Benchmark re-run (step 5) — DONE.** `python -m benchmark.run --force` re-ran all 14
  genomes (NCBI verified working; 14 OK / 0 errors). Fresh `predictions.csv` (no stale
  `SCL-co-MCL`), fresh `cache/`, results written to `benchmark/RESULTS_run02.md`.
  **Honest result (N=14, regression check — NOT a generalization claim):** recall 1.00
  (95% CI 0.70–1.00, n=9), 0 FP, **3 abstentions**, PHA-type 8/9; specificity meaningless (n=2).
  - The decontamination **worked as intended**: *Synechocystis* flipped SCL → "belirsiz"
    (its old SCL call was leakage-driven), turning a memorized true-positive into an honest
    abstention. *S. meliloti* stayed SCL (robust, not leakage). The PHBV/phaJ fix turned the
    two old `SCL-co-MCL` errors (*P. sacchari*, *B. thuringiensis*) into correct SCL.
- **`pipeline_eski/` — DONE.** Removed the 4.3 MB full-duplicate backup from the working tree
  (88 tracked files + untracked pyc); retained in git history.

### Raw-assembly support (Prodigal) — DONE (2026-06-30, third pass)

- New `phascout/input/gene_caller.py`: auto-detects nucleotide vs protein input
  (alphabet fraction), runs Prodigal (`-p single`, falling back to `-p meta` for
  <100 kb inputs), and extracts per-gene coordinates **from Prodigal's FASTA
  headers** (`# start # end # strand`) into `gff_data` — so raw genomes get the
  same operon/synteny analysis as annotated proteomes.
- Wired into `pipeline.py` `fasta_file`/`fasta_text` branches (auto-detect; no new
  flag needed). README + METHODOLOGY updated.
- Verified end-to-end on a synthetic phaCAB nucleotide contig: nucleotide detected
  → 3 genes called → phaC Class I functional → operon synteny (phaA 193 bp, phaB
  1564 bp from phaC) → SCL potential. New tests `tests/test_gene_caller.py` (5).
  Full suite: **17 passed**.

### Still outstanding (larger roadmap — not started)

- **Single unified prediction layer** (collapse `pathway_engine` + `pha_type`).
- **Expanded negative set** — n=2 cannot support any specificity claim.
- **Archaeal / monomer-route HMM coverage** (Haloferax, Synechocystis abstain on phaA/phaB).
- **Fresh held-out validation** — the current 14 genomes are now burned; the fixes must be
  confirmed on never-touched genomes.
- **ML scorer circularity** (`hmm_score` feature vs. same-source HMMs) — document/redesign.

---

## 0. How this was verified

- Tools confirmed on PATH in conda env `bioinfo` (Python 3.11.15): `python, blastn, blastp, hmmsearch, hmmbuild, bakta, prodigal, datasets, snakemake, diamond, mmseqs, mafft, seqkit` — all present.
- `python -m pytest tests` → **9 passed** (test_pha_type.py, test_triad.py).
- Pipeline executed end-to-end **offline** on a local FASTA (phaC+phaA+phaB) → correct Class I call, triad+box found, alpha pathway active. The core works.
- All `phascout/` modules, both entry points (`cli.py`, `app.py`), the `benchmark/` harness, `scripts/train_ml_scorer.py`, and the data tree were read.

---

## 1. What the project actually does

A 9-step proteome-screening pipeline (`phascout/pipeline.py`):

1. **Input** — NCBI accession (downloads proteome + GFF via Datasets API) or local FASTA.
2. **HMM scan (layer 1, broad net)** — `pyhmmer` vs local Pfam profiles (`hmm_scanner.py`).
3. **Double-layer BLOSUM62 filter (layer 2, specificity)** — for phaA/phaB/phaG/phaJ only; phaC/phaP/phaR/phaE pass through (`double_layer.py`).
4. **Gene vector** — boolean presence of phaC/A/B/J/G/P/R/E.
5. **PhaC classification + functional validation** — class HMMs + catalytic triad (Cys/Asp/His on locked HMM columns) + lipase-box, via `pyhmmer` alignment coordinates (`phac_validator.py`). Multi-candidate ranking favors functional + route-complete synthases.
6. **Operon/synteny evidence** — uses GFF coordinates to support phaA/phaB calls and rescue near-threshold ones (`operon_analyzer.py`).
7. **Subunit check** — Class III→PhaE, Class IV→PhaR (`subunit_checker.py`).
8. **Pathway engine (boolean)** — alpha/beta/gamma/delta/epsilon (`pathway_engine.py`).
9. **ML auxiliary confidence** — Random Forest on HMM bit-score + physicochemical features; **non-overriding** (`ml_scorer.py`).
10. **PHA-type potential** — honest "potential" layer (SCL/MCL/none/belirsiz) (`pha_type.py`).
11. **Report** — JSON + text (`report_generator.py`); Streamlit UI (`app.py`); CLI (`cli.py`).

The codebase is bilingual: code/comments are in Turkish, public docs in English.

---

## 2. What is scientifically strong

- **Functional PhaC validation is genuine.** The catalytic triad is checked at specific, derived HMM match-state columns (`config.CATALYTIC_HMM_COLUMNS`, tolerance ±2), with the lipase-box Cys pinned to the catalytic Cys. This is well above "BLAST said synthase." Pseudogenes/non-functional homologs are correctly demoted.
- **The deterministic decision is correctly separated from ML.** ML adjusts confidence and produces warnings; it never flips the call. This is the right design and is documented in code.
- **`pha_type.py` is careful and correct.** It uses "potential" language, abstains ("belirsiz") when a synthase has no monomer-supplying route, and — importantly — **does not over-claim SCL-co-MCL from phaJ alone** on SCL-class synthases. This is biochemically sound and is regression-tested.
- **The honest benchmark harness (`benchmark/`) is the best part of the repo.** It separates *detection* from *type accuracy*, treats "belirsiz" as abstention (not a forced positive), reports Wilson 95% CIs, breaks down by phylum/class, lists every FP/FN, **forbids annotation-only evidence for positive labels** (anti-circularity), and explicitly warns that the negative set is too weak to trust specificity/precision. The `benchmark/README.md` "iron rules" are genuinely good scientific hygiene.
- **Operon/synteny rescue is biologically motivated** (PhaB↔FabG and PhaA↔other-thiolase ambiguity has a sequence-only ceiling; synteny is the right tiebreaker) and correctly does nothing when no operon/GFF exists.

---

## 3. What is weak, missing, or misleading

### 3A. CRITICAL — Documentation/claims are overstated (the "fake" layer)

These are the items that make the project look dishonest even though the engine is honest. They are the **first thing to fix**.

1. **README headline numbers are not supported by the honest benchmark.**
   - README claims *"93.3% Accuracy, 100% Precision (0 False Positives), 94.7% F1 on 90 Genomes."* Git history: *"100% Accuracy"*, *"Flawless pipeline."*
   - The only methodologically defensible benchmark in the repo (`benchmark/`) is **n=16** (15 evaluated), and its own metrics code prints an explicit warning that the negative set is too small to trust precision/specificity. The "90-genome / 100% precision" numbers come from `scripts/run_massive_benchmark.py` using **NCBI-annotation labels** — which the project's *own iron rules* forbid as circular. **These numbers should be removed or relabeled as non-rigorous.**

2. **Genome → "producer" claims (the exact thing to avoid).**
   - README Phase 3: *"Identified complete, functional P3HB pathways in 3 'No Data' species … proving them to be novel producers."* You cannot *prove* production from a genome. This must become "genomic potential, requires experimental validation."
   - `app.py` Streamlit UI prints **"✅ PHA Üreticisi (Pozitif)" / "❌ PHA Üreticisi Değil"** ("IS / IS NOT a PHA producer"). This is a production verdict, not a potential statement. (The backend `pha_potential` layer is honest; the UI label is not.)
   - `app.py` subtitle: *"…tahmin eder yüksek doğrulukla"* ("predicts with high accuracy") and *"Yapay Zeka Destekli"* ("AI-powered") — marketing language over a 16-genome tool.

3. **METHODOLOGY.md describes behavior the code no longer has.**
   - §5 claims an **SCL-co-MCL** output path; `pha_type.py` can only emit SCL / MCL / none / belirsiz. SCL-co-MCL is **dead in the current prediction path**.
   - §2 lists **phaZ** among scanned genes; `config.py` explicitly removed phaZ.
   - §4 states a **BLOSUM 0.71 filter "for PhaC"**; there is **no BLOSUM filter on phaC at all** (phaC is passthrough). 0.7161 is the *phaA* threshold, misattributed.
   - §2 says **E-value < 1e-3**; code uses **1e-5** (`HMM_EVALUE_THRESHOLD`).
   - Net effect: the methodology oversells and misdescribes. Anyone reproducing from the doc will be confused.

### 3B. CRITICAL — Data leakage into the "held-out" benchmark

The benchmark README claims strict held-out discipline, but the reference/training data contain benchmark organisms:

| Benchmark genome | Leaks into | Evidence |
|---|---|---|
| *Allochromatium vinosum* DSM 180 (`GCF_000025485.1`) | phaA + phaB BLOSUM refs **and** phaC class-HMM training | `P45369`/`P45375` in `positive/phaa,phab.fasta`; 5 hits in phaC training fastas |
| *Synechocystis* PCC 6803 (`GCF_000009725.1`) | phaA + phaB BLOSUM refs | `P73825`/`P73826` (SYNY3) in `positive/phaa,phab.fasta` |
| *Sinorhizobium meliloti* 1021 (`GCF_000006965.1`) | phaA + phaB BLOSUM refs | `P50174`/`P50205` (RHIME) in `positive/phaa,phab.fasta` |

Consequence: phaA/phaB detection (and, for *A. vinosum*, phaC classification) on these genomes is partly memorization, not generalization. The benchmark's detection numbers are inflated for at least 3 of 15 genomes. **The reference sets must exclude any organism that appears in the benchmark** (or those genomes must be dropped from the benchmark).

> Note: the ML scorer has a related circularity — its `hmm_score` feature is computed against the same class HMMs built from the same `phac_classes_clean` positives it trains on. It is auxiliary/non-overriding, so impact is limited, but the reported CV ROC-AUC/MCC are optimistic and should not be quoted as generalization performance.

### 3C. HIGH — PHBV / 3HV logic is an over-claim in the pathway layer

The user's specific concern is real and currently violated:

- `pathway_engine.py` marks the **delta** pathway *active → "Urun egilimi: P(3HB-co-3HV)"* whenever **phaC + phaA + phaB** are present (phaJ is optional). The end-to-end test run confirmed this: a bare phaC/phaA/phaB input prints **"delta … AKTIF → P(3HB-co-3HV)."**
- **3HV requires a propionyl-CoA → 3-ketovaleryl-CoA route** (e.g. a C5-accepting β-ketothiolase such as BktB, plus odd-chain VFA feeding). phaA/phaB alone (C4 acetoacetyl route) do **not** confer 3HV capability. The pipeline checks for **no** 3HV-precursor gene.
- The report is **internally inconsistent**: the honest `pha_type` layer says "SCL" for the same organism while the pathways table says "P(3HB-co-3HV) active." Two layers, two answers, in one report.
- There is **no dedicated PHBV module** and no separate handling for 3HV precursor support. PHBV is currently "phaC present + we listed propionate as a carbon source," which is not evidence of PHBV potential.

### 3D. MEDIUM — Reproducibility / packaging

1. **`requirements.txt` is incomplete — a fresh env cannot run the project.** Imports used but missing from requirements: **`click` (cli.py), `streamlit` + `pandas` (app.py), `joblib` + `pandas` (ml_scorer.py), `numpy`, `reportlab` (scripts)**. Listed `pylint` is a dev tool, not a runtime dep. As written, `pip install -r requirements.txt` then `python cli.py` fails on `import click`.
2. **`__pycache__` is committed and compiled with cpython-3.14**, but the runtime env is 3.11. The committed `benchmark/predictions.csv` (see below) was produced by that other interpreter/older code.
3. **`predictions.csv` is stale.** It contains `pred_potential = SCL-co-MCL` for `GCF_000785435.2` and `GCF_000092165.1`, a value the current `pha_type.py` can never emit. The committed benchmark output does not match the committed code → re-run required before any number is quoted.
4. **`pipeline_eski/` is a full duplicate of the project** (second copy of `phascout/`, `data/`, etc.) committed as a backup. It roughly doubles repo size and will confuse readers/graders. Should be removed from the working tree (history retains it).

### 3E. LOW — Smaller correctness / doc nits

- `subunit_checker.py` docstring says PhaE = **PF08333**; `config.py` uses **PF09712**. (Code path uses config; doc is wrong.)
- `benchmark/metrics.py`: `--dataset` override sets module-global `DATASET` but `PREDICTIONS` stays at the default path — fine for default use, surprising for smoke runs.
- `hmm_scanner.download_pfam_hmm` reaches InterPro at runtime if a profile is missing. All 11 Pfam HMMs are present locally, so offline runs work today, but a missing file silently triggers a network call — should be explicit/optional for a reproducible build.
- Heavy reliance on a single `numbered` pathway taxonomy (alpha…epsilon) duplicated between `config.PATHWAYS` and the `pha_type` logic; the two can (and now do) disagree.

---

## 4. Tooling reality check (installed vs integrated)

| Tool | Installed | Integrated in pipeline? |
|---|---|---|
| **HMMER / pyhmmer** | ✅ | ✅ core (HMM scan + triad alignment) |
| **BLOSUM62 (Biopython)** | ✅ | ✅ double-layer filter |
| **Biopython** | ✅ | ✅ parsing, ProtParam, alignment |
| **scikit-learn / joblib** | ✅ | ✅ RF auxiliary scorer |
| **NCBI Datasets (REST)** | ✅ | ✅ input (via `requests`, not the `datasets` CLI) |
| **BLAST** | ✅ | ❌ not used by pipeline (referenced only in prose) |
| **Bakta** | ✅ | ❌ not used (proteomes come pre-annotated from NCBI) |
| **Prodigal** | ✅ | ✅ integrated (raw nucleotide genome input → gene calling + coords; `gene_caller.py`) |
| **DIAMOND / MMseqs2** | ✅ | ❌ not used in pipeline (MMseqs only in a clustering helper script) |
| **MAFFT** | ✅ | ❌ used only offline to build HMMs (`scripts/`), not at runtime |
| **Snakemake** | ✅ | ❌ no workflow file exists |
| **hmmbuild** | ✅ | ⚠️ used offline in `scripts/` to build class HMMs, not at runtime |

**Takeaway:** the pipeline genuinely integrates HMMER + BLOSUM + Biopython + sklearn + NCBI. The other 7 tools are installed but unused by the runtime path. That's fine — but the docs should not imply a BLAST/Bakta/Prodigal/Snakemake pipeline that doesn't exist. The biggest missing capability is **raw-genome → gene-calling (Prodigal/Bakta)**: today the tool requires a pre-annotated proteome, so it cannot screen unannotated assemblies despite the README implying "raw genomic sequences."

---

## 5. PHA gene logic review (per the requested gene list)

| Gene | Status | Note |
|---|---|---|
| **phaC** | ✅ strong | class HMM + triad (locked columns) + box; multi-candidate ranking; passthrough on BLOSUM (validated structurally instead) |
| **phaA** | ✅ + leakage | PF00108/PF02803 + BLOSUM 0.716 + length + synteny rescue. Refs leak benchmark organisms (§3B). |
| **phaB** | ✅ + leakage | PF13561 + BLOSUM 0.49 + synteny. FabG confusion handled by synteny. Refs leak (§3B). |
| **phaP** | ⚠️ thin | PF09361, HMM-only passthrough, not used in decision |
| **phaR** | ⚠️ dual role | PF07879; used as Class IV subunit. Also nominally a regulator — only the subunit role is wired in |
| **phaE** | ✅ | PF09712; Class III subunit check |
| **phaZ** | ❌ removed | dropped from config (PF05898 returned no HMM); still advertised in METHODOLOGY |
| **phaJ** | ✅ corrected | PF01575 (MaoC, correct (R)-specific hotdog fold). Old PF00767/crotonase was wrong and caught nothing. phaJ no longer over-claims co-polymer — good. |
| **phaG** | ✅ corrected | PF00561 only (PF07167 removed because it caught the synthase). Diversity-lowered threshold 0.35. |

phaC/A/B/J/G logic is sound and the recent fixes (phaJ profile, phaG profile, phaJ co-polymer over-claim) are correct, biochemically reasoned, and tested.

---

## 6. Brutally honest status

**What already works**
- Functional PhaC detection with real catalytic-triad verification.
- Class I–IV classification on canonical inputs.
- The honest "potential + abstention" framing in `pha_type.py` and the report summary.
- The benchmark harness methodology (detection vs type, abstention, Wilson CI, anti-circularity rules).
- End-to-end run offline; tests pass.

**What is fake / overclaimed**
- README "93.3%/100% precision/90 genomes" and "100% accuracy / flawless" — not supported; partly from annotation-labeled (circular) runs.
- "Proving novel producers" and the Streamlit "IS / IS NOT a PHA producer" verdict — genome→producer claims you explicitly want avoided.
- METHODOLOGY describing SCL-co-MCL output, phaZ scanning, and a phaC BLOSUM filter — none exist in current code.
- "Held-out" benchmark — leaks ≥3 organisms into the reference/training data.
- `predictions.csv` — stale, doesn't match current code.

**What must be fixed first (blockers for honesty/credibility)**
1. Strip/relabel the unsupported accuracy claims (README, METHODOLOGY, PDF/Table3).
2. Replace producer verdicts with potential language in `app.py` UI (backend is already honest).
3. Remove benchmark organisms from reference/training fastas (or drop them from the benchmark) and **re-run** `benchmark/run.py` + `metrics.py`; commit fresh `predictions.csv`.
4. Fix the PHBV/3HV over-claim: the delta pathway must not assert P(3HB-co-3HV) without a 3HV-precursor signal; make the pathway layer consistent with `pha_type`.
5. Fix `requirements.txt` so the project installs and runs.

**What can wait**
- Removing `pipeline_eski/` and committed `__pycache__` (hygiene).
- Adding Prodigal/Bakta for raw-genome input (real feature, not a correctness bug).
- Snakemake workflow (nice-to-have; currently unused).
- ML scorer circularity (it's auxiliary; document the caveat now, redesign later).
- Doc nits (subunit PF id, E-value, metrics `--dataset`).

**What a publishable / portfolio-level version needs**
- A genuinely held-out, wet-lab-labeled benchmark of meaningful size (the current iron-rule framework is the right basis — just needs more genomes and zero reference overlap), with the negative set expanded (the current n≈3 negatives cannot support any specificity claim).
- One consistent prediction layer (collapse `pathway_engine` + `pha_type` so the report never self-contradicts).
- Explicit PHBV handling tied to 3HV-precursor evidence.
- Honest, reproducible docs whose every claim maps to a re-runnable command and a committed result.
- Optional but valuable: raw-assembly support (Prodigal), and a Snakemake/CI wrapper so reviewers can reproduce in one command.

---

## 7. Minimal, safe next-action plan (in order; small diffs)

> Principle: fix honesty and reproducibility before features. Each step is independently committable.

1. **`requirements.txt`** — add `click`, `streamlit`, `pandas`, `numpy`, `joblib`, `reportlab`; move `pylint` to a dev section/comment. *(2-line-ish change, unblocks install.)*
2. **`app.py` labels only** — change "PHA Üreticisi (Pozitif/Değil)" → "PHA genomik potansiyeli: VAR / BELİRSİZ / YOK" and soften the subtitle ("predicts with high accuracy" → "estimates genomic potential; requires experimental validation"). No logic change.
3. **README + METHODOLOGY** — remove the 90-genome/100%-precision/"flawless"/"proving novel producers" claims; replace with the honest n=16 framing + caveats; fix the phaZ / SCL-co-MCL / phaC-BLOSUM / E-value descriptions to match code.
4. **Decontaminate references** — drop *A. vinosum*, *Synechocystis* 6803, *S. meliloti* 1021 entries from `positive/phaa.fasta`, `positive/phab.fasta`, and the phaC class training fastas (or remove those 3 genomes from `benchmark/dataset.csv`). Document the choice.
5. **Re-run the honest benchmark** and commit fresh `predictions.csv` + a short results note; delete the stale CSV.
6. **PHBV consistency fix (smallest correct version)** — in `pathway_engine.py`, stop emitting `delta` as an active P(3HB-co-3HV) pathway on phaA+phaB alone; either gate it behind a 3HV-precursor signal or downgrade it to a clearly-labeled conditional ("P(3HB); P(3HB-co-3HV) only with odd-chain VFA co-substrate + C5-accepting thiolase — not detected"). Align wording with `pha_type`. Add a regression test mirroring `test_pha_type.py`.
7. **Hygiene** — `git rm -r --cached` the `__pycache__` dirs, add to `.gitignore`; decide on `pipeline_eski/` (recommend removing from working tree).

Steps 1–3 are pure-honesty wins with near-zero risk and should go first. Steps 4–6 are the scientific-integrity fixes. Step 7 is cleanup.
