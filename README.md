# smallRNAseq-DREAM

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Zenodo DOI: pending](https://img.shields.io/badge/Zenodo%20DOI-pending-lightgrey.svg)](data/MANIFEST.md)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Data Source: NASA OSDR](https://img.shields.io/badge/Data%20Source-NASA%20OSDR-red.svg)](https://osdr.nasa.gov)

**A FAIR, cross-species small RNA-seq (miRNA) analysis framework, with an inventory of NASA OSDR plant RNA-seq datasets.**

🌐 **Project site:** https://dr-richard-barker.github.io/smallRNAseq-DREAM/

---

## Overview

**smallRNAseq-DREAM** consolidates small RNA-seq analysis work: four analysis engines, a
benchmarking harness, NASA OSDR retrieval tools, and an exploratory scan of OSDR plant RNA-seq
tables for miRNA-related signal.

> **Status of the plant miRNA scan (evidence check, 2026-10-28): its outputs are not interpretable yet.**
>
> - **OSDR holds no plant small RNA-seq** ([`osdr/README.md`](osdr/README.md)). OSD-437 is
>   standard poly(A) RNA-seq: its OSDR record lists one RNA-Seq assay (19 samples, library
>   selection "polyA enrichment", NEB Next Ultra II RNA Library Prep Kit), not small RNA-seq.
> - **The locus map mislabels most IDs.** `MIRNA_LOCI_MAP` in `ncrna_quantifier.py` is hand-typed.
>   Of the 29 IDs that reach `mirna_spaceflight_profiles.csv`, 3 are the named MIR gene; 16 are
>   protein-coding (e.g. AT1G53160 is *SPL4*, labelled "MIR156b"), 2 are lncRNA, 1 is snoRNA and
>   7 are a different MIR gene (Ensembl Plants REST lookup). The "miRNA" rows are therefore mostly
>   mRNA values.
> - **The log2FC values are not flight vs ground.** The scripts take the first `Log2fc_` column of
>   each GeneLab differential-expression table, which is:
>   OSD-437 `(4 day & 1G by centrifugation & Space Flight)v(4 day & uG & Space Flight)` (both arms in flight);
>   OSD-120 `(Col-0 & Ground Control & Dark Treatment)v(Col-0 & Ground Control & Light Treatment)`;
>   OSD-218 `(4 day & Ground Control & Col-0)v(4 day & Ground Control & WS)`;
>   OSD-217 `(leaf & Ground Control)v(leaf & Space Flight)` (the reverse orientation);
>   OSD-565 `(Plant Leaves & Treatment 1 & Not Applicable)v(Plant Leaves & Treatment 2 & Not Applicable)`;
>   OSD-782 `(cesium-137 gamma radiation & 10 centigray & 1 hour)v(cesium-137 gamma radiation & 10 centigray & 24 hour)`.
> - **TODO:** rebuild the locus map from Ensembl Plants genes with biotype `miRNA`; select the
>   flight-vs-ground contrast in each table by column name; then regenerate every output and figure.
>   Until then, no miRNA result from this scan should be cited.

What the repository contains:

1. Code that indexes 12 NASA OSDR plant RNA-seq studies (table below, read from the OSDR study
   records): 9 spaceflight, 1 suborbital (high-altitude) and 2 ground studies, across
   *Arabidopsis thaliana*, *Brachypodium distachyon*, *Solanum lycopersicum* and *Brassica rapa*.
2. A scan of GeneLab differential-expression tables for IDs in the typed MIR map (see status above),
   written to `mirna_spaceflight_profiles.csv` (122 rows = 29 IDs across 6 studies).
3. 24 hand-curated literature miRNA–target pairs (typed in `mirna_target_network.py`; no target
   database was queried), joined to the same tables in five studies to give 115 study-level pairs.
   In 45 of the 115, the "miRNA" row and the target have opposite-sign log2FC. No correlation is claimed.
4. An **illustrative** Fisher test per pathway. Its background is typed, not measured (25,000 genes,
   3,500 "DE" genes), and each target is counted once per study, so the p-values are not an
   enrichment result.
5. Seven figures regenerated from those tables. Figures 5 and 7 are schematics.

---

## The Decision & Analysis Engines

From hands-on testing of four small RNA-seq analysis engines (`nf-core/smrnaseq`, `sRNAtoolbox`,
`miRDeep2`, `miRDeep-P2`), the author's recommended route for cross-species spaceflight analysis is
below ([`docs/00_pipeline_decision.md`](docs/00_pipeline_decision.md)). The quantitative four-engine
benchmark has not been run yet ([`benchmarking/README.md`](benchmarking/README.md)).

```
SRA  →  mirnaQC  →  sRNAbench  →  sRNAde  →  mirNOVO
```

### Four Engine Comparison

| Engine | Kingdom | Path |
|---|---|---|
| `nf-core/smrnaseq` | Both | [`pipelines/nfcore_smrnaseq/`](pipelines/nfcore_smrnaseq/) |
| `sRNAtoolbox` (`sRNAbench` / `sRNAde` / `mirNOVO`) | Both | [`pipelines/srnatoolbox/`](pipelines/srnatoolbox/) |
| `miRDeep2` | Animal | [`pipelines/mirdeep2_animal/`](pipelines/mirdeep2_animal/) |
| `miRDeep-P2` (`miRDP2`) | Plant | [`pipelines/mirdp2_plant/`](pipelines/mirdp2_plant/) |

---

## NASA OSDR Plant Dataset Inventory

Read from each study's OSDR record (`https://osdr.nasa.gov/osdr/data/osd/meta/<n>`) on 2026-10-28.
"RNA-seq samples" counts the samples in the study's RNA-Seq assay. A blank cell means the OSDR record
has no value for that field.

| Study | Species | OSDR project type | Mission | Platform | RNA-seq samples | Library selection (ISA) |
| :--- | :--- | :--- | :--- | :--- | :---: | :--- |
| **OSD-437** | *Arabidopsis thaliana* | Spaceflight Study | SpaceX-9 | European Modular Cultivation System (EMCS), Experimental Reference Module (ERM2) | 19 | polyA enrichment |
| **OSD-120** | *Arabidopsis thaliana* | Spaceflight Study | SpaceX-3 |  | 36 | ribo-depletion |
| **OSD-218** | *Arabidopsis thaliana* | Spaceflight Study | SpaceX-5 | Vegetable Production System (VPS/Veggie) | 32 | polyA enrichment |
| **OSD-217** | *Arabidopsis thaliana* | Spaceflight Study | SpaceX-5 | Vegetable Production System (VPS/Veggie) | 12 | cDNA |
| **OSD-375** | *Brachypodium distachyon* | Spaceflight Study | SpaceX-14 | Vegetable Production System (VEGGIE), APEX Growth Unit | 48 | Ribo-depletion (32), cDNA (16) |
| **OSD-480** | *Arabidopsis thaliana* | Spaceflight Study | SpaceX-13 | European Modular Cultivation System (EMCS) | 400 |  |
| **OSD-522** | *Arabidopsis thaliana* | Spaceflight Study | SpaceX-13, SpaceX-14 | BRIC-LED | 12 | polyA enrichment |
| **OSD-565** | *Arabidopsis thaliana* | Ground Study |  | KSC Fixation Tubes (KFTs) | 45 | cDNA |
| **OSD-624** | *Arabidopsis thaliana* | High Altitude Study | Unity 22 | Kennedy Space Center Fixation Tubes (KFTs) | 36 | polyA enrichment |
| **OSD-767** | *Solanum lycopersicum* | Spaceflight Study | SpaceX-26 | Vegetable Production System (Veggie) | 36 | polyA enrichment |
| **OSD-782** | *Arabidopsis thaliana* | Ground Study |  |  | 36 | polyA enrichment |
| **OSD-59** | *Brassica rapa* | Spaceflight Study | Expedition 22 | Lada | 2 | mRNA enrichment |

---

## Repository Structure

```
smallRNAseq-DREAM/
├── run_pipeline.py                    # Driver: runs the four stages below in order
├── osdr_api_client.py                 # Stage 1: OSDR study inventory + GeneLab table URLs
├── ncrna_quantifier.py                # Stage 2: scan of DE tables for IDs in the typed MIR map (see Status)
├── mirna_target_network.py            # Stage 3: curated miRNA-target pairing + illustrative Fisher test
├── build_dream_manuscript_figures.py  # Stage 4: Figures 1-7
├── plant_osdr_inventory.csv / .json   # Stage 1 output
├── mirna_spaceflight_profiles.csv     # Stage 2 output
├── mirna_mrna_regulatory_network.csv  # Stage 3 output
├── target_gene_pathway_enrichment.csv # Stage 3 output (illustrative)
├── figures/                           # Figures 1-7 (Figures 5 and 7 are schematics)
├── osdr/                              # OSDR search/download tool + small-RNA survey
├── meta-analysis/                     # Cross-study miRNA count-matrix tools + cohort table
├── smallrna_from_rnaseq/              # Precursor-level signal from standard plant RNA-seq
├── pipelines/                         # Code and notes for each analysis engine
├── target_prediction/                 # miRNA target prediction (BLAST) draft
├── benchmarking/                      # Synthetic ground-truth accuracy harness (engines not yet run)
├── docs/                              # Pipeline decision, methods tables, manuscript outline
├── data/                              # MANIFEST.md only; OSDR downloads land in data/OSD-*/ (gitignored)
├── tables/                            # 7 CSVs not referenced by any code or doc here (TODO: remove or label)
├── third_party/                       # Notes on third-party tool licences
├── requirements.txt / environment.yml
├── .zenodo.json / CITATION.cff
└── LICENSE                            # MIT License
```

---

## Quickstart & Execution Guide

### Prerequisites
Python 3.10+ is required. Set up a virtual environment and install dependencies:

```bash
git clone https://github.com/dr-richard-barker/smallRNAseq-DREAM.git
cd smallRNAseq-DREAM
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Running the Complete Pipeline
To execute the complete analysis end-to-end:

```bash
python run_pipeline.py
```

Or execute individual pipeline stages:

```bash
# Stage 1: Mine OSDR API metadata & data matrices
python osdr_api_client.py

# Stage 2: Scan DE tables for IDs in the typed MIR map (not interpretable yet; see Status)
python ncrna_quantifier.py

# Stage 3: Curated miRNA-target pairing & illustrative Fisher test
python mirna_target_network.py

# Stage 4: Render publication figures (Figures 1-7)
python build_dream_manuscript_figures.py
```

---

## Analytical Methodology

### 1. Data retrieval
Study metadata and file listings come from the OSDR API (`https://osdr.nasa.gov/osdr/data/osd/meta/<n>`
and `https://osdr.nasa.gov/osdr/data/osd/files/<n>`). The GeneLab-processed differential-expression,
count and sample tables are downloaded into `data/OSD-*/` (gitignored). No read alignment or
normalisation is done in this repository: every log2FC and p-value is GeneLab's, taken from the
first `Log2fc_` / `P.value_` column of each table (see Status for what those contrasts are).

### 2. Locus scan and target pairing
`ncrna_quantifier.py` keeps rows whose ID is in the typed `MIRNA_LOCI_MAP` (mostly mislabelled; see
Status). `mirna_target_network.py` joins 24 typed miRNA–target pairs to the target genes' rows in five
studies (115 pairs), using the mean of the matching "miRNA" rows as the miRNA-side value. A pair is
marked "Inverse" when the two values have opposite signs (45 of 115).

### 3. Illustrative Fisher test
For each pathway, `scipy.stats.fisher_exact` compares the pairs whose target has P ≤ 0.05 against a
typed background (25,000 genes, 3,500 "DE" genes, constants in `mirna_target_network.py`). Because the
background is assumed and each target is counted once per study, the p-values in
`target_gene_pathway_enrichment.csv` and Figure 4 are illustrative, not an enrichment result.
TODO: use each study's measured gene universe and DE set.

---

## Provenance & Consolidation

This framework consolidates work from five earlier repositories by the author, which are no longer
publicly available (each returned HTTP 404 to a public request on 2026-10-28):
`smallRNAseq_mirDeep2_NFcore_DRB`, `sRNA_toolbox`, `mirDeep2_accuracy_synthetic_microRNA`,
`miRDeep-P2_pipeline` (a fork of [`TF-Chan-Lab/miRDeep-P2_pipeline`](https://github.com/TF-Chan-Lab/miRDeep-P2_pipeline))
and `smallRNA_targets`.

---

## Data Availability & Zenodo Archiving

All input data are public NASA Open Science Data Repository (OSDR) datasets. A Zenodo deposit of the
benchmark bundle is prepared but not yet made: **Zenodo DOI: pending** (see
[`data/MANIFEST.md`](data/MANIFEST.md)).

## License
Distributed under the **MIT License**. See `LICENSE` for more information.
