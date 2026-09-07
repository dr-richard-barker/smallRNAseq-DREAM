# smallRNAseq-DREAM

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.smallRNAseq--DREAM-blue.svg)](https://zenodo.org)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Data Source: NASA OSDR](https://img.shields.io/badge/Data%20Source-NASA%20OSDR-red.svg)](https://osdr.nasa.gov)

**A FAIR, cross-species small RNA-seq (miRNA) analysis framework and NASA OSDR spaceflight audit for spaceflight omics.**

🌐 **Project site:** https://dr-richard-barker.github.io/smallRNAseq-DREAM/

---

## Overview & Abstract

The **smallRNAseq-DREAM** project presents a comprehensive computational framework and empirical audit of microRNAs (miRNAs) and non-coding RNAs (ncRNAs) across plant spaceflight biological data from the **NASA Open Science Data Repository (OSDR)**.

While transcriptomic profiling in spaceflight has historically focused on protein-coding mRNAs, small non-coding RNAs play vital post-transcriptional regulatory roles in plant acclimation to microgravity, space radiation, and spaceflight hardware stressors. This repository provides an end-to-end reproducible bioinformatic workflow that:
1. Programmatically mines dataset inventories and metadata from 12 core NASA OSDR plant spaceflight studies (*Arabidopsis thaliana*, *Solanum lycopersicum*, and *Brachypodium distachyon*).
2. Profiles 122 microRNA/ncRNA loci fold-changes and abundance metrics across dedicated small RNA-Seq (**OSD-437**) and bulk spaceflight transcriptomes.
3. Maps 115 empirically verified miRNA-target gene interaction pairs with anti-correlated log fold-change coupling ($r = -0.842$).
4. Evaluates functional pathway enrichment using rigorous SciPy Fisher Exact Tests.
5. Renders seven 300-DPI publication-grade figures, including multimodal knowledge graphs, Sankey flow diagrams, and spatial pathway coupling maps.

---

## The Decision & Analysis Engines

After benchmarking four small RNA-seq analysis engines (`nf-core/smrnaseq`, `sRNAtoolbox`, `miRDeep2`, `miRDeep-P2`), the optimal route for cross-species spaceflight analysis is:

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

## NASA OSDR Plant Spaceflight Dataset Inventory

| Study Accession | Species | Spaceflight Mission / Hardware | Sample Count | Primary Assay Modality |
| :--- | :--- | :--- | :---: | :--- |
| **OSD-437** | *Arabidopsis thaliana* | SpaceX-14 (Advanced Plant Habitat / APH) | 12 | **Small RNA-Seq** (Dedicated miRNA Profiling) |
| **OSD-120** | *Arabidopsis thaliana* | SpaceX-3 (BRIC-18 / Petri dish) | 12 | Poly(A)-Selected Bulk RNA-Seq |
| **OSD-218** | *Arabidopsis thaliana* | SpaceX-10 (APEX-04 / Petri dish) | 12 | Poly(A)-Selected Bulk RNA-Seq |
| **OSD-217** | *Arabidopsis thaliana* | ISS (EMCS / Tropi-2) | 16 | Total RNA-Seq with rRNA Depletion |
| **OSD-375** | *Arabidopsis thaliana* | ISS (APH / Plant Habitat-01) | 18 | Total RNA-Seq with rRNA Depletion |
| **OSD-480** | *Arabidopsis thaliana* | SpaceX-21 (BRIC-24 / Lighted) | 12 | Poly(A)-Selected Bulk RNA-Seq |
| **OSD-522** | *Arabidopsis thaliana* | ISS (Veggie / VEG-03D) | 16 | Poly(A)-Selected Bulk RNA-Seq |
| **OSD-565** | *Arabidopsis thaliana* | SpaceX-23 (APEX-08) | 12 | Total RNA-Seq with rRNA Depletion |
| **OSD-624** | *Arabidopsis thaliana* | ISS (Veggie / VEG-03E) | 16 | Poly(A)-Selected Bulk RNA-Seq |
| **OSD-767** | *Arabidopsis thaliana* | ISS (EMCS / MULTIGEN-1) | 14 | Poly(A)-Selected Bulk RNA-Seq |
| **OSD-782** | *Solanum lycopersicum* | ISS (Veggie / VEG-05 Tomato) | 12 | Multi-Omics / Bulk RNA-Seq |
| **OSD-59** | *Brachypodium distachyon* | Space Flight (BRIC-16) | 8 | Bulk Transcriptomics |

---

## Repository Structure

```
smallRNAseq-DREAM/
├── run_pipeline.py                 # Master execution driver (runs steps 1-4)
├── osdr_api_client.py              # Stage 1: OSDR REST API dataset mining & inventorying
├── ncrna_quantifier.py             # Stage 2: Small RNA & microRNA spaceflight profiling
├── mirna_target_network.py         # Stage 3: miRNA-target pairing & Fisher exact enrichment
├── build_dream_manuscript_figures.py # Stage 4: 300-DPI publication figure generator
├── docs/                           # Pipeline decision, methods tables, engine comparison
├── pipelines/                      # Code for each analysis engine
├── target_prediction/              # miRNA target prediction & BLAST matching
├── benchmarking/                   # Synthetic / accuracy testing harness
├── data/                           # Processed count matrices & dataset inventory
│   ├── plant_osdr_inventory.csv
│   ├── mirna_spaceflight_profiles.csv
│   ├── mirna_mrna_regulatory_network.csv
│   └── target_gene_pathway_enrichment.csv
├── figures/                        # Publication figures (Figures 1-7)
│   ├── Figure1_OSDR_Plant_Dataset_Landscape.png
│   ├── Figure2_microRNA_Differential_Expression.png
│   ├── Figure3_miRNA_mRNA_Regulatory_Coupling.png
│   ├── Figure4_miRNA_Functional_Pathway_Enrichment.png
│   ├── Figure5_miRNA_Target_Ontology_Sankey.png
│   ├── Figure6_Gene_Abundance_Knowledge_Graph.png
│   └── Figure7_KEGG_Pathway_RNA_Coupling.png
├── requirements.txt                # Python dependencies
├── README.md                       # Repository documentation & execution guide
├── .zenodo.json                    # Zenodo DOI metadata
├── CITATION.cff                    # Citation file format
└── LICENSE                         # MIT License
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

# Stage 2: Quantify microRNA/ncRNA spaceflight abundance
python ncrna_quantifier.py

# Stage 3: Map miRNA-mRNA target networks & run Fisher exact enrichment
python mirna_target_network.py

# Stage 4: Render publication figures (Figures 1-7)
python build_dream_manuscript_figures.py
```

---

## Analytical Methodology

### 1. Data Mining & Normalization
Dataset metadata were programmatically pulled from the NASA OSDR REST API (`https://osdr.nasa.gov/osdr/data/study/meta/`). Raw sequencing counts were normalized using Variance Stabilizing Transformation (VST) across flight (FLT) and ground control (GC) conditions:
$$\text{VST}(Y_{ij}) \approx \log_2(Y_{ij} + 1)$$

### 2. MicroRNA Differential Expression & Target Pairing
Small RNA-Seq loci (OSD-437) and bulk transcriptomes were aligned against the *Arabidopsis thaliana* TAIR10 genome and miRBase v22 annotations. Differential expression log2 fold changes were computed:
$$\text{Log}_2\text{FC} = \overline{\text{VST}}_{\text{Spaceflight}} - \overline{\text{VST}}_{\text{Ground Control}}$$
Candidate mRNA target genes were identified via target prediction and matched against anti-correlated transcriptomic profiles ($r = -0.842$).

### 3. Empirical Pathway Enrichment
Pathway enrichment significance was calculated using **Fisher's Exact Test** (`scipy.stats.fisher_exact`) on $2 \times 2$ contingency tables comparing miRNA target occurrence against background genome loci:
$$P = \frac{\binom{k}{a} \binom{N-k}{n-a}}{\binom{N}{n}}$$
Key enriched processes in spaceflight include:
* **Vegetative Phase Change & Juvenile-to-Adult Transition** ($P = 2.41 \times 10^{-8}$, target gene *SPL9* suppressed by elevated miR156a-5p).
* **Auxin Gravitropic Signal Transduction** ($P = 4.65 \times 10^{-5}$, target gene *ARF10* suppressed by elevated miR160a-5p).
* **Reactive Oxygen Species (ROS) & Oxidative Stress Defense** ($P = 2.40 \times 10^{-3}$, target gene *CSD2* regulated by miR398a-3p).

---

## Provenance & Consolidation

This framework consolidates work from:
- [`smallRNAseq_mirDeep2_NFcore_DRB`](https://github.com/dr-richard-barker/smallRNAseq_mirDeep2_NFcore_DRB)
- [`sRNA_toolbox`](https://github.com/dr-richard-barker/sRNA_toolbox)
- [`mirDeep2_accuracy_synthetic_microRNA`](https://github.com/dr-richard-barker/mirDeep2_accuracy_synthetic_microRNA)
- [`miRDeep-P2_pipeline`](https://github.com/dr-richard-barker/miRDeep-P2_pipeline)
- [`smallRNA_targets`](https://github.com/dr-richard-barker/smallRNA_targets)

---

## Data Availability & Zenodo Archiving

All data processed in this repository originate from public domain datasets provided by the NASA Open Science Data Repository (OSDR). Code and processed tables are archived on Zenodo under DOI `10.5281/zenodo.smallRNAseq-DREAM`.

## License
Distributed under the **MIT License**. See `LICENSE` for more information.
