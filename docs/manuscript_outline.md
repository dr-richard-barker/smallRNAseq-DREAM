# Manuscript outline

*Working outline for the smallRNAseq-DREAM methods/benchmarking paper. Grounded in the
material already in this repository. Everything under "Results" that depends on numbers is
marked **[NEEDS BENCHMARK]** — do not write results until `benchmarking/` has been run for
real. No numbers are invented here.*

---

## Working titles

1. *smallRNAseq-DREAM: a FAIR, cross-species small RNA-seq analysis framework benchmarked
   across four miRNA pipelines for spaceflight omics*
2. *Choosing a small RNA-seq pipeline for cross-kingdom spaceflight biology: a benchmark of
   nf-core/smrnaseq, sRNAtoolbox, miRDeep2 and miRDeep-P2*
3. *A reproducible small RNA-seq workflow for NASA OSDR reanalysis: benchmarking known-miRNA
   quantification and novel-miRNA discovery across plants and animals*

## Target journals (suggested, fit noted)

- **GigaScience / GigaByte** — tools + FAIR data + reproducibility; strong fit for a
  benchmarked, containerised, openly-deposited framework.
- **NAR Genomics and Bioinformatics** — methods/benchmark audience.
- **BMC Bioinformatics** — pipeline + benchmark, established format.
- **npj Microgravity** or **Life (MDPI)** — if the spaceflight application is the headline
  rather than the method.
- **Bioinformatics (Application Note)** — if kept short (2 pages + the tool).

---

## Abstract (structured, ~250 words)

- **Background** — miRNAs regulate stress responses; spaceflight small RNA-seq data in NASA
  OSDR span plants and animals; no single agreed pipeline for cross-species reanalysis.
- **Results** — we deconstructed nf-core/smrnaseq and benchmarked four engines against a
  synthetic ground-truth dataset for known-miRNA quantification and novel-miRNA discovery,
  split by kingdom. **[NEEDS BENCHMARK: headline recall/precision/correlation numbers]**
  We recommend the route SRA → mirnaQC → sRNAbench → sRNAde → mirNOVO for cross-species work.
- **Conclusions** — a FAIR, containerised framework (code + synthetic benchmark + data DOI)
  that makes spaceflight small RNA-seq reanalysis reproducible.
- **Availability** — GitHub + Zenodo DOI.

## Keywords

microRNA; small RNA-seq; spaceflight; NASA OSDR/GeneLab; benchmarking; FAIR; miRDeep2;
sRNAtoolbox; cross-species; isomiR.

---

## 1. Introduction

1.1 miRNA biology and why small RNA-seq is hard (short reads, adapter handling, isomiRs,
    novel-miRNA calling, plant vs animal biogenesis differences).
1.2 The spaceflight context — NASA OSDR/GeneLab, cross-species and cross-kingdom datasets,
    the reproducibility/FAIR imperative. *(Cite the author's related OSDR/GeneLab work.)*
1.3 The problem: many tools, no consensus route; plant and animal predictors differ
    (miRDeep2 vs miRDeep-P2); browser tools vs command-line/high-throughput.
1.4 Contribution: (i) a deconstruction of nf-core/smrnaseq, (ii) a four-engine benchmark on
    synthetic ground truth, (iii) a recommended cross-species route, (iv) a FAIR,
    containerised, openly-deposited framework.

## 2. Materials and Methods

2.1 **Pipelines compared** — Table 1 (tool/version/role/kingdom/citation; already drafted in
    [`methods_tables.md`](methods_tables.md)):
    - nf-core/smrnaseq (Plans A/B; the deconstruction) — [`../pipelines/nfcore_smrnaseq/`](../pipelines/nfcore_smrnaseq/)
    - sRNAtoolbox: sRNAbench / sRNAde / mirNOVO (+ mirnaQC) — [`../pipelines/srnatoolbox/`](../pipelines/srnatoolbox/)
    - miRDeep2 (animal) — [`../pipelines/mirdeep2_animal/`](../pipelines/mirdeep2_animal/)
    - miRDeep-P2 / miRDP2 (plant; fork of TF-Chan-Lab) — [`../pipelines/mirdp2_plant/`](../pipelines/mirdp2_plant/)
2.2 **Synthetic benchmark design** — reads simulated from miRBase miRNAs at defined
    abundances; a held-out subset as ground-truth novels; isomiR/error/adapter/decoy model;
    mature vs genome-anchored modes. Reference [`../benchmarking/`](../benchmarking/) and the
    seeded `config.yaml` for exact parameters.
2.3 **Real datasets** — SRR950892–95 (and any OSDR small RNA-seq accessions used); GeneLab
    API / sra-tools acquisition; adapter and QC handling.
2.4 **Scoring** — known-miRNA recall/precision/F1; quantification correlation (Spearman/
    Pearson, log); novel precision/recall (sequence match ≤2 edits); cross-engine concordance
    (Jaccard). Implemented in `score.py`.
2.5 **Target prediction** — sRNAtoolbox miRNAconsTargets (animal & plant) and the BLAST
    approach in [`../target_prediction/`](../target_prediction/).
2.6 **Deployment / reproducibility** — Singularity/Docker, Nextflow, SLURM + OSDR
    integration; conda `environment.yml`; FAIR data deposition.

## 3. Results

3.1 **Deconstructing nf-core/smrnaseq** — the step map (Plan A → tidied Plan B); Fig 1.
3.2 **Known-miRNA quantification accuracy** — per-engine recall/precision and abundance
    correlation on synthetic data. **[NEEDS BENCHMARK]** Fig 2, Table 3.
3.3 **Novel-miRNA discovery** — precision/recall on held-out miRNAs; where each engine
    over-/under-calls. **[NEEDS BENCHMARK]** (use `genome` mode for the fair test.)
3.4 **Cross-kingdom comparison** — animal (miRDeep2) vs plant (miRDeep-P2); does the sRNAtoolbox
    route hold across kingdoms without retuning? **[NEEDS BENCHMARK]** Fig 3.
3.5 **Concordance on real spaceflight data** — engine agreement on SRR/OSDR runs (truth-free).
    **[NEEDS BENCHMARK]** Fig 4.
3.6 **Recommended route** — evidence for SRA → mirnaQC → sRNAbench → sRNAde → mirNOVO as the
    cross-species default; throughput/deployment notes.

## 4. Discussion

- Why the recommended route wins for cross-species (accuracy vs generality vs throughput).
- Plant vs animal caveats (biogenesis, hairpin criteria, miRDP2 variant rule).
- Limitations: synthetic-data realism, miRBase completeness, novel-call validation needs
  wet-lab confirmation (qPCR), DB-version sensitivity of parsers.
- Fit into the broader OSDR/GeneLab reanalysis and cross-species meta-analysis effort.

## 5. Conclusions

A reproducible, FAIR, cross-species small RNA-seq framework with an evidence-based pipeline
recommendation and an openly-deposited synthetic benchmark others can reuse.

---

## Figures & Tables

| # | Item | Status | Source |
|---|---|---|---|
| Fig 1 | nf-core/smrnaseq deconstruction (Plan A/B schematic) | **have** | migrate PNGs to `docs/schematics/` |
| Fig 2 | Known-miRNA recall/precision + quantification scatter | NEEDS BENCHMARK | `benchmarking/results/` |
| Fig 3 | Cross-kingdom (animal vs plant) comparison | NEEDS BENCHMARK | `benchmarking/results/` |
| Fig 4 | Cross-engine concordance on real OSDR/SRR data | NEEDS BENCHMARK | `benchmarking/results/concordance_heatmap.png` |
| Fig 5 | (optional) miRNA→target network / example | optional | `target_prediction/` |
| Table 1 | Tools compared (version/role/kingdom/citation) | **have (draft)** | `docs/methods_tables.md` |
| Table 2 | Datasets (synthetic params + SRR/OSDR accessions) | partial | `benchmarking/config.yaml` |
| Table 3 | Accuracy metrics per engine | NEEDS BENCHMARK | `benchmarking/results/metrics_summary.tsv` |

## Data & code availability

- Code: `https://github.com/dr-richard-barker/smallRNAseq-DREAM` (**TODO: push**).
- Synthetic benchmark data + truth tables: **Zenodo DOI TODO** (see [`../data/MANIFEST.md`](../data/MANIFEST.md)).
- Real data: NCBI-SRA (SRR950892–95) and NASA OSDR accessions **[list]**.

## Author contributions / funding / competing interests / acknowledgements

- **TODO** — author list, ORCID, funding (NASA/GeneLab program acknowledgement?), competing
  interests. Retain attribution for TF-Chan-Lab (miRDeep-P2) and upstream tool authors.

---

## What's still needed before submission (checklist)

- [ ] Run `benchmarking/` for real on all four engines → fills every **[NEEDS BENCHMARK]**.
- [ ] Verify the `CONFIRM:` parser lines against one real run of each engine.
- [ ] Migrate the Plan A/B schematic PNGs into `docs/schematics/` (Fig 1).
- [ ] Finalise Table 1 citations (full author lists, years, DOIs).
- [ ] List the OSDR spaceflight accessions used (Table 2).
- [ ] Deposit synthetic data to Zenodo; insert the DOI.
- [ ] Decide target journal → match length/format (full paper vs Application Note).
