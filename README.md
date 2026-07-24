# smallRNAseq-DREAM

**A FAIR, cross-species small RNA-seq (miRNA) analysis framework for spaceflight omics.**

🌐 **Project site:** https://dr-richard-barker.github.io/smallRNAseq-DREAM/

This repository consolidates work previously split across five repositories into a
single pipeline framework. It benchmarks four small RNA-seq analysis engines against
one another and recommends a single cross-species route for reanalysing NASA OSDR /
GeneLab data.

> **Status:** consolidation in progress. Content migrated from the repositories listed
> under [Provenance](#provenance). Sections marked `TODO` are placeholders to be filled
> by the author — no results have been invented.

## The decision, in one line

After deconstructing the `nf-core/smrnaseq` pipeline and testing the alternatives, the
route found optimal for **cross-species** analysis is:

```
SRA  →  mirnaQC  →  sRNAbench  →  sRNAde  →  mirNOVO
```

The full reasoning (Plans A/B/C/D) is in [`docs/00_pipeline_decision.md`](docs/00_pipeline_decision.md).

## Why four engines?

The same workflow — QC → align/annotate → count → differential expression → target
prediction — was implemented with four different engines, split by kingdom:

| Engine | Kingdom | Where |
|---|---|---|
| nf-core/smrnaseq | Both | [`pipelines/nfcore_smrnaseq/`](pipelines/nfcore_smrnaseq/) |
| sRNAtoolbox (sRNAbench / sRNAde / mirNOVO) | Both | [`pipelines/srnatoolbox/`](pipelines/srnatoolbox/) |
| miRDeep2 | Animal | [`pipelines/mirdeep2_animal/`](pipelines/mirdeep2_animal/) |
| miRDeep-P2 (miRDP2) | Plant | [`pipelines/mirdp2_plant/`](pipelines/mirdp2_plant/) |

Downstream, [`target_prediction/`](target_prediction/) covers miRNA target search, and
[`benchmarking/`](benchmarking/) holds the synthetic-data accuracy harness.

## Repository layout

```
smallRNAseq-DREAM/
├── docs/                 # the pipeline decision, methods tables, engine comparison
├── pipelines/            # one directory per analysis engine
├── target_prediction/    # miRNA target search (BLAST)
├── benchmarking/         # synthetic / SRR accuracy tests (code only; data on Zenodo)
├── third_party/          # notes on vendored / forked upstream tools
└── data/                 # data manifest — large files live on Zenodo / Releases
```

## Data

Large reference and test files (genomes, miRBase dumps, synthetic FASTQ) are **not**
stored in git. See [`data/MANIFEST.md`](data/MANIFEST.md) for what to download and from
where.

## Provenance

Consolidated from:

- [`smallRNAseq_mirDeep2_NFcore_DRB`](https://github.com/dr-richard-barker/smallRNAseq_mirDeep2_NFcore_DRB)
- [`sRNA_toolbox`](https://github.com/dr-richard-barker/sRNA_toolbox)
- [`mirDeep2_accuracy_synthetic_microRNA`](https://github.com/dr-richard-barker/mirDeep2_accuracy_synthetic_microRNA)
- [`miRDeep-P2_pipeline`](https://github.com/dr-richard-barker/miRDeep-P2_pipeline) (forked from TF-Chan-Lab)
- [`smallRNA_targets`](https://github.com/dr-richard-barker/smallRNA_targets)

## Licensing

See [`NOTICE.md`](NOTICE.md). This framework wraps several third-party tools that carry
their own licenses (retained in `third_party/` and the relevant pipeline directories).

## Author

Richard John Barker. Contributions and co-development welcome.
