# The pipeline decision: Plans A, B, C and D

*Source: migrated and lightly edited from the `smallRNAseq_mirDeep2_NFcore_DRB` README.
This is the narrative of how the recommended route was chosen. It is intended to become
the backbone of the manuscript's Methods rationale.*

## How we got here

As we explored different small RNA-seq analysis pipelines we started by breaking down
the `nf-core/smrnaseq` pipeline to understand how it works. Writing up these early stages
is good practice while developing bioinformatics pipelines even if they don't make it
into the final code.

Initially the `nf-core/smrnaseq` pipeline was assessed and then broken down into its
components (see `docs/schematics/`). Example commands were created to test its effect on
accuracy using synthetic small RNA-seq data from SRA.

## The four plans

- **Plan A** — the full `nf-core/smrnaseq` "schematic explosion": every step exploded out
  (FastQC → adapter trim → collapse → Bowtie vs mature/hairpin → edgeR → mirtop →
  miRDeep2 → miRTrace → MultiQC). See [`pipelines/nfcore_smrnaseq/`](../pipelines/nfcore_smrnaseq/).
- **Plan B** — the tidied schematic of the same pipeline after deconstruction.
- **Plan C** — the sRNAtoolbox route (sRNAbench / sRNAde / mirNOVO), run from the command
  line inside a Singularity instance for high throughput.
  See [`pipelines/srnatoolbox/`](../pipelines/srnatoolbox/).
- **Plan D (Option D)** — the **miRDeep2-centric** branch, split by kingdom:
  miRDeep2 for animals ([`pipelines/mirdeep2_animal/`](../pipelines/mirdeep2_animal/)) and
  miRDeep-P2 for plants ([`pipelines/mirdp2_plant/`](../pipelines/mirdp2_plant/)).

## The conclusion

> After testing these tools I've come to believe this pipeline
> **SRA → mirnaQC → sRNAbench → sRNAde → mirNOVO** is optimal for cross-species analysis.

Remaining work noted at consolidation time:

- A Snakemake (or equivalent) wrapper to parse data through the pipeline for high
  throughput. The full v3 method is written up in [`methods_tables.md`](methods_tables.md).
- Continued work on the miRDeep2 / miRDeep-P2 plant-specific branch as part of OSDR.
- **TODO (manuscript gap):** quantified accuracy numbers from the synthetic benchmark
  (known-miRNA recovery, novel-prediction precision). See [`../benchmarking/`](../benchmarking/).

## OSDR / SLURM integration

A SLURM + OSDR integration plan exists for running the nf-core route on cluster
infrastructure — migrated to [`pipelines/nfcore_smrnaseq/`](../pipelines/nfcore_smrnaseq/).
