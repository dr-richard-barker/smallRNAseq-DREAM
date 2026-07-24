# nf-core/smrnaseq (Plans A & B)

*Source: migrated from `smallRNAseq_mirDeep2_NFcore_DRB` (`NF_core_smallRNA_pipeline`,
`smallRNAseq_nf_cor_slurm_v2_for_OSDR`).*

The starting point: the community `nf-core/smrnaseq` best-practice pipeline, deconstructed
step-by-step (Plan A = full "schematic explosion", Plan B = the tidied schematic). Adapted
from https://hub.docker.com/r/nfcore/smrnaseq/.

## Quick start

```bash
# Requires Nextflow >=20.04.0 and one of Docker / Singularity / Podman / Shifter / Charliecloud
nextflow run nf-core/smrnaseq -profile test,<container>
nextflow run nf-core/smrnaseq -profile <container> \
  --input '*_R{1,2}.fastq.gz' --genome GRCh37
```

## Pipeline steps (as deconstructed)

FastQC → FastP (adapter trim) → Bowtie2 (contaminant filtering) → Bowtie (align vs miRBase
mature + hairpin) → SAMtools (feature counting) → edgeR (normalisation, MDS, heatmap) →
Bowtie (align vs reference genome, QC) → mirtop (miRNA + isomiR annotation) → **miRDeep2**
(known + novel discovery) → miRTrace (QC) → MultiQC (aggregate report).

## Source files (archived offline)

The original supporting files lived in the `smallRNAseq_mirDeep2_NFcore_DRB` repository,
which was **deleted from GitHub** during consolidation before these were copied out:

- `NF_core_smallRNA_pipeline` — the annotated pipeline breakdown + example commands for
  testing accuracy on synthetic SRA data. *(Its step-by-step content is summarised above.)*
- `smallRNAseq_nf_cor_slurm_v2_for_OSDR` — SLURM + OSDR cluster integration plan.
  **Not recovered** — supply from a local backup if available.
- Pipeline schematics — see [`../../docs/schematics/`](../../docs/schematics/) for status.

> **Note:** the step list above preserves the pipeline logic. The full original text of
> `NF_core_smallRNA_pipeline` and the SLURM/OSDR script were not captured before the source
> repo was deleted; restore them here if a copy surfaces.
