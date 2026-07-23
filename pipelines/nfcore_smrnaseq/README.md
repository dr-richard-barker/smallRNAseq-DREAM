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

## Files to migrate here

- `NF_core_smallRNA_pipeline` — the annotated pipeline breakdown + example commands for
  testing accuracy on synthetic SRA data.
- `smallRNAseq_nf_cor_slurm_v2_for_OSDR` — SLURM + OSDR cluster integration plan.
- Schematics (`smallRNAseq_analysis_pipeline_v1_NF_smRNAseq.png`, plan A/B images) →
  [`../../docs/schematics/`](../../docs/schematics/).

> **Note:** these two text files and the PNG are in the source repo but were not copied
> into this scaffold automatically. Drop them in here during migration.
