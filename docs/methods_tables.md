# Methods tables — the DREAM small RNA-seq pipeline (v3)

*Source: migrated from `sRNA_toolbox/drb-dream-machine-learning-smallrna-detection-pipeline.md`.
This is the manuscript-ready methods write-up for the recommended route.*

## Abstract

A workflow for analyzing small RNA sequencing (sRNA-seq) data using readily available
tools to address quality control, annotation, differential expression, and de novo miRNA
discovery. By integrating mirnaQC, sRNAbench, sRNAde and mirNOVO, raw NCBI-SRA data is
transformed into insights about miRNA function and regulation.

## Table 1 — tools used

| Tool | Inputs | Data products | Function | Plants/Animals | Reference |
|---|---|---|---|---|---|
| mirnaQC | FASTQ | BAM, QC reports | Pre-processing & QC for miRNA-seq | Both | Alexandrov et al. (2010) |
| sRNAbench | FASTQ | FASTQ, BAM, count tables | miRNA annotation, DE analysis | Both | Langenberger et al. (2013) |
| sRNAblast | FASTQ | BLAST report, alignments | miRNA target prediction | Both | Wei et al. (2014) |
| mirNOVO | FASTQ | Novel miRNA candidates, FASTA | De novo miRNA discovery | Both | Bonnet et al. (2010) |
| sRNAde | counts.txt, factors.txt | DE results, plots | DE analysis for small RNAs | Both | Anders et al. (2010) |

## Workflow stages

1. **Data acquisition & preprocessing** — identify datasets in NCBI-SRA (GEO accessions
   or keywords); download with sra-tools (`fastq-dump`) or the GeneLab API; run **mirnaQC**
   for adapter trimming, quality analysis and contaminant identification.
2. **Annotation & mapping** — **sRNAbench** annotates known miRNAs against reference
   genomes and miRBase; Bowtie aligns filtered reads; count tables are generated per miRNA.
3. **Differential expression** — count tables into **sRNAde** (DESeq2 / edgeR); visualise
   with heatmaps, volcano plots, expression profiles.
4. **De novo miRNA discovery** — **mirNOVO** on unmapped reads; analyse hairpin secondary
   structure, flanking regions and conservation; validate candidates experimentally
   (e.g. qPCR).

## Tool installation notes

Installation notes for mirnaQC, sRNAbench and mirNOVO (Docker and Python routes, plus the
UGR web tools) are preserved in the source document and in
[`../pipelines/srnatoolbox/README.md`](../pipelines/srnatoolbox/README.md).

> **TODO:** complete the citation list with full author lists, years, journals and DOIs
> before submission. The links in the source document are a starting point.
