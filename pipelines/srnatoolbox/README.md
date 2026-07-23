# sRNAtoolbox pipeline (the recommended cross-species route)

*Source: migrated from `sRNA_toolbox` (`README (1).md` + `sRNA_toolbox_Code.md`).*

This module implements the **SRA → mirnaQC → sRNAbench → sRNAde → mirNOVO** route judged
optimal for cross-species analysis. The goal is to run the sRNAtoolbox from the command
line inside a Singularity instance for high throughput (the browser GUI was shown to be
accurate for known small RNAs and capable of predicting novel miRNAs and their targets).

## Software

| Program | Version | Link |
|---|---|---|
| Bowtie | 1.3.1 | https://bowtie-bio.sourceforge.net/tutorial.shtml |
| Vienna RNA | 2.6.4 | https://www.tbi.univie.ac.at/RNA/ |
| sRNAtoolboxDB | — | https://bioinfo2.ugr.es/srnatoolbox/standalone/ |

## Processing overview

### 1. Set up the sRNA database
```bash
# 1a. Build a bowtie index of the genome and add it to the sRNAtoolbox index folder
bowtie-build genome.fa genome
mv genome.*.ebwt /opt/sRNAtoolboxDB/index/

# 1b. Build the prepared genome sequence object and add to seqOBJ folder
java -jar makeSeqObj.jar genome.fa
mv genome.zip /opt/sRNAtoolboxDB/seqOBJ/genome.zip
```
Genome preprocessing may need the NCBI parser: https://arn.ugr.es/srnatoolbox/helper/ncbiparser/

### 2. Preprocessing / adapter removal
```bash
java -jar /opt/sRNAtoolboxDB/exec/sRNAbench.jar \
  input=/path/to/*.fq output=/opt/sRNAtoolboxDB/out/pre/ \
  adapterMinLength=6 adapter=TCGTATGCCG
```

### 3. microRNA profiling (genome mapping mode)
```bash
java -jar /opt/sRNAtoolboxDB/exec/sRNAbench.jar \
  input=/opt/sRNAtoolboxDB/out/pre/reads_orig.fa \
  output=/opt/sRNAtoolboxDB/out/miR microRNA=hsa species=genome
```

### 4. Novel microRNA prediction
```bash
java -jar /opt/sRNAtoolboxDB/exec/sRNAbench.jar \
  input=/opt/sRNAtoolboxDB/out/pre/reads_orig.fa \
  output=/opt/sRNAtoolboxDB/out/prediction \
  microRNA=hsa species=genome predict=true minReadLength=19 maxReadLength=25
# For plants add: kingdom=plants
```

### 5. isomiR / isoRNA detection
```bash
java -jar /opt/sRNAtoolboxDB/exec/sRNAbench.jar \
  input=/opt/sRNAtoolboxDB/out/pre/reads_orig.fa \
  output=/opt/sRNAtoolboxDB/out/libs microRNA=hsa isoMiR=true
```

### 6. Differential expression (sRNAde)
```bash
java -jar /opt/sRNAtoolboxDB/exec/sRNAde.jar \
  input=/opt/sRNAtoolboxDB/out/ output=/path/to/DE/out/ \
  grpString=<f1_1:f2_1#f1_2:f2_2> diffExpr=true
# Outputs edgeR / DESeq / NOISeq matrices, TMM-normalised tables and heatmaps.
```

### 7. Consensus miRNA target prediction
```bash
# From https://github.com/bioinfoUGR/sRNAtoolbox/tree/master/exec
java -jar miRNAconsTargets.jar          # animals
java -jar miRNAconsTargets_plants.py    # plants
```

> Full per-step input/output field descriptions are preserved in the source document.
> Installation notes (Docker / Python / UGR web tools) for mirnaQC, sRNAbench and mirNOVO
> are in [`../../docs/methods_tables.md`](../../docs/methods_tables.md).

## TODO
- Wrap steps 2–6 in a Snakemake / Nextflow driver for high-throughput batch runs.
- Provide the Singularity definition file used for deployment.
