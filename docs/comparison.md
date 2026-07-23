# Engine comparison

*This is the results/discussion scaffold for the manuscript. It lists what each engine
does and the axes on which they were (or should be) compared. Accuracy numbers are
`TODO` — to be filled from the benchmarking harness, not invented.*

## Engines compared

| Axis | nf-core/smrnaseq | sRNAtoolbox | miRDeep2 | miRDeep-P2 |
|---|---|---|---|---|
| Kingdom | Both | Both | Animal | Plant |
| Known-miRNA quantification | Bowtie vs miRBase mature/hairpin | sRNAbench | miRDeep2.pl | Bowtie vs predicted |
| Novel prediction | miRDeep2 module | sRNAbench `predict=true` / mirNOVO | miRDeep2.pl | miRDP2 |
| isomiR detection | mirtop | sRNAbench `isoMiR=true` | — | — |
| Differential expression | edgeR | sRNAde (edgeR/DESeq/NOISeq) | external | DESeq2 |
| Target prediction | — | miRNAconsTargets (animal & plant) | — | external BLAST |
| Deployment | Nextflow (Docker/Singularity/SLURM) | JAR + Singularity | conda/mamba | Perl + BioPerl |
| High-throughput ready | Yes (Nextflow) | via Singularity wrapper | needs loop/wrapper | needs wrapper |

## Comparison axes for the paper

- **Accuracy on known miRNAs** — recovery rate against a truth set (synthetic + SRR950892–95).
- **Novel-prediction precision/recall** — on the synthetic set where ground truth is known.
- **Cross-species robustness** — same workflow across kingdoms without per-species retuning.
- **Throughput / reproducibility** — containerisation, cluster deployment, FAIR outputs.

## Headline finding (author's stated conclusion)

The **SRA → mirnaQC → sRNAbench → sRNAde → mirNOVO** route was judged optimal for
cross-species analysis. The kingdom-specific miRDeep2 / miRDeep-P2 branch remains the
comparator for prediction accuracy.

> **TODO:** populate the accuracy table from [`../benchmarking/`](../benchmarking/) results.
