# miRDeep-P2 (miRDP2) pipeline (plant)

> **Attribution:** this module is a **fork of** [`TF-Chan-Lab/miRDeep-P2_pipeline`](https://github.com/TF-Chan-Lab/miRDeep-P2_pipeline).
> The upstream LICENSE and the Perl parsing scripts belong to their original authors and
> must be retained. First described in https://doi.org/10.1002/tpg2.20103. miRDeep-P2
> tool: https://academic.oup.com/bioinformatics/article/35/14/2521/5232218

*Source: migrated from `miRDeep-P2_pipeline`.*

A pipeline for plant miRNA prediction and differential expression from sRNA-seq data.
miRDeep-P2 predicts miRNAs; predictions are parsed, annotated against known miRNAs
(BLAST), then quantified (Bowtie) and tested for DE (DESeq2).

## Environment

Linux CLI · Perl + BioPerl (Bio::SeqIO, Bio::Seq, Bio::SearchIO) · miRDeep-P2 · BLAST ·
Bowtie · samtools · R + DESeq2.

## Scripts (retained from upstream, under `script/`)

`fasta_U2T.pl` · `unique_fasta_v1.3.pl` · `parse_miRDP2_prediction.pl` ·
`general_blast_parser.pl` · `parse_parsed_blast_known_plants.pl` ·
`filter_lines_by_key_words_list.pl` · `bam2ref_counts.pl` · `combine_htseq_counts.pl`

> These Perl scripts are **not** copied into this scaffold — pull them from the upstream
> repo (or the existing fork) so provenance and license stay intact:
> `git clone https://github.com/TF-Chan-Lab/miRDeep-P2_pipeline.git`

## Workflow (summary)

1. **Prediction** — `miRDP2-v1.1.4_pipeline.bash -g ref.fa -x ref.fa -q -b fastq_list.txt -o result`,
   then `parse_miRDP2_prediction.pl` → `miRDP2_mature.fa`, `miRDP2_arms.txt`.
2. **Annotation** — `blastn` predicted vs miRBase (U→T converted, collapsed), parse into
   *known* / *variant* (≤2 mismatches, no indels) / *novel*.
3. **Differential expression** — Bowtie map each sample to `miRNA_ref.fa`, `bam2ref_counts.pl`
   → per-sample counts, `combine_htseq_counts.pl` → count table, then DESeq2.

The full command sequence is preserved in the upstream README.
