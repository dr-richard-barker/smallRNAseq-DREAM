#!/usr/bin/env bash
# Templated engine runs for the benchmark. Each block runs one engine on the SAME synthetic
# reads and reference, writing into the run_dir the config expects. Uncomment/adjust the
# blocks for the engines you have installed. Adapter used by make_synthetic.py:
#   TGGAATTCTCGGGTGCCAAGG   (Illumina small RNA 3' adapter)
set -euo pipefail

READS=synthetic/synthetic_reads.fastq.gz
REF=synthetic/reference_mature.fa          # known miRNAs only (novels are held out)
GENOME=refs/genome.fa                      # whitespace-stripped headers
HAIRPIN=refs/hairpin.fa

# --- miRDeep2 (animal) ---------------------------------------------------------------
# mkdir -p runs/mirdeep2 && cd runs/mirdeep2
# seqkit fq2fa ../../$READS | sed 's/ /_/g' > reads.fa
# mapper.pl reads.fa -c -j -k TGGAATTCTCGGGTGCCAAGG -l 18 -m -p genome_idx \
#     -s reads_collapsed.fa -t reads_vs_genome.arf -v -n
# miRDeep2.pl reads_collapsed.fa ../../$GENOME reads_vs_genome.arf \
#     ../../$REF none ../../$HAIRPIN -t Human 2>report.log
# quantifier.pl -p ../../$HAIRPIN -m ../../$REF -r reads_collapsed.fa   # -> miRNAs_expressed_all_samples*.csv
# cd ../..

# --- sRNAtoolbox (sRNAbench + sRNAde) -----------------------------------------------
# mkdir -p runs/srnatoolbox
# java -jar $SRNATOOLBOX/sRNAbench.jar input=$READS output=runs/srnatoolbox/ \
#     microRNA=hsa species=genome predict=true isoMiR=true
#   # -> mature_sense.grouped, novel_mature.fa

# --- miRDeep-P2 (plant) --------------------------------------------------------------
# mkdir -p runs/mirdp2
# miRDP2-v1.1.4_pipeline.bash -g $GENOME -x $GENOME -q -b fastq_list.txt -o runs/mirdp2
# perl parse_miRDP2_prediction.pl miRDP2_predictions_list.txt runs/mirdp2/miRDP2
#   # then Bowtie map + bam2ref_counts.pl + combine_htseq_counts.pl -> runs/mirdp2/count_table.txt

# --- nf-core/smrnaseq ----------------------------------------------------------------
# nextflow run nf-core/smrnaseq -profile singularity \
#     --input samplesheet.csv --genome GRCh38 --outdir runs/nfcore \
#     --three_prime_adapter TGGAATTCTCGGGTGCCAAGG

echo "Edit this file: uncomment the engine block(s) you have installed, then re-run:"
echo "  python run_benchmark.py --config config.yaml --step score"
