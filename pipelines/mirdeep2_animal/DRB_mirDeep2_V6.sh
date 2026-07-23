#!/usr/bin/env bash
## DRB miRDeep2 pipeline (V6)
## Source: mirDeep2_accuracy_synthetic_microRNA/DRB_mirDeep2_V6.sh
##
## Inputs : genome.fa, mature.fa, hairpin.fa, smallRNAseq.fastq.gz
## Outputs: quantification of known and novel miRNAs / hairpins as text, .bed and .html
##
## One-time setup (uncomment as needed) -----------------------------------------------
## Install miniconda + mirdeep2 via mamba:
#   curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
#   conda install -y -n base -c conda-forge mamba
#   mamba create -y -n mirdeep2 -c conda-forge -c bioconda -c defaults mirdeep2
#   conda activate mirdeep2
##
## Download miRBase references (once):
#   wget -O hairpin.fa https://mirbase.org/download/CURRENT/hairpin.fa
#   wget -O mature.fa  https://mirbase.org/download/CURRENT/mature.fa
##
## Split references by species (human example):
#   grep -A 1 --no-group-separator "^>hsa" mature_ut.fa > mature_hsa_ut.fa
#   grep -v "^>hsa" mature_ut.fa > mature_ut_No-HSA.fa
##
## Prepare genome (once per genome):
#   sed 's/ /_/g' Homo_sapiens.GRCh38.dna.primary_assembly.fna > ..._no_space.fna
#   bowtie-build Homo_sapiens.GRCh38...no_space_pl.fa Homo_sapiens.GRCh38...no_space_pl
## -------------------------------------------------------------------------------------

# Loop over samples (edit the accession list / adapter as needed)
for sample in SRR950892 SRR950893 SRR950894 SRR950895
do
    printf "\n\n    Working on: ${sample}\n\n"

    # Pre-processing (uncomment for a fresh run):
    #   seqkit fq2fa ${sample}_trimmed.fq.gz -o ${sample}.fa
    #   remove_white_space_in_id.pl ${sample}.fa > ${sample}_no_whitespace.fa
    #   sed 's/ /_/g' ${sample}_trimmed.fa > ${sample}_trimmed_no_spaces.fa

    # Map reads to the indexed genome
    mapper.pl ${sample}_trimmed_no_spaces.fa -c -j -k TCGTATGCCGTCTTCTGCTTGT -l 18 -m \
        -p Homo_sapiens.GRCh38.dna.primary_assembly_no_space_pl \
        -s ${sample}_trimmed_collapsed.fa \
        -t ${sample}_trimmed_collapsed_vs_genome_no_space_pl.arf -v -n

    # Quantify known + predict novel miRNAs
    miRDeep2.pl ${sample}_trimmed_collapsed.fa \
        Homo_sapiens.GRCh38.dna.primary_assembly_no_space_pl.fa \
        ${sample}_trimmed_collapsed_vs_genome_no_space_pl.arf \
        mature_ut.part_hsa_no_whitespace.fasta \
        mature_ut_No-HSA_NOWHITESPACE.fa \
        hairpin_ut.part_hsa_no_whitespace.fasta \
        -t Human 2>report.log
done

## Next: extract counts and summary stats from results.html / result tables.
