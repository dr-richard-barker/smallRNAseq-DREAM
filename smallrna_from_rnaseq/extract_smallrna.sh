#!/usr/bin/env bash
# Quantify miRNA signal from STANDARD (not small-RNA) RNA-seq — the plant route.
#
# Why this exists: OSDR has no plant small RNA-seq (see ../osdr/README.md), but it has
# Arabidopsis RNA-seq (OSD-208, OSD-437). Standard RNA-seq reads are fragmented to ~200 bp+
# and sequenced at 50-150 bp, so they CANNOT be mature ~21 nt miRNAs. However, plant
# pri-miRNAs are Pol II transcripts (capped + polyadenylated), so poly-A / total RNA-seq DOES
# carry MIR-locus / precursor signal. This script therefore quantifies PRECURSOR-level miRNA
# gene expression, and runs a read-length diagnostic so you can see empirically whether any
# small-RNA-length reads survived the prep.
#
# Output counts.tsv uses layer = precursor and feeds ../meta-analysis/ (kept distinct from the
# animal mature-miRNA counts — different measurement units; do not pool naively).
#
# Requirements (your compute): python3, seqkit, bowtie2, samtools.
# Usage:
#   OSD=OSD-437 N=2 ./extract_smallrna.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
OSD="${OSD:-OSD-437}"
N="${N:-2}"
OUT="${OUT:-out/$OSD}"
SPECIES="${SPECIES:-ath}"          # miRBase prefix (ath = Arabidopsis thaliana)

need(){ command -v "$1" >/dev/null 2>&1 || { echo "MISSING TOOL: $1"; MISSING=1; }; }
MISSING=0; for t in python3 seqkit bowtie2 samtools; do need "$t"; done
[ "$MISSING" = 1 ] && { echo "Install the toolchain first (see ../environment.yml)."; exit 1; }
mkdir -p "$OUT"; cd "$OUT"

echo ">> 1/5  Fetch $N RNA-seq FASTQ from $OSD"
python3 "$HERE/../osdr/osdr_fetch.py" download "$OSD" --raw --limit "$N" -o reads

echo ">> 2/5  Read-length diagnostic (does any small-RNA-length signal exist?)"
{
  echo -e "sample\tmin_len\tmed_len\tmax_len\treads_le_30nt\treads_total\tpct_le_30nt"
  for fq in reads/*.f*q.gz; do
    [ -e "$fq" ] || continue
    s=$(basename "$fq" | sed 's/\.\(fastq\|fq\)\.gz$//')
    # length distribution + fraction of reads <=30 nt (miRNA-length window)
    read min med max <<<"$(seqkit stats -T "$fq" | awk 'NR==2{print $6, $7, $8}')" || true
    le30=$(seqkit seq -M 30 "$fq" 2>/dev/null | seqkit stats -T 2>/dev/null | awk 'NR==2{print $4}'); le30=${le30:-0}
    tot=$(seqkit stats -T "$fq" | awk 'NR==2{print $4}')
    pct=$(python3 -c "print(f'{100*$le30/max(1,$tot):.3f}')")
    echo -e "$s\t$min\t$med\t$max\t$le30\t$tot\t$pct"
  done
} | tee read_length_diagnostic.tsv
echo "   -> if pct_le_30nt is ~0, mature miRNAs were removed by the prep (expected); precursor signal below is the informative layer."

echo ">> 3/5  Get $SPECIES miRBase precursors (hairpins) + mature, U->T"
if [ ! -s "${SPECIES}_hairpin.fa" ]; then
  curl -sSL https://www.mirbase.org/download/CURRENT/hairpin.fa \
    | seqkit grep -nrp "^${SPECIES}-" | seqkit seq --rna2dna > "${SPECIES}_hairpin.fa"
  curl -sSL https://www.mirbase.org/download/CURRENT/mature.fa \
    | seqkit grep -nrp "^${SPECIES}-" | seqkit seq --rna2dna > "${SPECIES}_mature.fa"
fi
bowtie2-build -q "${SPECIES}_hairpin.fa" "${SPECIES}_hairpin" >/dev/null

echo ">> 4/5  Align RNA-seq reads to precursors (bowtie2 --local) + count per precursor"
: > counts.long.tsv
for fq in reads/*.f*q.gz; do
  [ -e "$fq" ] || continue
  s=$(basename "$fq" | sed 's/\.\(fastq\|fq\)\.gz$//')
  bowtie2 --local --no-unal -p 4 -x "${SPECIES}_hairpin" -U "$fq" 2>"${s}.bt2.log" \
    | samtools sort -o "${s}.bam" - && samtools index "${s}.bam"
  # reads per precursor (mirna_id = precursor)
  samtools idxstats "${s}.bam" | awk -v s="$s" '$1!="*" && $3>0{print $1"\t"s"\t"$3}' >> counts.long.tsv
done

echo ">> 5/5  Write counts.tsv (common format; layer=precursor)"
python3 - <<'PY'
import pandas as pd, os
df = pd.read_csv("counts.long.tsv", sep="\t", names=["mirna_id","sample","count"])
df.to_csv("counts.tsv", sep="\t", index=False)
print(f"  wrote counts.tsv: {df['mirna_id'].nunique()} precursors x {df['sample'].nunique()} samples")
PY

echo "DONE -> $OUT/counts.tsv   (feed to ../meta-analysis with layer=precursor)"
echo "Review read_length_diagnostic.tsv to see whether any mature-length small RNAs survived."
