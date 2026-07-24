#!/usr/bin/env bash
# End-to-end demo: OSDR small RNA-seq -> miRDeep2 -> counts.tsv (common format).
#
# Demonstrates the animal (mature-miRNA) route on real spaceflight data: OSD-483
# (astronaut sEV small RNA-seq). The output counts.tsv feeds ../meta-analysis/.
#
# This runs the real toolchain on YOUR compute — it is not run in the dev environment.
# Requirements (install via ../environment.yml or conda): python3, bowtie, seqkit, mirdeep2
# (provides mapper.pl, miRDeep2.pl, quantifier.pl). References: human miRBase mature+hairpin
# (from the Zenodo bundle, or miRBase directly) and a bowtie index of the human genome.
#
# Usage:
#   ./run_demo.sh                 # defaults: OSD-483, 2 samples, ./demo_out
#   OSD=OSD-334 N=4 ./run_demo.sh # a mouse HZE-radiation miRNA study instead
set -euo pipefail

OSD="${OSD:-OSD-483}"
N="${N:-2}"                       # number of raw FASTQ to pull for the demo
OUT="${OUT:-demo_out/$OSD}"
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"

# reference files — point these at your copies (Zenodo bundle has the human set)
GENOME="${GENOME:-refs/genome_no_space.fa}"          # bowtie-indexed, whitespace-free headers
GENOME_IDX="${GENOME_IDX:-refs/genome_no_space}"     # bowtie index prefix
MATURE="${MATURE:-refs/mature_ut.part_hsa_no_whitespace.fasta}"
MATURE_OTHER="${MATURE_OTHER:-refs/mature_ut_No-HSA_NOWHITESPACE.fa}"
HAIRPIN="${HAIRPIN:-refs/hairpin_ut.part_hsa_no_whitespace.fasta}"
ADAPTER="${ADAPTER:-TGGAATTCTCGGGTGCCAAGG}"           # confirm per study from OSDR metadata

need(){ command -v "$1" >/dev/null 2>&1 || { echo "MISSING TOOL: $1"; MISSING=1; }; }
MISSING=0; for t in python3 bowtie seqkit mapper.pl miRDeep2.pl; do need "$t"; done
[ "$MISSING" = 1 ] && { echo "Install the toolchain first (see environment.yml)."; exit 1; }

mkdir -p "$OUT"; cd "$OUT"

echo ">> 1/4  Fetch $N raw FASTQ from $OSD"
python3 "$HERE/osdr_fetch.py" download "$OSD" --raw --limit "$N" -o reads

echo ">> 2/4  Collapse + map reads with the miRDeep2 mapper"
: > all_counts_input.txt
for fq in reads/*.fastq.gz reads/*.fq.gz; do
    [ -e "$fq" ] || continue
    s=$(basename "$fq" | sed 's/\.\(fastq\|fq\)\.gz$//')
    seqkit fq2fa "$fq" | sed 's/ /_/g' > "$s.fa"
    mapper.pl "$s.fa" -c -j -k "$ADAPTER" -l 18 -m -p "$REPO/$GENOME_IDX" \
        -s "${s}_collapsed.fa" -t "${s}_vs_genome.arf" -v -n
    echo ">> 3/4  Quantify known + predict novel miRNAs for $s"
    miRDeep2.pl "${s}_collapsed.fa" "$REPO/$GENOME" "${s}_vs_genome.arf" \
        "$REPO/$MATURE" "$REPO/$MATURE_OTHER" "$REPO/$HAIRPIN" -t none 2>"${s}_report.log" || true
    quantifier.pl -p "$REPO/$HAIRPIN" -m "$REPO/$MATURE" -r "${s}_collapsed.fa" -y "$s" || true
done

echo ">> 4/4  Parse miRDeep2 output -> counts.tsv (common format)"
PYTHONPATH="$REPO/benchmarking" python3 - <<'PY'
import os
from parse_outputs import parse_mirdeep2
counts, preds = parse_mirdeep2(os.getcwd())
counts.to_csv("counts.tsv", sep="\t", index=False)
preds.to_csv("predictions.tsv", sep="\t", index=False)
print(f"  wrote counts.tsv ({len(counts)} rows) and predictions.tsv ({len(preds)} novel)")
PY

echo "DONE -> $OUT/counts.tsv"
echo "Next: repeat for more studies, then  cd ../meta-analysis && python build_count_matrix.py --runs <dir>"
