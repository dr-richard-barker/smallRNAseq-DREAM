#!/usr/bin/env bash
# Fetch + process the studies in studies.tsv, dispatching each to the right pipeline by its
# `route` column, producing runs/<OSD-id>/counts.tsv, then build the combined matrix + plots.
#
#   route=mirdeep2   -> animal small RNA-seq (mature)     via ../osdr/run_demo.sh
#   route=precursor  -> plant standard RNA-seq (precursor) via ../smallrna_from_rnaseq/extract_smallrna.sh
#
# Runs the real toolchain on YOUR compute. Start small with N (samples/study).
# Optional filter:  KINGDOM=plant ./fetch_studies.sh   (or animal); default = all.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
N="${N:-2}"
RUNS="${RUNS:-runs}"
KINGDOM="${KINGDOM:-all}"
TSV="$HERE/studies.tsv"

mkdir -p "$RUNS"
# iterate studies.tsv (skip header), respecting the KINGDOM filter
tail -n +2 "$TSV" | while IFS=$'\t' read -r osd organism kingdom layer route paired factor raw size title; do
    [ -z "$osd" ] && continue
    if [ "$KINGDOM" != "all" ] && [ "$KINGDOM" != "$kingdom" ]; then continue; fi
    echo "=================  $osd  ($kingdom / $route)  ================="
    case "$route" in
        mirdeep2)
            OSD="$osd" N="$N" OUT="$RUNS/$osd" "$HERE/../osdr/run_demo.sh" ;;
        precursor)
            OSD="$osd" N="$N" OUT="$RUNS/$osd" PAIRED="$paired" \
                "$HERE/../smallrna_from_rnaseq/extract_smallrna.sh" ;;
        *) echo "  unknown route '$route' for $osd — skipping" ;;
    esac
done

echo "=================  meta-analysis + plots  ================="
python3 "$HERE/build_count_matrix.py" --runs "$RUNS" --studies "$TSV" --out combined
python3 "$HERE/plot_meta.py" --combined combined --out combined/figures || \
    echo "  (plotting skipped — install matplotlib)"
echo "Combined matrix + figures in ./combined/"
