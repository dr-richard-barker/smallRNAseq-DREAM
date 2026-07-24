#!/usr/bin/env bash
# Fetch + process the animal small-RNA studies listed in studies.tsv, producing one
# runs/<OSD-id>/counts.tsv per study, then build the combined matrix.
#
# Runs the real toolchain on YOUR compute (see ../osdr/run_demo.sh requirements).
# For a quick trial, set N low (few samples per study) and start with one study.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
N="${N:-2}"                       # samples per study for a trial run
RUNS="${RUNS:-runs}"

# animal (mature) studies verified to hold real small RNA-seq raw reads:
STUDIES=(OSD-334 OSD-335 OSD-336 OSD-337 OSD-483)

mkdir -p "$RUNS"
for osd in "${STUDIES[@]}"; do
    echo "=================  $osd  ================="
    OSD="$osd" N="$N" OUT="$RUNS/$osd" "$HERE/../osdr/run_demo.sh"
done

echo "=================  meta-analysis  ================="
python3 "$HERE/build_count_matrix.py" --runs "$RUNS" --studies "$HERE/studies.tsv" --out combined
echo "Combined matrix in ./combined/ — see family_prevalence.tsv for cross-species hits."
