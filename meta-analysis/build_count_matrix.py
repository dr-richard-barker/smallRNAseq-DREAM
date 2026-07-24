#!/usr/bin/env python3
"""Build a combined miRNA count matrix across OSDR studies (meta-analysis scaffold).

Each study is expected to have produced a `counts.tsv` in the common pipeline format:

    mirna_id <tab> sample <tab> count

(that is what pipelines/* and osdr/run_demo.sh emit). This script merges them into:

    combined_long.tsv     mirna_id, sample, count, cpm, study, organism, kingdom, layer, factor
    combined_wide.tsv     miRNA x sample matrix of raw counts
    combined_wide_cpm.tsv miRNA x sample matrix of CPM (within-sample normalised)
    mirna_prevalence.tsv  per-miRNA: n_samples, n_studies, mean_cpm, kingdoms, layers

Layer honesty
-------------
`layer` (from studies.tsv) records what the counts actually measure:
  mature     = mature miRNA reads (dedicated small RNA-seq; e.g. the mouse/human studies)
  precursor  = MIR-locus / pri-miRNA reads recovered from standard RNA-seq (the plant route;
               see ../smallrna_from_rnaseq). These are NOT the same measurement unit as
               mature counts and must not be pooled naively. This script keeps `layer` on
               every row and in the prevalence summary so downstream analysis can stratify.

Normalisation
-------------
CPM is computed *within sample* (count / sample_total * 1e6) only to make samples roughly
comparable. It does NOT correct cross-study batch effects, differing library preps, or the
mature-vs-precursor layer difference. Proper meta-analysis needs batch correction
(e.g. ComBat-seq) and layer stratification — flagged, not done here.

Usage
-----
    python build_count_matrix.py --runs runs/ --studies studies.tsv --out combined/
        # expects runs/OSD-334/counts.tsv, runs/OSD-483/counts.tsv, ...
"""

import argparse
import os
import sys
from glob import glob

import pandas as pd


import re


def norm_id(x):
    return str(x).strip().lower().lstrip(">")


def mirna_family(mirna_id):
    """Strip the species prefix so homologous miRNAs align across species.
    e.g. 'mmu-miR-21-5p' and 'hsa-miR-21-5p' -> 'mir-21-5p'; 'ath-MIR159a' -> 'mir159a'.
    This is a heuristic for cross-species pooling — verify against miRBase families for
    anything load-bearing (seed-based families can group differently)."""
    s = norm_id(mirna_id)
    s = re.sub(r"^[a-z]{3,4}-", "", s)      # drop 3-4 letter organism code + dash
    s = re.sub(r"^micrornas?", "mir", s)    # 'microRNA159' -> 'mir159'
    s = re.sub(r"^mir[-_]?", "mir-", s)     # normalise 'mir21'/'miR-21'/'MIR21' -> 'mir-21'
    return s


def load_studies(path):
    if not path or not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep="\t")
    return {r["osd_id"]: r.to_dict() for _, r in df.iterrows()}


def find_count_files(runs_dir):
    """Return {osd_id: path} for every runs/<OSD-id>/counts.tsv found."""
    out = {}
    for p in sorted(glob(os.path.join(runs_dir, "*", "counts.tsv"))):
        osd = os.path.basename(os.path.dirname(p))
        out[osd] = p
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default="runs", help="dir containing <OSD-id>/counts.tsv")
    ap.add_argument("--studies", default="studies.tsv", help="study metadata table")
    ap.add_argument("--out", default="combined", help="output dir")
    args = ap.parse_args()

    meta = load_studies(args.studies)
    files = find_count_files(args.runs)
    if not files:
        sys.exit(f"No <OSD-id>/counts.tsv found under {args.runs}/. "
                 f"Run the pipeline per study first (see osdr/run_demo.sh).")
    os.makedirs(args.out, exist_ok=True)

    frames = []
    for osd, path in files.items():
        df = pd.read_csv(path, sep="\t")
        if not {"mirna_id", "sample", "count"} <= set(df.columns):
            print(f"[skip] {osd}: {path} missing required columns", file=sys.stderr)
            continue
        df["mirna_id"] = df["mirna_id"].map(norm_id)
        df["family"] = df["mirna_id"].map(mirna_family)
        # prefix sample with study so samples stay unique across studies
        df["sample"] = osd + ":" + df["sample"].astype(str)
        # within-sample CPM
        totals = df.groupby("sample")["count"].transform("sum").replace(0, pd.NA)
        df["cpm"] = df["count"] / totals * 1e6
        m = meta.get(osd, {})
        for col in ("organism", "kingdom", "layer", "factor"):
            df[col] = m.get(col, "unknown")
        df["study"] = osd
        frames.append(df)
        print(f"[ok] {osd}: {df['mirna_id'].nunique()} miRNAs x "
              f"{df['sample'].nunique()} samples  (layer={m.get('layer','?')})")

    long = pd.concat(frames, ignore_index=True)
    long.to_csv(os.path.join(args.out, "combined_long.tsv"), sep="\t", index=False)

    wide = long.pivot_table(index="mirna_id", columns="sample", values="count",
                            aggfunc="sum", fill_value=0)
    wide.to_csv(os.path.join(args.out, "combined_wide.tsv"), sep="\t")
    wide_cpm = long.pivot_table(index="mirna_id", columns="sample", values="cpm",
                                aggfunc="sum", fill_value=0)
    wide_cpm.to_csv(os.path.join(args.out, "combined_wide_cpm.tsv"), sep="\t")

    # prevalence: how widely is each miRNA seen, and in which kingdoms/layers
    present = long[long["count"] > 0]
    prev = present.groupby("mirna_id").agg(
        n_samples=("sample", "nunique"),
        n_studies=("study", "nunique"),
        mean_cpm=("cpm", "mean"),
        kingdoms=("kingdom", lambda s: ",".join(sorted(set(s)))),
        layers=("layer", lambda s: ",".join(sorted(set(s)))),
    ).sort_values(["n_studies", "n_samples"], ascending=False)
    prev.to_csv(os.path.join(args.out, "mirna_prevalence.tsv"), sep="\t")

    # family-level prevalence — collapses species prefixes so homologues pool across studies
    fam = present.groupby("family").agg(
        n_samples=("sample", "nunique"),
        n_studies=("study", "nunique"),
        n_ids=("mirna_id", "nunique"),
        mean_cpm=("cpm", "mean"),
        kingdoms=("kingdom", lambda s: ",".join(sorted(set(s)))),
        layers=("layer", lambda s: ",".join(sorted(set(s)))),
    ).sort_values(["n_studies", "n_samples"], ascending=False)
    fam.to_csv(os.path.join(args.out, "family_prevalence.tsv"), sep="\t")

    n_studies = long["study"].nunique()
    print(f"\nCombined: {wide.shape[0]} miRNAs x {wide.shape[1]} samples "
          f"across {n_studies} studies -> {args.out}/")
    core = prev[prev["n_studies"] == n_studies]
    fam_core = fam[fam["n_studies"] == n_studies]
    print(f"miRNA IDs detected in ALL {n_studies} studies: {len(core)}  "
          f"| miRNA FAMILIES in all {n_studies}: {len(fam_core)} (cross-species)")
    if (long["layer"].nunique() > 1):
        print("NOTE: mixed layers present "
              f"({', '.join(sorted(long['layer'].unique()))}) — stratify before pooling "
              "(mature vs precursor are different measurement units).")


if __name__ == "__main__":
    main()
