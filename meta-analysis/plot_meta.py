#!/usr/bin/env python3
"""Manuscript figures from the combined meta-analysis matrix.

Reads the outputs of build_count_matrix.py and writes:
  heatmap_top_mirnas.png     top-variance miRNAs x samples (log CPM), samples colour-barred
                             by kingdom+layer so mature (animal) vs precursor (plant) is visible
  family_prevalence.png      miRNA families by number of studies, cross-species families
                             (seen in >1 study) highlighted

Usage:
  python plot_meta.py --combined combined --out combined/figures --top 30
"""

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

KING_COLOR = {"animal": "#4f46e5", "plant": "#047857", "microbe": "#b45309", "unknown": "#94a3b8"}


def sample_annot(combined):
    """sample -> (kingdom, layer, study) from combined_long.tsv."""
    long = pd.read_csv(os.path.join(combined, "combined_long.tsv"), sep="\t")
    a = long.drop_duplicates("sample").set_index("sample")
    return a[["kingdom", "layer", "study"]]


def heatmap(combined, out, top):
    cpm = pd.read_csv(os.path.join(combined, "combined_wide_cpm.tsv"), sep="\t", index_col=0)
    if cpm.empty:
        print("  [heatmap] no data"); return
    annot = sample_annot(combined).reindex(cpm.columns)
    L = np.log1p(cpm)
    # top miRNAs by variance across samples
    keep = L.var(axis=1).sort_values(ascending=False).head(top).index
    M = L.loc[keep]

    fig_h = max(4, 0.28 * len(keep) + 1.5)
    fig_w = max(6, 0.42 * M.shape[1] + 3)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(M.values, aspect="auto", cmap="magma")
    ax.set_yticks(range(len(keep))); ax.set_yticklabels(keep, fontsize=7)
    ax.set_xticks(range(M.shape[1]))
    ax.set_xticklabels([c.split(":")[-1] for c in M.columns], rotation=90, fontsize=6)
    # kingdom/layer colour bar above the columns
    for j, s in enumerate(M.columns):
        k = annot.loc[s, "kingdom"] if s in annot.index else "unknown"
        ax.add_patch(plt.Rectangle((j - 0.5, -1.4), 1, 0.9, color=KING_COLOR.get(k, "#94a3b8"),
                                   clip_on=False))
    ax.set_ylim(len(keep) - 0.5, -1.6)
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02); cbar.set_label("log(1+CPM)", fontsize=8)
    handles = [plt.Line2D([0], [0], marker="s", ls="", color=c, label=k)
               for k, c in KING_COLOR.items() if k in set(annot["kingdom"])]
    ax.legend(handles=handles, title="kingdom", loc="upper left", bbox_to_anchor=(1.12, 1),
              fontsize=7, title_fontsize=8, frameon=False)
    ax.set_title(f"Top {len(keep)} variable miRNAs across OSDR studies", fontsize=11)
    fig.tight_layout()
    p = os.path.join(out, "heatmap_top_mirnas.png"); fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig); print(f"  wrote {p}")


def family_bar(combined, out, top):
    fp = os.path.join(combined, "family_prevalence.tsv")
    if not os.path.exists(fp):
        print("  [family] family_prevalence.tsv missing"); return
    fam = pd.read_csv(fp, sep="\t")
    fam = fam.sort_values(["n_studies", "mean_cpm"], ascending=False).head(top)
    def col(row):
        ks = str(row["kingdoms"])
        if "," in ks: return "#e11d48"                 # cross-KINGDOM family (rare, notable)
        return KING_COLOR.get(ks, "#94a3b8")
    colors = [col(r) for _, r in fam.iterrows()]
    fig, ax = plt.subplots(figsize=(7, max(3, 0.3 * len(fam) + 1)))
    y = range(len(fam))
    ax.barh(list(y), fam["n_studies"], color=colors)
    ax.set_yticks(list(y)); ax.set_yticklabels(fam["family"], fontsize=8)
    ax.invert_yaxis(); ax.set_xlabel("number of studies detecting the family")
    ax.set_title(f"Cross-study miRNA family prevalence (top {len(fam)})", fontsize=11)
    for i, (_, r) in enumerate(fam.iterrows()):
        ax.text(r["n_studies"] + 0.02, i, f"{int(r['n_ids'])} id(s)", va="center", fontsize=6, color="#475569")
    handles = [plt.Line2D([0], [0], marker="s", ls="", color="#4f46e5", label="animal"),
               plt.Line2D([0], [0], marker="s", ls="", color="#047857", label="plant"),
               plt.Line2D([0], [0], marker="s", ls="", color="#e11d48", label="cross-kingdom")]
    ax.legend(handles=handles, fontsize=7, frameon=False, loc="lower right")
    fig.tight_layout()
    p = os.path.join(out, "family_prevalence.png"); fig.savefig(p, dpi=150)
    plt.close(fig); print(f"  wrote {p}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--combined", default="combined")
    ap.add_argument("--out", default="combined/figures")
    ap.add_argument("--top", type=int, default=30)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    heatmap(args.combined, args.out, args.top)
    family_bar(args.combined, args.out, args.top)


if __name__ == "__main__":
    main()
