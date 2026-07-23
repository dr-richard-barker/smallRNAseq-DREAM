#!/usr/bin/env python3
"""Score engine outputs against the synthetic truth set.

Produces the numbers the manuscript needs:

  known-miRNA recovery      recall / precision / F1 of KNOWN miRNA detection
  quantification accuracy   Spearman & Pearson of true vs detected abundance (log)
  novel prediction          precision / recall / F1 vs the held-out truth novels
  cross-engine concordance  pairwise Jaccard of detected miRNA sets (truth-free)

Outputs (into scoring.out_dir):
  metrics_summary.tsv
  quantification_<engine>.png
  recovery_bar.png
  concordance_heatmap.png

Usage: called by run_benchmark.py, or standalone:
  python score.py --config config.yaml
"""

import argparse
import itertools
import os

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from parse_outputs import parse_engine, _norm_id, read_fasta


# --------------------------------------------------------------------------------------
# Sequence matching for novel predictions
# --------------------------------------------------------------------------------------
def bounded_edit_distance(a, b, cap):
    """Edit distance with early exit once it exceeds `cap` (returns cap+1)."""
    if abs(len(a) - len(b)) > cap:
        return cap + 1
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i] + [0] * len(b)
        best = cur[0]
        for j, cb in enumerate(b, 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb))
            best = min(best, cur[j])
        if best > cap:                       # whole row already worse than cap
            return cap + 1
        prev = cur
    return prev[-1]


def seq_matches(pred, truth_seqs, cap):
    """True if `pred` is within `cap` edits of any truth sequence (allowing a mature
    sequence to sit inside a slightly longer prediction and vice-versa)."""
    for t in truth_seqs:
        short, long = sorted((pred, t), key=len)
        if bounded_edit_distance(pred, t, cap) <= cap:
            return True
        # sliding window: short sequence against windows of the long one
        w = len(short)
        for k in range(0, len(long) - w + 1):
            if bounded_edit_distance(short, long[k:k + w], cap) <= cap:
                return True
    return False


# --------------------------------------------------------------------------------------
# Metrics
# --------------------------------------------------------------------------------------
def collapse_counts(counts_df):
    """Sum over samples -> Series indexed by normalised mirna_id."""
    if counts_df.empty:
        return pd.Series(dtype=float)
    c = counts_df.copy()
    c["mirna_id"] = c["mirna_id"].map(_norm_id)
    return c.groupby("mirna_id")["count"].sum()


def known_recovery(detected, truth, min_count):
    known = {_norm_id(i) for i in truth.loc[truth.category == "known", "mirna_id"]}
    novel = {_norm_id(i) for i in truth.loc[truth.category == "novel", "mirna_id"]}
    hits = set(detected[detected >= min_count].index)
    tp = len(hits & known)
    fn = len(known - hits)
    # false positives = detected things that are neither a known truth miRNA nor a held-out
    # novel (i.e. decoy noise or spurious ids)
    fp = len(hits - known - novel)
    recall = tp / (tp + fn) if (tp + fn) else float("nan")
    precision = tp / (tp + fp) if (tp + fp) else float("nan")
    f1 = (2 * precision * recall / (precision + recall)
          if precision and recall and not np.isnan(precision) and not np.isnan(recall) else float("nan"))
    return dict(known_tp=tp, known_fn=fn, known_fp=fp,
                known_recall=recall, known_precision=precision, known_f1=f1)


def quant_accuracy(detected, truth, min_count):
    known = truth[truth.category == "known"].copy()
    known["mirna_id"] = known["mirna_id"].map(_norm_id)
    known = known.set_index("mirna_id")
    rows = []
    for mid, tab in known["true_abundance"].items():
        d = detected.get(mid, 0.0)
        if d >= min_count:
            rows.append((tab, d))
    if len(rows) < 3:
        return dict(quant_spearman=float("nan"), quant_pearson=float("nan"), quant_n=len(rows))
    arr = np.array(rows, dtype=float)
    t, d = np.log1p(arr[:, 0]), np.log1p(arr[:, 1])
    # Spearman via rank-Pearson (avoids a scipy hard dependency)
    def pearson(x, y):
        x, y = x - x.mean(), y - y.mean()
        denom = np.sqrt((x * x).sum() * (y * y).sum())
        return float((x * y).sum() / denom) if denom else float("nan")
    def rank(v):
        order = v.argsort()
        r = np.empty_like(order, dtype=float)
        r[order] = np.arange(len(v))
        return r
    return dict(quant_spearman=pearson(rank(t), rank(d)),
                quant_pearson=pearson(t, d), quant_n=len(rows))


def novel_scores(preds_df, novel_truth, cap):
    truth_seqs = list(novel_truth.values())
    if preds_df.empty:
        return dict(novel_tp=0, novel_fp=0, novel_fn=len(truth_seqs),
                    novel_precision=float("nan"), novel_recall=0.0, novel_f1=float("nan"))
    pred_seqs = [s.upper().replace("U", "T") for s in preds_df["sequence"].dropna()]
    tp_preds = sum(1 for p in pred_seqs if seq_matches(p, truth_seqs, cap))
    fp = len(pred_seqs) - tp_preds
    matched_truth = sum(1 for t in truth_seqs if seq_matches(t, pred_seqs, cap))
    fn = len(truth_seqs) - matched_truth
    precision = tp_preds / len(pred_seqs) if pred_seqs else float("nan")
    recall = matched_truth / len(truth_seqs) if truth_seqs else float("nan")
    f1 = (2 * precision * recall / (precision + recall)
          if precision and recall else float("nan"))
    return dict(novel_tp=tp_preds, novel_fp=fp, novel_fn=fn,
                novel_precision=precision, novel_recall=recall, novel_f1=f1)


# --------------------------------------------------------------------------------------
# Plots
# --------------------------------------------------------------------------------------
def plot_quant(detected, truth, engine, out_dir, min_count):
    known = truth[truth.category == "known"].copy()
    known["mirna_id"] = known["mirna_id"].map(_norm_id)
    xs, ys = [], []
    for _, r in known.iterrows():
        d = detected.get(r["mirna_id"], 0.0)
        if d >= min_count:
            xs.append(r["true_abundance"]); ys.append(d)
    if len(xs) < 3:
        return
    plt.figure(figsize=(5, 5))
    plt.scatter(np.log1p(xs), np.log1p(ys), s=12, alpha=0.6)
    lim = [0, max(np.log1p(xs + ys))]
    plt.plot(lim, lim, "--", color="grey", lw=1)
    plt.xlabel("log(1 + true abundance)")
    plt.ylabel("log(1 + detected count)")
    plt.title(f"Quantification accuracy — {engine}")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"quantification_{engine}.png"), dpi=150)
    plt.close()


def plot_recovery(summary, out_dir):
    metrics = ["known_recall", "known_precision", "novel_recall", "novel_precision"]
    engines = summary["engine"].tolist()
    x = np.arange(len(engines))
    w = 0.2
    plt.figure(figsize=(1.6 * len(engines) + 3, 4.5))
    for i, m in enumerate(metrics):
        plt.bar(x + i * w, summary[m].fillna(0), w, label=m.replace("_", " "))
    plt.xticks(x + 1.5 * w, engines)
    plt.ylim(0, 1)
    plt.ylabel("score")
    plt.title("Detection recall / precision")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "recovery_bar.png"), dpi=150)
    plt.close()


def plot_concordance(detected_sets, out_dir):
    engines = list(detected_sets)
    n = len(engines)
    if n < 2:
        return
    mat = np.ones((n, n))
    for i, j in itertools.product(range(n), range(n)):
        a, b = detected_sets[engines[i]], detected_sets[engines[j]]
        mat[i, j] = len(a & b) / len(a | b) if (a | b) else 0.0
    plt.figure(figsize=(1.1 * n + 2, 1.1 * n + 2))
    plt.imshow(mat, vmin=0, vmax=1, cmap="viridis")
    plt.colorbar(label="Jaccard")
    plt.xticks(range(n), engines, rotation=45, ha="right")
    plt.yticks(range(n), engines)
    for i, j in itertools.product(range(n), range(n)):
        plt.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center",
                 color="white" if mat[i, j] < 0.6 else "black", fontsize=8)
    plt.title("Cross-engine concordance")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "concordance_heatmap.png"), dpi=150)
    plt.close()


# --------------------------------------------------------------------------------------
def main():
    import yaml
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    args = ap.parse_args()
    cfg = yaml.safe_load(open(args.config))
    sc = cfg["scoring"]
    out_dir = sc["out_dir"]
    os.makedirs(out_dir, exist_ok=True)

    truth = pd.read_csv(sc["truth_counts"], sep="\t")
    novel_truth = read_fasta(sc["novel_truth_fa"]) if os.path.exists(sc["novel_truth_fa"]) else {}

    summary_rows = []
    detected_sets = {}
    for eng in cfg["engines"]:
        if not eng.get("enabled", True):
            continue
        name = eng["name"]
        try:
            counts_df, preds_df = parse_engine(eng["parser"], eng["run_dir"])
        except Exception as e:                      # noqa: BLE001 - report and skip a missing engine
            print(f"[skip] {name}: {e}")
            continue

        detected = collapse_counts(counts_df)
        detected_sets[name] = set(detected[detected >= sc["min_count"]].index)

        row = {"engine": name}
        row.update(known_recovery(detected, truth, sc["min_count"]))
        row.update(quant_accuracy(detected, truth, sc["min_count"]))
        row.update(novel_scores(preds_df, novel_truth, sc["novel_max_edit"]))
        summary_rows.append(row)
        plot_quant(detected, truth, name, out_dir, sc["min_count"])
        print(f"[ok] {name}: "
              f"known recall={row['known_recall']:.2f} "
              f"quant rho={row['quant_spearman']:.2f} "
              f"novel P/R={row['novel_precision']:.2f}/{row['novel_recall']:.2f}")

    if not summary_rows:
        print("No engines scored. Check run_dir paths in the config.")
        return

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(os.path.join(out_dir, "metrics_summary.tsv"), sep="\t", index=False)
    plot_recovery(summary, out_dir)
    plot_concordance(detected_sets, out_dir)
    print(f"\nWrote {out_dir}/metrics_summary.tsv and plots.")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
