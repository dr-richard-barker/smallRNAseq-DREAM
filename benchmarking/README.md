# Benchmarking — synthetic / SRR accuracy harness

The accuracy work lives here: it produces the numbers the manuscript needs — known-miRNA
recovery, quantification accuracy, and novel-prediction precision/recall — by comparing
each engine's output against a **synthetic dataset with known ground truth**. Code only;
test data lives on Zenodo / Releases (see [`../data/MANIFEST.md`](../data/MANIFEST.md)).

## Why synthetic data

Real SRA runs (SRR950892–95) have no absolute truth, so they can only give cross-engine
*concordance*. To get precision/recall you need a dataset where the right answer is known.
[`make_synthetic.py`](make_synthetic.py) simulates reads from real miRBase miRNAs at
defined abundances and **holds a subset out of the reference** — those held-out miRNAs are
the ground truth for de-novo (novel) detection.

## Files

| File | Role |
|---|---|
| `config.yaml` | all paths and parameters |
| `make_synthetic.py` | simulate reads + write truth tables (`mature` or `genome` mode) |
| `run_engines.sh` | templated commands to run each engine on the synthetic reads |
| `parse_outputs.py` | normalise each engine's output into a common format |
| `score.py` | compute metrics + plots |
| `run_benchmark.py` | orchestrator (synth → engines → score) |

## Run it

```bash
pip install numpy pandas pyyaml matplotlib biopython   # or: conda env create -f ../environment.yml
cd benchmarking

# 1. Simulate the truth set (writes synthetic/)
python run_benchmark.py --config config.yaml --step synth

# 2. Run the engines on synthetic/synthetic_reads.fastq.gz + synthetic/reference_mature.fa
#    (edit run_engines.sh for the engines you have installed)
bash run_engines.sh

# 3. Score everything -> results/
python run_benchmark.py --config config.yaml --step score
```

Each engine's output is read by its parser in `parse_outputs.py`. If you'd rather not rely
on the built-in parsers, drop a `counts.tsv` (`mirna_id, sample, count`) and/or
`predictions.tsv` (`predicted_id, sequence, score`) into an engine's `run_dir` and the
`generic` parser picks them up.

## Metrics produced (`results/metrics_summary.tsv`)

| Metric | Meaning |
|---|---|
| `known_recall` / `known_precision` / `known_f1` | detection of the known truth miRNAs |
| `quant_spearman` / `quant_pearson` | true vs detected abundance (log scale) |
| `novel_precision` / `novel_recall` / `novel_f1` | held-out miRNAs recovered de novo (sequence match ≤ `novel_max_edit` edits) |
| concordance heatmap | truth-free Jaccard between engines (also works on the SRR runs) |

Plots: `quantification_<engine>.png`, `recovery_bar.png`, `concordance_heatmap.png`.

## Two simulation modes

- **`mature`** (default, self-contained) — reads drawn from mature sequences. Solid for
  known recovery + quantification; novel scores are *conservative* because there's no
  genomic hairpin for the folding-based predictors.
- **`genome`** (needs `genome.fa` + miRBase GFF3) — reads drawn from precursor loci, so
  held-out miRNAs can be rediscovered from hairpin structure. Use this for the rigorous
  novel-prediction test that goes in the paper.

## Status

The harness is tested end-to-end on mock data (simulation + scoring verified). The
per-engine parsers target each tool's documented output layout; lines marked `CONFIRM:` in
`parse_outputs.py` should be checked against one real run of each engine, since exact
filenames/columns vary with tool version and miRBase release.

## Truth / real data

- Public SRR runs for concordance: **SRR950892, SRR950893, SRR950894, SRR950895**.
- Synthetic set is regenerated from `config.yaml` (seeded, reproducible).
