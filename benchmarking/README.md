# Benchmarking — synthetic / SRR accuracy harness

The accuracy work lives here: comparing engine outputs against a truth set to produce the
numbers the manuscript needs. Code only — test data lives on Zenodo / Releases
(see [`../data/MANIFEST.md`](../data/MANIFEST.md)).

## Truth sets used

- Synthetic small RNA-seq reads (miRNAs with known identity/abundance).
- Public SRR runs used during miRDeep2 testing: **SRR950892, SRR950893, SRR950894, SRR950895**.

## What to measure (manuscript gap — currently unquantified)

| Metric | Definition |
|---|---|
| Known-miRNA recovery | fraction of truth-set miRNAs correctly quantified |
| Novel-prediction precision | predicted novel miRNAs that are real (on synthetic set) |
| Novel-prediction recall | real novel miRNAs recovered |
| Cross-species stability | metric spread across kingdoms without retuning |

Run each engine (`../pipelines/*`) on the same inputs, collect counts/predictions, and
score against the truth set. Populate the table in
[`../docs/comparison.md`](../docs/comparison.md).

> This is the last experiment before the paper is submittable: the pipelines "work", but
> the accuracy numbers that justify the recommendation are not yet computed.
