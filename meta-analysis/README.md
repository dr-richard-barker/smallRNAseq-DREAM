# Meta-analysis — combined miRNA count matrix across OSDR studies

Pulls several spaceflight/radiation small RNA-seq studies from OSDR, runs each through the
pipeline, and merges the per-study counts into one matrix for cross-study (and cross-species)
meta-analysis.

## The cohort ([`studies.tsv`](studies.tsv))

Verified against the live OSDR API (raw FASTQ present):

| OSD | Organism | Layer | Factor | Raw FASTQ |
|---|---|---|---|---|
| OSD-334/335/336/337 | *Mus musculus* | mature | HZE radiation | ~78–80 each |
| OSD-483 | *Homo sapiens* | mature | Spaceflight (STS), astronaut sEV | 42 |
| OSD-208, OSD-437 | *Arabidopsis thaliana* | **precursor** | Spaceflight / microgravity | — (see below) |

The mouse + human studies are genuine small RNA-seq (mature miRNA). The Arabidopsis studies
are **standard RNA-seq, not small RNA-seq** — handled by the precursor-level route in
[`../smallrna_from_rnaseq/`](../smallrna_from_rnaseq/).

## Run it

```bash
# 1. process the animal studies (real toolchain, your compute) -> runs/<OSD>/counts.tsv
N=2 ./fetch_studies.sh            # N samples/study for a trial; drop N for all samples

# 2. or build the matrix from whatever runs/ you already have:
python build_count_matrix.py --runs runs/ --studies studies.tsv --out combined/
```

## Outputs (`combined/`)

| File | What |
|---|---|
| `combined_long.tsv` | tidy: mirna_id, family, sample, count, cpm, study, organism, kingdom, layer, factor |
| `combined_wide.tsv` | miRNA × sample raw-count matrix |
| `combined_wide_cpm.tsv` | miRNA × sample CPM matrix |
| `mirna_prevalence.tsv` | per-ID: how many samples/studies detect it |
| `family_prevalence.tsv` | **per-family** (species prefix stripped) — cross-species pooling |

## Two things this scaffold does NOT do (on purpose)

- **No batch correction.** CPM is within-sample only. Cross-study comparison needs proper
  batch/library-prep correction (e.g. ComBat-seq) — add it before drawing conclusions.
- **No naive layer pooling.** Every row keeps a `layer` tag: `mature` (animal small RNA-seq)
  vs `precursor` (plant, from RNA-seq). These are different measurement units — the builder
  warns when both are present. Stratify by layer; treat cross-kingdom comparisons as
  precursor-vs-mature, not like-for-like.

Family collapsing (`family_prevalence.tsv`) is a prefix-strip heuristic (`mmu-miR-21-5p` and
`hsa-miR-21-5p` → `mir-21-5p`); verify against miRBase families for anything load-bearing.
