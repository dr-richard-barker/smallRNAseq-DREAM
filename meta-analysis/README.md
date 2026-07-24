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

`fetch_studies.sh` reads `studies.tsv` and dispatches each study by its **`route`** column:
`mirdeep2` (animal small RNA-seq → `../osdr/run_demo.sh`) or `precursor` (plant RNA-seq →
`../smallrna_from_rnaseq/extract_smallrna.sh`, paired-end aware). Then it builds the matrix
and figures.

```bash
N=2 ./fetch_studies.sh                 # all studies, N samples each (trial)
KINGDOM=animal N=2 ./fetch_studies.sh  # just the mouse/human mature route
KINGDOM=plant  N=2 ./fetch_studies.sh  # just the Arabidopsis precursor route

# or, if you already have runs/<OSD>/counts.tsv:
python build_count_matrix.py --runs runs/ --studies studies.tsv --out combined/
python plot_meta.py --combined combined/ --out combined/figures --top 30
```

## Outputs (`combined/`)

| File | What |
|---|---|
| `combined_long.tsv` | tidy: mirna_id, family, sample, count, cpm, study, organism, kingdom, layer, factor |
| `combined_wide.tsv` | miRNA × sample raw-count matrix |
| `combined_wide_cpm.tsv` | miRNA × sample CPM matrix |
| `mirna_prevalence.tsv` | per-ID: how many samples/studies detect it |
| `family_prevalence.tsv` | **per-family** (species prefix stripped) — cross-species pooling |
| `figures/heatmap_top_mirnas.png` | top-variance miRNAs × samples (log CPM), samples colour-barred by kingdom |
| `figures/family_prevalence.png` | families by number of studies; cross-kingdom families highlighted |

Figures are manuscript-ready (150 dpi PNG); regenerate any time with `plot_meta.py`.

## Two things this scaffold does NOT do (on purpose)

- **No batch correction.** CPM is within-sample only. Cross-study comparison needs proper
  batch/library-prep correction (e.g. ComBat-seq) — add it before drawing conclusions.
- **No naive layer pooling.** Every row keeps a `layer` tag: `mature` (animal small RNA-seq)
  vs `precursor` (plant, from RNA-seq). These are different measurement units — the builder
  warns when both are present. Stratify by layer; treat cross-kingdom comparisons as
  precursor-vs-mature, not like-for-like.

Family collapsing (`family_prevalence.tsv`) is a prefix-strip heuristic (`mmu-miR-21-5p` and
`hsa-miR-21-5p` → `mir-21-5p`); verify against miRBase families for anything load-bearing.
