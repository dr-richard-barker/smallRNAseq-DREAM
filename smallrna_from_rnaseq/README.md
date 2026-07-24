# Extracting miRNA signal from standard RNA-seq (the plant route)

OSDR has **no plant small RNA-seq** ([../osdr/README.md](../osdr/README.md)), but it has
Arabidopsis **RNA-seq** (OSD-208, OSD-437). This module recovers what miRNA signal it can
from those standard libraries — with the biology made explicit so the numbers aren't
over-interpreted.

## What is and isn't recoverable — the biology

| | Mature miRNA (~21 nt) | Precursor / MIR locus (pri-/pre-miRNA) |
|---|---|---|
| In standard RNA-seq? | **No** — reads are fragmented to ~200 bp+ and sequenced at 50–150 bp, so a read cannot *be* a 21 nt mature miRNA | **Yes** — plant pri-miRNAs are Pol II transcripts (capped + polyadenylated), so poly-A / total RNA-seq carries them |
| What this module reports | a **diagnostic** only (read-length distribution; fraction ≤30 nt) | **counts per precursor** — the informative layer |

So "small RNAs detected in RNA-seq" in practice means **MIR-locus / precursor expression**,
a proxy for miRNA *gene activity* — not mature miRNA abundance. The library preps confirm
the constraint:

- **OSD-437** — NEBNext Ultra II RNA (standard mRNA-seq; rRNA only measured, not depleted).
- **OSD-208** — Ovation Pico WTA + Nextera, total-RNA input.

Poly-A selection + fragmentation + size selection all bias against the mature ~21 nt product,
while retaining the polyadenylated precursor. The read-length diagnostic (`step 2`) lets you
**see this empirically** per sample rather than assume it.

## Run it

```bash
OSD=OSD-437 N=2 ./extract_smallrna.sh     # your compute; needs seqkit, bowtie2, samtools
```

Steps: fetch RNA-seq → read-length diagnostic → pull *ath* miRBase hairpins+mature (U→T) →
align (`bowtie2 --local`) to precursors → per-precursor counts → `counts.tsv`
(`mirna_id, sample, count`), tagged **layer = precursor** for the meta-analysis.

## Feeding meta-analysis

`counts.tsv` drops straight into [`../meta-analysis/`](../meta-analysis/). Keep it under a
`runs/OSD-437/` dir and the builder tags it `layer=precursor` from `studies.tsv`. The builder
**warns when mature and precursor layers are mixed** and never pools them naively — cross-kingdom
comparison is precursor-vs-mature, not like-for-like.

## Caveats (don't skip)

- **Precursor ≠ mature.** A highly expressed MIR locus does not guarantee an abundant mature
  miRNA (processing/turnover intervene). Treat these as gene-level miRNA activity.
- **A more rigorous route** aligns to the genome (HISAT2/STAR) and counts over miRBase MIR
  GFF loci with `featureCounts` — better for multi-exon / overlapping loci than aligning to
  hairpin sequences directly. The bowtie2-to-hairpin approach here is a fast scaffold.
- **Validate.** If precursor signal looks interesting, confirm with dedicated plant small
  RNA-seq (the gap this framework is built to fill once such data exists).
