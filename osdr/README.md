# OSDR integration — pulling spaceflight small RNA-seq data

[`osdr_fetch.py`](osdr_fetch.py) queries the **NASA Open Science Data Repository (OSDR)**
public API to find and download small RNA-seq data for the smallRNAseq-DREAM pipelines.
Standard-library only; validated against the live API.

```bash
python osdr_fetch.py search                       # find + classify small-RNA studies
python osdr_fetch.py files OSD-483                 # list a study's files
python osdr_fetch.py meta  OSD-483                 # metadata summary
python osdr_fetch.py download OSD-483 --raw -o data/OSD-483   # download raw FASTQ
python osdr_fetch.py download OSD-483 --raw --dry-run         # preview sizes first
```

Endpoints used (all public):

| Purpose | Endpoint |
|---|---|
| Search | `https://osdr.nasa.gov/osdr/data/search?term=<t>&type=cgene&size=<n>` |
| Files | `https://osdr.nasa.gov/osdr/data/osd/files/<n>` |
| Download | `https://osdr.nasa.gov<remote_url>` (from the files listing) |

## Finding: plant small RNA-seq is a gap in OSDR (surveyed 2026-07)

A survey with this tool ([`osdr_smallrna_survey.tsv`](osdr_smallrna_survey.tsv)) returned
**176 studies** across OSDR's small-RNA search terms:

| Kingdom | Studies returned |
|---|---|
| Animal | 156 |
| Plant | 12 |
| Microbe | 7 |
| Other | 1 |

**None of the 12 "plant" studies are actually small RNA-seq** — file-level inspection
(`osdr_fetch.py files ...`) shows every one is DNA microarray or bulk mRNA-seq; the
small-RNA terms matched incidental text, not the assay. So **OSDR currently contains no
plant spaceflight small RNA-seq dataset.**

Existing spaceflight small RNA-seq is animal — mouse-dominated, plus astronaut
extracellular-vesicle miRNA. The confirmed demonstration dataset with raw reads is:

- **OSD-483** — *small-RNA sequencing of sEV isolated from plasma of astronauts*
  (*Homo sapiens*, ISS/STS): 42 raw FASTQ files (~37 GB total) downloadable via the API.

## What this means for the framework

1. **The tool can be validated now** on real OSDR animal small RNA-seq (OSD-483 + mouse
   studies) — the animal (miRDeep2) arm of the framework.
2. **The plant (miRDeep-P2) arm is positioned to fill a documented gap**: when plant
   spaceflight small RNA-seq is deposited (e.g. from VEGGIE/APH-type missions, or the
   author's own datasets), smallRNAseq-DREAM is ready to process it cross-kingdom.
3. **Meta-analysis angle**: multiple animal small-RNA studies across missions/timepoints
   can be pulled and integrated now; plant is the forward-looking target.

> The heavy pipeline execution (alignment, miRNA calling) runs on your compute — this
> module handles discovery and retrieval, then hands FASTQ to `pipelines/`.
