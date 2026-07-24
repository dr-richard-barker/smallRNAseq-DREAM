# Data manifest

Large files are **not** stored in git. This manifest lists what each pipeline needs, where
to get it, and (for the fixed test set) where the archived copy will live.

## Re-downloadable references (do not commit)

| File | Source |
|---|---|
| miRBase mature miRNAs | https://mirbase.org/download/CURRENT/mature.fa |
| miRBase hairpin miRNAs | https://mirbase.org/download/CURRENT/hairpin.fa |
| Human genome (GRCh38) | https://www.ncbi.nlm.nih.gov/genome/guide/human/ |
| Arabidopsis genome (TAIR/GCF_000001735.3) | https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001735.3/ |
| SRR runs (benchmark) | `prefetch SRR950892 SRR950893 SRR950894 SRR950895` (sra-tools) |

## Fixed benchmark set → Zenodo / GitHub Release

The synthetic reads and curated/derived FASTAs that made the accuracy test reproducible
have been **assembled into a Zenodo-ready bundle** (rescued from
`mirDeep2_accuracy_synthetic_microRNA` @ `51e15be` before archiving):

> **Bundle location (not in git):** `../../smallRNAseq-DREAM-zenodo/`
> — 7 files, ~10 MB, with `MANIFEST.tsv`, `CHECKSUMS.sha256`, `.zenodo.json` and
> `LICENSE_NOTES.md`. See that folder's `UPLOAD_INSTRUCTIONS.md` (web UI + scripted API
> route via `zenodo_upload.py`) to deposit it and mint the DOI.

| Item | File(s) | In bundle |
|---|---|---|
| Trimmed reads (SRR950892, 228,549 reads) | `SRR950892_trimmed.fq.gz` (3.8 MB) | `reads/` |
| miRBase mature, all species, U→T | `mature_ut.fa` (48,885 seqs) | `reference/` |
| miRBase human / non-human / hairpin subsets | `mature*_hsa*`, `mature*No-HSA*`, `hairpin_ut*` | `reference/` |

Once deposited, add:

```
Benchmark data: Zenodo DOI 10.5281/zenodo.XXXXXXX
```

to this file and to [`../benchmarking/README.md`](../benchmarking/README.md).

> Depositing here also gives the consolidated work a citable snapshot for the manuscript,
> and lets the five original repositories be archived without losing the test data.
