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

The synthetic reads and any curated/derived FASTAs that made the accuracy test
reproducible should be deposited once and referenced by DOI:

| Item | Previously in git as | New home |
|---|---|---|
| Synthetic / trimmed reads | `SRR950892_trimmed.fq.gz` (~3.8 MB) | Zenodo record — **TODO: add DOI** |
| Curated miRBase subsets | `mature_ut*.fa`, `hairpin_ut*.fasta` (~5 MB) | Zenodo record — **TODO: add DOI** |

Once deposited, add:

```
Benchmark data: Zenodo DOI 10.5281/zenodo.XXXXXXX
```

to this file and to [`../benchmarking/README.md`](../benchmarking/README.md).

> Depositing here also gives the consolidated work a citable snapshot for the manuscript,
> and lets the five original repositories be archived without losing the test data.
