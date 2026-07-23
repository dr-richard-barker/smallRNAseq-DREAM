# miRDeep2 pipeline (animal)

*Source: migrated from `mirDeep2_accuracy_synthetic_microRNA`.*

Install-and-run wrapper for miRDeep2, used both as an **animal** miRNA quantification
route and as the **accuracy benchmark** against synthetic / SRR data. The driver script
is [`DRB_mirDeep2_V6.sh`](DRB_mirDeep2_V6.sh).

miRDeep2 tutorial: https://drmirdeep.github.io/mirdeep2_tutorial.html

## Preliminary files required

| File | Description |
|---|---|
| `species.fa` | reference genome assembly (whitespace removed from headers) |
| `mature_ref_this_species.fa` | miRBase mature miRNAs for the species |
| `mature_ref_other_species.fa` | miRBase mature miRNAs for related species |
| `precursors_ref_this_species.fa` | miRBase precursor (hairpin) miRNAs for the species |
| `reads.fa` | the deep-sequencing reads (FASTA) |

## Notes carried over

- miRDeep2 requires **no whitespace** in FASTA headers — hence the `remove_white_space_in_id.pl`
  and `sed 's/ /_/g'` steps.
- Recommended pre-processing not scripted here: TrimGalore, post-trim FastX QC, and
  Kraken2 for contamination assessment.
- Install miRDeep2 via `conda/mamba` (see the script header). **Do not vendor** the
  upstream miRDeep2 patch scripts — see [`../../third_party/README.md`](../../third_party/README.md).

## Data

The original repo carried large reference/test files in git (miRBase dumps, `SRR950892_trimmed.fq.gz`,
hairpin/mature FASTAs, ~5.5 MB). These are **not** migrated into git — see
[`../../data/MANIFEST.md`](../../data/MANIFEST.md).

## Citation

Friedländer M.R., Mackowiak S.D., Li N., Chen W., Rajewsky N. *miRDeep2 accurately
identifies known and hundreds of novel microRNA genes in seven animal clades.*
Nucleic Acids Research, 40:37–52, 2012.
