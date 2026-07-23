# miRNA target prediction

*Source: migrated from `smallRNA_targets` (two near-identical BioPython drafts merged).*

[`find_mirna_targets.py`](find_mirna_targets.py) reads a CSV of miRNA names + sequences,
BLASTs each against NCBI, and writes candidate hit sites to one CSV.

```bash
pip install biopython pandas
python find_mirna_targets.py -i input.csv -o output.csv
```

`input.csv`:

```csv
name,sequence
miR-159,TTTGGATTGAAGGGAGCTCTA
```

## Status & caveats

- The original scripts were marked "To be tested…". Treat as a **draft**.
- `NCBIWWW.qblast` runs over the network and is rate-limited — fine for a few sequences,
  unworkable at scale. For real runs use **local blastn** against a downloaded database,
  or the sRNAtoolbox consensus target predictors (`miRNAconsTargets`, animal & plant),
  which are the recommended route — see [`../pipelines/srnatoolbox/README.md`](../pipelines/srnatoolbox/README.md).

## TODO
- Add a local-blastn mode (no network dependency).
- Validate against a known miRNA→target set before manuscript use.
