# Notice — licensing

This repository consolidates the author's own wrapper code, documentation and analysis
plans with references to several third-party bioinformatics tools.

- **Author's own material** (documentation, wrapper scripts, the merged target-prediction
  script, benchmarking harness): choose a license and record it in `LICENSE`. A permissive
  license (MIT / BSD-3) is suggested for maximum reuse — **TODO: author to confirm.**

- **Third-party tools** (miRDeep2, miRDeep-P2, sRNAtoolbox, mirnaQC, Bowtie, etc.) are
  **not** relicensed here. They are installed from upstream and retain their own licenses.
  See [`third_party/README.md`](third_party/README.md).

- **`pipelines/mirdp2_plant/`** derives from `TF-Chan-Lab/miRDeep-P2_pipeline`. If any of
  its Perl scripts are copied into this repo, its upstream LICENSE **must** be retained in
  that directory.

> **TODO:** add a top-level `LICENSE` file and confirm compatibility with any vendored
> third-party code before making the repository public.
