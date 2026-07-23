# Third-party tools

This framework wraps several external tools. **Prefer installing them via conda/mamba or
their own installers** rather than vendoring their source into this repo — that keeps
provenance and licensing clean and avoids carrying GPL code into your own tree.

| Tool | Role | Install | Upstream |
|---|---|---|---|
| miRDeep2 | animal miRNA prediction | `mamba create -n mirdeep2 -c bioconda mirdeep2` | https://github.com/rajewsky-lab/mirdeep2 |
| miRDeep2 patch scripts | performance-improved miRDeep2 scripts | `git clone https://github.com/Drmirdeep/mirdeep2_patch.git` | Drmirdeep |
| miRDeep-P2 (miRDP2) | plant miRNA prediction | https://sourceforge.net/projects/mirdp2/ | — |
| miRDeep-P2_pipeline | plant DE pipeline (fork) | `git clone https://github.com/TF-Chan-Lab/miRDeep-P2_pipeline.git` | TF-Chan-Lab |
| sRNAtoolbox (sRNAbench/sRNAde) | annotation, DE, prediction | https://bioinfo2.ugr.es/srnatoolbox/standalone/ | UGR |
| mirnaQC | QC | https://arn.ugr.es/mirnaqc/ | UGR / Hackenberg |
| Bowtie, ViennaRNA, samtools, BLAST | core deps | `environment.yml` | — |

## Why not vendor?

The original `mirDeep2_accuracy_synthetic_microRNA` repo committed upstream miRDeep2 patch
scripts (`install.pl`, `make_html2.pl`, `miRDeep2_patch.pl`, etc.). These carry the
upstream project's license (miRDeep2 is GPL). Bundling them into a differently-licensed
repo is a licensing hazard. Install them at runtime instead; if you must keep a copy for
reproducibility, place it here **with its original LICENSE file untouched**.

> **TODO:** verify and retain each upstream tool's own license terms before publishing.
