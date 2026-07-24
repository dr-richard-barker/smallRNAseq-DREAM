# Manuscript outline

*Working outline for the smallRNAseq-DREAM methods/benchmarking paper. Grounded in the
material already in this repository. Everything under "Results" that depends on numbers is
marked **[NEEDS BENCHMARK]** — do not write results until `benchmarking/` has been run for
real. No numbers are invented here.*

---

## Remaining fill-ins before submission

Drafted in full prose: Abstract, §1 Introduction, §2 Methods (2.1–2.8), §4 Discussion,
§5 Conclusions, and captions for Fig 1–5 + S1–S2. What is left, grouped by what unblocks it:

- [ ] **Run the benchmark** (`benchmarking/` on all four engines) — the single biggest gate.
  Unblocks: **§3 Results** (3.2–3.5 prose), **Fig 2, Fig 3, Fig S1**, **Table 3**, the Abstract
  and Discussion §4.1 headline numbers, and the `[N]` sample/study counts in Fig 4/5.
  All `**[NEEDS BENCHMARK]**` tags resolve here. First verify the `CONFIRM:` parser lines in
  `benchmarking/parse_outputs.py` against one real run per engine.
- [ ] **Deposit the Zenodo bundle** → replace `[TODO]` DOI in the Abstract, Availability,
  `data/MANIFEST.md`, `CITATION.cff`, and the archive-notice banner. Cite the **concept DOI**.
- [ ] **Citations** — fill every `[ref]` (Intro §1.1–1.4) and `[Barker et al., ref]` (author's
  OSDR/GeneLab work), plus Table 1 tool citations (full authors/years/DOIs).
- [ ] **Tool versions / DB release** — fill each `[version]` (nf-core, sRNAtoolbox, miRDeep2,
  miRDeep-P2, Bowtie 2, samtools, seqkit, sra-tools, Nextflow) and `[release]` (miRBase) from
  the run environment. `environment.yml` pins Bowtie 1.3.1 / ViennaRNA 2.6.4.
- [ ] **Survey count `[N]`** — the OSDR small-RNA study total (176 at survey time; refresh from
  `osdr/osdr_smallrna_survey.tsv`) in §2.7 and §3.7.
- [ ] **Fig 1** — swap in the original schematic if a backup surfaces (current panel is a
  reproducible remake; original lost with the deleted hub repo).
- [ ] **Fig S2** — optional; include only if a worked miRNA→target example is run.
- [ ] **Front/back matter** — author list + ORCID, funding/acknowledgements, competing
  interests; **decide target journal** → set length/format (full paper vs Application Note).

---

## Working titles

1. *smallRNAseq-DREAM: a FAIR, cross-species small RNA-seq analysis framework benchmarked
   across four miRNA pipelines for spaceflight omics*
2. *Choosing a small RNA-seq pipeline for cross-kingdom spaceflight biology: a benchmark of
   nf-core/smrnaseq, sRNAtoolbox, miRDeep2 and miRDeep-P2*
3. *A reproducible small RNA-seq workflow for NASA OSDR reanalysis: benchmarking known-miRNA
   quantification and novel-miRNA discovery across plants and animals*

## Target journals (suggested, fit noted)

- **GigaScience / GigaByte** — tools + FAIR data + reproducibility; strong fit for a
  benchmarked, containerised, openly-deposited framework.
- **NAR Genomics and Bioinformatics** — methods/benchmark audience.
- **BMC Bioinformatics** — pipeline + benchmark, established format.
- **npj Microgravity** or **Life (MDPI)** — if the spaceflight application is the headline
  rather than the method.
- **Bioinformatics (Application Note)** — if kept short (2 pages + the tool).

---

## Abstract (structured, ~250 words) *(draft prose)*

**Background.** Small RNA sequencing (small RNA-seq) is the standard method for profiling
microRNAs (miRNAs), but its analysis is unusually tool-dependent, and because miRNA biogenesis
differs between animals and plants no single pipeline is agreed for cross-species work.
Spaceflight omics in NASA's Open Science Data Repository (OSDR) span both kingdoms, yet FAIR,
reproducible reanalysis of their small-RNA data remains difficult.

**Results.** We consolidated five overlapping analysis efforts into one framework,
deconstructed the community nf-core/smrnaseq pipeline, and benchmarked four engines
(nf-core/smrnaseq, sRNAtoolbox, miRDeep2 and miRDeep-P2) against synthetic data with absolute
ground truth, scoring known-miRNA recovery, quantification accuracy and novel-miRNA discovery
in both animal and plant modes [NEEDS BENCHMARK: headline recall/precision/correlation]. On
this basis we recommend the cross-species route SRA → mirnaQC → sRNAbench → sRNAde → mirNOVO.
Applying the framework to OSDR, we find it currently contains no plant small RNA-seq; we
therefore demonstrate a cross-study meta-analysis of the available animal spaceflight and
radiation miRNA datasets, and recover precursor-level miRNA signal from standard *Arabidopsis*
RNA-seq, keeping mature and precursor measurements as distinct layers.

**Conclusions.** smallRNAseq-DREAM is a FAIR, containerised small RNA-seq framework that
provides an evidence-based cross-species pipeline recommendation, a reusable synthetic
benchmark, and NASA OSDR integration — reproducible today on animal spaceflight data and
positioned to process plant spaceflight small RNA-seq as it becomes available.

**Availability.** Code: https://github.com/dr-richard-barker/smallRNAseq-DREAM (project site
https://dr-richard-barker.github.io/smallRNAseq-DREAM/); reference and benchmark data on Zenodo
(DOI [TODO]).

## Keywords

microRNA; small RNA-seq; spaceflight; NASA OSDR/GeneLab; benchmarking; FAIR; miRDeep2;
sRNAtoolbox; cross-species; isomiR.

---

## 1. Introduction

*Draft prose. Literature citations are marked `[ref]`; the author's own related work is marked
`[Barker et al., ref]`. Fill from the reference manager before submission.*

### 1.1 microRNAs and the challenges of small RNA-seq

MicroRNAs (miRNAs) are ~21–22-nucleotide endogenous non-coding RNAs that regulate gene
expression post-transcriptionally and are central to development and stress responses across
eukaryotes [ref]. Measuring them by sequencing (small RNA-seq) carries analytical challenges
that set it apart from bulk RNA-seq. Because the mature molecule is shorter than the sequenced
fragment, the 3′ adapter must be detected and removed before reads can be mapped, and identical
reads are collapsed to unique sequences for quantification. Mature miRNAs must then be
distinguished both from other small-RNA classes (small interfering RNAs, piwi-interacting RNAs,
and tRNA- and rRNA-derived fragments) and from the length and sequence variants of a single
miRNA (isomiRs) [ref]. Identifying *novel* miRNAs is harder still: it requires evaluating the
predicted secondary (hairpin) structure of the candidate precursor rather than relying on
alignment alone [ref]. As a result, numerous specialised tools have been developed, and
reported miRNA repertoires and abundances can depend materially on the choices made at each
analytical step.

### 1.2 Kingdom-specific miRNA biology

miRNA biogenesis differs fundamentally between animals and plants, and this difference has
direct consequences for analysis. In animals, primary transcripts are processed sequentially
by the Drosha/DGCR8 complex and Dicer, and mature miRNAs typically repress targets through
seed-region base-pairing, most often in 3′ untranslated regions [ref]. In plants, processing is
carried out largely by a single Dicer-like enzyme (DCL1), mature miRNAs are 2′-O-methylated by
HEN1, and they generally direct cleavage of highly complementary target transcripts [ref].
These distinctions are embedded in the prediction software: animal-oriented callers such as
miRDeep2 and plant-specific callers such as miRDeep-P2 apply different precursor length,
structure and conservation criteria [ref]. Any framework intended to operate across kingdoms
must therefore accommodate both models rather than assume one.

### 1.3 Spaceflight biology and the reproducibility imperative

Spaceflight exposes organisms to microgravity and ionising radiation, provoking stress and
adaptive responses in which miRNAs have been implicated in both animal and plant systems [ref].
NASA's Open Science Data Repository (OSDR; formerly the GeneLab Data System) makes spaceflight
omics data openly available across a wide range of organisms, creating an opportunity for
cross-species and even cross-kingdom reanalysis [ref]. Realising that opportunity depends on
FAIR (Findable, Accessible, Interoperable, Reusable) and reproducible analysis pipelines
[Barker et al., ref]. For small RNA-seq specifically, that goal is undercut by the tool
heterogeneity described above and by the practical gap between convenient browser-based
services and the command-line, containerised workflows needed for reproducible, high-throughput
reanalysis of repository-scale data.

### 1.4 This work

We address these needs with a consolidated, cross-species small RNA-seq analysis framework. We
(i) deconstruct the community nf-core/smrnaseq pipeline to expose its per-step structure and
tool choices; (ii) benchmark four analysis engines (nf-core/smrnaseq, sRNAtoolbox, miRDeep2 and
miRDeep-P2) against synthetic data with absolute ground truth, quantifying known-miRNA recovery,
quantification accuracy and novel-miRNA discovery across kingdoms; (iii) recommend a single
cross-species route, SRA → mirnaQC → sRNAbench → sRNAde → mirNOVO; and (iv) package the result
as a FAIR, containerised, openly deposited framework and demonstrate it on NASA OSDR
spaceflight and radiation datasets. In doing so we document that OSDR currently contains no
plant small RNA-seq, and show how the framework both integrates the available animal
small RNA-seq into a cross-study meta-analysis and recovers precursor-level miRNA signal from
existing plant RNA-seq — positioning it to process dedicated plant small RNA-seq as such data
become available.

## 2. Materials and Methods

### 2.1 Pipelines compared *(draft prose)*

Four small RNA-seq analysis engines were evaluated within a single framework, each
implemented as an independent module ([`../pipelines/`](../pipelines/)). The community
nf-core/smrnaseq pipeline [version] was used as the baseline and deconstructed into its
component steps (Fig 1). Three further engines were run under a common interface: the
sRNAtoolbox suite (mirnaQC, sRNAbench, sRNAde and mirNOVO) [version]; miRDeep2 [version] for
animal miRNA prediction; and miRDeep-P2 (miRDP2) [version] for plant miRNA prediction, the
latter derived from the TF-Chan-Lab pipeline (retaining upstream attribution and licence). On
the basis of this comparison we defined a recommended cross-species route,
SRA → mirnaQC → sRNAbench → sRNAde → mirNOVO. Tool versions, roles, kingdom applicability and
citations are given in Table 1 ([`methods_tables.md`](methods_tables.md)).

### 2.2 Synthetic benchmark design *(draft prose)*

Because real sequencing runs lack a known truth set, engine accuracy was assessed on
synthetic small RNA-seq data with defined ground truth (`benchmarking/make_synthetic.py`).
Reads were simulated from miRBase [release] mature sequences: a fixed number of miRNAs were
retained in the reference supplied to each engine (known set) while a disjoint subset was
withheld to serve as ground-truth novel miRNAs. Per-miRNA abundances were drawn from a
log-normal distribution and reads generated to a target depth, incorporating isomiR variation
(1–2-nt 5′/3′ trimming or templated/non-templated addition), per-base sequencing error, 3′
adapter read-through, and a defined fraction of non-miRNA decoy reads. Two simulation modes
were provided: a self-contained *mature* mode (reads drawn from mature sequences; used for
known-miRNA recovery and quantification) and a *genome-anchored* mode (reads drawn from
precursor loci via a miRBase GFF3, so that folding-based predictors can rediscover held-out
miRNAs from hairpin structure; used for the novel-discovery comparison, Fig 3). All parameters
are fixed and version-controlled in a seeded configuration file (`benchmarking/config.yaml`;
defaults: 150 known and 30 held-out miRNAs, 2 × 10⁶ reads, 20 % decoy fraction, 50-nt reads,
0.1 % error rate, 30 % isomiR rate), making each dataset exactly reproducible.

### 2.3 Real datasets *(draft prose)*

Two sources of real data were used. Public human small RNA-seq runs SRR950892–SRR950895 were
obtained from NCBI-SRA and used during miRDeep2 benchmarking; and spaceflight/radiation
datasets were retrieved from NASA OSDR as described in §2.7 (animal small RNA-seq OSD-334–337
and OSD-483; *Arabidopsis* RNA-seq OSD-208 and OSD-437). Reads were acquired with sra-tools
[version] (`prefetch`/`fastq-dump`) or the OSDR API, converted to FASTA where required
(seqkit [version]), and adapter-trimmed; the small-RNA 3′ adapter and read-length windows were
set per study from the associated metadata (benchmark default adapter TGGAATTCTCGGGTGCCAAGG,
minimum read length 18 nt). FASTA headers were stripped of whitespace as required by miRDeep2.

### 2.4 Scoring *(draft prose)*

Each engine's native output was normalised to a common format — a per-miRNA count table
(`mirna_id, sample, count`) and a table of predicted novel miRNAs with sequences — using
engine-specific parsers (`benchmarking/parse_outputs.py`). Four metrics were then computed
(`benchmarking/score.py`). Known-miRNA detection was scored as recall, precision and F1, with
a miRNA counted as detected above a fixed read threshold (default 5). Quantification accuracy
was assessed as Spearman and Pearson correlations between simulated abundance and detected
count on a log(1 + x) scale. Novel-miRNA discovery was scored as precision and recall against
the held-out truth set, a prediction counted as correct when its mature sequence lay within a
bounded edit distance (≤2, allowing a sliding offset) of a held-out sequence — mirroring the
miRDeep-P2 "variant" criterion. Finally, cross-engine concordance was computed as the pairwise
Jaccard index of detected-miRNA sets, a truth-free measure applicable to real runs.

### 2.5 miRNA target prediction *(draft prose)*

Predicted and known miRNAs were screened for targets using the sRNAtoolbox consensus target
predictors (miRNAconsTargets), which provide separate animal and plant models, complemented by
a BLAST-based search of candidate miRNA sequences against a nucleotide database
([`../target_prediction/`](../target_prediction/)). The BLAST route is provided as an
exploratory utility; for production use a local BLAST database or the consensus predictors are
recommended over networked queries, and predicted targets require experimental validation.

### 2.6 Deployment and reproducibility *(draft prose)*

The framework is distributed as version-controlled code with a conda environment
specification (`environment.yml`; e.g. Bowtie 1.3.1, ViennaRNA 2.6.4). Engines are run in
containers for portability — Singularity for the sRNAtoolbox Java tools and Docker images
where available — and the nf-core route is executed under Nextflow [version] with
Docker/Singularity profiles; a SLURM configuration supports cluster and OSDR execution. In
keeping with FAIR principles, all code is openly available on GitHub, large reference and
test-data files are deposited to Zenodo (DOI [TODO]) rather than held in version control, and
the synthetic benchmark is fully regenerable from its seeded configuration.

### 2.7 OSDR data retrieval and cross-study integration *(draft prose)*

Spaceflight and radiation datasets were retrieved from the NASA Open Science Data Repository
(OSDR) through its public REST API using a purpose-written client (`osdr/osdr_fetch.py`;
Python standard library only). Candidate studies were identified with the search endpoint
(`/osdr/data/search`, `type=cgene`) across the query terms *small RNA*, *microRNA*, *sRNA-seq*
and *miRNA*, and each returned study was classified by organism into kingdom (animal, plant,
microbe). This survey returned [N] studies (see Results); per-study file listings were then
obtained from the files endpoint (`/osdr/data/osd/files/{id}`) and raw sequencing reads
downloaded from the URLs it provides. For the meta-analysis we used the animal small RNA-seq
cohort OSD-334, OSD-335, OSD-336 and OSD-337 (*Mus musculus*, high-charge-and-energy [HZE]
radiation) and OSD-483 (*Homo sapiens*, astronaut plasma extracellular-vesicle small RNA-seq),
together with the *Arabidopsis thaliana* RNA-seq studies OSD-208 and OSD-437 (Table 4). Each
study was processed to a per-miRNA count table in a common long format
(`mirna_id, sample, count`); tables were merged across studies with `meta-analysis/
build_count_matrix.py`, which computes within-sample counts-per-million (CPM) and, to allow
homologous miRNAs to be compared across species, collapses organism prefixes to a family label
(e.g. mmu-miR-21-5p and hsa-miR-21-5p → miR-21-5p). Because library preparation, sequencing
depth and organism differ between studies, CPM is reported only for within-sample scaling and
no cross-study batch correction was applied; comparisons are therefore treated as
hypothesis-generating. Counts derived from mature small RNA-seq and from RNA-seq precursor
recovery (§2.8) are tagged with a `layer` field and never pooled (see §2.8).

### 2.8 Recovering miRNA precursor signal from standard RNA-seq (plant route) *(draft prose)*

No plant small RNA-seq dataset was present in OSDR at the time of analysis (§3.7); the
available *Arabidopsis* spaceflight data are standard RNA-seq (OSD-208, Ovation Pico WTA;
OSD-437, NEBNext Ultra II RNA). Because such libraries are fragmented and sequenced at read
lengths well above the ~21-nt mature-miRNA size, mature miRNAs cannot be observed directly;
however, plant primary miRNAs (pri-miRNAs) are RNA-polymerase-II transcripts that are capped
and polyadenylated, so poly-A and total-RNA libraries retain miRNA-locus (precursor) signal.
We therefore quantified precursor-level miRNA expression (`smallrna_from_rnaseq/
extract_smallrna.sh`). Raw reads were retrieved as above; per sample, a read-length
distribution was recorded (seqkit [version]) and the fraction of reads ≤30 nt reported as a
diagnostic of whether any mature-length small RNAs survived library preparation. *Arabidopsis*
(ath) mature and hairpin (precursor) sequences were obtained from miRBase [release] and
converted from RNA to DNA. Reads were aligned to the hairpin precursors with Bowtie 2 [version]
in local-alignment mode (`--local --no-unal`); paired-end libraries were aligned as read pairs
(`--no-mixed --no-discordant`) and quantified as fragments rather than mates to avoid
double-counting. Per-precursor read counts (samtools [version] `idxstats`) were written to the
same common long format with `layer = precursor`, and integrated with the animal cohort as in
§2.7. These values represent miRNA-gene (precursor) activity rather than mature-miRNA
abundance and are analysed as a distinct measurement layer; a more rigorous alternative aligns
reads to the genome and counts over miRBase MIR loci with a feature-counting step, which is
preferable for multi-exon or overlapping loci.

## 3. Results

3.1 **Deconstructing nf-core/smrnaseq** — the step map (Plan A → tidied Plan B); Fig 1.
3.2 **Known-miRNA quantification accuracy** — per-engine recall/precision and abundance
    correlation on synthetic data. **[NEEDS BENCHMARK]** Fig 2, Table 3.
3.3 **Novel-miRNA discovery** — precision/recall on held-out miRNAs; where each engine
    over-/under-calls. **[NEEDS BENCHMARK]** (use `genome` mode for the fair test.)
3.4 **Cross-kingdom comparison** — animal (miRDeep2) vs plant (miRDeep-P2); does the sRNAtoolbox
    route hold across kingdoms without retuning? **[NEEDS BENCHMARK]** Fig 3.
3.5 **Concordance on real spaceflight data** — engine agreement on SRR/OSDR runs (truth-free).
    **[NEEDS BENCHMARK]** Fig S1.
3.6 **Recommended route** — evidence for SRA → mirnaQC → sRNAbench → sRNAde → mirNOVO as the
    cross-species default; throughput/deployment notes.
3.7 **Cross-study, cross-species meta-analysis (OSDR demonstration)** — a documented finding
    that OSDR holds no plant small RNA-seq (176 small-RNA studies surveyed: 156 animal, 12
    plant — all microarray/mRNA-seq, 0 true small RNA-seq); the framework instead (i) integrates
    a real animal spaceflight/radiation miRNA cohort (mouse HZE OSD-334–337 + astronaut sEV
    OSD-483) into one matrix with cross-species family pooling (**Fig 4, Fig 5, Table 4**), and
    (ii) recovers plant miRNA *precursor* signal from standard Arabidopsis RNA-seq (OSD-208/437),
    kept as a separate measurement layer. **[FIGURES: tool tested; panels from the full run]**

## 4. Discussion

*Draft prose. Claims that depend on benchmark values are written conditionally and marked
**[NEEDS BENCHMARK]**; design- and survey-based statements stand as written.*

### 4.1 An integrated route for cross-species reanalysis

Deconstructing nf-core/smrnaseq (Fig 1) makes clear that a small RNA-seq workflow is a chain of
separable decisions — QC, adapter handling, alignment, quantification, differential expression
and novel-miRNA discovery — each with several tool options. The route we recommend
(SRA → mirnaQC → sRNAbench → sRNAde → mirNOVO) collapses most of these phases into a few
interoperable tools, reducing the integration and glue-code burden of assembling the full
multi-step pipeline while retaining known-miRNA quantification, isomiR annotation, differential
expression and de-novo discovery. The case for this route rests on three axes — accuracy,
cross-species generality and throughput. Its throughput and generality advantages follow from
its design and containerised deployment; whether it matches or exceeds the alternatives on
accuracy is the question the benchmark answers directly, and should be stated from those
results [NEEDS BENCHMARK: relative known-miRNA recovery, quantification correlation and
novel-discovery precision/recall across engines].

### 4.2 Analysis must respect kingdom-specific biology

Because miRNA biogenesis differs between animals and plants (§1.2), no single predictor is
appropriate for both, and the framework keeps miRDeep2 (animal) and miRDeep-P2 (plant) as
parallel arms rather than forcing one model across kingdoms. This choice is reflected even in
evaluation: our scoring counts a predicted novel miRNA as correct when its mature sequence lies
within two edits of a held-out truth sequence, mirroring the miRDeep-P2 "variant" criterion, so
that the same yardstick is defensible in both kingdoms. A framework that ignores these
differences risks systematically mis-calling novel miRNAs in whichever kingdom it was not tuned
for.

### 4.3 A gap in the spaceflight record: no plant small RNA-seq

Surveying OSDR returned no plant small RNA-seq dataset: every plant study among the small-RNA
search hits is a microarray or bulk mRNA-seq experiment, whereas the small RNA-seq record is
animal-dominated. Given the breadth of plant spaceflight biology, this is a substantive gap —
miRNA-level regulation of plant responses to microgravity and radiation is effectively
uncharacterised by dedicated sequencing in the public spaceflight record. It argues for
generating and depositing plant small RNA-seq from missions where plant material is already
flown. In the interim, recovering miRNA-precursor signal from existing plant RNA-seq offers a
partial, hypothesis-generating view: because plant pri-miRNAs are polyadenylated Pol II
transcripts, standard libraries retain MIR-locus expression even though the mature ~21-nt
product is lost to fragmentation and size selection. We stress that these precursor counts
measure miRNA-gene activity, not mature-miRNA abundance — precursor level and mature level can
be decoupled by processing and turnover — and we therefore keep the two as distinct measurement
layers and never pool them. The per-sample read-length diagnostic makes the underlying
limitation visible rather than assumed.

### 4.4 Limitations

Several limitations bound the interpretation. The benchmark uses simulated reads: while these
provide the absolute ground truth real runs lack, they may not reproduce every library artefact,
and the self-contained *mature* simulation mode is deliberately conservative for novel discovery
(the genome-anchored mode is the fair test). Results depend on miRBase completeness and release,
and the per-engine output parsers are sensitive to tool version and output layout (flagged in
the code) and should be re-checked against one real run of each engine. In the meta-analysis,
counts are CPM-normalised within sample only; no cross-study batch correction is applied, so
cross-study patterns are hypothesis-generating, and the family-level pooling used for
cross-species comparison is a prefix-strip heuristic that should be verified against curated
miRBase families for any load-bearing claim. Novel-miRNA predictions, in either kingdom, require
experimental validation (e.g. qPCR or dedicated small RNA-seq). Finally, the OSDR demonstration
is necessarily on animal data, because the plant small RNA-seq it would ideally use does not yet
exist.

### 4.5 Fit and outlook

The framework is intended to slot into the broader OSDR/GeneLab reanalysis and cross-species
meta-analysis effort: studies are added by listing an accession, and the same common count
format carries animal mature-miRNA and plant precursor data through to a joint, layer-aware
matrix. Its immediate utility is reproducible reanalysis of the existing animal spaceflight and
radiation miRNA datasets; its forward value is that it is ready to process dedicated plant
spaceflight small RNA-seq the moment such data are deposited — turning the gap identified here
into a concrete target for the community.

## 5. Conclusions *(draft prose)*

smallRNAseq-DREAM consolidates five overlapping efforts into a single, FAIR, containerised
framework for cross-species small RNA-seq analysis. It deconstructs the nf-core/smrnaseq
pipeline, benchmarks four engines against synthetic data with absolute ground truth, and — on
that evidence [NEEDS BENCHMARK] — recommends the route SRA → mirnaQC → sRNAbench → sRNAde →
mirNOVO for reanalysis across kingdoms. Applied to NASA OSDR, it integrates the available
animal spaceflight and radiation miRNA datasets into a cross-study, layer-aware meta-analysis
and recovers precursor-level miRNA signal from existing plant RNA-seq, while documenting that
dedicated plant spaceflight small RNA-seq is still absent. With its code, a reusable synthetic
benchmark and deposited reference data openly available, the framework is reproducible on
animal spaceflight data today and positioned to process plant spaceflight small RNA-seq as soon
as it is generated.

---

## Figures & Tables

| # | Item | Status | Source |
|---|---|---|---|
| Fig 1 | nf-core/smrnaseq deconstruction vs. DREAM route | **have (remake)** — regenerated as `fig1_pipeline_overview.svg`; original PNG/PDF lost with the deleted hub repo | `docs/schematics/` |
| Fig 2 | Known-miRNA recall/precision + quantification scatter | NEEDS BENCHMARK | `benchmarking/results/` |
| Fig 3 | Cross-kingdom (animal vs plant) comparison | NEEDS BENCHMARK | `benchmarking/results/` |
| Fig 4 | Cross-study miRNA abundance across OSDR datasets (heatmap) | **have (tool; demo tested)** — from real run | `meta-analysis/plot_meta.py` → `combined/figures/heatmap_top_mirnas.png` |
| Fig 5 | Cross-species miRNA family prevalence (bar) | **have (tool; demo tested)** — from real run | `meta-analysis/plot_meta.py` → `combined/figures/family_prevalence.png` |
| Fig S1 | Cross-engine concordance on real OSDR/SRR data | NEEDS BENCHMARK | `benchmarking/results/concordance_heatmap.png` |
| Fig S2 | (optional) miRNA→target network / example | optional | `target_prediction/` |
| Table 1 | Tools compared (version/role/kingdom/citation) | **have (draft)** | `docs/methods_tables.md` |
| Table 2 | Datasets (synthetic params + SRR/OSDR accessions) | partial | `benchmarking/config.yaml`, `meta-analysis/studies.tsv` |
| Table 3 | Accuracy metrics per engine | NEEDS BENCHMARK | `benchmarking/results/metrics_summary.tsv` |
| Table 4 | OSDR meta-analysis cohort (OSD id, organism, kingdom, layer, factor) | **have** | `meta-analysis/studies.tsv` |

## Figure captions (draft)

*Captions marked **[NEEDS BENCHMARK]** describe the tested figure layout; replace with the
panel from the full run and confirm the numbers. Supplementary captions (S1–S2) follow Fig 5.*

**Figure 1. The nf-core/smrnaseq pipeline deconstructed alongside the recommended
cross-species (DREAM) route.** Each small RNA-seq analysis phase — QC and preprocessing;
alignment, annotation and counts; differential expression; novel-miRNA discovery; reporting —
is shown as a horizontal band. Left (A): the community nf-core/smrnaseq pipeline broken into
its component steps and tools (FastQC, FastP, Bowtie2, Bowtie, SAMtools, edgeR, mirtop,
miRDeep2, miRTrace, MultiQC). Right (B): the recommended DREAM route, in which each phase
collapses to a single integrated tool (mirnaQC → sRNAbench → sRNAde → mirNOVO). The lower band
notes the kingdom-specific novel-miRNA predictors used downstream (miRDeep2 for animals,
miRDeep-P2 for plants). Workflow schematic only; no data are shown. This panel is a
reconstruction drawn from the documented pipeline steps — the original figure was lost with a
since-deleted repository — and is generated reproducibly by `docs/schematics/make_fig1.py`.

**Figure 2. Known-miRNA detection and quantification accuracy on synthetic ground truth.**
(A) Detection performance of each engine (nf-core/smrnaseq, sRNAtoolbox, miRDeep2, miRDeep-P2)
on a synthetic small RNA-seq dataset simulated from miRBase miRNAs at defined abundances:
recall and precision for known miRNAs (grouped bars). (B) Quantification accuracy: for each
recovered known miRNA, detected count against simulated abundance on a log(1 + x) scale with
the identity line; per-engine Spearman and Pearson correlations are reported in Table 3. A
miRNA is scored as detected above a fixed read-count threshold (Methods). Synthetic data
provide the absolute ground truth that real SRA runs cannot. Generated by
`benchmarking/score.py` (harness verified on mock data). **[NEEDS BENCHMARK — panels from the
full four-engine run.]**

**Figure 3. Cross-kingdom comparison of novel-miRNA discovery.** Precision and recall of
de-novo (novel) miRNA prediction for the animal (miRDeep2) and plant (miRDeep-P2) engines,
scored against a set of miRBase miRNAs withheld from the reference supplied to each engine
(genome-anchored synthetic mode, so held-out miRNAs can be rediscovered from hairpin
structure; Methods). A predicted novel miRNA is counted as correct when its mature sequence
matches a held-out truth sequence within ≤2 edits (the miRDeep-P2 "variant" criterion). The
comparison shows whether the recommended workflow generalises across kingdoms without
per-species retuning, and where each kingdom-specific predictor over- or under-calls.
Generated by `benchmarking/score.py`, run separately per kingdom (hsa/animal, ath/plant).
**[NEEDS BENCHMARK.]**

**Figure 4. Cross-study miRNA abundance across NASA OSDR spaceflight datasets.** Heatmap of
the top-variance miRNAs (rows) across all samples (columns) from the combined OSDR
meta-analysis, coloured by log(1 + counts-per-million). The colour bar above the columns
annotates each sample's kingdom: animal *mature*-miRNA counts (mouse HZE-radiation studies
OSD-334–337; astronaut extracellular-vesicle small RNA-seq OSD-483) and plant *precursor*-level
counts (Arabidopsis OSD-208/OSD-437, recovered from standard RNA-seq; see Methods). Mature and
precursor layers are displayed together but never pooled — the block structure reflects
study- and kingdom-specific detection rather than a single shared abundance scale. Counts are
within-sample CPM-normalised only; no cross-study batch correction is applied. Generated by
`meta-analysis/plot_meta.py`. *n* studies = [N], *n* samples = [N].

**Figure 5. Cross-species miRNA family prevalence.** miRNA families ranked by the number of
OSDR studies in which they are detected, after collapsing organism prefixes so homologues
pool across species (e.g. mmu-miR-21-5p + hsa-miR-21-5p → miR-21-5p). Bars are coloured by
kingdom; families detected in more than one kingdom are highlighted, and the number of
distinct member IDs per family is annotated. This surfaces conserved, spaceflight/
radiation-associated miRNA families across the animal cohort while keeping the plant
precursor layer visually distinct. Family collapsing is a prefix-strip heuristic — verify
against miRBase families for load-bearing claims. Generated by `meta-analysis/plot_meta.py`.

**Figure S1. Cross-engine concordance on real data.** Pairwise agreement between the four
engines (nf-core/smrnaseq, sRNAtoolbox, miRDeep2, miRDeep-P2) run on the same real input
(human benchmark runs SRR950892–95 and/or an OSDR small RNA-seq study), shown as a heatmap of
the Jaccard index between the sets of miRNAs each engine detects above the count threshold.
Because it requires no ground truth, this measure applies to real data where the true miRNA
complement is unknown, and complements the synthetic-data accuracy metrics (Fig 2, Table 3) by
showing where engines agree or diverge on identical input — high concordance indicating that
the choice of engine is not critical for a given dataset, low concordance flagging
engine-sensitive calls that warrant scrutiny. Generated by `benchmarking/score.py`
(`concordance_heatmap.png`). **[NEEDS BENCHMARK — from the four-engine run on real data.]**

**Figure S2. Example miRNA–target relationships (optional/illustrative).** Predicted targets
for a selected set of known and/or novel miRNAs, obtained with the sRNAtoolbox consensus target
predictor (miRNAconsTargets, using the plant or animal model as appropriate) and/or the
BLAST-based utility, displayed as a bipartite miRNA→target network (or table). Included to
illustrate the downstream target-prediction step of the framework rather than to make
biological claims: all interactions shown are computational predictions and require
experimental validation. Generated from [`../target_prediction/`](../target_prediction/).
**[Optional — include only if a worked example is run.]**

## Data & code availability

- Code: https://github.com/dr-richard-barker/smallRNAseq-DREAM (live); project site
  https://dr-richard-barker.github.io/smallRNAseq-DREAM/ .
- Synthetic benchmark data + truth tables: **Zenodo DOI TODO** (bundle prepared; see [`../data/MANIFEST.md`](../data/MANIFEST.md)).
- Real data (NASA OSDR): animal small RNA-seq OSD-334, OSD-335, OSD-336, OSD-337 (mouse, HZE
  radiation), OSD-483 (human, astronaut sEV); plant RNA-seq OSD-208, OSD-437 (*Arabidopsis*).
  Benchmark reads: NCBI-SRA SRR950892–95.

## Author contributions / funding / competing interests / acknowledgements

- **TODO** — author list, ORCID, funding (NASA/GeneLab program acknowledgement?), competing
  interests. Retain attribution for TF-Chan-Lab (miRDeep-P2) and upstream tool authors.

---

## Progress log (done)

*Open items live in the "Remaining fill-ins before submission" checklist at the top.*

- [x] All narrative sections drafted in prose: Abstract, §1, §2.1–2.8, §4, §5.
- [x] Captions drafted for Fig 1–5 and Fig S1–S2.
- [x] Fig 1 regenerated as a reproducible remake (`docs/schematics/fig1_pipeline_overview.svg`).
- [x] Fig 4 & Fig 5 wired (`meta-analysis/plot_meta.py`), captions drafted.
- [x] OSDR cohort verified + tabulated (Table 4 / Availability): OSD-334–337, 483, 208, 437.
- [x] OSDR fetch, meta-analysis, and plant precursor modules built and tested on mock data.
- [x] Zenodo bundle assembled with checksums + metadata (deposit pending → DOI).
