#!/usr/bin/env python3
"""Generate a synthetic small RNA-seq dataset with known ground truth.

Purpose
-------
The manuscript needs quantified accuracy numbers. Real SRA runs have no absolute truth,
so this script simulates reads from known miRBase miRNAs at defined abundances and writes
the truth table alongside. A subset of miRNAs is HELD OUT of the reference the engines are
given — those are the ground-truth "novel" miRNAs used to score de-novo prediction.

Two modes
---------
mature  (self-contained): reads are drawn from miRBase mature sequences directly. Good for
        known-miRNA recovery and quantification accuracy. Novel-prediction scoring is
        weaker here because there is no genomic hairpin context for the folding-based
        predictors (miRDeep2 / miRDP2) — use `genome` mode for a fair novel test.
genome  (needs genome.fa + miRBase GFF3): reads are drawn from the genomic precursor loci,
        so held-out miRNAs can be rediscovered from hairpin structure. This is the rigorous
        setup for novel-prediction precision/recall.

Outputs (into <out_dir>/)
-------------------------
  synthetic_reads.fastq.gz   simulated reads (feed this to every engine)
  reference_mature.fa        ONLY the `known` miRNAs — give this to the engines as the reference
  truth_counts.tsv           mirna_id, sequence, true_abundance, category(known|novel)
  novel_truth.fa             the held-out miRNAs (ground truth for novel detection)

Usage
-----
  python make_synthetic.py --config config.yaml
"""

import argparse
import gzip
import os
import sys

import numpy as np
import yaml
from Bio import SeqIO


def u2t(seq):
    return str(seq).upper().replace("U", "T")


def load_mature(mature_fa, species_prefix):
    """Return {id: dna_seq} for one species from a miRBase mature.fa."""
    out = {}
    for rec in SeqIO.parse(mature_fa, "fasta"):
        if species_prefix and not rec.id.lower().startswith(species_prefix.lower()):
            continue
        seq = u2t(rec.seq)
        if 18 <= len(seq) <= 26:            # sane mature-miRNA length window
            out[rec.id] = seq
    if not out:
        sys.exit(f"No sequences with prefix '{species_prefix}' in {mature_fa}")
    return out


def mutate(seq, error_rate, rng):
    """Apply per-base substitutions at error_rate."""
    bases = "ACGT"
    s = list(seq)
    for i, b in enumerate(s):
        if rng.random() < error_rate:
            s[i] = rng.choice([x for x in bases if x != b])
    return "".join(s)


def make_isomir(seq, rng):
    """Trim or extend 1–2 nt at the 5'/3' end to mimic isomiRs."""
    end = rng.choice(["5", "3"])
    op = rng.choice(["trim", "add"])
    n = rng.integers(1, 3)  # 1 or 2
    if op == "trim":
        return seq[n:] if end == "5" else seq[:-n]
    tail = "".join(rng.choice(list("ACGT")) for _ in range(n))
    return tail + seq if end == "5" else seq + tail


def build_read(small_rna, adapter, read_length, error_rate, isomir_rate, rng):
    frag = small_rna
    if rng.random() < isomir_rate:
        frag = make_isomir(frag, rng)
    read = (frag + adapter)[:read_length]
    if len(read) < read_length:                      # pad short reads with adapter-ish noise
        read += "".join(rng.choice(list("ACGT")) for _ in range(read_length - len(read)))
    return mutate(read, error_rate, rng)


def assign_abundances(ids, depth, sigma, rng):
    """Log-normal abundance per miRNA, scaled so the total equals `depth`."""
    weights = rng.lognormal(mean=0.0, sigma=sigma, size=len(ids))
    weights = weights / weights.sum()
    counts = np.round(weights * depth).astype(int)
    counts[counts < 1] = 1
    return dict(zip(ids, counts))


def write_fastq(path, reads, rng):
    with gzip.open(path, "wt") as fh:
        for i, seq in enumerate(reads):
            qual = "I" * len(seq)                     # flat high quality; good enough for benchmarking
            fh.write(f"@synth_read_{i}\n{seq}\n+\n{qual}\n")


def gather_genome_loci(genome_fa, gff3, species_prefix):
    """genome mode: return {precursor_id: genomic_dna} from a miRBase GFF3."""
    genome = SeqIO.to_dict(SeqIO.parse(genome_fa, "fasta"))
    loci = {}
    with open(gff3) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "miRNA_primary_transcript":
                continue
            chrom, start, end, strand, attrs = f[0], int(f[3]), int(f[4]), f[6], f[8]
            name = dict(kv.split("=", 1) for kv in attrs.split(";") if "=" in kv).get("Name", "")
            if species_prefix and not name.lower().startswith(species_prefix.lower()):
                continue
            if chrom not in genome:
                continue
            sub = genome[chrom].seq[start - 1:end]
            if strand == "-":
                sub = sub.reverse_complement()
            loci[name] = u2t(sub)
    if not loci:
        sys.exit(f"No precursor loci matched prefix '{species_prefix}' in {gff3}")
    return loci


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="config.yaml")
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    s = cfg["synthetic"]
    rng = np.random.default_rng(cfg.get("seed", 42))
    os.makedirs(s["out_dir"], exist_ok=True)

    # 1. Load the miRNA pool ------------------------------------------------------------
    mature = load_mature(s["mature_fa"], s["species_prefix"])
    ids = list(mature)
    rng.shuffle(ids)
    n_known, n_novel = s["n_known"], s["n_novel"]
    if len(ids) < n_known + n_novel:
        sys.exit(f"Only {len(ids)} miRNAs available; need n_known+n_novel={n_known+n_novel}")
    known_ids = ids[:n_known]
    novel_ids = ids[n_known:n_known + n_novel]

    # 2. Abundances for every simulated miRNA (known + novel) ---------------------------
    mirna_depth = int(s["depth_total"] * (1 - s["decoy_fraction"]))
    abundances = assign_abundances(known_ids + novel_ids, mirna_depth,
                                   s["abundance_lognorm_sigma"], rng)

    # 3. Source sequence for read simulation -------------------------------------------
    if s["mode"] == "genome":
        loci = gather_genome_loci(s["genome_fa"], s["mirbase_gff3"], s["species_prefix"])
        # reads span the precursor around the mature; fall back to mature if no locus
        def source_seq(mid):
            return loci.get(mid, mature[mid])
    else:
        def source_seq(mid):
            return mature[mid]

    # 4. Simulate reads -----------------------------------------------------------------
    reads = []
    for mid, n in abundances.items():
        base = mature[mid]                # the biologically real mature sequence for that read
        for _ in range(n):
            reads.append(build_read(base, s["adapter"], s["read_length"],
                                    s["error_rate"], s["isomir_rate"], rng))

    # 5. Decoy noise reads (non-miRNA) --------------------------------------------------
    n_decoy = s["depth_total"] - len(reads)
    for _ in range(max(0, n_decoy)):
        length = int(rng.integers(18, 27))
        frag = "".join(rng.choice(list("ACGT")) for _ in range(length))
        reads.append(build_read(frag, s["adapter"], s["read_length"],
                                s["error_rate"], 0.0, rng))
    rng.shuffle(reads)

    # 6. Write outputs ------------------------------------------------------------------
    out = s["out_dir"]
    write_fastq(os.path.join(out, "synthetic_reads.fastq.gz"), reads, rng)

    with open(os.path.join(out, "reference_mature.fa"), "w") as fh:
        for mid in known_ids:
            fh.write(f">{mid}\n{mature[mid]}\n")

    with open(os.path.join(out, "novel_truth.fa"), "w") as fh:
        for mid in novel_ids:
            fh.write(f">{mid}\n{mature[mid]}\n")

    with open(os.path.join(out, "truth_counts.tsv"), "w") as fh:
        fh.write("mirna_id\tsequence\ttrue_abundance\tcategory\n")
        for mid in known_ids:
            fh.write(f"{mid}\t{mature[mid]}\t{abundances[mid]}\tknown\n")
        for mid in novel_ids:
            fh.write(f"{mid}\t{mature[mid]}\t{abundances[mid]}\tnovel\n")

    print(f"Wrote {len(reads):,} reads to {out}/synthetic_reads.fastq.gz")
    print(f"  known miRNAs (in reference): {len(known_ids)}")
    print(f"  novel miRNAs (held out):     {len(novel_ids)}")
    print(f"  mode: {s['mode']}")
    if s["mode"] == "mature":
        print("  NOTE: novel-prediction scores are conservative in 'mature' mode "
              "(no genomic hairpin context). Use mode: genome for the rigorous novel test.")


if __name__ == "__main__":
    main()
