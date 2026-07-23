#!/usr/bin/env python3
"""Normalise each engine's output into a common format for scoring.

Every parser returns a tuple (counts_df, predictions_df):

  counts_df       columns: mirna_id, sample, count
                  -> quantification of KNOWN miRNAs (matched to the reference)
  predictions_df  columns: predicted_id, sequence, score
                  -> de-novo NOVEL predictions with their mature sequence

The four engines write very different files. The parsers below target each engine's
documented output layout, but exact filenames/columns depend on tool version and DB.
Where that is the case it is flagged `CONFIRM:` — verify against one real run and adjust
the glob/column names. If a parser can't find its file it raises a clear error rather than
returning wrong data.

A universal escape hatch is also provided: if a run_dir contains `counts.tsv` and/or
`predictions.tsv` already in the common format, `parse_generic` uses them directly.
"""

import glob
import os

import pandas as pd


def _norm_id(x):
    """Canonicalise a miRNA id for cross-engine matching."""
    return str(x).strip().lower().lstrip(">")


def _find(run_dir, *patterns):
    for pat in patterns:
        hits = sorted(glob.glob(os.path.join(run_dir, "**", pat), recursive=True))
        if hits:
            return hits[0]
    return None


def read_fasta(path):
    """Minimal FASTA reader -> dict{id: seq} (no Biopython dependency here)."""
    out, cur, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur is not None:
                    out[cur] = "".join(buf)
                cur, buf = line[1:].split()[0], []
            elif line:
                buf.append(line.upper().replace("U", "T"))
    if cur is not None:
        out[cur] = "".join(buf)
    return out


# --------------------------------------------------------------------------------------
# Common-format escape hatch
# --------------------------------------------------------------------------------------
def parse_generic(run_dir):
    counts = pd.DataFrame(columns=["mirna_id", "sample", "count"])
    preds = pd.DataFrame(columns=["predicted_id", "sequence", "score"])
    cpath = os.path.join(run_dir, "counts.tsv")
    ppath = os.path.join(run_dir, "predictions.tsv")
    if os.path.exists(cpath):
        counts = pd.read_csv(cpath, sep="\t")
    if os.path.exists(ppath):
        preds = pd.read_csv(ppath, sep="\t")
    return counts, preds


# --------------------------------------------------------------------------------------
# miRDeep2  (mirDeep2_accuracy / pipelines/mirdeep2_animal)
# --------------------------------------------------------------------------------------
def parse_mirdeep2(run_dir):
    # Known counts: miRDeep2 quantifier writes miRNAs_expressed_all_samples_<ts>.csv (TSV)
    # CONFIRM: columns are '#miRNA', 'read_count', 'precursor', 'total', <per-sample>...
    kpath = _find(run_dir, "miRNAs_expressed_all_samples*.csv", "miRNAs_expressed*.csv")
    counts_rows = []
    if kpath:
        df = pd.read_csv(kpath, sep="\t")
        mcol = df.columns[0]                                   # '#miRNA'
        # per-sample columns come after the fixed header block; fall back to 'read_count'
        sample_cols = [c for c in df.columns if c not in
                       (mcol, "precursor", "total", "read_count")]
        if not sample_cols:
            sample_cols = ["read_count"]
        for _, r in df.iterrows():
            for sc in sample_cols:
                counts_rows.append({"mirna_id": _norm_id(r[mcol]),
                                    "sample": sc, "count": float(r[sc])})
    else:
        raise FileNotFoundError(f"miRDeep2 expression file not found under {run_dir}")

    # Novel predictions: result_<ts>.csv has a 'novel miRNAs predicted by miRDeep2' block.
    # CONFIRM: cols 'provisional id', 'miRDeep2 score', 'consensus mature sequence'.
    preds = []
    rpath = _find(run_dir, "result_*.csv", "result*.csv")
    if rpath:
        raw = open(rpath).read()
        # the novel block is a TSV section; parse lines with a score + a mature sequence
        for line in raw.splitlines():
            parts = line.split("\t")
            if len(parts) < 14:
                continue
            pid, score, mature = parts[0].strip(), parts[1].strip(), parts[13].strip()
            try:
                score = float(score)
            except ValueError:
                continue
            seq = mature.upper().replace("U", "T")
            if set(seq) <= set("ACGTN") and 15 <= len(seq) <= 30:
                preds.append({"predicted_id": pid, "sequence": seq, "score": score})

    return (pd.DataFrame(counts_rows),
            pd.DataFrame(preds, columns=["predicted_id", "sequence", "score"]))


# --------------------------------------------------------------------------------------
# sRNAtoolbox (sRNAbench / sRNAde)  (pipelines/srnatoolbox)
# --------------------------------------------------------------------------------------
def parse_srnatoolbox(run_dir):
    # Known counts: sRNAbench mature_sense.grouped OR the sRNAde adjusted count matrix (*.mat)
    # CONFIRM: 'mature_sense.grouped' is TSV with a 'name' col + count col per sample.
    counts_rows = []
    kpath = _find(run_dir, "*.mat", "mature_sense.grouped", "*grouped*")
    if kpath:
        df = pd.read_csv(kpath, sep="\t")
        idcol = df.columns[0]
        num_cols = [c for c in df.columns[1:] if pd.api.types.is_numeric_dtype(df[c])]
        for _, r in df.iterrows():
            for sc in (num_cols or df.columns[1:2]):
                counts_rows.append({"mirna_id": _norm_id(r[idcol]),
                                    "sample": sc, "count": float(r[sc])})
    else:
        raise FileNotFoundError(f"sRNAtoolbox count matrix not found under {run_dir}")

    # Novel predictions: novel.txt summary + novel_mature.fa sequences
    preds = []
    fa = _find(run_dir, "novel_mature.fa", "*novel*mature*.fa")
    if fa:
        for pid, seq in read_fasta(fa).items():
            preds.append({"predicted_id": pid, "sequence": seq, "score": float("nan")})

    return (pd.DataFrame(counts_rows),
            pd.DataFrame(preds, columns=["predicted_id", "sequence", "score"]))


# --------------------------------------------------------------------------------------
# miRDeep-P2 (plant)  (pipelines/mirdp2_plant)
# --------------------------------------------------------------------------------------
def parse_mirdp2(run_dir):
    # Known counts: combine_htseq_counts.pl -> count_table.txt
    # (last 4 cols are mean/median/variance/CV -> dropped here)
    counts_rows = []
    kpath = _find(run_dir, "count_table.txt", "*count_table*")
    if kpath:
        df = pd.read_csv(kpath, sep="\t")
        idcol = df.columns[0]
        drop = {"mean", "median", "variance", "cv", "coefficient_of_variation"}
        sample_cols = [c for c in df.columns[1:] if c.strip().lower() not in drop]
        for _, r in df.iterrows():
            for sc in sample_cols:
                if pd.api.types.is_number(r[sc]):
                    counts_rows.append({"mirna_id": _norm_id(r[idcol]),
                                        "sample": sc, "count": float(r[sc])})
    else:
        raise FileNotFoundError(f"miRDP2 count_table.txt not found under {run_dir}")

    # Novel predictions: parse_miRDP2_prediction.pl -> miRDP2_mature.fa
    preds = []
    fa = _find(run_dir, "miRDP2_mature.fa", "*miRDP2*mature*.fa")
    if fa:
        for pid, seq in read_fasta(fa).items():
            preds.append({"predicted_id": pid, "sequence": seq, "score": float("nan")})

    return (pd.DataFrame(counts_rows),
            pd.DataFrame(preds, columns=["predicted_id", "sequence", "score"]))


# --------------------------------------------------------------------------------------
# nf-core/smrnaseq  (pipelines/nfcore_smrnaseq)
# --------------------------------------------------------------------------------------
def parse_nfcore(run_dir):
    # Known counts: edgeR / mirtop table. CONFIRM against your results/ layout;
    # common candidates are mirtop.tsv or the mature-miRNA count matrix.
    counts_rows = []
    kpath = _find(run_dir, "mirtop.tsv", "mirna.tsv", "*mature*counts*.tsv", "*.mirna.tsv")
    if kpath:
        df = pd.read_csv(kpath, sep="\t")
        idcol = df.columns[0]
        num_cols = [c for c in df.columns[1:] if pd.api.types.is_numeric_dtype(df[c])]
        for _, r in df.iterrows():
            for sc in (num_cols or df.columns[1:2]):
                counts_rows.append({"mirna_id": _norm_id(r[idcol]),
                                    "sample": sc, "count": float(r[sc])})
    else:
        raise FileNotFoundError(f"nf-core mature-count table not found under {run_dir}")

    # Novel predictions: miRDeep2 module output within the nf-core run
    preds = []
    fa = _find(run_dir, "novel_mature*.fa", "*novel*.fa")
    if fa:
        for pid, seq in read_fasta(fa).items():
            preds.append({"predicted_id": pid, "sequence": seq, "score": float("nan")})

    return (pd.DataFrame(counts_rows),
            pd.DataFrame(preds, columns=["predicted_id", "sequence", "score"]))


PARSERS = {
    "mirdeep2": parse_mirdeep2,
    "srnatoolbox": parse_srnatoolbox,
    "mirdp2": parse_mirdp2,
    "nfcore": parse_nfcore,
    "generic": parse_generic,
}


def parse_engine(name, run_dir):
    """Try the named parser; fall back to the common-format files if present."""
    parser = PARSERS.get(name, parse_generic)
    try:
        return parser(run_dir)
    except FileNotFoundError:
        counts, preds = parse_generic(run_dir)
        if counts.empty and preds.empty:
            raise
        return counts, preds
