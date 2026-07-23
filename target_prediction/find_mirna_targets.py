#!/usr/bin/env python3
"""Search for microRNA targets by BLAST.

Consolidated from smallRNA_targets: python_Find_microRNA_Target_1.0 and
BioPython_script2_BING (two near-identical drafts merged into one script).

Reads a CSV of miRNA names + sequences, runs blastn against NCBI 'nt' over the network,
and writes one combined CSV of candidate hit sites.

Input CSV columns : name, sequence   (header row required)
Output CSV columns: name, sequence, hit_id, hit_def, hit_start, hit_end, e_value

Note: NCBIWWW.qblast submits over the network and is slow / rate-limited. For anything
beyond a handful of sequences, run blastn locally against a downloaded database instead.
"""

import argparse
import csv
import tempfile

import pandas as pd
from Bio.Blast import NCBIWWW, NCBIXML


def blast_sequence(sequence, database="nt", program="blastn"):
    """Run one BLAST query and yield parsed records."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa") as tmp:
        tmp.write(f">query\n{sequence}\n")
        tmp.flush()
        handle = NCBIWWW.qblast(program, database, tmp.name)
        yield from NCBIXML.parse(handle)


def main():
    ap = argparse.ArgumentParser(description="BLAST miRNA sequences to find target sites.")
    ap.add_argument("-i", "--input", default="input.csv", help="CSV with name,sequence columns")
    ap.add_argument("-o", "--output", default="output.csv", help="combined results CSV")
    ap.add_argument("-d", "--database", default="nt", help="BLAST database (default: nt)")
    args = ap.parse_args()

    data = pd.read_csv(args.input)
    rows = []
    for _, row in data.iterrows():
        name, sequence = row["name"], row["sequence"]
        print(f"Processing {name}...")
        for record in blast_sequence(sequence, database=args.database):
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    rows.append({
                        "name": name,
                        "sequence": sequence,
                        "hit_id": alignment.hit_id,
                        "hit_def": alignment.hit_def,
                        "hit_start": hsp.sbjct_start,
                        "hit_end": hsp.sbjct_end,
                        "e_value": hsp.expect,
                    })

    pd.DataFrame(rows).to_csv(args.output, index=False, quoting=csv.QUOTE_MINIMAL)
    print(f"Results saved to {args.output}")


if __name__ == "__main__":
    main()
