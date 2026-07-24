#!/usr/bin/env python3
"""Find and fetch small RNA-seq data from the NASA OSDR (Open Science Data Repository) API.

Feeds real spaceflight small RNA-seq data into the smallRNAseq-DREAM pipelines. Uses only
public OSDR endpoints (validated against the live API):

  search : https://osdr.nasa.gov/osdr/data/search?term=<t>&type=cgene&size=<n>
  files  : https://osdr.nasa.gov/osdr/data/osd/files/<n>
  file   : https://osdr.nasa.gov<remote_url>   (from the files listing)

Commands
--------
  osdr_fetch.py search                 # find small-RNA studies, classify by kingdom
  osdr_fetch.py files OSD-483          # list a study's files (with sizes)
  osdr_fetch.py download OSD-483 -o data/OSD-483 --raw   # download raw FASTQ
  osdr_fetch.py meta OSD-483           # print study metadata summary

No dependencies beyond the standard library.

Note (finding, 2026-07): a survey with this tool found that OSDR currently holds NO plant
small RNA-seq — every plant study is microarray or bulk mRNA-seq. Existing spaceflight
small RNA-seq is animal (mouse-dominated + astronaut sEV miRNA, e.g. OSD-483). See README.md.
"""

import argparse
import json
import os
import sys
import urllib.parse
import urllib.request

API = "https://osdr.nasa.gov"
SMALLRNA_TERMS = ["small RNA", "microRNA", "sRNA-seq", "miRNA", "small RNA sequencing"]

PLANT = ("arabidopsis", "brassica", "oryza", "zea", "glycine", "medicago", "triticum",
         "solanum", "populus", "physcomitr", "ceratopteris", "lolium", "moss", "fern", "plant")
ANIMAL = ("homo sapiens", "mus musculus", "rattus", "drosophila", "danio", "caenorhabditis",
          "oryzias", "helix", "leptopilina")


def _get(url, timeout=60):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return json.load(r)


def classify_kingdom(organism):
    o = (organism or "").lower()
    if any(p in o for p in PLANT):
        return "plant"
    if any(a in o for a in ANIMAL):
        return "animal"
    if any(m in o for m in ("microbiota", "aspergillus", "saccharomyces", "serratia", "staphylococcus")):
        return "microbe"
    return "other"


def search(terms=SMALLRNA_TERMS, size=200):
    """Return {accession: metadata} across the given search terms (deduplicated)."""
    out = {}
    for term in terms:
        url = f"{API}/osdr/data/search?" + urllib.parse.urlencode(
            {"term": term, "type": "cgene", "size": size})
        try:
            d = _get(url)
        except Exception as e:                       # noqa: BLE001
            print(f"  [warn] search '{term}' failed: {e}", file=sys.stderr)
            continue
        for h in d.get("hits", {}).get("hits", []):
            s = h["_source"]
            acc = s.get("Accession")
            if acc and acc not in out:
                out[acc] = {
                    "accession": acc,
                    "organism": s.get("organism"),
                    "kingdom": classify_kingdom(s.get("organism")),
                    "title": (s.get("Study Title") or "").strip(),
                    "assay": s.get("Study Assay Technology Type"),
                    "gse": s.get("Data Source Accession"),
                }
    return out


def list_files(osd_id):
    n = str(osd_id).split("-")[-1]
    d = _get(f"{API}/osdr/data/osd/files/{n}")
    study = d.get("studies", {}).get(f"OSD-{n}", {})
    return study.get("study_files", [])


def is_fastq(name):
    return name.lower().endswith((".fastq.gz", ".fq.gz", ".fastq", ".fq"))


def download(osd_id, out_dir, raw_only=False, pattern=None, dry_run=False, limit=None):
    files = list_files(osd_id)
    picks = []
    for f in files:
        name = f["file_name"]
        if raw_only and not is_fastq(name):
            continue
        if raw_only and "raw" not in name.lower() and not is_fastq(name):
            continue
        if pattern and pattern.lower() not in name.lower():
            continue
        picks.append(f)
    if raw_only:
        picks = [f for f in picks if is_fastq(f["file_name"])]
    if limit:
        picks = picks[:limit]
    total = sum(f["file_size"] for f in picks)
    print(f"{osd_id}: {len(picks)} file(s), {total/1e6:.1f} MB total")
    if dry_run:
        for f in picks:
            print(f"  [dry] {f['file_name']}  ({f['file_size']/1e6:.1f} MB)")
        return
    os.makedirs(out_dir, exist_ok=True)
    for f in picks:
        dest = os.path.join(out_dir, f["file_name"])
        url = API + f["remote_url"]
        print(f"  -> {f['file_name']}  ({f['file_size']/1e6:.1f} MB)")
        urllib.request.urlretrieve(url, dest)
    print(f"done: {out_dir}")


def cmd_search(args):
    studies = search(size=args.size)
    from collections import Counter
    by_king = Counter(v["kingdom"] for v in studies.values())
    print(f"# {len(studies)} studies returned by OSDR small-RNA search terms")
    print(f"# by kingdom: " + ", ".join(f"{k}={n}" for k, n in by_king.most_common()))
    print("accession\tkingdom\torganism\tassay\tgse\ttitle")
    for v in sorted(studies.values(), key=lambda x: (x["kingdom"], x["accession"])):
        if args.kingdom and v["kingdom"] != args.kingdom:
            continue
        print(f"{v['accession']}\t{v['kingdom']}\t{v['organism']}\t{v['assay']}\t{v['gse']}\t{v['title'][:70]}")


def cmd_files(args):
    for f in list_files(args.osd):
        print(f"{f['file_size']/1e6:8.1f} MB  {f['category']:32.32s}  {f['file_name']}")


def cmd_meta(args):
    studies = search(size=200)
    v = studies.get(args.osd if args.osd.startswith("OSD-") else f"OSD-{args.osd}")
    if not v:
        print(f"{args.osd} not found among small-RNA search hits."); return
    for k, val in v.items():
        print(f"{k:10s}: {val}")


def cmd_download(args):
    download(args.osd, args.out, raw_only=args.raw, pattern=args.pattern,
             dry_run=args.dry_run, limit=args.limit)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("search", help="find + classify small-RNA studies")
    s.add_argument("--size", type=int, default=200)
    s.add_argument("--kingdom", choices=["plant", "animal", "microbe", "other"])
    s.set_defaults(func=cmd_search)

    f = sub.add_parser("files", help="list a study's files")
    f.add_argument("osd")
    f.set_defaults(func=cmd_files)

    m = sub.add_parser("meta", help="study metadata summary")
    m.add_argument("osd")
    m.set_defaults(func=cmd_meta)

    d = sub.add_parser("download", help="download study files")
    d.add_argument("osd")
    d.add_argument("-o", "--out", default="data")
    d.add_argument("--raw", action="store_true", help="raw FASTQ only")
    d.add_argument("--pattern", help="only files containing this substring")
    d.add_argument("--limit", type=int, help="cap number of files")
    d.add_argument("--dry-run", action="store_true")
    d.set_defaults(func=cmd_download)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
