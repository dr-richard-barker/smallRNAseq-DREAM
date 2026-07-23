#!/usr/bin/env python3
"""Orchestrate the benchmark: generate synthetic data -> (run engines) -> score.

  python run_benchmark.py --config config.yaml            # full flow
  python run_benchmark.py --config config.yaml --step synth
  python run_benchmark.py --config config.yaml --step score

Running the engines themselves is deliberately left to `run_engines.sh` (or your cluster
scheduler), because each engine needs its own toolchain and compute. This driver handles
the two pure-Python bookends — simulation and scoring — and tells you what to run in
between.
"""

import argparse
import subprocess
import sys

import yaml


def run(cmd):
    print(f"$ {' '.join(cmd)}")
    subprocess.check_call(cmd)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--step", choices=["synth", "engines", "score", "all"], default="all")
    args = ap.parse_args()
    cfg = yaml.safe_load(open(args.config))
    py = sys.executable

    if args.step in ("synth", "all"):
        run([py, "make_synthetic.py", "--config", args.config])

    if args.step in ("engines", "all"):
        s = cfg["synthetic"]
        print("\n--- Run each engine on the synthetic reads, then re-run scoring ---")
        print(f"  reads     : {s['out_dir']}/synthetic_reads.fastq.gz")
        print(f"  reference : {s['out_dir']}/reference_mature.fa  (known miRNAs only)")
        for eng in cfg["engines"]:
            if eng.get("enabled", True):
                print(f"  -> {eng['name']:12s} write output to: {eng['run_dir']}")
        print("  See run_engines.sh for templated commands.\n")

    if args.step in ("score", "all"):
        run([py, "score.py", "--config", args.config])


if __name__ == "__main__":
    main()
