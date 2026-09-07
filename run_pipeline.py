#!/usr/bin/env python3
"""
Master Execution Driver for smallRNAseq-DREAM
Comprehensive microRNA & Non-Coding RNA Audit across NASA OSDR Plant Spaceflight Datasets

This driver executes the four pipeline stages sequentially:
  1. osdr_api_client.py          - Fetch metadata and data assets from NASA OSDR REST API
  2. ncrna_quantifier.py         - Quantify miRNA/ncRNA abundance and differential expression
  3. mirna_target_network.py     - Construct miRNA-mRNA inverse regulatory target network & perform pathway enrichment
  4. build_dream_manuscript_figures.py - Generate 300-DPI publication figures (Figures 1-7)

Author: Dr. Richard Barker / Antigravity AI Research Team
Repository: https://github.com/dr-richard-barker/smallRNAseq-DREAM
License: MIT
"""

import sys
import subprocess
import time
from pathlib import Path

def run_step(script_name, description):
    print("\n" + "="*80)
    print(f"STAGE: {description}")
    print(f"Executing: python {script_name}")
    print("="*80)
    
    start_time = time.time()
    result = subprocess.run([sys.executable, script_name], check=False)
    elapsed = time.time() - start_time
    
    if result.returncode != 0:
        print(f"[-] ERROR: Step {script_name} failed with exit code {result.returncode}.")
        sys.exit(result.returncode)
    
    print(f"[+] Completed {script_name} in {elapsed:.2f} seconds.")

def main():
    print("================================================================================")
    print("      smallRNAseq-DREAM: NASA OSDR Plant Spaceflight microRNA Audit Pipeline")
    print("================================================================================")
    
    workspace = Path(__file__).parent.resolve()
    print(f"Workspace Directory: {workspace}")
    
    # 1. Fetch OSDR Datasets
    run_step("osdr_api_client.py", "1. NASA OSDR API Plant Dataset Mining & Asset Extraction")
    
    # 2. ncRNA Quantification
    run_step("ncrna_quantifier.py", "2. Small RNA & microRNA Spaceflight Abundance Profiling")
    
    # 3. miRNA Target Network & Pathway Enrichment
    run_step("mirna_target_network.py", "3. miRNA-mRNA Target Inverse Regulatory Network & Pathway Enrichment")
    
    # 4. Figure Rendering
    run_step("build_dream_manuscript_figures.py", "4. Publication & Advanced Visualization Figure Generation (Figures 1-7)")
    
    print("\n================================================================================")
    print("[SUCCESS] All pipeline stages completed successfully!")
    print("Generated data outputs:")
    print("  - plant_osdr_inventory.csv & .json")
    print("  - mirna_spaceflight_profiles.csv")
    print("  - mirna_mrna_regulatory_network.csv")
    print("  - target_gene_pathway_enrichment.csv")
    print("Generated figures in figures/:")
    print("  - Figure 1: NASA OSDR Plant Dataset Landscape & Assay Modality Stacked Bar")
    print("  - Figure 2: microRNA Differential Expression & Volcano Plot")
    print("  - Figure 3: miRNA-mRNA Anti-Correlated Regulatory Target Coupling Network")
    print("  - Figure 4: Fisher Exact Test Functional Pathway Enrichment")
    print("  - Figure 5: miRNA -> Target Gene -> Ontology Sankey Flow Diagram")
    print("  - Figure 6: Multimodal Genetic Feature & Abundance Knowledge Graph")
    print("  - Figure 7: Spatial Pathway-Level RNA Abundance Coupling Map")
    print("================================================================================")

if __name__ == "__main__":
    main()
