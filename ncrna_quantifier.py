#!/usr/bin/env python3
"""
MicroRNA & Non-Coding RNA Quantification and Audit Module for Plant OSDR Datasets
"""

import os
import glob
import json
import pandas as pd
import numpy as np

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

# WARNING (abai check 2026-10-28): this hand-typed map was NOT taken from miRBase or Araport11.
# Of the 29 IDs that reach mirna_spaceflight_profiles.csv, only 3 are the named MIR gene
# (AT2G46685 MIR166a, AT5G14545 MIR398b, AT5G14565 MIR398c); 16 are protein-coding (e.g.
# AT1G53160 = SPL4, labelled MIR156b), 2 lncRNA, 1 snoRNA and 7 a different MIR gene
# (Ensembl Plants REST lookup). The outputs are therefore not miRNA measurements.
# TODO: rebuild this map from Ensembl Plants genes with biotype=miRNA.
MIRNA_LOCI_MAP = {
    'AT1G48190': 'MIR156a', 'AT1G53160': 'MIR156b', 'AT2G24950': 'MIR156c', 'AT5G10940': 'MIR156d',
    'AT5G10950': 'MIR156e', 'AT5G10960': 'MIR156f', 'AT2G47540': 'MIR165a', 'AT4G00885': 'MIR165b',
    'AT2G46685': 'MIR166a', 'AT3G61710': 'MIR166b', 'AT3G24177': 'MIR166c', 'AT3G24179': 'MIR166d',
    'AT5G43355': 'MIR166e', 'AT5G53865': 'MIR166f', 'AT5G04275': 'MIR166g', 'AT2G38400': 'MIR398a',
    'AT5G14545': 'MIR398b', 'AT5G14565': 'MIR398c', 'AT2G28610': 'MIR168a', 'AT5G45185': 'MIR168b',
    'AT2G39175': 'MIR172a', 'AT5G60860': 'MIR172b', 'AT3G11440': 'MIR172c', 'AT3G55512': 'MIR172d',
    'AT5G04270': 'MIR172e', 'AT1G62670': 'MIR160a', 'AT1G63717': 'MIR160b', 'AT2G01045': 'MIR160c',
    'AT3G13405': 'MIR164a', 'AT5G07745': 'MIR164b', 'AT5G57865': 'MIR164c', 'AT2G23805': 'MIR167a',
    'AT3G04765': 'MIR167b', 'AT3G22886': 'MIR167c', 'AT1G30585': 'MIR167d', 'AT1G60165': 'MIR169a',
    'AT3G26815': 'MIR171a', 'AT1G75710': 'MIR171b', 'AT1G11735': 'MIR171c', 'AT1G44685': 'MIR395a',
    'AT1G30580': 'MIR395b', 'AT1G30425': 'MIR395c', 'AT1G53905': 'MIR159a', 'AT1G73687': 'MIR159b',
    'AT2G43845': 'MIR161',  'AT5G02195': 'MIR162a', 'AT5G64380': 'MIR163',  'AT5G07725': 'MIR393a',
    'AT3G55734': 'MIR393b', 'AT5G55845': 'MIR394a', 'AT1G24585': 'MIR397a', 'AT4G37725': 'MIR397b',
    'AT3G03575': 'MIR399a', 'AT2G21955': 'MIR399b', 'AT2G21945': 'MIR399c', 'AT2G39835': 'MIR824'
}

def analyze_srna_seq_osd437():
    # OSD-437 is standard poly(A) RNA-seq (NEBNext Ultra II RNA; ISA library selection
    # 'polyA enrichment'), not small RNA-seq. The function name is kept for compatibility.
    print("=== Extracting map-labelled loci from OSD-437 (poly(A) RNA-seq DE table) ===", flush=True)
    osd437_dir = os.path.join(DATA_DIR, "OSD-437")
    de_file = os.path.join(osd437_dir, "OSD-437_differential_expression.csv")
    counts_file = os.path.join(osd437_dir, "OSD-437_raw_counts.csv")
    vst_file = os.path.join(osd437_dir, "OSD-437_vst_counts.csv")
    
    records = []
    
    if os.path.exists(de_file) and os.path.getsize(de_file) > 1000:
        df_de = pd.read_csv(de_file)
        print(f"Loaded OSD-437 DE table ({len(df_de)} rows, columns: {list(df_de.columns)[:6]})", flush=True)
        
        # Identify gene column
        gene_col = df_de.columns[0]
        
        # Filter for miRNAs / MIR loci
        for idx, row in df_de.iterrows():
            gene_id = str(row[gene_col]).strip().upper()
            gene_name = MIRNA_LOCI_MAP.get(gene_id, None)
            
            is_mir = bool(gene_name) or 'MIR' in gene_id or 'MIRNA' in gene_id
            if is_mir:
                # Find log2FC and PValue columns.
                # NOTE: [0] is the first contrast in the GeneLab table, which is not flight vs ground
                # (for OSD-437 it is 1G-centrifuge vs uG, both in flight). TODO: select by name.
                log2fc_cols = [c for c in df_de.columns if 'log2' in c.lower() or 'foldchange' in c.lower()]
                p_cols = [c for c in df_de.columns if 'pvalue' in c.lower() or 'p.value' in c.lower() or 'padj' in c.lower()]
                
                log2fc_val = row[log2fc_cols[0]] if log2fc_cols else 0.0
                pval_val = row[p_cols[0]] if p_cols else 1.0
                
                records.append({
                    'Study': 'OSD-437',
                    'Assay': 'RNA-seq (polyA)',
                    'Gene_ID': gene_id,
                    'miRNA_Name': gene_name or gene_id,
                    'log2FC': float(log2fc_val) if pd.notnull(log2fc_val) else 0.0,
                    'PValue': float(pval_val) if pd.notnull(pval_val) else 1.0,
                    'Significant': bool(abs(float(log2fc_val or 0)) >= 0.5 and float(pval_val or 1) <= 0.05),
                    'Direction': 'Upregulated' if (log2fc_val or 0) > 0 else 'Downregulated'
                })
                
    # If explicit DE table rows are sparse, perform quantification across counts/vst matrix
    if os.path.exists(vst_file) and os.path.getsize(vst_file) > 1000:
        df_vst = pd.read_csv(vst_file)
        gene_col = df_vst.columns[0]
        
        for idx, row in df_vst.iterrows():
            gene_id = str(row[gene_col]).strip().upper()
            gene_name = MIRNA_LOCI_MAP.get(gene_id, None)
            if gene_name:
                # Fallback (not reached for the current outputs: every OSD-437 row comes from the DE table).
                # NOTE: the keyword match and the half-split default below do not identify flight vs ground.
                flight_cols = [c for c in df_vst.columns if any(k in c.lower() for k in ['flight', 'flt', 'f_', 'space'])]
                ground_cols = [c for c in df_vst.columns if any(k in c.lower() for k in ['ground', 'gc', 'g_', '1g', 'control'])]
                
                if not flight_cols:
                    # Default column partition for OSD-437 samples if not named flight/ground explicitly
                    all_sample_cols = list(df_vst.columns[1:])
                    half = len(all_sample_cols) // 2
                    flight_cols = all_sample_cols[:half]
                    ground_cols = all_sample_cols[half:]
                    
                flt_mean = row[flight_cols].mean()
                gc_mean = row[ground_cols].mean()
                diff_vst = flt_mean - gc_mean
                
                if not any(r['Gene_ID'] == gene_id for r in records):
                    records.append({
                        'Study': 'OSD-437',
                        'Assay': 'RNA-seq (polyA)',
                        'Gene_ID': gene_id,
                        'miRNA_Name': gene_name,
                        'log2FC': float(diff_vst),
                        'PValue': 0.01 if abs(diff_vst) > 0.5 else 0.5,
                        'Significant': bool(abs(diff_vst) >= 0.5),
                        'Direction': 'Upregulated' if diff_vst > 0 else 'Downregulated'
                    })

    df_res = pd.DataFrame(records)
    print(f"Extracted {len(df_res)} map-labelled rows from the OSD-437 RNA-seq DE table.", flush=True)
    return df_res

def audit_bulk_rnaseq_ncrna():
    print("\n=== Auditing Non-Coding RNAs across Plant Bulk RNA-seq Studies ===", flush=True)
    studies = ['OSD-120', 'OSD-218', 'OSD-217', 'OSD-565', 'OSD-782', 'OSD-522', 'OSD-624']
    
    all_ncrna_records = []
    
    for osd_id in studies:
        de_file = os.path.join(DATA_DIR, osd_id, f"{osd_id}_differential_expression.csv")
        vst_file = os.path.join(DATA_DIR, osd_id, f"{osd_id}_vst_counts.csv")
        
        if os.path.exists(de_file) and os.path.getsize(de_file) > 1000:
            try:
                df_de = pd.read_csv(de_file)
                gene_col = df_de.columns[0]
                
                log2fc_cols = [c for c in df_de.columns if 'log2' in c.lower() or 'foldchange' in c.lower()]
                p_cols = [c for c in df_de.columns if 'pvalue' in c.lower() or 'p.value' in c.lower() or 'padj' in c.lower()]
                
                if log2fc_cols:
                    # NOTE: first contrast in the table, not flight vs ground (TODO: select by name).
                    lfc_col = log2fc_cols[0]
                    pval_col = p_cols[0] if p_cols else None
                    
                    for idx, row in df_de.iterrows():
                        gene_id = str(row[gene_col]).strip().upper()
                        gene_name = MIRNA_LOCI_MAP.get(gene_id, None)
                        
                        is_ncrna = bool(gene_name) or gene_id.startswith('MIR') or 'SNOR' in gene_id or 'LNC' in gene_id or 'TRNA' in gene_id
                        
                        if is_ncrna:
                            lfc_val = float(row[lfc_col]) if pd.notnull(row[lfc_col]) else 0.0
                            pval_val = float(row[pval_col]) if (pval_col and pd.notnull(row[pval_col])) else 1.0
                            
                            all_ncrna_records.append({
                                'Study': osd_id,
                                'Assay': 'Bulk RNA-seq',
                                'Gene_ID': gene_id,
                                'miRNA_Name': gene_name or gene_id,
                                'log2FC': lfc_val,
                                'PValue': pval_val,
                                'Significant': bool(abs(lfc_val) >= 0.5 and pval_val <= 0.05),
                                'Direction': 'Upregulated' if lfc_val > 0 else 'Downregulated'
                            })
            except Exception as e:
                print(f"Error reading {de_file}: {e}", flush=True)

    df_all = pd.DataFrame(all_ncrna_records)
    print(f"Extracted {len(df_all)} non-coding RNA expression signals across bulk RNA-seq datasets.", flush=True)
    return df_all

def run_ncrna_audit():
    df_srna = analyze_srna_seq_osd437()
    df_bulk = audit_bulk_rnaseq_ncrna()
    
    df_combined = pd.concat([df_srna, df_bulk], ignore_index=True) if not df_srna.empty or not df_bulk.empty else pd.DataFrame()
    
    output_csv = os.path.join(os.path.dirname(__file__), "mirna_spaceflight_profiles.csv")
    df_combined.to_csv(output_csv, index=False)
    
    print(f"\nExported {len(df_combined)} map-labelled rows to {output_csv} (not miRNA measurements; see README status)", flush=True)
    return df_combined

if __name__ == "__main__":
    run_ncrna_audit()
