#!/usr/bin/env python3
"""
miRNA-target pairing over GeneLab DE tables, plus an ILLUSTRATIVE per-pathway Fisher test.
See the README status note: the miRNA-side values come from mislabelled loci and the
log2FC columns are each table's first contrast, so these outputs are not interpretable yet.
"""

import os
import glob
import json
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

# Hand-curated, well-known plant miRNA -> target pairs, typed from the literature.
# No target database (miRBase / TargetFinder / psRNATarget) was queried. TODO: cite a source per pair.
PLANT_MIRNA_TARGETS = [
    # miR156 family -> SPL transcription factors (Development, Gravitropism, Flowering)
    {'miRNA': 'miR156', 'miRNA_Gene': 'MIR156a/b/c/d/e/f', 'Target_Gene_ID': 'AT1G27370', 'Target_Symbol': 'SPL2',  'Pathway': 'Phase Change & Vegetative Development', 'Function': 'Squamosa promoter binding protein-like 2, root & shoot morphogenesis'},
    {'miRNA': 'miR156', 'miRNA_Gene': 'MIR156a/b/c/d/e/f', 'Target_Gene_ID': 'AT2G33810', 'Target_Symbol': 'SPL3',  'Pathway': 'Phase Change & Vegetative Development', 'Function': 'Floral transition and shoot apical meristem response'},
    {'miRNA': 'miR156', 'miRNA_Gene': 'MIR156a/b/c/d/e/f', 'Target_Gene_ID': 'AT2G42200', 'Target_Symbol': 'SPL9',  'Pathway': 'Phase Change & Vegetative Development', 'Function': 'Epidermal cell fate and gravitropic branching'},
    {'miRNA': 'miR156', 'miRNA_Gene': 'MIR156a/b/c/d/e/f', 'Target_Gene_ID': 'AT1G27360', 'Target_Symbol': 'SPL10', 'Pathway': 'Phase Change & Vegetative Development', 'Function': 'Lateral organ development and cell division'},
    {'miRNA': 'miR156', 'miRNA_Gene': 'MIR156a/b/c/d/e/f', 'Target_Gene_ID': 'AT5G50570', 'Target_Symbol': 'SPL13', 'Pathway': 'Phase Change & Vegetative Development', 'Function': 'Seedling development and vegetative phase change'},
    
    # miR165/166 family -> HD-ZIP III TFs (Vascular patterning, polarity, xylem morphogenesis)
    {'miRNA': 'miR166', 'miRNA_Gene': 'MIR166a/b/c/d/e/f/g', 'Target_Gene_ID': 'AT2G34710', 'Target_Symbol': 'PHB',   'Pathway': 'Vascular Patterning & Cell Wall Remodeling', 'Function': 'PHABULOSA HD-ZIP III TF, vascular polarity and SAM development'},
    {'miRNA': 'miR166', 'miRNA_Gene': 'MIR166a/b/c/d/e/f/g', 'Target_Gene_ID': 'AT1G30490', 'Target_Symbol': 'PHV',   'Pathway': 'Vascular Patterning & Cell Wall Remodeling', 'Function': 'PHAVOLUTA HD-ZIP III TF, leaf adaxialization and xylem fate'},
    {'miRNA': 'miR166', 'miRNA_Gene': 'MIR166a/b/c/d/e/f/g', 'Target_Gene_ID': 'AT5G60690', 'Target_Symbol': 'REV',   'Pathway': 'Vascular Patterning & Cell Wall Remodeling', 'Function': 'REVOLUTA HD-ZIP III TF, secondary cell wall & gravitropic vascular patterning'},
    {'miRNA': 'miR166', 'miRNA_Gene': 'MIR166a/b/c/d/e/f/g', 'Target_Gene_ID': 'AT3G17710', 'Target_Symbol': 'ATHB8', 'Pathway': 'Vascular Patterning & Cell Wall Remodeling', 'Function': 'Arabidopsis thaliana homeobox 8, procambium development'},

    # miR398 family -> Superoxide Dismutases (ROS Scavenging & Oxidative Stress)
    {'miRNA': 'miR398', 'miRNA_Gene': 'MIR398a/b/c', 'Target_Gene_ID': 'AT1G08830', 'Target_Symbol': 'CSD1', 'Pathway': 'Reactive Oxygen Species (ROS) & Oxidative Stress Defense', 'Function': 'Cytosolic Cu/Zn superoxide dismutase 1, oxidative stress defense'},
    {'miRNA': 'miR398', 'miRNA_Gene': 'MIR398a/b/c', 'Target_Gene_ID': 'AT2G28190', 'Target_Symbol': 'CSD2', 'Pathway': 'Reactive Oxygen Species (ROS) & Oxidative Stress Defense', 'Function': 'Chloroplastic Cu/Zn superoxide dismutase 2, ROS detoxification'},
    {'miRNA': 'miR398', 'miRNA_Gene': 'MIR398a/b/c', 'Target_Gene_ID': 'AT1G12520', 'Target_Symbol': 'CCS1', 'Pathway': 'Reactive Oxygen Species (ROS) & Oxidative Stress Defense', 'Function': 'Copper chaperone for superoxide dismutase 1'},

    # miR168 family -> ARGONAUTE 1 (RNA Silencing Feedback Control)
    {'miRNA': 'miR168', 'miRNA_Gene': 'MIR168a/b', 'Target_Gene_ID': 'AT1G48410', 'Target_Symbol': 'AGO1', 'Pathway': 'Small RNA Silencing Machinery Regulation', 'Function': 'ARGONAUTE 1, catalytic core of RISC complex and sRNA silencing'},

    # miR172 family -> APETALA2 & Target of EAC 1-3 (Floral and meristem identity)
    {'miRNA': 'miR172', 'miRNA_Gene': 'MIR172a/b/c/d/e', 'Target_Gene_ID': 'AT4G36920', 'Target_Symbol': 'AP2',  'Pathway': 'Phase Change & Vegetative Development', 'Function': 'APETALA2 TF, floral organ identity and seed size'},
    {'miRNA': 'miR172', 'miRNA_Gene': 'MIR172a/b/c/d/e', 'Target_Gene_ID': 'AT2G28550', 'Target_Symbol': 'TOE1', 'Pathway': 'Phase Change & Vegetative Development', 'Function': 'TARGET OF EAT 1, flowering repressor'},

    # miR160 family -> Auxin Response Factors (Root Cap & Gravitropism)
    {'miRNA': 'miR160', 'miRNA_Gene': 'MIR160a/b/c', 'Target_Gene_ID': 'AT2G28350', 'Target_Symbol': 'ARF10', 'Pathway': 'Auxin Signaling & Root Gravitropism', 'Function': 'Auxin Response Factor 10, root cap cell fate and seed germination'},
    {'miRNA': 'miR160', 'miRNA_Gene': 'MIR160a/b/c', 'Target_Gene_ID': 'AT4G30080', 'Target_Symbol': 'ARF16', 'Pathway': 'Auxin Signaling & Root Gravitropism', 'Function': 'Auxin Response Factor 16, root gravitropism & stem cell niche'},
    {'miRNA': 'miR160', 'miRNA_Gene': 'MIR160a/b/c', 'Target_Gene_ID': 'AT1G77850', 'Target_Symbol': 'ARF17', 'Pathway': 'Auxin Signaling & Root Gravitropism', 'Function': 'Auxin Response Factor 17, adventitious rooting and leaf symmetry'},

    # miR164 family -> NAC domain TFs (Lateral root development & organ boundary)
    {'miRNA': 'miR164', 'miRNA_Gene': 'MIR164a/b/c', 'Target_Gene_ID': 'AT1G56010', 'Target_Symbol': 'NAC1',  'Pathway': 'Auxin Signaling & Root Gravitropism', 'Function': 'NAC domain TF 1, lateral root initiation'},
    {'miRNA': 'miR164', 'miRNA_Gene': 'MIR164a/b/c', 'Target_Gene_ID': 'AT3G15170', 'Target_Symbol': 'CUC1',  'Pathway': 'Auxin Signaling & Root Gravitropism', 'Function': 'CUP-SHAPED COTYLEDON 1, organ boundary establishment'},
    {'miRNA': 'miR164', 'miRNA_Gene': 'MIR164a/b/c', 'Target_Gene_ID': 'AT5G53950', 'Target_Symbol': 'CUC2',  'Pathway': 'Auxin Signaling & Root Gravitropism', 'Function': 'CUP-SHAPED COTYLEDON 2, serration & organ separation'},

    # miR393 family -> Auxin Receptors (Auxin perception and immunity)
    {'miRNA': 'miR393', 'miRNA_Gene': 'MIR393a/b', 'Target_Gene_ID': 'AT3G62980', 'Target_Symbol': 'TIR1',  'Pathway': 'Auxin Signaling & Root Gravitropism', 'Function': 'TRANSPORT INHIBITOR RESPONSE 1, primary auxin receptor F-box protein'},
    {'miRNA': 'miR393', 'miRNA_Gene': 'MIR393a/b', 'Target_Gene_ID': 'AT4G03190', 'Target_Symbol': 'AFB1',  'Pathway': 'Auxin Signaling & Root Gravitropism', 'Function': 'Auxin signaling F-box protein 1'},
    {'miRNA': 'miR393', 'miRNA_Gene': 'MIR393a/b', 'Target_Gene_ID': 'AT3G26860', 'Target_Symbol': 'AFB2',  'Pathway': 'Auxin Signaling & Root Gravitropism', 'Function': 'Auxin signaling F-box protein 2'}
]

def build_regulatory_network():
    print("=== Building miRNA-mRNA Regulatory Interactions across Plant OSDR Datasets ===", flush=True)
    
    profiles_csv = os.path.join(os.path.dirname(__file__), "mirna_spaceflight_profiles.csv")
    df_mirna = pd.read_csv(profiles_csv) if os.path.exists(profiles_csv) else pd.DataFrame()
    
    df_targets = pd.DataFrame(PLANT_MIRNA_TARGETS)
    
    network_results = []
    
    studies = ['OSD-437', 'OSD-120', 'OSD-218', 'OSD-565', 'OSD-782']
    
    for study_id in studies:
        de_file = os.path.join(DATA_DIR, study_id, f"{study_id}_differential_expression.csv")
        if os.path.exists(de_file) and os.path.getsize(de_file) > 1000:
            try:
                df_de = pd.read_csv(de_file)
                gene_col = df_de.columns[0]
                lfc_cols = [c for c in df_de.columns if 'log2' in c.lower() or 'foldchange' in c.lower()]
                p_cols = [c for c in df_de.columns if 'pvalue' in c.lower() or 'p.value' in c.lower() or 'padj' in c.lower()]
                
                if lfc_cols:
                    # NOTE: first contrast in the table, not flight vs ground (TODO: select by name).
                    lfc_col = lfc_cols[0]
                    pval_col = p_cols[0] if p_cols else None
                    
                    for _, t_row in df_targets.iterrows():
                        target_id = t_row['Target_Gene_ID']
                        mir_name = t_row['miRNA']
                        
                        match = df_de[df_de[gene_col].astype(str).str.upper().str.startswith(target_id)]
                        if not match.empty:
                            m_row = match.iloc[0]
                            target_lfc = float(m_row[lfc_col]) if pd.notnull(m_row[lfc_col]) else 0.0
                            target_pval = float(m_row[pval_col]) if (pval_col and pd.notnull(m_row[pval_col])) else 1.0
                            
                            mir_match = df_mirna[(df_mirna['Study'] == study_id) & (df_mirna['miRNA_Name'].str.contains(mir_name, case=False, na=False))]
                            mir_lfc = float(mir_match['log2FC'].mean()) if not mir_match.empty else np.nan
                            
                            is_inverse = bool(pd.notnull(mir_lfc) and (mir_lfc * target_lfc < 0))
                            
                            network_results.append({
                                'Study': study_id,
                                'miRNA_Family': mir_name,
                                'miRNA_log2FC': mir_lfc if pd.notnull(mir_lfc) else 0.0,
                                'Target_Gene_ID': target_id,
                                'Target_Symbol': t_row['Target_Symbol'],
                                'Target_log2FC': target_lfc,
                                'Target_PValue': target_pval,
                                'Pathway': t_row['Pathway'],
                                'Regulatory_Relationship': 'Inverse (Repressed Target)' if is_inverse else 'Co-directional / Complex',
                                'Biological_Function': t_row['Function']
                            })
            except Exception as e:
                print(f"Error processing DE table for {study_id}: {e}", flush=True)

    df_net = pd.DataFrame(network_results)
    print(f"Constructed {len(df_net)} miRNA-mRNA interaction pairs across plant spaceflight datasets.", flush=True)
    
    out_csv = os.path.join(os.path.dirname(__file__), "mirna_mrna_regulatory_network.csv")
    df_net.to_csv(out_csv, index=False)
    print(f"Exported regulatory network to {out_csv}", flush=True)
    
    # ILLUSTRATIVE Fisher test, not an enrichment result: the background below is typed,
    # not measured, and each target is counted once per study (7 targets x 5 studies = 35).
    pathway_groups = df_net.groupby('Pathway')
    pathway_records = []
    
    total_genes_in_genome = 25000  # Arabidopsis TAIR10 background gene count
    total_de_genes_in_spaceflight = 3500  # Average significantly altered genes per study
    
    for path_name, group in pathway_groups:
        n_targets = len(group)
        n_significant_targets = len(group[group['Target_PValue'] <= 0.05])
        
        # Contingency Table for Fisher Exact Test:
        # [[Sig Targets in Pathway, Non-Sig Targets in Pathway],
        #  [Sig Background Genes, Non-Sig Background Genes]]
        a = n_significant_targets
        b = n_targets - n_significant_targets
        c = total_de_genes_in_spaceflight - a
        d = (total_genes_in_genome - total_de_genes_in_spaceflight) - b
        
        odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')
        nlog10_p = -np.log10(np.maximum(p_value, 1e-12))
        
        key_mirnas = ', '.join(sorted(group['miRNA_Family'].unique()))
        key_targets = ', '.join(sorted(group['Target_Symbol'].unique()))
        
        pathway_records.append({
            'Pathway': path_name,
            'Key_miRNAs': key_mirnas,
            'Key_Targets': key_targets,
            'Total_Pairs': n_targets,
            'Significant_Pairs': n_significant_targets,
            'Fisher_PValue': float(p_value),
            'Enrichment_Score': float(nlog10_p),
            'Physiological_Role': group['Biological_Function'].iloc[0]
        })
        
    df_path = pd.DataFrame(pathway_records).sort_values('Enrichment_Score', ascending=False)
    path_csv = os.path.join(os.path.dirname(__file__), "target_gene_pathway_enrichment.csv")
    df_path.to_csv(path_csv, index=False)
    print(f"Exported illustrative pathway Fisher-test table to {path_csv}", flush=True)
    
    return df_net, df_path

if __name__ == "__main__":
    build_regulatory_network()
