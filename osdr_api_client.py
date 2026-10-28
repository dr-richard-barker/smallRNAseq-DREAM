#!/usr/bin/env python3
"""
OSDR API Data Ingestion & Metadata Refinement Client for Plant Spaceflight Studies
"""

import os
import json
import urllib3
import requests
import pandas as pd

urllib3.disable_warnings()

BASE_META_URL = "https://osdr.nasa.gov/osdr/data/osd/meta"
BASE_FILES_URL = "https://osdr.nasa.gov/osdr/data/osd/files"
BASE_DOWNLOAD_URL = "https://osdr.nasa.gov"

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

# Species and assay labels are typed here. 'assay_type' is the RNA-seq assay's ISA
# 'Library Selection' as read from each OSDR study record on 2026-10-28 (abai check);
# OSDR holds no plant small RNA-seq, and OSD-437 is standard poly(A) RNA-seq.
PLANT_OSD_MAP = [
    {'id': 437, 'species': 'Arabidopsis thaliana', 'assay_type': 'polyA enrichment', 'desc': 'Seedling development in microgravity (poly(A) RNA-seq)'},
    {'id': 120, 'species': 'Arabidopsis thaliana', 'assay_type': 'ribo-depletion', 'desc': 'Arabidopsis spaceflight transcriptome (SpaceX-3)'},
    {'id': 218, 'species': 'Arabidopsis thaliana', 'assay_type': 'polyA enrichment', 'desc': 'ISS Flight Root/Shoot Transcriptomics'},
    {'id': 217, 'species': 'Arabidopsis thaliana', 'assay_type': 'cDNA', 'desc': 'Extraterrestrial Epigenetics & RNA-Seq'},
    {'id': 375, 'species': 'Brachypodium distachyon', 'assay_type': 'Ribo-depletion (32), cDNA (16)', 'desc': 'ISS Brachypodium Root/Shoot Transcriptomics'},
    {'id': 480, 'species': 'Arabidopsis thaliana', 'assay_type': 'not recorded in OSDR', 'desc': 'Plant Gravity Perception (PGP) Flight Study'},
    {'id': 522, 'species': 'Arabidopsis thaliana', 'assay_type': 'polyA enrichment', 'desc': 'Integrative Transcriptomics & Proteomics'},
    {'id': 565, 'species': 'Arabidopsis thaliana', 'assay_type': 'cDNA', 'desc': 'KSC Fixation Tube study (OSDR project type: Ground Study)'},
    {'id': 624, 'species': 'Arabidopsis thaliana', 'assay_type': 'polyA enrichment', 'desc': 'Virgin Galactic Suborbital Spaceflight'},
    {'id': 767, 'species': 'Solanum lycopersicum', 'assay_type': 'polyA enrichment', 'desc': 'ISS Tomato Crop Transcriptomics'},
    {'id': 782, 'species': 'Arabidopsis thaliana', 'assay_type': 'polyA enrichment', 'desc': 'Low-Dose Ionizing Radiation Response'},
    {'id': 59,  'species': 'Brassica rapa', 'assay_type': 'mRNA enrichment', 'desc': 'ROS Gene Network in Microgravity'}
]

def fetch_and_refine_inventory():
    os.makedirs(DATA_DIR, exist_ok=True)
    inventory_rows = []
    
    print("Fetching and refining OSDR plant spaceflight dataset metadata...", flush=True)
    
    for item in PLANT_OSD_MAP:
        num_id = item['id']
        osd_id = f"OSD-{num_id}"
        meta_url = f"{BASE_META_URL}/{num_id}"
        files_url = f"{BASE_FILES_URL}/{num_id}"
        
        try:
            r_meta = requests.get(meta_url, verify=False, timeout=15).json()
            r_files = requests.get(files_url, verify=False, timeout=15).json()
            
            study_meta = r_meta.get('study', {}).get(osd_id, {})
            study_files = r_files.get('studies', {}).get(osd_id, {}).get('study_files', [])
            
            s_list = study_meta.get('studies', [])
            title = s_list[0].get('title', item['desc']) if s_list else item['desc']
            
            diff_expr_url = None
            unnorm_counts_url = None
            vst_counts_url = None
            sample_table_url = None
            bam_count = 0
            
            for f in study_files:
                fname = f.get('file_name', '')
                remote_url = f.get('remote_url', '')
                if remote_url and not remote_url.startswith('http'):
                    remote_url = f"{BASE_DOWNLOAD_URL}{remote_url}"
                    
                if 'differential_expression' in fname.lower() and fname.endswith('.csv'):
                    diff_expr_url = remote_url
                elif ('unnormalized_counts' in fname.lower() or 'star_unnormalized' in fname.lower() or 'rsem_unnormalized' in fname.lower()) and fname.endswith('.csv'):
                    if not unnorm_counts_url or 'rsem' in fname.lower():
                        unnorm_counts_url = remote_url
                elif ('vst_counts' in fname.lower() or 'normalized_counts' in fname.lower()) and fname.endswith('.csv'):
                    if not vst_counts_url or 'vst' in fname.lower():
                        vst_counts_url = remote_url
                elif 'sampletable' in fname.lower() and fname.endswith('.csv'):
                    sample_table_url = remote_url
                    
                if fname.endswith('.bam'):
                    bam_count += 1
                    
            row = {
                'Accession': osd_id,
                'Numeric_ID': num_id,
                'Title': title,
                'Organism': item['species'],
                'Assay': 'RNA Sequencing (RNA-Seq)',
                'Assay_Subtype': item['assay_type'],
                'Total_Files': len(study_files),
                'BAM_Files_Count': bam_count,
                'Has_Diff_Expr_CSV': bool(diff_expr_url),
                'Has_Unnorm_Counts_CSV': bool(unnorm_counts_url),
                'Diff_Expr_URL': diff_expr_url or '',
                'Unnorm_Counts_URL': unnorm_counts_url or '',
                'VST_Counts_URL': vst_counts_url or '',
                'Sample_Table_URL': sample_table_url or ''
            }
            inventory_rows.append(row)
            print(f"Indexed {osd_id}: {item['species']} | Subtype: {item['assay_type']} | Files: {len(study_files)}", flush=True)

        except Exception as e:
            print(f"Error indexing {osd_id}: {e}", flush=True)

    df_inv = pd.DataFrame(inventory_rows).sort_values('Numeric_ID')
    df_inv.to_csv(os.path.join(os.path.dirname(__file__), 'plant_osdr_inventory.csv'), index=False)
    with open(os.path.join(os.path.dirname(__file__), 'plant_osdr_inventory.json'), 'w') as f:
        json.dump(inventory_rows, f, indent=2)
        
    print(f"\nInventory Metadata Refinement Complete! Saved to plant_osdr_inventory.csv ({len(df_inv)} datasets).", flush=True)
    return df_inv

if __name__ == "__main__":
    fetch_and_refine_inventory()
