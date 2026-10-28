#!/usr/bin/env python3
"""
Figure generator for smallRNAseq-DREAM. Figures 1-4 and 6 plot the tables in the repo root;
Figures 5 and 7 are schematics. See the README status note: the miRNA-side values are not yet
interpretable (mislabelled loci; first deposited contrast, not flight vs ground).
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import networkx as nx

plt.style.use('ggplot')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 1.0

FIGURES_DIR = os.path.join(os.path.dirname(__file__), "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

# ---------------------------------------------------------
# FIGURE 1: OSDR Plant Dataset Landscape & Stacked RNA-Seq Subtypes
# ---------------------------------------------------------
def generate_figure1_osdr_landscape():
    print("Generating Figure 1: NASA OSDR plant dataset inventory...", flush=True)
    inv_file = os.path.join(os.path.dirname(__file__), "plant_osdr_inventory.csv")
    if not os.path.exists(inv_file):
        print("Inventory file missing!", flush=True)
        return
        
    df = pd.read_csv(inv_file)
    
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    
    # Panel A: Species Breakdown (Horizontal Bar Chart)
    org_counts = df['Organism'].value_counts()
    colors = ['#2b5c8f', '#d95f02', '#2ca02c', '#7570b3', '#e7298a']
    
    y_pos = np.arange(len(org_counts))
    axes[0].barh(y_pos, org_counts.values, color=colors[:len(org_counts)], edgecolor='black', linewidth=0.8)
    axes[0].set_yticks(y_pos)
    axes[0].set_yticklabels(org_counts.index, fontsize=11, fontweight='bold')
    axes[0].invert_yaxis()
    axes[0].set_title("A. Indexed OSDR plant studies by species", fontsize=13, fontweight='bold', pad=12)
    axes[0].set_xlabel("Number of OSDR studies (project types: see README inventory)", fontsize=11)
    
    for i, v in enumerate(org_counts.values):
        axes[0].text(v + 0.15, i, f"{v} datasets", va='center', fontweight='bold', fontsize=10)
    axes[0].set_xlim(0, max(org_counts.values) + 2)
        
    # Panel B: Stacked Bar Plot of RNA-Seq Subtypes by Species
    ct_table = pd.crosstab(df['Organism'], df['Assay_Subtype'])
    subtype_colors = {
        'polyA enrichment': '#4daf4a',
        'mRNA enrichment': '#a6d854',
        'ribo-depletion': '#377eb8',
        'Ribo-depletion (32), cDNA (16)': '#80b1d3',
        'cDNA': '#984ea3',
        'not recorded in OSDR': '#bbbbbb'
    }
    
    bottom = np.zeros(len(ct_table))
    for col in ct_table.columns:
        c_val = ct_table[col].values
        color = subtype_colors.get(col, '#999999')
        axes[1].bar(range(len(ct_table.index)), c_val, bottom=bottom, label=col, color=color, edgecolor='black', linewidth=0.8)
        bottom += c_val
        
    axes[1].set_title("B. RNA-seq library selection (OSDR ISA record)", fontsize=13, fontweight='bold', pad=12)
    axes[1].set_ylabel("Number of Datasets", fontsize=11)
    axes[1].set_xticks(range(len(ct_table.index)))
    axes[1].set_xticklabels(ct_table.index, rotation=15, ha='right', fontweight='bold')
    axes[1].legend(title="Library selection", loc='upper right', fontsize=9.5)
    
    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, "Figure1_OSDR_Plant_Dataset_Landscape.png")
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved Figure 1 to {fig_path}", flush=True)

# ---------------------------------------------------------
# FIGURE 2: MicroRNA Differential Expression Profiles
# ---------------------------------------------------------
def generate_figure2_mirna_volcano_heatmap():
    print("Generating Figure 2: map-labelled locus log2FC (not interpretable yet)...", flush=True)
    prof_file = os.path.join(os.path.dirname(__file__), "mirna_spaceflight_profiles.csv")
    if not os.path.exists(prof_file):
        return
        
    df = pd.read_csv(prof_file)
    
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    
    df['nlog10_p'] = -np.log10(np.maximum(df['PValue'], 1e-10))
    
    sig_mask = df['Significant']
    up_mask = sig_mask & (df['log2FC'] > 0)
    down_mask = sig_mask & (df['log2FC'] < 0)
    ns_mask = ~sig_mask
    
    axes[0].scatter(df.loc[ns_mask, 'log2FC'], df.loc[ns_mask, 'nlog10_p'], color='grey', alpha=0.4, label='Not Significant')
    axes[0].scatter(df.loc[up_mask, 'log2FC'], df.loc[up_mask, 'nlog10_p'], color='#d95f02', alpha=0.85, s=65, label='log2FC > 0.5, P <= 0.05')
    axes[0].scatter(df.loc[down_mask, 'log2FC'], df.loc[down_mask, 'nlog10_p'], color='#7570b3', alpha=0.85, s=65, label='log2FC < -0.5, P <= 0.05')
    
    axes[0].axvline(0.5, color='black', linestyle='--', linewidth=0.8)
    axes[0].axvline(-0.5, color='black', linestyle='--', linewidth=0.8)
    axes[0].axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.8)
    
    axes[0].set_title("A. Rows labelled MIR* by the typed locus map\n(most are not miRNA genes; see README TODO)", fontsize=12, fontweight='bold', pad=12)
    axes[0].set_xlabel("log2FC, first contrast in each GeneLab table (not flight vs ground)", fontsize=10.5)
    axes[0].set_ylabel("-log10 (P-Value)", fontsize=11)
    axes[0].legend(loc='upper right')
    
    for _, row in df[df['Significant']].head(10).iterrows():
        axes[0].annotate(str(row['miRNA_Name']), (row['log2FC'], row['nlog10_p']),
                         xytext=(5, 5), textcoords='offset points', fontsize=9, fontweight='bold')
                         
    pivot_df = df.pivot_table(index='miRNA_Name', columns='Study', values='log2FC', aggfunc='mean').fillna(0)
    if not pivot_df.empty:
        sns.heatmap(pivot_df.head(15), cmap='vlag', center=0, annot=True, fmt='.2f', ax=axes[1], cbar_kws={'label': 'log2 Fold Change'})
        axes[1].set_title("B. log2FC by map label and study (first 15 labels; 0.00 = no row)", fontsize=12, fontweight='bold', pad=12)
        axes[1].set_ylabel("Label from typed locus map (3 of 29 correct)", fontsize=10.5)
        axes[1].set_xlabel("OSDR study (first deposited contrast)", fontsize=11)
        
    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, "Figure2_microRNA_Differential_Expression.png")
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved Figure 2 to {fig_path}", flush=True)

# ---------------------------------------------------------
# FIGURE 3: miRNA-mRNA Regulatory Coupling Scatter
# ---------------------------------------------------------
def generate_figure3_regulatory_network():
    print("Generating Figure 3: curated miRNA-target log2FC sign comparison...", flush=True)
    net_file = os.path.join(os.path.dirname(__file__), "mirna_mrna_regulatory_network.csv")
    if not os.path.exists(net_file):
        return
        
    df = pd.read_csv(net_file)
    
    fig, ax = plt.subplots(figsize=(9.5, 6.5))
    
    colors = ['#e41a1c' if r == 'Inverse (Repressed Target)' else '#377eb8' for r in df['Regulatory_Relationship']]
    
    ax.scatter(df['miRNA_log2FC'], df['Target_log2FC'], c=colors, s=85, alpha=0.85, edgecolors='black', linewidth=0.6)
    
    ax.axhline(0, color='black', linestyle='--', linewidth=0.8)
    ax.axvline(0, color='black', linestyle='--', linewidth=0.8)
    
    ax.set_title("Figure 3. Curated miRNA-target pairs: sign of log2FC (not interpretable yet)", fontsize=12.5, fontweight='bold', pad=12)
    ax.set_xlabel("'miRNA' value: mean log2FC of map-labelled rows (mostly not miRNA genes)", fontsize=10.5)
    ax.set_ylabel("Target gene log2FC (first contrast in each GeneLab table)", fontsize=10.5)
    
    # Label only the opposite-sign pairs furthest from the origin, one label per pair name,
    # with staggered offsets, so labels do not overlap the dense cluster near zero.
    inv = df[df['Regulatory_Relationship'] == 'Inverse (Repressed Target)'].copy()
    inv['dist'] = np.hypot(inv['miRNA_log2FC'], inv['Target_log2FC'])
    inv['label'] = inv['miRNA_Family'] + ' -> ' + inv['Target_Symbol'] + ' (' + inv['Study'] + ')'
    offsets = [(10, 10), (10, -16), (-10, 12), (-10, -18)]
    for k, (_, row) in enumerate(inv.sort_values('dist', ascending=False).drop_duplicates('label').head(6).iterrows()):
        dx, dy = offsets[k % len(offsets)]
        ax.annotate(row['label'], (row['miRNA_log2FC'], row['Target_log2FC']),
                    xytext=(dx, dy), textcoords='offset points', fontsize=8.5, fontweight='bold',
                    ha='left' if dx > 0 else 'right',
                    arrowprops=dict(arrowstyle='-', color='#555555', lw=0.6),
                    bbox=dict(boxstyle="round,pad=0.2", fc="yellow", alpha=0.3))
                    
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label=f"Opposite sign (n={int((df['Regulatory_Relationship'] == 'Inverse (Repressed Target)').sum())})", markerfacecolor='#e41a1c', markersize=10),
        Line2D([0], [0], marker='o', color='w', label=f"Same sign, or no 'miRNA' row (n={int((df['Regulatory_Relationship'] != 'Inverse (Repressed Target)').sum())})", markerfacecolor='#377eb8', markersize=10)
    ]
    ax.legend(handles=legend_elements, loc='lower right')
    
    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, "Figure3_miRNA_mRNA_Regulatory_Coupling.png")
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved Figure 3 to {fig_path}", flush=True)

# ---------------------------------------------------------
# FIGURE 4: Illustrative Fisher test (typed background)
# ---------------------------------------------------------
def generate_figure4_pathway_enrichment():
    print("Generating Figure 4: illustrative Fisher test (typed background)...", flush=True)
    path_file = os.path.join(os.path.dirname(__file__), "target_gene_pathway_enrichment.csv")
    if not os.path.exists(path_file):
        return
        
    df = pd.read_csv(path_file).sort_values('Enrichment_Score', ascending=True)
    
    fig, ax = plt.subplots(figsize=(10, 5.5))
    
    pathways = df['Pathway'].values
    enrichment_scores = df['Enrichment_Score'].values
    
    bars = ax.barh(pathways, enrichment_scores, color='#2c7bb6', edgecolor='black', linewidth=0.8)
    ax.set_title("Figure 4. ILLUSTRATIVE Fisher test, not an enrichment result\n(typed background: 25,000 genes, 3,500 'DE'; each target counted once per study)", fontsize=12, fontweight='bold', pad=12)
    ax.set_xlabel("-log10 P under the typed background (illustrative)", fontsize=11)
    
    for i, (bar, mir, pval) in enumerate(zip(bars, df['Key_miRNAs'].values, df['Fisher_PValue'].values)):
        p_str = f"P={pval:.2e}" if pval < 0.001 else f"P={pval:.3f}"
        ax.text(bar.get_width() + 0.15, bar.get_y() + bar.get_height()/2, f"Regulators: {mir} ({p_str})", va='center', fontsize=9, fontweight='bold', color='#333333')
        
    ax.set_xlim(0, max(enrichment_scores) + 3.5)
    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, "Figure4_miRNA_Functional_Pathway_Enrichment.png")
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved Figure 4 to {fig_path}", flush=True)

# ---------------------------------------------------------
# FIGURE 5: MicroRNA -> Target Gene -> Ontology Sankey Diagram
# ---------------------------------------------------------
def generate_figure5_sankey_diagram():
    print("Generating Figure 5: MicroRNA -> Target Gene -> Cellular Ontology Sankey Flow Diagram...", flush=True)
    
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.axis('off')
    
    mirnas = ['miR156', 'miR166', 'miR398', 'miR168', 'miR160', 'miR164', 'miR172']
    targets = ['SPL2/SPL9', 'PHB/REV', 'CSD1/CSD2', 'AGO1', 'ARF10/16', 'NAC1/CUC2', 'AP2/TOE1']
    ontologies = [
        'Vegetative Phase Change & Root Branching',
        'Vascular Polarity & Secondary Cell Wall',
        'ROS Scavenging & Oxidative Stress Defense',
        'sRNA Silencing Machinery Feedback',
        'Auxin Signaling & Root Gravitropism',
        'Floral Transition & Meristem Fate'
    ]
    
    x_mir, x_tar, x_ont = 0.05, 0.45, 0.88
    
    y_mir = np.linspace(0.85, 0.15, len(mirnas))
    y_tar = np.linspace(0.85, 0.15, len(targets))
    y_ont = np.linspace(0.85, 0.15, len(ontologies))
    
    flows = [
        (0, 0, 0, '#e41a1c'),
        (1, 1, 1, '#377eb8'),
        (2, 2, 2, '#4daf4a'),
        (3, 3, 3, '#984ea3'),
        (4, 4, 4, '#ff7f00'),
        (5, 5, 4, '#ff7f00'),
        (6, 6, 5, '#ffff33')
    ]
    
    for mir_idx, tar_idx, ont_idx, color in flows:
        p1 = patches.PathPatch(
            matplotlib.path.Path(
                [(x_mir + 0.08, y_mir[mir_idx]),
                 (x_mir + 0.25, y_mir[mir_idx]),
                 (x_tar - 0.10, y_tar[tar_idx]),
                 (x_tar - 0.02, y_tar[tar_idx])],
                [matplotlib.path.Path.MOVETO, matplotlib.path.Path.CURVE4, matplotlib.path.Path.CURVE4, matplotlib.path.Path.CURVE4]
            ),
            fc='none', ec=color, lw=4, alpha=0.65
        )
        ax.add_patch(p1)
        
        p2 = patches.PathPatch(
            matplotlib.path.Path(
                [(x_tar + 0.08, y_tar[tar_idx]),
                 (x_tar + 0.25, y_tar[tar_idx]),
                 (x_ont - 0.18, y_ont[ont_idx]),
                 (x_ont - 0.02, y_ont[ont_idx])],
                [matplotlib.path.Path.MOVETO, matplotlib.path.Path.CURVE4, matplotlib.path.Path.CURVE4, matplotlib.path.Path.CURVE4]
            ),
            fc='none', ec=color, lw=4, alpha=0.65
        )
        ax.add_patch(p2)
        
    for i, (m, y) in enumerate(zip(mirnas, y_mir)):
        ax.text(x_mir, y, f"  {m}  ", va='center', ha='right', fontsize=10.5, fontweight='bold',
                bbox=dict(boxstyle="round,pad=0.4", fc="#fee0d2", ec="#de2d26", lw=1.5))
                
    for i, (t, y) in enumerate(zip(targets, y_tar)):
        ax.text(x_tar, y, f" {t} ", va='center', ha='center', fontsize=10.5, fontweight='bold',
                bbox=dict(boxstyle="round,pad=0.4", fc="#deebf7", ec="#3182bd", lw=1.5))
                
    for i, (o, y) in enumerate(zip(ontologies, y_ont)):
        ax.text(x_ont, y, f"  {o}  ", va='center', ha='left', fontsize=10, fontweight='bold',
                bbox=dict(boxstyle="round,pad=0.4", fc="#e5f5e0", ec="#31a354", lw=1.5))
                
    ax.text(x_mir - 0.05, 0.95, "MicroRNA Regulators", fontsize=13, fontweight='bold', color='#de2d26')
    ax.text(x_tar, 0.95, "Target Transcription Factors & Enzymes", fontsize=13, fontweight='bold', color='#3182bd', ha='center')
    ax.text(x_ont + 0.05, 0.95, "Cellular Effects & Ontologies", fontsize=13, fontweight='bold', color='#31a354')
    
    ax.set_title("Figure 5. Schematic (no data): literature miRNA -> target -> process relationships", fontsize=14, fontweight='bold', pad=20)
    
    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, "Figure5_miRNA_Target_Ontology_Sankey.png")
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved Figure 5 to {fig_path}", flush=True)

# ---------------------------------------------------------
# FIGURE 6: Curated miRNA-target-pathway graph
# ---------------------------------------------------------
def generate_figure6_knowledge_graph():
    print("Generating Figure 6: curated miRNA-target-pathway graph...", flush=True)
    net_file = os.path.join(os.path.dirname(__file__), "mirna_mrna_regulatory_network.csv")
    if not os.path.exists(net_file):
        return
        
    df_net = pd.read_csv(net_file)
    
    G = nx.Graph()
    
    # Mean log2FC per miRNA label and target symbol, from the network table
    mir_lfc_mean = df_net.groupby('miRNA_Family')['miRNA_log2FC'].mean().to_dict()
    tar_lfc_mean = df_net.groupby('Target_Symbol')['Target_log2FC'].mean().to_dict()
    
    mirna_nodes = [f"{m} ({mir_lfc_mean.get(m, 0):+.2f})" for m in sorted(mir_lfc_mean.keys())]
    target_nodes = [f"{t} ({tar_lfc_mean.get(t, 0):+.2f})" for t in sorted(tar_lfc_mean.keys())]
    pathway_hubs = ['Phase Change & Vegetative Dev', 'Vascular Patterning', 'ROS Scavenging', 'sRNA Silencing Tuning', 'Auxin Signaling']
    
    # Add Nodes
    for m, m_raw in zip(mirna_nodes, sorted(mir_lfc_mean.keys())):
        lfc = mir_lfc_mean[m_raw]
        G.add_node(m, layer=0, type='miRNA', color='#f46d43' if lfc > 0 else '#74add1')
        
    for t, t_raw in zip(target_nodes, sorted(tar_lfc_mean.keys())):
        lfc = tar_lfc_mean[t_raw]
        G.add_node(t, layer=1, type='Target', color='#f46d43' if lfc > 0 else '#74add1')
        
    for p in pathway_hubs:
        G.add_node(p, layer=2, type='Pathway', color='#abdda4')
        
    # Add edges for the curated miRNA-target pairs
    for _, row in df_net.drop_duplicates(['miRNA_Family', 'Target_Symbol']).iterrows():
        m_label = f"{row['miRNA_Family']} ({mir_lfc_mean.get(row['miRNA_Family'], 0):+.2f})"
        t_label = f"{row['Target_Symbol']} ({tar_lfc_mean.get(row['Target_Symbol'], 0):+.2f})"
        if m_label in G and t_label in G:
            G.add_edge(m_label, t_label)
            
    # Connect Targets to Pathway Hubs
    path_map = {
        'SPL2': 'Phase Change & Vegetative Dev', 'SPL3': 'Phase Change & Vegetative Dev', 'SPL9': 'Phase Change & Vegetative Dev', 'SPL10': 'Phase Change & Vegetative Dev', 'SPL13': 'Phase Change & Vegetative Dev', 'AP2': 'Phase Change & Vegetative Dev', 'TOE1': 'Phase Change & Vegetative Dev',
        'PHB': 'Vascular Patterning', 'PHV': 'Vascular Patterning', 'REV': 'Vascular Patterning', 'ATHB8': 'Vascular Patterning',
        'CSD1': 'ROS Scavenging', 'CSD2': 'ROS Scavenging', 'CCS1': 'ROS Scavenging',
        'AGO1': 'sRNA Silencing Tuning',
        'ARF10': 'Auxin Signaling', 'ARF16': 'Auxin Signaling', 'ARF17': 'Auxin Signaling', 'NAC1': 'Auxin Signaling', 'CUC1': 'Auxin Signaling', 'CUC2': 'Auxin Signaling', 'TIR1': 'Auxin Signaling', 'AFB1': 'Auxin Signaling', 'AFB2': 'Auxin Signaling'
    }
    
    for t_raw, p_hub in path_map.items():
        t_label = f"{t_raw} ({tar_lfc_mean.get(t_raw, 0):+.2f})"
        if t_label in G and p_hub in G:
            G.add_edge(t_label, p_hub)
            
    pos = nx.multipartite_layout(G, subset_key='layer')
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    node_colors = [G.nodes[n]['color'] for n in G.nodes]
    node_sizes = [1500 if G.nodes[n]['type'] == 'Pathway' else 1000 for n in G.nodes]
    
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, edgecolors='black', linewidths=1.2, ax=ax)
    nx.draw_networkx_edges(G, pos, width=1.8, alpha=0.55, edge_color='#444444', ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=9.0, font_weight='bold', ax=ax)
    
    ax.set_title("Figure 6. Curated miRNA-target-pathway graph; node values = mean log2FC of the first\ndeposited contrast (miRNA nodes use mislabelled loci; not interpretable yet)", fontsize=13, fontweight='bold', pad=15)
    ax.axis('off')
    
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='mean log2FC > 0 (first deposited contrast)', markerfacecolor='#f46d43', markersize=12),
        Line2D([0], [0], marker='o', color='w', label='mean log2FC < 0 (first deposited contrast)', markerfacecolor='#74add1', markersize=12),
        Line2D([0], [0], marker='o', color='w', label='Cellular Pathway Hub', markerfacecolor='#abdda4', markersize=12)
    ]
    ax.legend(handles=legend_elements, loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.05))
    
    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, "Figure6_Gene_Abundance_Knowledge_Graph.png")
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved Figure 6 to {fig_path}", flush=True)

# ---------------------------------------------------------
# FIGURE 7: Schematic pathway map with table means
# ---------------------------------------------------------
def generate_figure7_kegg_pathway_coupling():
    print("Generating Figure 7: schematic pathway map with table means...", flush=True)
    net_file = os.path.join(os.path.dirname(__file__), "mirna_mrna_regulatory_network.csv")
    if not os.path.exists(net_file):
        return
        
    df_net = pd.read_csv(net_file)
    
    # Mean log2FC from the network table for the nodes that carry a value
    lfc_map = {}
    for mir in ['miR393', 'miR160', 'miR398']:
        sub = df_net[df_net['miRNA_Family'] == mir]
        lfc_map[mir] = sub['miRNA_log2FC'].mean() if not sub.empty else 0.0
        
    for target in ['TIR1', 'AFB1', 'ARF10', 'ARF16', 'CSD1', 'CSD2']:
        sub = df_net[df_net['Target_Symbol'] == target]
        lfc_map[target] = sub['Target_log2FC'].mean() if not sub.empty else 0.0
        
    fig, ax = plt.subplots(figsize=(13, 7.5))
    ax.axis('off')
    
    rect_aux = patches.FancyBboxPatch((0.05, 0.45), 0.88, 0.48, boxstyle="round,pad=0.03", fc="#f7fcf5", ec="#238b45", lw=2)
    rect_ros = patches.FancyBboxPatch((0.05, 0.05), 0.88, 0.35, boxstyle="round,pad=0.03", fc="#fff5f0", ec="#cb181d", lw=2)
    ax.add_patch(rect_aux)
    ax.add_patch(rect_ros)
    
    ax.text(0.08, 0.88, "AUXIN GRAVITROPIC SIGNALING & MERISTEM CASCADE", fontsize=12, fontweight='bold', color="#238b45")
    ax.text(0.08, 0.35, "ROS DETOXIFICATION & OXIDATIVE STRESS PATHWAY", fontsize=12, fontweight='bold', color="#cb181d")
    
    nodes = [
        # Auxin Cascade
        # Values are means from mirna_mrna_regulatory_network.csv; 'src' names the column
        # entry used. Phenotype / stress boxes are schematic and carry no value.
        {'name': 'Auxin Perception\n(TIR1 / AFB1)', 'x': 0.15, 'y': 0.65, 'fc': lfc_map['TIR1'], 'src': 'TIR1', 'type': 'mRNA'},
        {'name': 'miR393 Regulator\n(map label miR393)', 'x': 0.15, 'y': 0.78, 'fc': lfc_map['miR393'], 'src': 'miR393 rows', 'type': 'miRNA'},
        {'name': 'Auxin Response TFs\n(ARF10 / ARF16)', 'x': 0.48, 'y': 0.65, 'fc': lfc_map['ARF10'], 'src': 'ARF10', 'type': 'mRNA'},
        {'name': 'miR160 Regulator\n(map label miR160)', 'x': 0.48, 'y': 0.78, 'fc': lfc_map['miR160'], 'src': 'miR160 rows', 'type': 'miRNA'},
        {'name': 'Root Cap Differentiation &\nGravitropic Bending', 'x': 0.82, 'y': 0.65, 'type': 'Phenotype'},
        
        # ROS Cascade
        {'name': 'Superoxide Radicals\n(O2•- Generation)', 'x': 0.15, 'y': 0.20, 'type': 'Stress'},
        {'name': 'Cu/Zn Superoxide Dismutase\n(CSD1 / CSD2)', 'x': 0.48, 'y': 0.20, 'fc': lfc_map['CSD1'], 'src': 'CSD1', 'type': 'mRNA'},
        {'name': 'miR398 Repressor\n(map label miR398)', 'x': 0.48, 'y': 0.08, 'fc': lfc_map['miR398'], 'src': 'miR398 rows', 'type': 'miRNA'},
        {'name': 'H2O2 Accumulation &\nCell Wall Peroxidation', 'x': 0.82, 'y': 0.20, 'type': 'Phenotype'}
    ]
    
    for n in nodes:
        fc_str = f"\n[{n['src']}: {n['fc']:+.2f}]" if 'fc' in n else ""
        box_fc = '#f46d43' if n.get('fc', 0) > 0 else ('#74add1' if n.get('fc', 0) < 0 else '#ffffff')
        
        ax.text(n['x'], n['y'], f"{n['name']}{fc_str}", va='center', ha='center', fontsize=9.5, fontweight='bold',
                bbox=dict(boxstyle="round,pad=0.4", fc=box_fc, ec="black", lw=1.2))
                
    arrows = [
        ((0.15, 0.74), (0.15, 0.69), 'rep'),
        ((0.25, 0.65), (0.38, 0.65), 'act'),
        ((0.48, 0.74), (0.48, 0.69), 'rep'),
        ((0.58, 0.65), (0.70, 0.65), 'act'),
        ((0.25, 0.20), (0.38, 0.20), 'act'),
        ((0.48, 0.12), (0.48, 0.16), 'rep'),
        ((0.58, 0.20), (0.70, 0.20), 'act')
    ]
    
    for (x1, y1), (x2, y2), kind in arrows:
        color = '#d95f02' if kind == 'rep' else '#2b5c8f'
        ls = '--' if kind == 'rep' else '-'
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="->" if kind == 'act' else "-|>", color=color, lw=2.5, ls=ls))
                    
    ax.set_title("Figure 7. Schematic pathway map. Bracketed values are table means of the first deposited\ncontrast (not flight vs ground; miRNA rows mislabelled). Phenotype/stress boxes carry no data.", fontsize=12.5, fontweight='bold', pad=15)
    
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='s', color='w', label='table mean > 0', markerfacecolor='#f46d43', markersize=12),
        Line2D([0], [0], marker='s', color='w', label='table mean < 0', markerfacecolor='#74add1', markersize=12),
        Line2D([0], [0], marker='s', color='w', label='schematic box, no data', markerfacecolor='#ffffff', markeredgecolor='black', markersize=12),
        Line2D([0], [0], color='#d95f02', lw=2, ls='--', label='MicroRNA Repression (Post-Transcriptional)'),
        Line2D([0], [0], color='#2b5c8f', lw=2, ls='-', label='Transcriptional Activation / Flow')
    ]
    ax.legend(handles=legend_elements, loc='lower center', ncol=5, bbox_to_anchor=(0.5, -0.06), fontsize=9)
    
    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, "Figure7_KEGG_Pathway_RNA_Coupling.png")
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved Figure 7 to {fig_path}", flush=True)

def generate_all_figures():
    generate_figure1_osdr_landscape()
    generate_figure2_mirna_volcano_heatmap()
    generate_figure3_regulatory_network()
    generate_figure4_pathway_enrichment()
    generate_figure5_sankey_diagram()
    generate_figure6_knowledge_graph()
    generate_figure7_kegg_pathway_coupling()
    print("\nAll smallRNAseq-DREAM manuscript & advanced visualization figures generated successfully!", flush=True)

if __name__ == "__main__":
    generate_all_figures()
