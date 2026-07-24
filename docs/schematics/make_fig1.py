#!/usr/bin/env python3
"""Generate Figure 1 (a NEW remake) for the smallRNAseq-DREAM manuscript.

The original Plan A/B schematics were lost with the deleted hub repo. This regenerates an
equivalent figure from the documented pipeline steps: the nf-core/smrnaseq pipeline
(deconstructed) beside the recommended cross-species "DREAM" route, with a shared set of
analysis phases showing how the many nf-core steps collapse into the integrated
sRNAtoolbox tools.

Output: fig1_pipeline_overview.svg  (scalable; export to PNG/TIFF for the manuscript with
e.g.  rsvg-convert -f png -o fig1.png fig1_pipeline_overview.svg).

Content is grounded in pipelines/nfcore_smrnaseq/README.md and docs/00_pipeline_decision.md.
No numbers or results are depicted — this is a workflow schematic only.
"""

from xml.sax.saxutils import escape

W = 1140
PAD = 40
COL_W = 300
L_CX = 300          # left column centre (nf-core)
R_CX = 840          # right column centre (DREAM)
BOX_W = 300
BOX_H = 40
STEP = 52           # vertical distance between stacked boxes
BAND_PAD = 18

# phase -> (label, band fill, accent stroke)
PHASES = [
    ("QC & preprocessing",        "#eef2ff", "#6366f1"),
    ("Alignment · annotation · counts", "#ecfdf5", "#10b981"),
    ("Differential expression",   "#fef9c3", "#ca8a04"),
    ("Novel miRNA discovery",     "#ffe4e6", "#f43f5e"),
    ("Reporting",                 "#f1f5f9", "#64748b"),
]

# per phase: left (nf-core) steps and right (DREAM) steps as (tool, role)
LEFT = [
    [("FastQC", "read QC"), ("FastP", "adapter trimming"),
     ("Bowtie2", "contaminant filtering"), ("miRTrace", "sRNA QC")],
    [("Bowtie", "align vs miRBase mature + hairpin"), ("SAMtools", "feature counting"),
     ("Bowtie", "align vs genome (QC)"), ("mirtop", "miRNA + isomiR annotation")],
    [("edgeR", "normalisation · MDS · heatmap")],
    [("miRDeep2", "known + novel discovery")],
    [("MultiQC", "aggregate report")],
]
RIGHT = [
    [("mirnaQC", "QC + preprocessing")],
    [("sRNAbench", "annotation · counts · isomiRs")],
    [("sRNAde", "DESeq2 / edgeR / NOISeq")],
    [("mirNOVO", "de novo novel miRNAs")],
    [("sRNAtoolbox reports", "HTML / tables")],
]

svg = []


def rect(x, y, w, h, rx, fill, stroke, sw=1.5):
    svg.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
               f'rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')


def text(x, y, s, size=15, weight="normal", fill="#0f172a", anchor="middle", style=""):
    svg.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" font-weight="{weight}" '
               f'fill="{fill}" text-anchor="{anchor}" font-family="Helvetica,Arial,sans-serif" '
               f'{style}>{escape(s)}</text>')


def box(cx, y, tool, role, accent):
    x = cx - BOX_W / 2
    rect(x, y, BOX_W, BOX_H, 8, "#ffffff", accent, 1.8)
    text(cx, y + 17, tool, size=15, weight="bold")
    text(cx, y + 32, role, size=11, fill="#475569")


def varrow(cx, y1, y2, color="#94a3b8"):
    svg.append(f'<line x1="{cx:.1f}" y1="{y1:.1f}" x2="{cx:.1f}" y2="{y2:.1f}" '
               f'stroke="{color}" stroke-width="1.8" marker-end="url(#ah)"/>')


# ---- compute band heights from the taller column in each phase --------------------
band_heights = []
for i in range(len(PHASES)):
    n = max(len(LEFT[i]), len(RIGHT[i]))
    band_heights.append(n * STEP - (STEP - BOX_H) + 2 * BAND_PAD)

top = 150
band_tops = []
y = top
for h in band_heights:
    band_tops.append(y)
    y += h + 12
H = y + 90

# ---- header --------------------------------------------------------------------------
svg.insert(0, f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" '
              f'width="{W}" height="{H}" font-family="Helvetica,Arial,sans-serif">')
svg.append('<defs><marker id="ah" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" '
           'markerHeight="7" orient="auto-start-reverse">'
           '<path d="M0,0 L10,5 L0,10 z" fill="#94a3b8"/></marker></defs>')
rect(0, 0, W, H, 0, "#ffffff", "#ffffff", 0)
text(W / 2, 42, "Figure 1  |  Small RNA-seq analysis: nf-core/smrnaseq deconstructed vs. the recommended cross-species (DREAM) route",
     size=17, weight="bold")
text(W / 2, 66, "Shared analysis phases (left) show how the many nf-core steps collapse into the integrated sRNAtoolbox tools of the DREAM route.",
     size=12, fill="#475569")

# column headers
rect(L_CX - COL_W / 2, 92, COL_W, 34, 8, "#0f172a", "#0f172a", 0)
text(L_CX, 114, "A   nf-core/smrnaseq  (deconstructed)", size=14, weight="bold", fill="#ffffff")
rect(R_CX - COL_W / 2, 92, COL_W, 34, 8, "#0e7490", "#0e7490", 0)
text(R_CX, 114, "B   Recommended  'DREAM'  route", size=14, weight="bold", fill="#ffffff")

# ---- bands + boxes -------------------------------------------------------------------
for i, (label, fill, accent) in enumerate(PHASES):
    bt = band_tops[i]
    bh = band_heights[i]
    rect(PAD, bt, W - 2 * PAD, bh, 12, fill, accent, 1.2)
    # phase label (vertical accent chip on the far left)
    rect(PAD, bt, 6, bh, 3, accent, accent, 0)
    text(PAD + 20, bt + 22, label.upper(), size=11, weight="bold", fill=accent, anchor="start")

    for col_cx, steps in ((L_CX, LEFT[i]), (R_CX, RIGHT[i])):
        y0 = bt + BAND_PAD + 14
        for j, (tool, role) in enumerate(steps):
            by = y0 + j * STEP
            box(col_cx, by, tool, role, accent)
            if j > 0:
                varrow(col_cx, by - (STEP - BOX_H) + 2, by - 2)
    # arrow crossing into next band
    if i < len(PHASES) - 1:
        nb = band_tops[i + 1]
        for col_cx in (L_CX, R_CX):
            varrow(col_cx, bt + bh + 1, nb - 1, color="#cbd5e1")

# ---- kingdom note under the novel-discovery phase -----------------------------------
kb = band_tops[3]
note_x = R_CX + COL_W / 2 + 20
# annotation near novel-discovery band, centered under the whole figure instead:
fy = band_tops[-1] + band_heights[-1] + 34
rect(PAD, fy, W - 2 * PAD, 40, 10, "#fff7ed", "#f59e0b", 1.2)
text(W / 2, fy + 17, "Kingdom-specific novel-miRNA prediction (Plan D):  miRDeep2 → animals    ·    miRDeep-P2 (miRDP2) → plants",
     size=12.5, weight="bold", fill="#9a3412")
text(W / 2, fy + 33, "New remake of the lost Plan A/B schematic — workflow only, generated by docs/schematics/make_fig1.py",
     size=10.5, fill="#92400e", style='font-style="italic"')

svg.append("</svg>")

out = "fig1_pipeline_overview.svg"
with open(out, "w") as fh:
    fh.write("\n".join(svg))
print(f"wrote {out}  ({W}x{H})")
