#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

SRC = Path("results/ncRF_17_26/figures/source_data")
OUT = Path("results/ncRF_17_26/figures/final")
OUT.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# Load frozen source data
# ------------------------------------------------------------

A = pd.read_csv(
    SRC/"Figure_ncRF_panelA_stage_remodeling.tsv",
    sep="\t"
)

B = pd.read_csv(
    SRC/"Figure_ncRF_panelB_Endosperm_classes.tsv",
    sep="\t"
)

C = pd.read_csv(
    SRC/"Figure_ncRF_panelC_LSU3prime_CANONICAL.tsv",
    sep="\t"
)

D = pd.read_csv(
    SRC/"Figure_ncRF_panelD_Endosperm_effects.tsv",
    sep="\t"
)

E = pd.read_csv(
    SRC/"Figure_ncRF_panelE_snoR21_targets.tsv",
    sep="\t"
)

F = pd.read_csv(
    SRC/"Figure_ncRF_panelF_fraction_negative_nulls.tsv",
    sep="\t"
)

FC = pd.read_csv(
    SRC/"Figure_ncRF_panelF_combined_effect_null.tsv",
    sep="\t"
)

# ------------------------------------------------------------
# Global appearance
# ------------------------------------------------------------

plt.rcParams.update({
    "font.size": 8,
    "axes.titlesize": 10,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

fig = plt.figure(figsize=(12, 14))

gs = fig.add_gridspec(
    3, 2,
    height_ratios=[0.85, 1.45, 1.15],
    hspace=0.42,
    wspace=0.32
)

axA = fig.add_subplot(gs[0,0])
axB = fig.add_subplot(gs[0,1])
axC = fig.add_subplot(gs[1,0])
axD = fig.add_subplot(gs[1,1])
axE = fig.add_subplot(gs[2,0])
axF = fig.add_subplot(gs[2,1])

def panel_label(ax, letter):
    ax.text(
        -0.12, 1.08, letter,
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        va="top"
    )

# ============================================================
# A — Stage specificity
# ============================================================

stage_col = "stage"

tested_col = next(
    c for c in A.columns
    if c.lower() in {"tested","n_tested"}
)

up_col = next(
    c for c in A.columns
    if c.lower() in {"up","n_up","strong_up"}
)

down_col = next(
    c for c in A.columns
    if c.lower() in {"down","n_down","strong_down"}
)

labels = [
    "Leaves",
    "Inflorescence\nmeristem",
    "Pollen",
    "Unfertilized\novules",
    "Fertilized\novules",
    "Embryo",
    "Endosperm",
    "Progeny"
]

x=np.arange(len(A))

axA.bar(
    x,
    A[up_col],
    label="UP"
)

axA.bar(
    x,
    -A[down_col],
    label="DOWN"
)

axA.axhline(0, linewidth=0.7)

axA.set_xticks(x)
axA.set_xticklabels(labels, rotation=40, ha="right")

axA.set_ylabel("Strong endpoint families")
axA.set_title("Stage-specific ncRF remodeling")
axA.legend(frameon=False)

for i,r in A.iterrows():
    total=int(r[up_col]+r[down_col])
    if total:
        y = r[up_col] if r[up_col] else -r[down_col]
        axA.text(
            i,
            y + (0.8 if y >= 0 else -0.8),
            str(total),
            ha="center",
            va="bottom" if y >= 0 else "top",
            fontsize=7
        )

panel_label(axA,"A")

# ============================================================
# B — Endosperm composition
# ============================================================

class_col = next(
    c for c in B.columns
    if c.lower() in {"ncrna_class","class"}
)

bup = next(
    c for c in B.columns
    if c.lower() in {"up","n_up","strong_up"}
)

bdown = next(
    c for c in B.columns
    if c.lower() in {"down","n_down","strong_down"}
)

classes=["tRNA","snoRNA","rRNA","snRNA"]

BB=B.set_index(class_col).loc[classes]

y=np.arange(len(BB))

axB.barh(
    y,
    BB[bup],
    label="UP"
)

axB.barh(
    y,
    -BB[bdown],
    label="DOWN"
)

axB.axvline(0,linewidth=0.7)
axB.set_yticks(y)
axB.set_yticklabels(classes)
axB.invert_yaxis()

axB.set_xlabel("Strong endpoint families")
axB.set_title("Endosperm: 31 UP / 2 DOWN")
axB.legend(frameon=False)

axB.text(
    0.98,0.05,
    "26/33 robust to\nalternative normalization",
    transform=axB.transAxes,
    ha="right",
    va="bottom",
    fontsize=7
)

panel_label(axB,"B")

# ============================================================
# C — canonical nested LSU 3' series
# ============================================================

CC=C.sort_values("length")

matrix=CC[["CT1","CT2","TR1","TR2"]].to_numpy(dtype=float)

# log10(count+1) only for visualization
im=axC.imshow(
    np.log10(matrix+1),
    aspect="auto"
)

axC.set_xticks(range(4))
axC.set_xticklabels(["CT1","CT2","TR1","TR2"])

axC.set_yticks(range(len(CC)))
axC.set_yticklabels(
    [
        f'{int(r.length)} nt  {r.sequence}'
        for _,r in CC.iterrows()
    ],
    fontfamily="monospace",
    fontsize=6.2
)

axC.set_xlabel("Endosperm library")
axC.set_title("Nested LSU 3′-terminal ncRFs")

cb=fig.colorbar(
    im,
    ax=axC,
    fraction=0.046,
    pad=0.03
)
cb.set_label("log10(raw reads + 1)")

axC.text(
    0.02,-0.14,
    "Canonical series: 99.7% of CT and 98.8% of TR LSU 3′ reads\n"
    "Family-level TMM: log2FC = +4.862; FDR = 2.76×10⁻⁵",
    transform=axC.transAxes,
    fontsize=7,
    va="top"
)

panel_label(axC,"C")

# ============================================================
# D — all 33 Endosperm effects
# ============================================================

DD=D.copy().sort_values("log2FC")

labelsD=[
    f"{r.family} {str(r.endpoint).replace('prime_CCA','′-CCA').replace('prime','′')}"
    for _,r in DD.iterrows()
]

y=np.arange(len(DD))

axD.scatter(
    DD["log2FC"],
    y,
    s=22
)

axD.axvline(0,linewidth=0.8)

axD.set_yticks(y)
axD.set_yticklabels(labelsD,fontsize=5.8)

axD.set_xlabel("log2 fold change (TMM)")
axD.set_title("33 strong Endosperm endpoint families")

# mark normalization-sensitive calls
robust = (
    DD["normalization_robust"]
    .astype(str)
    .str.lower()
    .isin(["yes","true","1","t","y"])
)

for yy,(xx,rr) in enumerate(
    zip(DD["log2FC"],robust)
):
    if not rr:
        axD.scatter(
            [xx],[yy],
            s=55,
            facecolors="none",
            edgecolors="black",
            linewidths=0.8
        )

axD.text(
    0.98,0.02,
    "Open ring = not strong\nunder total-library scaling",
    transform=axD.transAxes,
    ha="right",
    va="bottom",
    fontsize=6.5
)

panel_label(axD,"D")

# ============================================================
# E — snoR21 candidate-target RNA-seq response
# ============================================================

EE=E.sort_values("log2FC")
y=np.arange(len(EE))

axE.hlines(
    y,
    0,
    EE["log2FC"],
    linewidth=1
)

axE.scatter(
    EE["log2FC"],
    y,
    s=22
)

axE.axvline(0,linewidth=0.8)

axE.set_yticks([])
axE.set_xlabel("mRNA log2 fold change")
axE.set_ylabel("25 high-complementarity targets")
axE.set_title("snoR21 3′ ncRF candidate-target responses")

axE.text(
    0.03,0.95,
    "19/25 (76%) negative\n9 strong DOWN; 2 strong UP",
    transform=axE.transAxes,
    va="top",
    fontsize=8
)

panel_label(axE,"E")

# ============================================================
# F — null tests
# ============================================================

expr = F.loc[
    F["null_type"]=="expression_matched",
    "fraction_negative"
].to_numpy()

seq = F.loc[
    F["null_type"]=="dinucleotide_sequence",
    "fraction_negative"
].to_numpy()

effect = FC["sequence_minus_matched"].to_numpy()

data=[expr,seq,effect]

axF.boxplot(
    data,
    positions=[1,2,3],
    widths=0.55,
    showfliers=False
)

# jittered subsample for visibility
rng=np.random.default_rng(20260917)

for pos,arr in zip([1,2,3],data):
    if len(arr)>1500:
        arr=rng.choice(arr,1500,replace=False)

    jitter=rng.normal(pos,0.045,len(arr))

    axF.scatter(
        jitter,
        arr,
        s=3,
        alpha=0.15
    )

# authentic values
axF.scatter(
    [1,2],
    [0.76,0.76],
    marker="D",
    s=45,
    zorder=10,
    label="Authentic snoR21"
)

axF.scatter(
    [3],
    [0.250772],
    marker="D",
    s=45,
    zorder=10
)

axF.set_xticks([1,2,3])
axF.set_xticklabels(
    [
        "Expression-\nmatched",
        "Dinucleotide\nsequence",
        "Sequence −\nmatched"
    ]
)

axF.set_ylabel("Fraction negative / excess")
axF.set_title("snoR21 sequence-specific inverse association")

axF.text(
    1,0.80,
    "P=0.0079",
    ha="center",
    fontsize=7
)

axF.text(
    2,0.80,
    "P=0.0060",
    ha="center",
    fontsize=7
)

axF.text(
    3,0.30,
    "P=0.038",
    ha="center",
    fontsize=7
)

axF.legend(
    frameon=False,
    loc="lower left"
)

panel_label(axF,"F")

# ------------------------------------------------------------
# Finish
# ------------------------------------------------------------

fig.suptitle(
    "Endosperm-specific remodeling of ncRNA-derived terminal fragments",
    fontsize=14,
    y=0.995
)

pdf=OUT/"Figure1_ncRF_Endosperm_remodeling.pdf"
png=OUT/"Figure1_ncRF_Endosperm_remodeling.png"

fig.savefig(
    pdf,
    bbox_inches="tight"
)

fig.savefig(
    png,
    dpi=600,
    bbox_inches="tight"
)

plt.close(fig)

print("WROTE",pdf)
print("WROTE",png)
