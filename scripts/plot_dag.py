#!/usr/bin/env python3
"""
TEAPOT — nf-core subway style, polished.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from matplotlib.path import Path
import matplotlib.patches as mp
import numpy as np

FIG_W, FIG_H = 26, 10
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
ax.set_xlim(0, FIG_W)
ax.set_ylim(0, FIG_H)
ax.axis("off")
BG = "#e8e8e8"
fig.patch.set_facecolor(BG)
ax.set_facecolor(BG)

C = {
    "green":   "#1B7F2E",
    "blue":    "#1455B5",
    "orange":  "#D95E0A",
    "magenta": "#B0186E",
    "teal":    "#0A7A6E",
}
LW = 12

T = {
    "A": 8.0,
    "B": 6.6,
    "C": 5.2,
    "D": 3.8,
    "E": 2.6,
}

S = [0.8, 6.4, 13.2, 19.0, 25.6]

def seg(x0, x1, y, col, lw=LW, ls="-", z=2):
    ax.plot([x0, x1], [y, y], color=col, lw=lw, ls=ls,
            solid_capstyle="round", zorder=z)

def vseg(x, y0, y1, col, lw=LW, ls="-", z=2):
    ax.plot([x, x], [y0, y1], color=col, lw=lw, ls=ls,
            solid_capstyle="round", zorder=z)

def elbow_right_down(x0, y_start, x1, y_end, col, lw=LW, z=2, r=0.6):
    """Horizontal then curves down/up into vertical."""
    if abs(y_start - y_end) < 0.05:
        seg(x0, x1, y_start, col, lw, z=z); return
    # horizontal to near x1
    seg(x0, x1 - r, y_start, col, lw, z=z)
    # smooth bezier elbow
    verts = [(x1-r, y_start), (x1, y_start), (x1, y_end)]
    codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
    patch = mp.PathPatch(Path(verts, codes), fc="none", ec=col,
                          lw=lw, capstyle="round", joinstyle="round", zorder=z)
    ax.add_patch(patch)

def elbow_down_right(x_start, y0, x1, y_end, col, lw=LW, z=2, r=0.6):
    """Vertical then curves into horizontal."""
    if abs(y0 - y_end) < 0.05:
        seg(x_start, x1, y0, col, lw, z=z); return
    vseg(x_start, y0, y_end + r * np.sign(y_end - y0), col, lw, z=z)
    verts = [(x_start, y_end + r*np.sign(y_end-y0)), (x_start, y_end), (x1, y_end)]
    codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
    patch = mp.PathPatch(Path(verts, codes), fc="none", ec=col,
                          lw=lw, capstyle="round", joinstyle="round", zorder=z)
    ax.add_patch(patch)

def stop(x, y, z=6):
    ax.add_patch(plt.Circle((x, y), 0.24, color="white", zorder=z, ec="#111", lw=2.5))

def lbl(x, y, text, above=True, fs=9, col="#111"):
    dy, va = (0.48, "bottom") if above else (-0.48, "top")
    ax.text(x, y+dy, text, ha="center", va=va, fontsize=fs,
            color=col, fontweight="bold", multialignment="center", zorder=8)

def file_icon(x, y, tag="", col="#333"):
    w, h = 0.44, 0.56
    ax.add_patch(FancyBboxPatch((x-w/2, y-h/2), w, h,
                  boxstyle="round,pad=0.04", fc="white", ec=col, lw=1.8, zorder=9))
    de = 0.13
    ax.fill([x+w/2-de, x+w/2, x+w/2-de],
            [y+h/2, y+h/2, y+h/2-de], color=col, zorder=10)
    if tag:
        ax.text(x, y-0.04, tag, ha="center", va="center",
                fontsize=6.5, color=col, fontweight="bold", zorder=11)

def sbox(xa, xb, num, bg):
    ax.add_patch(FancyBboxPatch((xa+0.1, 0.6), xb-xa-0.2, FIG_H-1.2,
                  boxstyle="round,pad=0.2", fc=bg, ec="#bbb",
                  lw=1.5, alpha=0.38, zorder=0))
    ax.text(xa+0.5, FIG_H-0.7, num, fontsize=22, fontweight="bold",
            color="#444", va="top", ha="left", zorder=1)

# ── stage boxes ───────────────────────────────────────────────────────────────
sbox(S[0], S[1], "1", "#c8dcff")
sbox(S[1], S[2], "2", "#c8ffd8")
sbox(S[2], S[3], "3", "#fff8c0")
sbox(S[3], S[4], "4", "#ecd0ff")

# ────────────────────────────────────────────────────────────────────────────
# STAGE 1
# ────────────────────────────────────────────────────────────────────────────
seg(1.2, S[1], T["A"], C["green"])
stop(2.0, T["A"]); lbl(2.0, T["A"], "ProteomEdge")
stop(4.4, T["A"]); lbl(4.4, T["A"], "Library")

seg(1.2, S[1], T["B"], C["blue"])
stop(2.0, T["B"]); lbl(2.0, T["B"], "ms-experiment", above=False)
stop(4.4, T["B"]); lbl(4.4, T["B"], "SDRF")

seg(1.2, S[1], T["C"], C["orange"])
stop(2.0, T["C"]); lbl(2.0, T["C"], "Spiked Conc.", above=False)
stop(3.8, T["C"]); lbl(3.8, T["C"], "MS Analysis", above=False)
stop(5.8, T["C"]); lbl(5.8, T["C"], ".raw files", above=False)
file_icon(5.7, T["C"]+1.05, "RAW", C["orange"])

# ────────────────────────────────────────────────────────────────────────────
# STAGE 2
# ────────────────────────────────────────────────────────────────────────────
seg(S[1], S[2], T["A"], C["green"])
stop(8.2,  T["A"]); lbl(8.2,  T["A"], "OpenSWATH")
stop(11.0, T["A"]); lbl(11.0, T["A"], "Skyline")

seg(S[1], S[2], T["B"], C["blue"])
stop(9.2, T["B"]); lbl(9.2, T["B"], "EncyclopeDIA")

seg(S[1], S[2], T["C"], C["orange"])
stop(9.2, T["C"]); lbl(9.2, T["C"], "DIA-NN", above=False)

# Library dashed into OpenSWATH (same track A, so just a label dashed)
ax.annotate("", xy=(7.0, T["A"]), xytext=(4.4, T["A"]),
            arrowprops=dict(arrowstyle="-", color=C["green"], lw=2, ls="--"), zorder=3)
ax.annotate("", xy=(7.2, T["B"]), xytext=(4.4, T["B"]),
            arrowprops=dict(arrowstyle="-", color=C["blue"], lw=2, ls="--"), zorder=3)

# ────────────────────────────────────────────────────────────────────────────
# STAGE 3
# ────────────────────────────────────────────────────────────────────────────
seg(S[2], S[3], T["A"], C["green"])
stop(14.5, T["A"]); lbl(14.5, T["A"], "Context")
stop(17.2, T["A"]); lbl(17.2, T["A"], "DIAthem")

seg(S[2], S[3], T["B"], C["blue"])
stop(14.5, T["B"]); lbl(14.5, T["B"], "Downstream A.")
stop(17.2, T["B"]); lbl(17.2, T["B"], "Batch corr.")

seg(S[2], S[3], T["C"], C["orange"])
stop(14.5, T["C"]); lbl(14.5, T["C"], "QC", above=False)
stop(17.2, T["C"]); lbl(17.2, T["C"], "File Conversion", above=False)

# pmultiQC branch (magenta): drops from A at x=13.5, runs horizontal at T["D"]
vseg(13.5, T["A"], T["D"], C["magenta"], lw=LW-3, z=2)
seg(13.5, S[3], T["D"], C["magenta"], lw=LW-3)
stop(16.0, T["D"]); lbl(16.0, T["D"], "pmultiQC", above=False)

# thin vertical tap from B to pmultiQC branch
vseg(13.5, T["B"], T["A"], C["magenta"], lw=4, z=3)

# ────────────────────────────────────────────────────────────────────────────
# STAGE 4 – all lines merge into kardemumma
# ────────────────────────────────────────────────────────────────────────────
KX = 21.0
RX = 24.0

# smooth elbows converge at (KX, T["E"])
for yt, col, lw_use in [
    (T["A"], C["green"],   LW),
    (T["B"], C["blue"],    LW),
    (T["C"], C["orange"],  LW),
    (T["D"], C["magenta"], LW-3),
]:
    elbow_right_down(S[3], yt, KX, T["E"], col, lw=lw_use, z=2, r=0.8)

# SDRF dashed annotation arc to kardemumma
ax.annotate("", xy=(KX-0.3, T["E"]+0.28), xytext=(4.4, T["B"]),
            arrowprops=dict(arrowstyle="-", color=C["blue"], lw=2, ls="--",
                            connectionstyle="arc3,rad=0.3"), zorder=3)

# kardemumma → result
seg(KX, RX+0.6, T["E"], C["teal"], lw=LW+3)
stop(KX, T["E"]); lbl(KX, T["E"], "kardemumma", fs=10.5)
stop(RX, T["E"]); lbl(RX, T["E"], "Result matrix", fs=10.5)
file_icon(RX+0.95, T["E"]+0.88, "TSV", C["teal"])

# ── title ──────────────────────────────────────────────────────────────────
ax.text(FIG_W/2, FIG_H-0.12, "TEAPOT — Pipeline DAG",
        ha="center", va="top", fontsize=17, fontweight="bold", color="#111")

# ── stage footers ──────────────────────────────────────────────────────────
for txt, xm in [
    ("Inputs & metadata",    (S[0]+S[1])/2),
    ("DIA search tools",     (S[1]+S[2])/2),
    ("Post-processing",      (S[2]+S[3])/2),
    ("Aggregation & output", (S[3]+S[4])/2),
]:
    ax.text(xm, 0.82, txt, ha="center", va="bottom",
            fontsize=9, color="#555", style="italic")

# ── legend ──────────────────────────────────────────────────────────────────
items = [
    mpatches.Patch(color=C["green"],   label="ProteomEdge / Library → OpenSWATH / Skyline / Context / DIAthem"),
    mpatches.Patch(color=C["blue"],    label="ms-experiment / SDRF → EncyclopeDIA → Downstream A. / Batch corr."),
    mpatches.Patch(color=C["orange"],  label="Spiked Conc. / MS Analysis → DIA-NN → QC / File Conversion"),
    mpatches.Patch(color=C["magenta"], label="pmultiQC branch"),
    mpatches.Patch(color=C["teal"],    label="kardemumma → Result matrix"),
]
leg = ax.legend(handles=items, loc="lower left", bbox_to_anchor=(0.005, 0.08),
                fontsize=8, framealpha=0.95, edgecolor="#ccc",
                title="METHOD", title_fontsize=8.5)

plt.tight_layout(pad=0.3)
for ext in ("png", "svg"):
    out = f"img/workflow/dag_teapot_nfcore.{ext}"
    plt.savefig(out, dpi=160, bbox_inches="tight", facecolor=BG)
    print(f"Saved {out}")
plt.close()