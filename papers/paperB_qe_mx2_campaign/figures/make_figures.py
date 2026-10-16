#!/usr/bin/env python3
"""Figures for the screening article (QE MX2 screening campaign), drawn at print size.

Every value is read from committed measurement files under
 or, where the committed record is a Markdown table
only, transcribed from that table with its path noted beside the numbers.
Widths: single column 88 mm, double column 183 mm (elsarticle 5p layout).
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[2]          # repository root
MM = 1 / 25.4
BLUE, ORANGE, AQUA, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#8a8984"
INK, INK2 = "#0b0b0b", "#52514e"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7.5,
    "axes.labelsize": 7.5, "axes.titlesize": 7.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.8,
    "axes.edgecolor": INK2, "axes.linewidth": 0.6,
    "xtick.color": INK2, "ytick.color": INK2,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "axes.labelcolor": INK, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.dpi": 600, "pdf.fonttype": 42,
})


def last_records(path):
    recs = {}
    for line in open(path):
        r = json.loads(line)
        recs[r["candidate"]] = r
    return recs


# ---------------------------------------------------------------- Fig. 1
def fig_gaps():
    """1T monolayers (first pass), 2H single layers and D3 AA' bilayers
    (corrected 2H structures, stacking_study/results_2088961)."""
    recs = last_records(PROJ / "dft_campaign/results/outcomes.jsonl")
    stack = last_records(PROJ / "stacking_study/results_2088961/outcomes.jsonl")
    metals = ["Ti", "Zr", "Hf", "Mo", "W"]
    chal = ["S", "Se", "Te"]
    comps = [f"{m}{x}2" for m in metals for x in chal]
    fig, ax = plt.subplots(figsize=(183 * MM, 62 * MM))
    ax.axhspan(1.0, 2.0, color="#e6eef9", zorder=0, lw=0)
    ax.text(len(comps) - 0.45, 1.95, "screening window 1.0–2.0 eV",
            ha="right", va="top", fontsize=6.8, color=INK2)
    series = (("1T", recs, "_1T", -0.22, ORANGE, "o"),
              ("2H", stack, "_mono2H", 0.0, BLUE, "s"),
              ("bi", stack, "_bi2H", 0.22, AQUA, "D"))
    for i, c in enumerate(comps):
        for _, src, suf, dx, col, mk in series:
            r = src[c + suf]
            f2 = r["stages"].get("F2", {})
            f3 = r["stages"].get("F3")
            if not f2.get("relax_converged") or not f3:
                ax.plot(i + dx, 0.03, marker="x", color=GREY, ms=4.5, mew=1.0, zorder=3)
                continue
            g = 0.0 if f3.get("metallic") else r["gap_ev"]
            ax.plot(i + dx, g, marker=mk, color=col, ms=4.5 if mk != "D" else 4,
                    mec="white", mew=0.6, zorder=3)
    for k in range(1, len(metals)):
        ax.axvline(3 * k - 0.5, color="#d9d8d4", lw=0.5, zorder=0)
    labels = [c.replace("2", "$_2$") for c in comps]
    ax.set_xticks(range(len(comps)))
    ax.set_xticklabels(labels, rotation=0, fontsize=6.4)
    ax.set_xlim(-0.6, len(comps) - 0.4)
    ax.set_ylim(-0.05, 2.0)
    ax.set_ylabel("PBE band gap (eV)")
    h = [plt.Line2D([], [], ls="", marker="o", color=ORANGE, mec="white", ms=5, label="1T monolayer"),
         plt.Line2D([], [], ls="", marker="s", color=BLUE, mec="white", ms=5, label="2H single layer"),
         plt.Line2D([], [], ls="", marker="D", color=AQUA, mec="white", ms=4.5, label="2H bilayer, D3 spacing"),
         plt.Line2D([], [], ls="", marker="x", color=GREY, ms=5, label="relaxation not converged")]
    ax.legend(handles=h, loc="upper left", frameon=False, ncol=4,
              bbox_to_anchor=(0.0, 1.13), handletextpad=0.3, columnspacing=1.2)
    fig.tight_layout(pad=0.3)
    fig.savefig(HERE / "fig1_gaps.pdf")
    plt.close(fig)


# ---------------------------------------------------------------- Fig. 2
def fig_cost():
    runs = json.load(open(PROJ / "scale_validation/results/scale_fit.json"))["runs"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(183 * MM, 62 * MM))
    for name, col, mk, lab in (("TiS2_1T", ORANGE, "o", "TiS$_2$ 1T"),
                               ("MoSe2_2H", BLUE, "s", "MoSe$_2$ 2H")):
        pts = sorted((r["nelec"], r["wall_s"]) for r in runs
                     if r["name"].startswith(name) and r["name"].endswith("_k6"))
        a.plot([p[0] for p in pts], [p[1] for p in pts], "-", color=col, lw=1.2,
               marker=mk, ms=4.5, mec="white", mew=0.5, label=f"{lab}, 6×6 k-grid")
    s4 = next(r for r in runs if r["name"] == "MoSe2_2H_s4_k3")
    a.plot(s4["nelec"], s4["wall_s"], ls="", marker="s", mfc="none", mec=BLUE, ms=4.5, mew=0.9,
           label="MoSe$_2$ 2H 4×4, 3×3 k-grid")
    # screening-campaign F1 SCFs (6x6 k, 16 ranks)
    recs = last_records(PROJ / "dft_campaign/results/outcomes.jsonl")
    xs = [r["stages"]["F1"]["nelec"] for r in recs.values()]
    ys = [r["stages"]["F1"]["wall_s"] for r in recs.values()]
    a.plot(xs, ys, ls="", marker="o", ms=2.6, color=GREY, alpha=0.8,
           label="campaign coarse SCFs (30)")
    a.set_xscale("log"); a.set_yscale("log")
    a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    a.set_xticks([20, 50, 100, 200, 500, 1000])
    a.set_xticklabels(["20", "50", "100", "200", "500", "1000"])
    a.set_xlabel("valence electrons per cell")
    a.set_ylabel("wall time per SCF (s)")
    a.legend(frameon=False, loc="upper left", handletextpad=0.3)
    a.text(-0.16, 1.02, "(a)", transform=a.transAxes, fontweight="bold")

    # memory: transcribed from memory_band_probe/results/RESULTS.md (re-parsed QE totals)
    probe = [  # candidate, electrons, lean GB, k-parallel GB
        ("TiS2 1T", 600, 18.2, 38.2), ("TiS2 1T", 864, 30.9, 75.5),
        ("MoSe2 2H", 1300, 23.3, 67.6), ("WTe2 2H", 1300, 30.8, 65.9),
        ("MoSe2 2H", 1872, 42.6, 137.2), ("WTe2 2H", 1872, 52.6, 133.4),
        ("MoSe2 2H", 3328, 112.0, 424.4)]
    b.axhline(380, color=INK2, lw=0.8, ls="--")
    b.text(560, 395, "node memory 380 GB", fontsize=6.8, color=INK2, va="bottom")
    b.plot([p[1] for p in probe], [p[2] for p in probe], ls="", marker="o", ms=4.5,
           color=AQUA, mec="white", mew=0.5, label="one k-point pool (32 ranks)")
    b.plot([p[1] for p in probe], [p[3] for p in probe], ls="", marker="D", ms=4.2,
           color=ORANGE, mec="white", mew=0.5, label="four k-point pools (32 ranks)")
    b.annotate("WTe$_2$", xy=(1872, 52.6), xytext=(2150, 62), fontsize=6.4, color=INK2,
               arrowprops=dict(arrowstyle="-", lw=0.5, color=INK2))
    b.annotate("MoSe$_2$", xy=(1872, 42.6), xytext=(2150, 33), fontsize=6.4, color=INK2,
               arrowprops=dict(arrowstyle="-", lw=0.5, color=INK2))
    b.set_xscale("log"); b.set_yscale("log")
    b.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    b.set_xticks([600, 1000, 2000, 3000])
    b.set_xticklabels(["600", "1000", "2000", "3000"])
    b.set_xlabel("valence electrons per cell")
    b.set_ylabel("estimated total memory (GB)")
    b.set_ylim(10, 700)
    b.legend(frameon=False, loc="lower right", handletextpad=0.3)
    b.text(-0.16, 1.02, "(b)", transform=b.transAxes, fontweight="bold")
    fig.tight_layout(pad=0.3, w_pad=2.0)
    fig.savefig(HERE / "fig2_cost_memory.pdf")
    plt.close(fig)


# ---------------------------------------------------------------- Fig. 3
def fig_gain():
    # Measured throughput ratio (sharing arm / no-sharing arm), values from:
    # paired_comparison_1/RESULTS.md (1.005), rescreen_3x3/results/GUEST_RESULTS.md
    # (cost ratio median 0.986, range 0.849-1.191 -> gain = 1/ratio),
    # cofit_probe/results/RESULTS.md (8645 s sequential / 8800 s concurrent),
    # paired_comparison_2/final/RESULTS.md (1.133 as run, 1.207 with F3 priced).
    rows = [
        ("3×3 cells, whole-node charging,\nrank/pool right-sizing", 1.005, None, None),
        ("3×3 cells, per-candidate\nclass-matched pools (29)", 1 / 0.986, 1 / 1.191, 1 / 0.849),
        ("6×6 cells, two jobs\nsharing one node", 8645 / 8800, None, None),
        ("8×8 cells, memory-reservation\ncharging, one memory bit", 1.133, None, None),
    ]
    fig, ax = plt.subplots(figsize=(88 * MM, 64 * MM))
    ax.axvspan(0.8, 1.05, color="#f0efec", lw=0, zorder=0)
    ax.axvline(1.0, color=INK2, lw=0.6, zorder=1)
    ax.axvline(1.05, color=INK2, lw=0.8, ls="--", zorder=1)
    ax.text(1.058, 1.0, "agreed\nthreshold 1.05", fontsize=6.4, color=INK2, va="center")
    for i, (lab, v, lo, hi) in enumerate(rows):
        y = len(rows) - 1 - i
        if lo is not None:
            ax.plot([lo, hi], [y, y], color=BLUE, lw=1.2, solid_capstyle="round", zorder=2)
        ax.plot(v, y, marker="o", ms=5.5, color=BLUE if v < 1.05 else AQUA,
                mec="white", mew=0.6, zorder=3)
        ax.text(v, y + (0.28 if lo is not None else 0.2), f"{v:.3f}", ha="center", va="bottom", fontsize=6.4, color=INK)
    ax.plot(1.207, 0, marker="o", ms=5.5, mfc="none", mec=AQUA, mew=1.0, zorder=3)
    ax.text(1.207, 0.2, "1.207", ha="center", va="bottom", fontsize=6.4)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows[::-1]], fontsize=6.4)
    ax.set_xlim(0.8, 1.25)
    ax.set_ylim(-0.5, len(rows) - 0.3)
    ax.set_xlabel("throughput ratio, sharing / no sharing")
    ax.tick_params(axis="y", length=0)
    fig.tight_layout(pad=0.3)
    fig.savefig(HERE / "fig3_gain.pdf")
    plt.close(fig)


if __name__ == "__main__":
    fig_gaps()
    fig_cost()
    fig_gain()
    print("figures written to", HERE)
