#!/usr/bin/env python3
"""Figures for the screening article (QE MX2 screening campaign), drawn at print size.

Every value is read from committed measurement files under
 or, where the committed record is a Markdown table
only, transcribed from that table with its path noted beside the numbers.
Widths: single column 88 mm, double column 183 mm (elsarticle 5p layout).

Figures (in order of appearance in the article):
  fig_method.pdf       workflow: fidelity ladder and paired comparisons (Methods)
  fig1_gaps.pdf        PBE and SOC gaps of the 30 candidates
  fig_stage_cost.pdf   primitive-cell cost per call by stage and cell
  fig2_cost_memory.pdf cost scaling with cell size and memory at large cells
  fig3_gain.pdf        every measured throughput ratio, with paired-bootstrap intervals
  fig_group_leak.pdf   research-group wins and timing leakage per comparison
"""
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
from matplotlib.patches import FancyBboxPatch

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[2]          # repository root
MM = 1 / 25.4
BLUE, ORANGE, AQUA, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#8a8984"
INK, INK2 = "#0b0b0b", "#52514e"
BROKER, PARTIAL = "#3b3a8f", "#9b9ac9"   # arrangements (kept apart from polytype colours)

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

B = 10000                       # bootstrap replicates
RNG_SEED = 0                    # same resampling seed as the companion disclosure article


def last_records(path):
    recs = {}
    for line in open(path):
        r = json.loads(line)
        recs[r["candidate"]] = r
    return recs


def chem(name):
    """'MoS2_2H' -> 'MoS$_2$-2H'."""
    base, _, poly = name.partition("_")
    return base.replace("2", "$_2$") + ("-" + poly if poly else "")


def letter(ax, s, x=-0.14, y=1.03):
    ax.text(x, y, s, transform=ax.transAxes, fontweight="bold", fontsize=7.5,
            va="bottom", ha="left")


def wilson(k, n, z=1.96):
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def boot_ratio(lab_s, ch_s, lab_d, ch_d, seed=RNG_SEED):
    """Paired percentile bootstrap (resampling seeds) of
    (sum labels_s / sum charge_s) / (sum labels_d / sum charge_d)."""
    rng = np.random.default_rng(seed)
    n = len(lab_s)
    idx = rng.integers(0, n, size=(B, n))
    r = (lab_s[idx].sum(1) / ch_s[idx].sum(1)) / (lab_d[idx].sum(1) / ch_d[idx].sum(1))
    return np.percentile(r, [2.5, 97.5])


# ---------------------------------------------------------------- workflow figure
def _box(ax, x, y, w, h, title, body, fc="#f4f3f0", ec=INK2, tsize=6.8, bsize=6.2):
    """Rounded box with a bold title and a body; offsets in points so the text
    position does not depend on the figure height."""
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.003,rounding_size=0.01",
                                fc=fc, ec=ec, lw=0.7))
    ax.annotate(title, (x + w / 2, y + h), xytext=(0, -3.5), textcoords="offset points",
                ha="center", va="top", fontsize=tsize, fontweight="bold", color=INK,
                linespacing=1.1)
    nt = title.count("\n") + 1
    ax.annotate(body, (x + w / 2, y + h), xytext=(0, -5.5 - 9.0 * nt), textcoords="offset points",
                ha="center", va="top", fontsize=bsize, color=INK2, linespacing=1.22)


def _arrow(ax, x0, y0, x1, y1, color=INK2, ls="-"):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", lw=0.8, color=color, mutation_scale=7, ls=ls))


def fig_method():
    fig = plt.figure(figsize=(183 * MM, 104 * MM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    # ---- (a) screening ladder
    ax.text(0.004, 0.992, "(a)", fontweight="bold", fontsize=7.5, va="top")
    ax.text(0.034, 0.992, "Screening ladder: every candidate climbs the same stages, and every call is timed",
            fontsize=7.0, va="top")
    y, h, w, gap = 0.705, 0.245, 0.152, 0.0148
    xs = [0.004 + i * (w + gap) for i in range(6)]
    boxes = [
        ("Candidates", "30 MX$_2$ layers:\nM = Ti, Zr, Hf,\nMo, W;\nX = S, Se, Te;\n1T and 2H", "#ffffff"),
        ("Coarse SCF", "idealized geometry,\n6×6 k-points\n\nmeasured:\nwall time, memory", "#f4f3f0"),
        ("Relaxation and\nproduction SCF", "lattice and ions free,\n12×12 k-points\n\nmeasured: lattice\nconstant, energy", "#f4f3f0"),
        ("Band calculation", "18×18 k-points\n\n\nmeasured:\nPBE gap against\n1.0–2.0 eV window", "#f4f3f0"),
        ("Γ phonons", "candidates in\nthe window\n\nmeasured:\nlowest frequency", "#f4f3f0"),
        ("Spin–orbit\ncoupling", "candidates in\nthe window\n\nmeasured:\ngap with SOC", "#f4f3f0"),
    ]
    for (t, b, fc), x in zip(boxes, xs):
        _box(ax, x, y, w, h, t, b, fc=fc)
    for i in range(5):
        _arrow(ax, xs[i] + w + 0.001, y + h / 2, xs[i + 1] - 0.001, y + h / 2)
    # follow-up branches, below the ladder
    yb2, hb2 = 0.565, 0.09
    bx1, bw1 = xs[1], 2 * w + gap
    bx2, bw2 = xs[3] + w / 2, 2.5 * w + 2 * gap - w / 2 + w / 2
    ax.add_patch(FancyBboxPatch((bx1, yb2), bw1, hb2, boxstyle="round,pad=0.003,rounding_size=0.01",
                                fc="#eef3fb", ec=INK2, lw=0.6))
    ax.text(bx1 + bw1 / 2, yb2 + hb2 / 2, "Corrected 2H structures, same ladder:\nsingle layers; "
            "bilayers with D3 (AA$'$, 3R)", ha="center", va="center", fontsize=6.2, color=INK2)
    ax.add_patch(FancyBboxPatch((bx2, yb2), bw2, hb2, boxstyle="round,pad=0.003,rounding_size=0.01",
                                fc="#eef3fb", ec=INK2, lw=0.6))
    ax.text(bx2 + bw2 / 2, yb2 + hb2 / 2, "Chalcogen-vacancy hops: climbing-image NEB\nin 3×3 cells "
            "(selected and Janus layers); measured: barrier", ha="center", va="center",
            fontsize=6.2, color=INK2)
    _arrow(ax, xs[2] + w / 2, y - 0.002, xs[2] + w / 2, yb2 + hb2 + 0.002, color=GREY, ls=(0, (3, 2)))
    _arrow(ax, xs[3] + w / 2 + 0.06, y - 0.002, xs[3] + w / 2 + 0.06, yb2 + hb2 + 0.002,
           color=GREY, ls=(0, (3, 2)))
    ax.text(0.004, 0.515, "Every call is recorded with wall time, MPI ranks, k-point pools and peak "
            "memory; node-hours = wall time × ranks / 64.\nSupercells (2×2 to 8×8) of the relaxed cells "
            "give the cost series.", fontsize=6.2, color=INK2, va="top", linespacing=1.3)

    # ---- (b) paired comparison
    ax.text(0.004, 0.445, "(b)", fontweight="bold", fontsize=7.5, va="top")
    ax.text(0.034, 0.445, "Paired comparison: each seed fixes one shuffled candidate order, run under "
            "both arrangements on the same physics", fontsize=7.0, va="top")
    yl, yh, hb = 0.035, 0.215, 0.165
    _box(ax, 0.004, 0.095, 0.15, 0.225, "Campaign seed",
         "shuffled order of\nthe relaxed\ncompounds; the\ngroup climbs the\nladder in that order",
         fc="#ffffff")
    _box(ax, 0.205, yh, 0.335, hb, "No sharing",
         "the centre sees the calculation class only\nand reserves worst-case memory\nfor every job",
         fc="#f4f3f0")
    _box(ax, 0.205, yl, 0.335, hb, "Broker",
         "the group discloses one bit per job (memory\nabove or below the pool median); low-memory\njobs share a node",
         fc="#ecebf6", ec=BROKER)
    _arrow(ax, 0.156, 0.21, 0.203, yh + hb / 2)
    _arrow(ax, 0.156, 0.21, 0.203, yl + hb / 2)
    _box(ax, 0.59, 0.035, 0.195, 0.345, "Measured outcomes",
         "computing centre:\nvalid calculations per\ncharged node-hour\n(criterion ≥ 1.05×)\n\n"
         "research group:\npaired win rate to first\nconfirmed candidate\n(criterion ≥ 0.60)",
         fc="#ffffff")
    _box(ax, 0.83, 0.035, 0.166, 0.345, "Leakage audit",
         "from opaque job IDs\nand completion times\nonly:\n\ncompound identity\naccuracy and\nmembership AUC",
         fc="#ffffff")
    _arrow(ax, 0.542, yh + hb / 2, 0.588, yh + hb / 2)
    _arrow(ax, 0.542, yl + hb / 2, 0.588, yl + hb / 2)
    _arrow(ax, 0.787, 0.21, 0.828, 0.21)
    fig.savefig(HERE / "fig_method.pdf")
    plt.close(fig)


# ---------------------------------------------------------------- Fig. 1
def fig_gaps():
    """1T monolayers (first pass), 2H single layers and D3 AA' bilayers
    (corrected 2H structures, stacking_study/results_2088961), the follow-up
    relaxation of the 2H TiS2 single layer (stacking_study/followup_2026-10-12),
    the first-pass PBE-spaced 2H bilayers, and the SOC gaps of the finalists."""
    recs = last_records(PROJ / "dft_campaign/results/outcomes.jsonl")
    stack = last_records(PROJ / "stacking_study/results_2088961/outcomes.jsonl")
    fu = last_records(PROJ / "stacking_study/followup_2026-10-12/results/relax/results_16208957/outcomes.jsonl")
    for k, r in fu.items():               # converged follow-up replaces the failed relaxation
        if r["stages"].get("F2", {}).get("relax_converged") and r["stages"].get("F3"):
            stack[k] = r
    soc1 = last_records(PROJ / "dft_campaign/soc/results/soc_outcomes.jsonl")
    soc2 = last_records(PROJ / "stacking_study/followup_2026-10-12/results/socph/results_16208958/soc_outcomes.jsonl")
    soc = {"_1T": soc1, "_mono2H": soc2}
    metals = ["Ti", "Zr", "Hf", "Mo", "W"]
    chal = ["S", "Se", "Te"]
    comps = [f"{m}{x}2" for m in metals for x in chal]
    METAL_Y = -0.22
    fig, ax = plt.subplots(figsize=(183 * MM, 66 * MM))
    ax.axhspan(1.0, 2.0, color="#e6eef9", zorder=0, lw=0)
    ax.axhspan(METAL_Y - 0.12, METAL_Y + 0.12, color="#f2f1ee", zorder=0, lw=0)
    ax.axhline(0, color="#d9d8d4", lw=0.5, zorder=0)
    ax.text(-0.5, 1.97, "screening window 1.0–2.0 eV",
            ha="left", va="top", fontsize=6.8, color=INK2)
    unrelaxed = False
    series = (("1T", recs, "_1T", -0.24, ORANGE, "o"),
              ("2H", stack, "_mono2H", 0.0, BLUE, "s"),
              ("bi", stack, "_bi2H", 0.24, AQUA, "D"))
    for i, c in enumerate(comps):
        for _, src, suf, dx, col, mk in series:
            r = src[c + suf]
            f2 = r["stages"].get("F2", {})
            f3 = r["stages"].get("F3")
            if not f2.get("relax_converged") or not f3:
                ax.plot(i + dx, METAL_Y, marker="x", color=GREY, ms=4.5, mew=1.0, zorder=3)
                unrelaxed = True
                continue
            y = METAL_Y if f3.get("metallic") else r["gap_ev"]
            ax.plot(i + dx, y, marker=mk, color=col, ms=4.6 if mk != "D" else 4.1,
                    mec="white", mew=0.6, zorder=3)
            s = soc.get(suf, {}).get(c + suf)
            if s is not None and not s.get("soc_metallic"):
                ax.plot(i + dx, s["soc_gap_ev"], marker=mk, mfc="none", mec=col,
                        ms=4.6, mew=0.9, zorder=4)
                ax.plot([i + dx, i + dx], [s["soc_gap_ev"], y], color=col, lw=0.5, zorder=2)
        # first-pass PBE-spaced 2H bilayer (6-atom cell, no dispersion)
        r = recs[c + "_2H"]
        if r["stages"].get("F2", {}).get("relax_converged") and r["stages"].get("F3"):
            g = METAL_Y if r["stages"]["F3"].get("metallic") else r["gap_ev"]
            ax.plot(i + 0.11, g, marker="_", color=INK2, ms=5, mew=1.0, zorder=3)
    for k in range(1, len(metals)):
        ax.axvline(3 * k - 0.5, color="#d9d8d4", lw=0.5, zorder=0)
    labels = [c.replace("2", "$_2$") for c in comps]
    ax.set_xticks(range(len(comps)))
    ax.set_xticklabels(labels, rotation=0, fontsize=6.8)
    ax.set_xlim(-0.6, len(comps) - 0.4)
    ax.set_ylim(METAL_Y - 0.13, 2.0)
    ax.set_yticks([METAL_Y, 0, 0.5, 1.0, 1.5, 2.0])
    ax.set_yticklabels(["metal", "0", "0.5", "1.0", "1.5", "2.0"])
    ax.set_ylabel("PBE band gap (eV)")
    h = [plt.Line2D([], [], ls="", marker="o", color=ORANGE, mec="white", ms=5, label="1T monolayer"),
         plt.Line2D([], [], ls="", marker="s", color=BLUE, mec="white", ms=5, label="2H single layer"),
         plt.Line2D([], [], ls="", marker="D", color=AQUA, mec="white", ms=4.5, label="2H bilayer, D3 spacing"),
         plt.Line2D([], [], ls="", marker="_", color=INK2, ms=6, mew=1.0, label="2H bilayer, plain PBE (first pass)"),
         plt.Line2D([], [], ls="", marker="s", mfc="none", mec=INK2, ms=5, mew=0.9, label="open: with spin–orbit")]
    if unrelaxed:
        h.append(plt.Line2D([], [], ls="", marker="x", color=GREY, ms=5, label="not relaxed"))
    ax.legend(handles=h, loc="upper left", frameon=False, ncol=len(h),
              bbox_to_anchor=(-0.01, 1.14), handletextpad=0.2, columnspacing=0.8, fontsize=6.5)
    fig.tight_layout(pad=0.3)
    fig.savefig(HERE / "fig1_gaps.pdf")
    plt.close(fig)


# ---------------------------------------------------------------- stage cost
def fig_stage_cost():
    """Per-call wall time by stage and polytype (primitive-cell campaign),
    and per-candidate cost against the PBE gap."""
    recs = last_records(PROJ / "dft_campaign/results/outcomes.jsonl")
    fig, (a, b) = plt.subplots(1, 2, figsize=(183 * MM, 58 * MM),
                               gridspec_kw={"width_ratios": [1.25, 1]})
    stages = [("F1", "coarse SCF"), ("F2", "relaxation\n+ SCF"), ("F3", "band\ncalculation"),
              ("F4", "Γ phonons")]
    rng = np.random.default_rng(3)
    for j, (st, _) in enumerate(stages):
        for poly, dx, col, mk in (("1T", -0.17, ORANGE, "o"), ("2H", 0.17, BLUE, "s")):
            w = [r["stages"][st]["wall_s"] for r in recs.values()
                 if r["phase"] == poly and st in r["stages"] and r["stages"][st].get("wall_s")]
            if not w:
                continue
            x = j + dx + rng.uniform(-0.07, 0.07, len(w))
            a.plot(x, w, ls="", marker=mk, ms=3.0, color=col, alpha=0.85, mec="white", mew=0.3)
            m = float(np.median(w))
            a.plot([j + dx - 0.13, j + dx + 0.13], [m, m], color=INK, lw=1.0)
            a.text(j + dx, min(w) * 0.8, f"n={len(w)}", ha="center", va="top", fontsize=6.2,
                   color=INK2)
    a.set_yscale("log")
    a.set_xticks(range(len(stages)))
    a.set_xticklabels([s[1] for s in stages])
    a.set_ylabel("wall time per call (s, 16 ranks)")
    a.set_ylim(1.5, 3000)
    a.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    h = [plt.Line2D([], [], ls="", marker="o", color=ORANGE, ms=4, label="1T (3 atoms)"),
         plt.Line2D([], [], ls="", marker="s", color=BLUE, ms=4, label="2H bilayer (6 atoms)"),
         plt.Line2D([], [], color=INK, lw=1.0, label="median")]
    a.legend(handles=h, frameon=False, loc="upper left", ncol=3, handletextpad=0.3,
             columnspacing=0.9)
    letter(a, "(a)", x=-0.12)

    for poly, col, mk in (("1T", ORANGE, "o"), ("2H", BLUE, "s")):
        xs, ys = [], []
        for r in recs.values():
            if r["phase"] != poly:
                continue
            nh = sum(s.get("node_hours", 0) or 0 for s in r["stages"].values())
            f3 = r["stages"].get("F3")
            if not f3:
                continue
            xs.append(0.0 if f3.get("metallic") else r["gap_ev"])
            ys.append(nh)
        b.plot(xs, ys, ls="", marker=mk, ms=4, color=col, mec="white", mew=0.5)
    b.axvspan(1.0, 2.0, color="#e6eef9", zorder=0, lw=0)
    b.set_yscale("log")
    b.set_xlabel("PBE gap (eV; 0 = metal)")
    b.set_ylabel("node-hours per candidate\n(all stages)")
    b.set_ylim(0.0025, 0.08)
    b.set_yticks([0.003, 0.01, 0.03])
    b.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    b.text(1.5, 0.07, "screening window", ha="center", va="top", fontsize=6.4, color=INK2)
    b.set_xlim(-0.08, 1.9)
    letter(b, "(b)", x=-0.17)
    fig.tight_layout(pad=0.3, w_pad=2.0)
    fig.savefig(HERE / "fig_stage_cost.pdf")
    plt.close(fig)


# ---------------------------------------------------------------- Fig. 2
def fig_cost():
    runs = json.load(open(PROJ / "scale_validation/results/scale_fit.json"))["runs"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(183 * MM, 66 * MM))
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
           label="campaign coarse SCFs (30, 16 ranks)")
    # S vacancy and heterobilayer: measured against the electron-count law
    # fitted to the Se- and Te-rich cells (contention_transfer_bench/results/bench_analysis.json)
    ba = json.load(open(PROJ / "contention_transfer_bench/results/bench_analysis.json"))["transfer"]
    br = {}
    for line in open(PROJ / "contention_transfer_bench/results/bench_results.jsonl"):
        r = json.loads(line)
        if r.get("kind") == "transfer":
            br[r["name"]] = r
    first = True
    for name, d in ba.items():
        if name not in br:
            continue
        ne, w = br[name]["nelec"], br[name]["wall_s"]
        a.plot([ne, ne], [w, d["mx2_law_pred_s"]], color=INK2, lw=0.6, ls=":")
        a.plot(ne, w, ls="", marker="^", ms=4.8, color=INK, mec="white", mew=0.4,
               label="S vacancy / MoS$_2$–WS$_2$ (measured)" if first else None)
        a.plot(ne, d["mx2_law_pred_s"], ls="", marker="^", ms=4.8, mfc="none", mec=INK,
               mew=0.7, label="same, size-only prediction" if first else None)
        first = False
    # guide lines N^2 and N^3 through the 3x3 MoSe2 point
    x0, y0 = 468, next(r["wall_s"] for r in runs if r["name"] == "MoSe2_2H_s3_k6")
    gx = np.array([150, 1000])
    for p, ls in ((2, (0, (4, 2))), (3, (0, (1, 1.5)))):
        a.plot(gx, y0 * (gx / x0) ** p, color="#b8b7b2", lw=0.8, ls=ls, zorder=0)
        a.text(gx[0] * 0.96, y0 * (gx[0] / x0) ** p, f"∝ N$^{p}$", fontsize=6.4, color=INK2,
               ha="right", va="center")
    a.set_xscale("log"); a.set_yscale("log")
    a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    a.set_xticks([20, 50, 100, 200, 500, 1000])
    a.set_xticklabels(["20", "50", "100", "200", "500", "1000"])
    a.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    a.set_ylim(1.5, 3000)
    a.set_xlabel("valence electrons per cell, N")
    a.set_ylabel("wall time per SCF (s, 32 ranks)")
    a.legend(frameon=False, loc="upper left", handletextpad=0.3, fontsize=6.3, labelspacing=0.3)
    letter(a, "(a)", x=-0.16)

    # memory: transcribed from memory_band_probe/results/RESULTS.md (re-parsed QE estimates)
    probe = [  # candidate, electrons, one pool GB, four pools GB
        ("TiS2", 600, 18.2, 38.2), ("TiS2", 864, 30.9, 75.5),
        ("MoSe2", 1300, 23.3, 67.6), ("WTe2", 1300, 30.8, 65.9),
        ("MoSe2", 1872, 42.6, 137.2), ("WTe2", 1872, 52.6, 133.4),
        ("MoSe2", 3328, 112.0, 424.4)]
    style = {"TiS2": ("o", ORANGE, "TiS$_2$ 1T"), "MoSe2": ("s", BLUE, "MoSe$_2$ 2H"),
             "WTe2": ("D", "#0f4c9c", "WTe$_2$ 2H")}
    b.axhline(380, color=INK2, lw=0.8, ls="--")
    b.text(560, 395, "node memory 380 GB", fontsize=6.8, color=INK2, va="bottom")
    shift = {"TiS2": 1.0, "MoSe2": 0.965, "WTe2": 1.035}   # separate equal electron counts
    for comp, (mk, col, lab) in style.items():
        pts = [p for p in probe if p[0] == comp]
        xx = [p[1] * shift[comp] for p in pts]
        b.plot(xx, [p[2] for p in pts], ls="", marker=mk, ms=4.6,
               color=col, mec="white", mew=0.5)
        b.plot(xx, [p[3] for p in pts], ls="", marker=mk, ms=4.6,
               mfc="white", mec=col, mew=0.9)
    hl = [plt.Line2D([], [], ls="", marker=style[c][0], color=style[c][1], ms=4.6,
                     label=style[c][2]) for c in style]
    hl += [plt.Line2D([], [], ls="", marker="o", color=INK2, ms=4.6, label="one k-point pool"),
           plt.Line2D([], [], ls="", marker="o", mfc="white", mec=INK2, mew=0.9, ms=4.6,
                      label="four k-point pools")]
    b.legend(handles=hl, frameon=False, loc="lower right", handletextpad=0.3, fontsize=6.3,
             labelspacing=0.3)
    b.set_xscale("log"); b.set_yscale("log")
    b.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    b.set_xticks([600, 1000, 2000, 3000])
    b.set_xticklabels(["600", "1000", "2000", "3000"])
    b.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    b.set_xlabel("valence electrons per cell (5×5 to 8×8 supercells)")
    b.set_ylabel("estimated total memory (GB, 32 ranks)")
    b.set_ylim(10, 700)
    letter(b, "(b)", x=-0.16)
    fig.tight_layout(pad=0.3, w_pad=2.0)
    fig.savefig(HERE / "fig2_cost_memory.pdf")
    plt.close(fig)


# ---------------------------------------------------------------- paired data
def no1_arrays():
    recs = [json.loads(l) for l in open(PROJ / "paired_comparison_1/results/campaigns.jsonl")]
    out = {}
    for arm in ("DTDE", "CEDT_P4"):
        rs = sorted([r for r in recs if r["arm"] == arm], key=lambda r: r["seed"])
        out[arm] = (np.array([r["labels_f2"] for r in rs], float),
                    np.array([r["nh"] for r in rs], float),
                    [r["nh_first_qual"] for r in rs])
    return out


def no2_arrays(key):
    rep = json.load(open(PROJ / "paired_comparison_2/f3_repair/results/repair_graded.json"))
    out = {}
    for arm in ("DTDE", "CEDT_P4"):
        rs = sorted([c for c in rep["campaigns"] if c["arm"] == arm], key=lambda c: c["seed"])
        out[arm] = (np.array([c["labels_f2"] for c in rs], float),
                    np.array([c[key] for c in rs], float))
    return out, rep


def no3_sets():
    """Complete collections of comparisons No. 3 and No. 3b (14 Oct 2026),
    read as in paired_comparison_3/final_grade_2026-10-14.py."""
    out = {}
    for key, folder in (("No3", PROJ / "paired_comparison_3"), ("No3b", PROJ / "paired_comparison_3b")):
        camps = {}
        for d in sorted((folder / "collected_2026-10-14").iterdir()):
            f = d / "campaigns.jsonl"
            if f.exists():
                rec = json.loads(f.read_text().splitlines()[-1])
                camps[(rec["arm"], rec["seed"])] = rec
        sites = json.load(open(folder / "pair_sites.json"))["campaigns"]
        out[key] = (camps, sites)
    return out


def paired(camps, sites, arm, same_cpu=False):
    seeds = [s for s in range(15) if (arm, s) in camps and ("DTDE", s) in camps]
    if same_cpu:
        seeds = [s for s in seeds if sites[f"{arm}_s{s}"]["cpu"] == sites[f"DTDE_s{s}"]["cpu"]]
    g = lambda a, k: np.array([camps[(a, s)][k] for s in seeds], float)
    return seeds, g(arm, "labels_f2"), g(arm, "billed"), g("DTDE", "labels_f2"), g("DTDE", "billed")


def wins(camps, arm, seeds):
    w = 0
    for s in seeds:
        x, y = camps[(arm, s)]["billed_first_qual"], camps[("DTDE", s)]["billed_first_qual"]
        if x is not None and (y is None or x < y):
            w += 1
    return w


# ---------------------------------------------------------------- Fig. 3
def fig_gain():
    rows = []   # (label, point, lo, hi, kind, n_text, extra)
    n1 = no1_arrays()
    (ls_, cs_, _), (ld, cd, _) = n1["CEDT_P4"], n1["DTDE"]
    r1 = (ls_.sum() / cs_.sum()) / (ld.sum() / cd.sum())
    rows.append(("3×3, settings matched to class\n(per-core node-hours)", r1, *boot_ratio(ls_, cs_, ld, cd),
                 "broker", "10 pairs", None))
    # 3x3 rescreen, per-compound cost ratio median 0.986, range 0.849-1.191 over 29
    # compounds (rescreen_3x3/results/GUEST_RESULTS.md); throughput = 1 / cost ratio
    rows.append(("3×3, settings matched\nper compound (rescreen)", 1 / 0.986, 1 / 1.191, 1 / 0.849,
                 "range", "29 compounds", None))
    # co-residency probe: 8645 s sequential / 8800 s concurrent (cofit_probe/results/RESULTS.md)
    rows.append(("6×6, two jobs\nsharing one node", 8645 / 8800, None, None, "single", "1 run", None))
    arr, rep = no2_arrays("billed_asrun")
    (ls_, cs_), (ld, cd) = arr["CEDT_P4"], arr["DTDE"]
    lo, hi = boot_ratio(ls_, cs_, ld, cd)
    arrc, _ = no2_arrays("billed_corrected")
    (lsc, csc), (ldc, cdc) = arrc["CEDT_P4"], arrc["DTDE"]
    rows.append(("8×8, memory-reservation charge,\nmemory class disclosed", rep["h_ratio_asrun"], lo, hi,
                 "broker", "15 pairs", (rep["h_ratio_corrected"], *boot_ratio(lsc, csc, ldc, cdc))))
    sets = no3_sets()
    for key, arm, lab in (("No3", "CEDT_P4", "8×8, guest nodes,\nmemory class disclosed"),
                          ("No3", "CTDE_P2", "8×8, guest nodes,\nshared logbook"),
                          ("No3", "MIX50", "8×8, guest nodes,\nhalf of jobs via broker"),
                          ("No3b", "CEDT_P4", "8×8, guest nodes, 2H single\nlayers, memory class disclosed")):
        camps, sites = sets[key]
        seeds, a1, b1, a0, b0 = paired(camps, sites, arm)
        r = (a1.sum() / b1.sum()) / (a0.sum() / b0.sum())
        lo, hi = boot_ratio(a1, b1, a0, b0)
        ss, c1, d1, c0, d0 = paired(camps, sites, arm, same_cpu=True)
        rs = (c1.sum() / d1.sum()) / (c0.sum() / d0.sum())
        rows.append((lab, r, lo, hi, "broker" if arm == "CEDT_P4" else "partial",
                     f"{len(seeds)} pairs ({len(ss)})", ("same", rs)))
    groups = [(0, 3, "dedicated\nnode"), (3, 4, ""), (4, 8, "guest\nnodes")]
    fig, ax = plt.subplots(figsize=(183 * MM, 84 * MM))
    ax.axvspan(0.6, 1.05, color="#f0efec", lw=0, zorder=0)
    ax.axvline(1.0, color=INK2, lw=0.6, zorder=1)
    ax.axvline(1.05, color=INK2, lw=0.8, ls="--", zorder=1)
    n = len(rows)
    ax.text(1.053, n - 0.35, "agreed threshold 1.05", fontsize=6.4, color=INK2, va="center")
    for i, (lab, v, lo, hi, kind, ntxt, extra) in enumerate(rows):
        y = n - 1 - i
        col = {"broker": BROKER, "partial": PARTIAL, "range": GREY, "single": BROKER}[kind]
        if kind == "range":
            ax.plot([lo, hi], [y, y], color=GREY, lw=0.7, zorder=2)
            ax.plot([lo, hi], [y, y], ls="", marker="|", color=GREY, ms=5, zorder=2)
        elif lo is not None:
            ax.plot([lo, hi], [y, y], color=col, lw=1.6, solid_capstyle="butt", zorder=2)
        ax.plot(v, y, marker="o", ms=5.5, color=col, mec="white", mew=0.6, zorder=3)
        ax.text(v, y + 0.2, f"{v:.3f}", ha="center", va="bottom", fontsize=6.3, color=INK)
        if extra and extra[0] == "same":
            ax.plot(extra[1], y, marker="^", ms=4.6, mfc="white", mec=col, mew=0.9, zorder=3)
        elif extra:
            ax.plot([extra[1], extra[2]], [y - 0.22, y - 0.22], color=col, lw=0.8, alpha=0.6)
            ax.plot(extra[0], y - 0.22, marker="o", ms=4.2, mfc="white", mec=col, mew=0.9, zorder=3)
            ax.text(extra[2] + 0.006, y - 0.22, f"{extra[0]:.3f}", ha="left", va="center",
                    fontsize=6.3, color=INK)
        ax.text(1.495, y, ntxt, ha="right", va="center", fontsize=6.3, color=INK2)
    ax.set_yticks(range(n))
    ax.set_yticklabels([r[0] for r in rows[::-1]], fontsize=6.6)
    ax.set_xlim(0.6, 1.50)
    ax.set_xticks(np.arange(0.6, 1.41, 0.1))
    ax.set_ylim(-0.6, n - 0.1)
    ax.set_xlabel("throughput ratio, sharing / no sharing (valid calculations per charged node-hour)")
    ax.tick_params(axis="y", length=0)
    ax.axhline(n - 3.5, color="#d9d8d4", lw=0.5)
    ax.axhline(n - 4.5, color="#d9d8d4", lw=0.5)
    h = [plt.Line2D([], [], color=BROKER, lw=1.6, marker="o", ms=5, mec="white",
                    label="broker, paired-bootstrap 95% interval"),
         plt.Line2D([], [], color=PARTIAL, lw=1.6, marker="o", ms=5, mec="white",
                    label="partial arrangement"),
         plt.Line2D([], [], color=GREY, lw=0.7, marker="|", ms=5, label="range over compounds"),
         plt.Line2D([], [], ls="", marker="^", mfc="white", mec=INK2, ms=4.6,
                    label="pairs on one processor model"),
         plt.Line2D([], [], ls="", marker="o", mfc="white", mec=INK2, ms=4.4,
                    label="band-structure stage also charged")]
    ax.legend(handles=h, frameon=False, loc="upper center", bbox_to_anchor=(0.4, -0.13),
              ncol=3, fontsize=6.3, handletextpad=0.3, columnspacing=1.2)
    ax.text(1.495, n - 0.35, "n", ha="right", va="center", fontsize=6.3, color=INK2,
            fontweight="bold")
    fig.tight_layout(pad=0.3)
    fig.savefig(HERE / "fig3_gain.pdf")
    plt.close(fig)
    return rows


# ---------------------------------------------------------------- group and leakage
def fig_group_leak():
    fig = plt.figure(figsize=(183 * MM, 118 * MM))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05], width_ratios=[1, 1],
                          hspace=0.45, wspace=0.42, left=0.175, right=0.985, top=0.95, bottom=0.115)
    a = fig.add_subplot(gs[0, 0]); b = fig.add_subplot(gs[0, 1])
    c = fig.add_subplot(gs[1, 0]); d = fig.add_subplot(gs[1, 1])
    # (a) research-group paired win rate
    n1 = no1_arrays()
    w1 = sum(1 for x, y in zip(n1["CEDT_P4"][2], n1["DTDE"][2])
             if x is not None and (y is None or x < y))
    rep = json.load(open(PROJ / "paired_comparison_2/f3_repair/results/repair_graded.json"))
    sets = no3_sets()
    items = [("No. 1, 3×3, broker", w1, len(n1["DTDE"][2]), BROKER),
             ("No. 2, 8×8, broker", rep["cedt_wins"], rep["n_pairs"], BROKER)]
    for key, arm, lab, col in (("No3", "CEDT_P4", "No. 3, broker", BROKER),
                               ("No3", "CTDE_P2", "No. 3, shared logbook", PARTIAL),
                               ("No3", "MIX50", "No. 3, half via broker", PARTIAL),
                               ("No3b", "CEDT_P4", "No. 3b, broker", BROKER)):
        camps, sites = sets[key]
        seeds = [s for s in range(15) if (arm, s) in camps and ("DTDE", s) in camps]
        items.append((lab, wins(camps, arm, seeds), len(seeds), col))
    y = np.arange(len(items))[::-1]
    for yy, (lab, k, nn, col) in zip(y, items):
        lo, hi = wilson(k, nn)
        a.plot([lo, hi], [yy, yy], color=col, lw=1.4)
        a.plot(k / nn, yy, "o", ms=5, color=col, mec="white", mew=0.5)
        a.text(1.0, yy, f"{k} of {nn}", ha="right", va="center", fontsize=6.3, color=INK2)
    a.axvline(0.6, color=INK2, ls="--", lw=0.8)
    a.axvline(0.5, color="#d9d8d4", lw=0.6)
    a.text(0.6, len(items) - 0.35, "criterion 0.60", fontsize=6.3, color=INK2, ha="center",
           va="bottom")
    a.set_yticks(y, [it[0] for it in items])
    a.set_xlim(0, 1.0); a.set_ylim(-0.6, len(items) - 0.1)
    a.set_xlabel("paired win rate of the sharing arm\n(Wilson 95% interval)")
    a.tick_params(axis="y", length=0)
    letter(a, "(a)", x=-0.5, y=1.02)

    # (b) comparison No. 2, cost to first confirmed qualifier per pair
    top = 11.0
    pts = []
    for p in rep["pairs"]:
        dd, cc = p["dtde_first"], p["cedt_first"]
        if dd is None and cc is None:
            continue
        name = p["cedt_first_qualifier"] or p["dtde_first_qualifier"]
        pts.append((top if dd is None else dd, top if cc is None else cc, name, p["cedt_win"]))
    for x0, y0, name, win in pts:
        b.plot(x0, y0, ls="", marker="o" if win else "s", ms=4.6,
               color=BROKER if win else GREY, mec="white", mew=0.5, zorder=3)
    # label placement: right of the point, except where two points crowd
    pts.sort(key=lambda t: (t[0], t[1]))
    for i, (x0, y0, name, win) in enumerate(pts):
        off = (5, -1)
        if x0 < 3:
            off = (6, 7) if i % 2 == 0 else (6, -7)
        if x0 >= top:
            off = (-4, 6)
        b.annotate(chem(name), (x0, y0), xytext=off, textcoords="offset points",
                   ha="right" if off[0] < 0 else "left", va="center", fontsize=6.0, color=INK2)
    b.plot([0, top], [0, top], color="#b8b7b2", lw=0.6)
    b.set_xlim(0, top + 0.6); b.set_ylim(0, top + 0.6)
    b.set_xticks([0, 2, 4, 6, 8, top], ["0", "2", "4", "6", "8", "none"])
    b.set_yticks([0, 2, 4, 6, 8, top], ["0", "2", "4", "6", "8", "none"])
    b.set_xlabel("no sharing: charged node-hours\nto first confirmed candidate")
    b.set_ylabel("broker: charged node-hours\nto first confirmed candidate")
    b.plot([], [], ls="", marker="o", color=BROKER, ms=4.6, label="broker lower")
    b.plot([], [], ls="", marker="s", color=GREY, ms=4.6, label="no sharing lower")
    b.legend(frameon=False, loc="upper left", handletextpad=0.2, fontsize=6.3)
    n_none = sum(1 for p in rep["pairs"] if p["dtde_first"] is None and p["cedt_first"] is None)
    b.text(top + 0.4, 4.6, f"{n_none} of {rep['n_pairs']} pairs: no\ncandidate in either arm",
           ha="right", va="bottom", fontsize=6.1, color=INK2)
    letter(b, "(b)", x=-0.2, y=1.02)

    # (c, d) timing leakage, seed-bootstrap intervals (bootstrap_ci/)
    s1 = json.load(open(PROJ / "bootstrap_ci/seed_bootstrap.json"))
    s3 = json.load(open(PROJ / "bootstrap_ci/seed_bootstrap3.json"))
    rows = [("No. 1, 3×3, no sharing", s1["comparison1"]["DTDE"], GREY),
            ("No. 1, 3×3, broker", s1["comparison1"]["CEDT_P4"], BROKER),
            ("No. 2, 8×8, no sharing", s1["comparison2"]["DTDE"], GREY),
            ("No. 2, 8×8, broker", s1["comparison2"]["CEDT_P4"], BROKER),
            ("No. 3, no sharing", s3["arms"]["DTDE"], GREY),
            ("No. 3, broker", s3["arms"]["CEDT_P4"], BROKER),
            ("No. 3, shared logbook", s3["arms"]["CTDE_P2"], PARTIAL),
            ("No. 3, half via broker", s3["arms"]["MIX50_routed"], PARTIAL)]
    yy = np.arange(len(rows))[::-1]
    for ax_, key, lim, lab, xlab in ((c, "identity", s3["identity_limit"], "(c)",
                                      "compound identity accuracy"),
                                     (d, "auc", s3["auc_limit"], "(d)",
                                      "campaign membership AUC")):
        for y0, (name, dd, col) in zip(yy, rows):
            q = dd[key]
            ax_.plot([q["lo"], q["hi"]], [y0, y0], color=col, lw=1.4)
            ax_.plot(q["point"], y0, ls="", marker="o", ms=4.6, color=col, mec="white", mew=0.5)
        ax_.axvline(lim, color=INK2, ls="--", lw=0.8)
        ax_.text(lim, len(rows) - 0.35, "limit 0.084" if key == "identity" else f"limit {lim:.2f}",
                 fontsize=6.3, color=INK2, ha="left" if key == "identity" else "center", va="bottom")
        ax_.set_yticks(yy, [r[0] for r in rows] if ax_ is c else [])
        ax_.set_xlim(0, 1.02); ax_.set_ylim(-0.6, len(rows) - 0.1)
        ax_.set_xlabel(xlab + "\n(seed-bootstrap 95% interval)")
        ax_.tick_params(axis="y", length=0)
        letter(ax_, lab, x=-0.5 if ax_ is c else -0.06, y=1.02)
    d.axvline(0.5, color="#d9d8d4", lw=0.6)
    fig.savefig(HERE / "fig_group_leak.pdf")
    plt.close(fig)
    return items


if __name__ == "__main__":
    fig_method()
    fig_gaps()
    fig_stage_cost()
    fig_cost()
    rows = fig_gain()
    for r in rows:
        print("gain", r[0].replace("\n", " "), round(r[1], 3), r[2], r[3], r[5], r[6])
    for it in fig_group_leak():
        print("wins", it[:3])
    print("figures written to", HERE)
