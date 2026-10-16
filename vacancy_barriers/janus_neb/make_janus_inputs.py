#!/usr/bin/env python3
"""Build the task manifest and 1x1 start cells for the Janus MXY vacancy-hop NEB runs.

Janus MXY monolayers (single 1H layer, metal M = Mo or W; chalcogen X on the top face,
Y on the bottom face). Start cell: a = mean of the relaxed single-layer 2H parents
(stacking_study, job 2088961), metal at z = 0.5, X and Y at the parents' M-X and M-Y
heights. The CHPC driver run_janus.py relaxes the 1x1 cell in plane, computes the PBE gap,
builds the 3x3 supercell and runs a climbing-image NEB for one chalcogen-vacancy hop on
the chosen face. The four MX2 parents run through the same driver as references.
Settings match ../dft_neb (PBE, GBRV 1.5 USPP, 40/320 Ry NEB, cold smearing, 5 images).
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
STACK = HERE.parents[1] / "stacking_study" / "results_2088961" / "stacking_summary.json"
C_ANG = 18.0
PSEUDO = {"Mo": ("95.95", "mo_pbe_v1.uspp.F.UPF"), "W": ("183.84", "w_pbe_v1.2.uspp.F.UPF"),
          "S": ("32.06", "s_pbe_v1.4.uspp.F.UPF"), "Se": ("78.971", "se_pbe_v1.uspp.F.UPF"),
          "Te": ("127.60", "te_pbe_v1.uspp.F.UPF")}
JANUS = [("Mo", "S", "Se"), ("Mo", "S", "Te"), ("Mo", "Se", "Te"),
         ("W", "S", "Se"), ("W", "S", "Te"), ("W", "Se", "Te")]
REFS = [("Mo", "S", "S"), ("W", "S", "S"), ("Mo", "Se", "Se"), ("W", "Se", "Se")]


def parents():
    out = {}
    for r in json.loads(STACK.read_text())["records"]:
        if r["candidate"].endswith("_mono2H") and r["metal"] in ("Mo", "W"):
            f2 = r["stages"]["F2"]
            pos = f2["final_positions_crystal"]
            dz = (pos[1][3] - pos[0][3]) * C_ANG          # metal-chalcogen height, A
            out[(r["metal"], r["chalcogen"])] = {"a": f2["a_relaxed_angstrom"], "dz": dz,
                                                 "gap": r["stages"]["F3"].get("gap_ev")}
    return out


def main():
    par = parents()
    tasks, cells = [], {}
    for m, x, y in JANUS + REFS:
        name = f"{m}{x}{y}" if x != y else f"{m}{x}2"
        a = 0.5 * (par[(m, x)]["a"] + par[(m, y)]["a"])
        cells[name] = {
            "metal": m, "top": x, "bottom": y, "a0": round(a, 4), "c": C_ANG,
            "species": sorted({m, x, y}, key=lambda s: (s != m, s)),
            "pseudo": {s: PSEUDO[s] for s in {m, x, y}},
            "frac": [[m, 1 / 3, 2 / 3, 0.5],
                     [x, 2 / 3, 1 / 3, 0.5 + par[(m, x)]["dz"] / C_ANG],
                     [y, 2 / 3, 1 / 3, 0.5 - par[(m, y)]["dz"] / C_ANG]],
            "parent_gaps_eV": [par[(m, x)]["gap"], par[(m, y)]["gap"]],
            "janus": x != y}
        sides = ("top", "bottom") if x != y else ("top",)
        for s in sides:
            tasks.append({"task": len(tasks), "compound": name, "face": s,
                          "vacancy": x if s == "top" else y})
    (HERE / "inputs" / "manifest.json").write_text(
        json.dumps({"cells": cells, "tasks": tasks}, indent=1))
    for t in tasks:
        print(t["task"], t["compound"], t["face"], t["vacancy"])


if __name__ == "__main__":
    main()
