#!/usr/bin/env python3
"""Grade the registered Janus predictions J0-J4 from results/*.json."""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
R = {}
for p in sorted((HERE / "results").glob("*.json")):
    if p.name == "grading.json":
        continue
    d = json.loads(p.read_text())
    R[(d["compound"], d["face"])] = d
E = {k: (v.get("E_hop_eV") if v.get("neb_converged") else None) for k, v in R.items()}
G = {k[0]: v.get("gap_eV") for k, v in R.items() if v.get("gap_eV") is not None}
JAN = ["MoSSe", "MoSTe", "MoSeTe", "WSSe", "WSTe", "WSeTe"]
PAR = {"S": {"Mo": "MoS2", "W": "WS2"}, "Se": {"Mo": "MoSe2", "W": "WSe2"}}
out = {"n_results": len(R), "barrier_eV": {f"{c}_{f}": e for (c, f), e in E.items()},
       "gap_eV": G, "verdicts": {}}


def spearman(x, y):
    rk = lambda v: [sorted(v).index(a) for a in v]
    rx, ry, n = rk(x), rk(y), len(x)
    return 1 - 6 * sum((a - b) ** 2 for a, b in zip(rx, ry)) / (n * (n * n - 1))


v = out["verdicts"]
a, b = E.get(("MoS2", "top")), E.get(("WS2", "top"))
v["J0"] = None if None in (a, b) else (abs(a - 2.14) <= 0.3 and abs(b - 2.66) <= 0.3)
pairs = [(E.get((c, "top")), E.get((c, "bottom"))) for c in JAN]
pairs = [p for p in pairs if None not in p]
v["J1"] = None if len(pairs) < 6 else sum(abs(t - u) >= 0.3 for t, u in pairs) >= 4
hits, n = 0, 0
for c in JAN:
    m = "Mo" if c.startswith("Mo") else "W"
    body = c[len(m):]
    x, y = (body[:2], body[2:]) if body.startswith("Se") else (body[:1], body[1:])
    for face, ch in (("top", x), ("bottom", y)):
        if ch in PAR and E.get((c, face)) is not None and E.get((PAR[ch][m], "top")) is not None:
            n += 1; hits += abs(E[(c, face)] - E[(PAR[ch][m], "top")]) <= 0.3
v["J2"] = None if n < 8 else hits >= 6
low = {c: min(E[(c, "top")], E[(c, "bottom")]) for c in JAN
       if E.get((c, "top")) is not None and E.get((c, "bottom")) is not None and c in G}
out["lower_face_barrier_eV"] = low
if len(low) == 6:
    rho = spearman([G[c] for c in low], list(low.values()))
    out["spearman_gap_vs_lower_barrier"] = rho
    v["J3"] = rho <= 0.4
    win = [c for c in low if 1.0 <= G[c] <= 2.0]
    share = sum(low[c] >= 1.5 for c in win) / len(win) if win else None
    out["window_share_clearing_1p5"] = share
    v["J4"] = None if share is None else share < 0.77
else:
    v["J3"] = v["J4"] = None
(HERE / "results" / "grading.json").write_text(json.dumps(out, indent=1))
print(json.dumps(v))
