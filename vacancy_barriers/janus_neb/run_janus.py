#!/usr/bin/env python3
"""One Janus (or parent) vacancy-hop task on a CHPC guest node; resume-safe, stdlib only.

Task = SLURM_ARRAY_TASK_ID into inputs/manifest.json["tasks"]. Stages, each skipped on a
requeue if its output finished:
  S1 1x1 in-plane relaxation (vc-relax, cell_dofree='2Dxy'), 12x12x1 k-points;
  S2 PBE gap: scf 12x12x1 then nscf 18x18x1 with dipole correction (non-fatal);
  S3 fixed-cell relax of the two vacancy end states in the 3x3 supercell (one chalcogen
     removed from the chosen face near the cell centre; its in-plane nearest neighbour on
     the same face moves into the site);
  S4 neb.x, 5 images, climbing image, path_thr 0.1 eV/A (as ../dft_neb).
Writes results/<compound>_<face>.json.
"""
import json
import math
import os
import re
import subprocess
import time
from pathlib import Path

ROOT = Path(os.environ.get("MPT_JANUS_ROOT", ".")).resolve()
MAN = json.loads((ROOT / "inputs" / "manifest.json").read_text())
T = MAN["tasks"][int(os.environ["SLURM_ARRAY_TASK_ID"])]
CELL = MAN["cells"][T["compound"]]
TAG = f"{T['compound']}_{T['face']}"
W = ROOT / "work" / TAG
W.mkdir(parents=True, exist_ok=True)
(ROOT / "results").mkdir(exist_ok=True)
NP = int(os.environ.get("SLURM_NTASKS", "32"))
NK = os.environ.get("MPT_NEB_NK", "2")
BOHR = 0.529177210903
log = {**T, "restart_count": int(os.environ.get("SLURM_RESTART_COUNT", "0")),
       "node": os.environ.get("SLURMD_NODENAME"), "ranks": NP, "npool": NK}
RES = ROOT / "results" / f"{TAG}.json"


def save():
    RES.write_text(json.dumps(log, indent=1))


def run(exe, inp, out, flag="-in"):
    t0 = time.time()
    rc = subprocess.call(f"mpirun -np {NP} {exe} -nk {NK} {flag} {inp} > {out} 2>&1",
                         shell=True, cwd=W)
    return rc, round(time.time() - t0, 1)


def done(out):
    p = W / out
    return p.exists() and "JOB DONE" in p.read_text(errors="ignore")


def cpu_s(txt, tag):
    m = re.findall(tag + r"\s*:\s*((?:[0-9.]+[hms]?\s*)+?)\s*CPU", txt)
    if not m:
        return None
    tok, sec, cur = re.sub(r"\s+", "", m[-1]), 0.0, ""
    for ch in tok:
        if ch in "hms":
            sec += float(cur) * {"h": 3600, "m": 60, "s": 1}[ch]; cur = ""
        else:
            cur += ch
    return sec + (float(cur) if cur else 0.0)


def species_block():
    return "".join(f"  {s}  {CELL['pseudo'][s][0]}  {CELL['pseudo'][s][1]}\n"
                   for s in CELL["species"])


def deck(prefix, calc, a, nat, pos, kpts, ecut=(40.0, 320.0), extra_sys="",
         extra_ctrl="", extra_el="", cell_ns=""):
    p = "\n".join(f"  {s}  {x:.10f}  {y:.10f}  {z:.10f}" for s, x, y, z in pos)
    return f"""&CONTROL
  calculation = '{calc}'
  prefix = '{prefix}'
  pseudo_dir = '../../pseudos'
  outdir = './out_{prefix}'
  disk_io = 'low'
  forc_conv_thr = 1.0d-3
  etot_conv_thr = 1.0d-5
  nstep = 150{extra_ctrl}
/
&SYSTEM
  ibrav = 4
  a = {a:.6f}
  c = {CELL['c']}
  nat = {nat}
  ntyp = {len(CELL['species'])}
  ecutwfc = {ecut[0]}
  ecutrho = {ecut[1]}
  occupations = 'smearing'
  smearing = 'cold'
  degauss = 0.015{extra_sys}
/
&ELECTRONS
  mixing_beta = 0.3
  conv_thr = 1.0d-8
  electron_maxstep = 150{extra_el}
/
&IONS
/
{cell_ns}ATOMIC_SPECIES
{species_block()}ATOMIC_POSITIONS crystal
{p}
K_POINTS automatic
  {kpts} {kpts} 1 0 0 0
"""


def final_coords(out, want_cell=False):
    t = (W / out).read_text(errors="ignore")
    if "Begin final coordinates" not in t:
        return None
    blk = t.split("Begin final coordinates")[1].split("End final coordinates")[0]
    a = None
    if want_cell:
        m = re.search(r"CELL_PARAMETERS \((alat=\s*([0-9.]+)|angstrom)\)\n(.*)\n", blk)
        v = [float(x) for x in m.group(3).split()]
        scale = float(m.group(2)) * BOHR if m.group(2) else 1.0
        a = scale * math.sqrt(v[0] ** 2 + v[1] ** 2)
    lines = blk.split("ATOMIC_POSITIONS (crystal)")[1].strip().splitlines()
    pos = [[l.split()[0]] + [float(x) for x in l.split()[1:4]] for l in lines if l.strip()]
    e = float(re.findall(r"!\s+total energy\s+=\s+(-?[0-9.]+) Ry", t)[-1])
    return a, pos, e, cpu_s(t, "PWSCF")


# ---- S1: 1x1 in-plane relaxation
if not done("s1.out"):
    (W / "s1.in").write_text(deck("s1", "vc-relax", CELL["a0"], 3, CELL["frac"], 12,
                                  ecut=(45.0, 360.0), cell_ns="&CELL\n  cell_dofree = '2Dxy'\n/\n"))
    log["s1_rc"], log["s1_wall_s"] = run("pw.x", "s1.in", "s1.out")
r = final_coords("s1.out", want_cell=True)
if r is None:
    log["error"] = "1x1 relaxation did not finish"; save(); raise SystemExit(1)
a1, pos1, log["s1_E_Ry"], log["s1_cpu_s"] = r
log["a_relaxed_A"] = round(a1, 4)
save()

# ---- S2: gap (scf + nscf, dipole correction across the vacuum); failure is recorded only
dip = ("\n  edir = 3\n  emaxpos = 0.95\n  eopreg = 0.04\n  eamp = 0.0", "\n  tefield = .true.\n  dipfield = .true.")
try:
    if not done("s2n.out"):
        (W / "s2.in").write_text(deck("s2", "scf", a1, 3, pos1, 12, ecut=(45.0, 360.0),
                                      extra_sys=dip[0], extra_ctrl=dip[1]))
        run("pw.x", "s2.in", "s2.out")
        nel = float(re.search(r"number of electrons\s+=\s+([0-9.]+)",
                              (W / "s2.out").read_text()).group(1))
        nocc = int(round(nel / 2))
        (W / "s2n.in").write_text(deck("s2", "nscf", a1, 3, pos1, 18, ecut=(45.0, 360.0),
                                       extra_sys=dip[0] + f"\n  nbnd = {nocc + 8}",
                                       extra_ctrl=dip[1] + "\n  verbosity = 'high'"))
        run("pw.x", "s2n.in", "s2n.out")
    t = (W / "s2n.out").read_text(errors="ignore")
    nel = float(re.search(r"number of electrons\s+=\s+([0-9.]+)", t).group(1))
    nocc = int(round(nel / 2))
    blocks = re.findall(r"bands \(ev\):\s*\n(.*?)(?:\n\s*\n\s*occupation|\n\s*\n\s*k =|\n\s*\n\s*the Fermi)",
                        t, re.S)
    vbm, cbm = -1e9, 1e9
    for b in blocks:
        ev = [float(x) for x in re.findall(r"-?\d+\.\d+", b)]
        vbm, cbm = max(vbm, ev[nocc - 1]), min(cbm, ev[nocc])
    log["gap_eV"] = round(cbm - vbm, 4) if blocks else None
    log["gap_nk"] = len(blocks)
except Exception as exc:                        # noqa: BLE001
    log["gap_error"] = repr(exc)
save()

# ---- S3: 3x3 supercell and vacancy end states
sup = [[s, (x + i) / 3, (y + j) / 3, z] for i in range(3) for j in range(3) for s, x, y, z in pos1]
A3 = 3 * a1


def cart(p):
    return (A3 * (p[1] - 0.5 * p[2]), A3 * math.sqrt(3) / 2 * p[2], CELL["c"] * p[3])


vac = T["vacancy"]
zsel = max if T["face"] == "top" else min
face_z = zsel(p[3] for p in sup if p[0] == vac and (p[3] > 0.5) == (T["face"] == "top"))
face = [k for k, p in enumerate(sup) if p[0] == vac and abs(p[3] - face_z) < 0.02]
ctr = cart(["", 0.5, 0.5, 0])
i = min(face, key=lambda k: math.dist(cart(sup[k])[:2], ctr[:2]))
j = min((k for k in face if k != i), key=lambda k: math.dist(cart(sup[k]), cart(sup[i])))
log["hop_distance_A"] = round(math.dist(cart(sup[i]), cart(sup[j])), 4)
keep = [k for k in range(len(sup)) if k != i]
ini = [sup[k][:] for k in keep]
fin = [p[:] for p in ini]
jj = keep.index(j)
fin[jj][1:] = sup[i][1:]
nat = len(ini)
for st, pos in (("ini", ini), ("fin", fin)):
    if not done(f"{st}.out"):
        (W / f"{st}.in").write_text(deck(st, "relax", A3, nat, pos, 3))
        log[f"{st}_rc"], log[f"{st}_wall_s"] = run("pw.x", f"{st}.in", f"{st}.out")
    rr = final_coords(f"{st}.out") if (W / f"{st}.out").exists() else None
    if rr is None:
        log["error"] = f"{st} relax did not finish"; save(); raise SystemExit(1)
    log[f"{st}_E_Ry"], log[f"{st}_cpu_s"] = rr[2], rr[3]
_, pi, _, _ = final_coords("ini.out")
_, pf, _, _ = final_coords("fin.out")
pf = [[q[0]] + [f - round(f - g) for f, g in zip(q[1:], p[1:])] for q, p in zip(pf, pi)]  # same image
save()

# ---- S4: NEB
fmt = lambda P: "\n".join(f"  {s}  {x:.10f}  {y:.10f}  {z:.10f}" for s, x, y, z in P)
engine = deck("neb", "scf", A3, nat, [], 3).replace("  calculation = 'scf'\n", "") \
    .replace("&IONS\n/\n", "").split("ATOMIC_POSITIONS crystal")[0]
restart = "restart" if (W / "neb.path").exists() else "from_scratch"
(W / "neb.in").write_text(f"""BEGIN
BEGIN_PATH_INPUT
&PATH
  restart_mode = '{restart}'
  string_method = 'neb'
  nstep_path = 120
  num_of_images = 5
  opt_scheme = 'broyden'
  CI_scheme = 'auto'
  path_thr = 0.1
/
END_PATH_INPUT
BEGIN_ENGINE_INPUT
{engine.replace("conv_thr = 1.0d-8", "conv_thr = 1.0d-7")}BEGIN_POSITIONS
FIRST_IMAGE
ATOMIC_POSITIONS crystal
{fmt(pi)}
LAST_IMAGE
ATOMIC_POSITIONS crystal
{fmt(pf)}
END_POSITIONS
K_POINTS automatic
  3 3 1 0 0 0
END_ENGINE_INPUT
END
""")
log["neb_restart_mode"] = restart
log["neb_rc"], log["neb_wall_s"] = run("neb.x", "neb.in", "neb.out", flag="-inp")
t = (W / "neb.out").read_text(errors="ignore")
fw = re.findall(r"activation energy \(->\)\s*=\s*(-?[0-9.]+)\s*eV", t)
bw = re.findall(r"activation energy \(<-\)\s*=\s*(-?[0-9.]+)\s*eV", t)
log["E_hop_eV"] = float(fw[-1]) if fw else None
log["E_hop_back_eV"] = float(bw[-1]) if bw else None
log["neb_converged"] = "neb: convergence achieved" in t
log["neb_iterations"] = len(re.findall(r"-{10,} iteration", t))
log["neb_cpu_s"] = cpu_s(t, "NEB")
if (W / "neb.dat").exists():
    log["profile_eV"] = [[float(x) for x in l.split()[:2]]
                         for l in (W / "neb.dat").read_text().splitlines() if l.strip()]
save()
print(json.dumps({k: log.get(k) for k in ("compound", "face", "gap_eV", "E_hop_eV", "neb_converged")}))
