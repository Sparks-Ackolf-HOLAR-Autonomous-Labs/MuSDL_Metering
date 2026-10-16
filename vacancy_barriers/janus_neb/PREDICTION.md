# Janus MXY vacancy-hop barriers: registered predictions

Registered 2026-10-12, before any Janus QE output exists. New NEB runs for the Janus MXY
family on a CHPC guest partition.

## Why this family

In the MX2 set the chalcogen-vacancy barrier follows the gap: about 77 % of compounds with
a PBE gap in 1.0-2.0 eV also clear the 1.5 eV barrier floor (vacancy goal, G3). A group's
own gap screen therefore already tells it most of what a peer's barrier result would, which
is why no exchange of results helped (`result_sharing/RESULTS.md`). A Janus layer has two
different chalcogen faces (MoSSe was grown by Lu et al. 2017 and Zhang et al. 2017), so a
vacancy has two barriers per compound and the gap is set by both faces together. If the
barriers come apart from the gap here, a barrier result carries information that a gap
screen does not, and the peer-exchange test can be repeated on a family where it might.

## Runs (array of 16 tasks, `inputs/manifest.json`)

- Janus 1H layers MoSSe, MoSTe, MoSeTe, WSSe, WSTe, WSeTe; one hop on each face (12 tasks).
- Single-layer references MoS2, WS2, MoSe2, WSe2, top face (4 tasks), through the same
  driver, as like-for-like anchors (the earlier 2H NEB values used AA' bilayer cells).
- Per task (`run_janus.py`): 1x1 in-plane relaxation (45/360 Ry, 12x12x1); PBE gap from
  scf 12x12x1 + nscf 18x18x1 with a dipole correction; 3x3 supercell, one chalcogen
  removed from the chosen face, its in-plane nearest neighbour hops into the site;
  fixed-cell end-state relaxations and a 5-image climbing-image NEB (40/320 Ry, 3x3x1,
  path_thr 0.1 eV/A), the settings of `../dft_neb`. NEB without dipole correction, as
  for the references. PBE, GBRV 1.5 USPP, non-spin-polarized, no SOC.

## Predictions

| Id | Prediction | Pass condition |
| --- | --- | --- |
| J0 | Single-layer anchors match the bilayer-cell values | MoS2 within 0.3 eV of 2.14 eV and WS2 within 0.3 eV of 2.66 eV |
| J1 | The two faces of a Janus layer have different barriers | top minus bottom at least 0.3 eV in magnitude in at least 4 of 6 Janus compounds |
| J2 | A face behaves like the parent with the same metal and chalcogen | within 0.3 eV of that parent's single-layer barrier in at least 6 of the 8 S and Se faces (no Te parent was run) |
| J3 | In the Janus family the barrier no longer follows the gap | Spearman rank correlation between the PBE gap and the lower of the two face barriers, over the 6 Janus compounds, at most 0.4 |
| J4 | The gap window no longer implies the barrier floor | among Janus compounds with a gap in 1.0-2.0 eV, the share whose lower face barrier is at least 1.5 eV is below 0.77 |

Decision rule: if J3 and J4 pass, the result-exchange simulation (`result_sharing/`) is
rerun with the Janus compounds added as a seventh family, with predictions registered
before that run. If either fails, the Janus family gives no new case for the disclosure article and the
result is reported only as a barrier data set.

Collection: `collect_janus.sh` (one Duo approval, or an open notchpeak ControlMaster),
then `analyze_janus.py` grades J0-J4.
