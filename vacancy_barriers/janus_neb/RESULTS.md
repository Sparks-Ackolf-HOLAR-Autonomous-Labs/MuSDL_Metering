# Janus MXY vacancy-hop barriers: results

Collected 2026-10-14 (array 16208159, 16/16 tasks COMPLETED on notchpeak-guest,
1.2 to 3.1 h each). Graded by `analyze_janus.py` against the predictions in
`PREDICTION.md` (registered 2026-10-12, before any output). Numbers:
`results/grading.json`; per-task records `results/*.json`; raw NEB outputs
`results/outputs/raw.tgz`.

| compound | PBE gap (eV) | barrier, top face (eV) | barrier, bottom face (eV) |
|---|---:|---:|---:|
| MoS2 (ref.) | 1.68 | 2.16 | |
| MoSe2 (ref.) | 1.46 | 1.91 | |
| WS2 (ref.) | 1.84 | 2.66 | |
| WSe2 (ref.) | 1.58 | 2.39 | |
| MoSSe | 1.58 | 2.04 | 1.95 |
| MoSTe | 0.98 | 1.70 | 1.60 |
| MoSeTe | 1.27 | 1.64 | 1.53 |
| WSSe | 1.72 | 2.56 | 2.44 |
| WSTe | 1.22 | 2.17 | 2.08 |
| WSeTe | 1.37 | 2.16 | 2.11 |

| Id | result | verdict |
|---|---|---|
| J0 | single-layer MoS2 2.16 eV (bilayer cell 2.14), WS2 2.66 eV (2.66) | pass |
| J1 | the two faces differ by only 0.04 to 0.12 eV (0 of 6 reach 0.3 eV) | fail |
| J2 | faces sit within 0.3 eV of the parent with the same metal and chalcogen | pass |
| J3 | Spearman rank correlation, gap vs lower face barrier = 0.60 (limit 0.4) | fail |
| J4 | all 5 Janus compounds in the gap window clear 1.5 eV (share 1.00, limit below 0.77) | fail |

## Reading

The two chalcogen faces of a Janus layer give almost the same hop barrier,
and the barrier is set mainly by the metal (W about 0.4 to 0.5 eV above Mo)
and by the heavier chalcogen, which also lowers the gap. The barrier
therefore still follows the gap, as in the MX2 set, and every Janus compound
inside the 1.0 to 2.0 eV window clears the 1.5 eV immobility floor. By the
registered decision rule the Janus family gives no new case for the
result-exchange study (`result_sharing/` is not rerun); the sixteen barriers
are reported as a data set (the screening article).
