# Compound-specific memory demand determines when a research group and a computing centre benefit from sharing information in two-dimensional semiconductor screening

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23413478.svg)](https://doi.org/10.5281/zenodo.23413478)

Code and data for the article "Compound-specific memory demand determines when a research group and a computing centre benefit from sharing information in two-dimensional semiconductor screening" (prepared for Computational Materials Science (Elsevier)).

Thirty transition-metal dichalcogenide layers (MX2; M = Ti, Zr, Hf, Mo, W; X = S, Se, Te) in 1T and 2H structures were screened with Quantum ESPRESSO and the PBE functional through a four-stage fidelity ladder, and every job was timed and its memory recorded. The article reports the screen (band gaps, spin-orbit gaps, Gamma-point phonons, the corrected 2H stacking, chalcogen-vacancy and Janus-layer migration barriers from nudged elastic band calculations) and the measured job costs from primitive cells to 8x8 supercells: randomized benchmarks, scale series, memory and co-residency probes, a multi-node benchmark and four paired comparisons of sharing arrangements between the research group and the computing centre.

The repository holds what the article reports and what is needed to rebuild
it: inputs, campaign drivers, registered predictions, raw records, result
summaries and the figure generator. Results are labelled **measured** (runs
on hardware) or **simulated** (the campaign simulator; simulated node-hours
are simulated units, not measured cost).

## Quick start

```
python3 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
cd papers/paperB_qe_mx2_campaign/figures && python3 make_figures.py
```

The figure generator regenerates every figure of the article from the
records in this repository; the manuscript and the rendered figures are not
included. `figure_test.json` records a run of the generator from a copy of
this repository with no access to anything else.

## Contents

- `papers/paperB_qe_mx2_campaign/figures/make_figures.py`: the figure generator.
- One folder per experiment:

| folder | content |
|---|---|
| `bootstrap_ci/` | cluster-bootstrap intervals over seeds for the measured timing attacks |
| `centralization_surface/` | simulated: gain for both parties over training and execution centralization, task scale, number of research parties and schedules |
| `cofit_probe/` | measured: two 6x6 jobs run side by side versus one after another |
| `contention_transfer_bench/` | measured: co-tenant contention and defect and heterostructure cost probes |
| `dft_campaign/` | measured: Quantum ESPRESSO screen of the 30 MX2 candidates (fidelity stages F1 to F4, spin-orbit finalists) |
| `fast_path/` | replay of recorded run times and the per-candidate run-time library |
| `memory_band_probe/` | measured: 5x5 to 8x8 supercells, configuration cost and memory |
| `multinode_bench/` | measured: 1, 2 and 4 node scaling |
| `paired_comparison_1/` | measured: comparison No. 1, 3x3 supercells, 480 QE calls |
| `paired_comparison_2/` | measured: comparison No. 2, 8x8 supercells, job packing and timing masking, band-structure repair pass |
| `paired_comparison_3/` | measured: comparison No. 3 on guest nodes (DTDE, CEDT, half-routed mixture and CTDE), timing attack |
| `paired_comparison_3b/` | measured: comparison No. 3b, single-layer 2H cells at 8x8 |
| `pilot/` | measured: randomized 81-job pilot |
| `prior_predictive/` | campaign simulator, prior specifications and the 30 candidate structures |
| `rescreen_3x3/` | measured: 3x3 rescreen on guest nodes |
| `scale_validation/` | measured: supercell, layer and strain cost series |
| `stacking_study/` | measured: 2H single layers and DFT-D3 bilayers (AA' and 3R), with relaxation, spin-orbit and phonon follow-up |
| `stage2_twin_fit/` | fit of the cost model to the measured QE calls |
| `timing_masking/` | measured completion times re-attacked under masking clauses |
| `vacancy_barriers/` | measured: DFT nudged elastic band vacancy barriers (7 compounds) and Janus-layer barriers, MACE-MP-0 estimates, kinetic Monte Carlo costs; simulated vacancy-lifetime objective |

- `prior_predictive/candidates/`: the 30 idealized MX2 structures and their
  metadata (CC0-1.0).
- `case_manifest.yaml`: the case definition (budgets, thresholds, leakage
  limits, fidelity ladder F0 to F4), fixed before the measured campaigns.
- Dated preregistration and results records use "frozen" for "fixed before
  the data were collected". They are kept as written, because their dates
  and content are what the preregistration states.

The folders `centralization_surface/`, `fast_path/`, `prior_predictive/`, `timing_masking/` are included because scripts in the reported folders read their files or stage them for a rerun.


Files the figure generator reads (`traced_inputs.json`):

- `bootstrap_ci/seed_bootstrap.json`
- `bootstrap_ci/seed_bootstrap3.json`
- `dft_campaign/results/outcomes.jsonl`
- `dft_campaign/soc/results/soc_outcomes.jsonl`
- `contention_transfer_bench/results/bench_analysis.json`
- `contention_transfer_bench/results/bench_results.jsonl`
- `paired_comparison_1/results/campaigns.jsonl`
- `paired_comparison_2/f3_repair/results/repair_graded.json`
- `paired_comparison_3/collected_2026-10-14/CEDT_P4_s0/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CEDT_P4_s1/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CEDT_P4_s10/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CEDT_P4_s11/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CEDT_P4_s12/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CEDT_P4_s13/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CEDT_P4_s14/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CEDT_P4_s2/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CEDT_P4_s3/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CEDT_P4_s4/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CEDT_P4_s5/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CEDT_P4_s6/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CEDT_P4_s7/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CEDT_P4_s8/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CEDT_P4_s9/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CTDE_P2_s0/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CTDE_P2_s1/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CTDE_P2_s10/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CTDE_P2_s11/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CTDE_P2_s12/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CTDE_P2_s13/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CTDE_P2_s14/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CTDE_P2_s2/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CTDE_P2_s3/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CTDE_P2_s4/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CTDE_P2_s5/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CTDE_P2_s6/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CTDE_P2_s7/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CTDE_P2_s8/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/CTDE_P2_s9/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/DTDE_s0/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/DTDE_s1/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/DTDE_s10/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/DTDE_s11/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/DTDE_s12/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/DTDE_s13/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/DTDE_s14/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/DTDE_s2/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/DTDE_s3/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/DTDE_s4/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/DTDE_s5/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/DTDE_s6/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/DTDE_s7/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/DTDE_s8/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/DTDE_s9/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/MIX50_s0/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/MIX50_s1/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/MIX50_s10/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/MIX50_s11/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/MIX50_s12/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/MIX50_s13/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/MIX50_s14/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/MIX50_s2/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/MIX50_s3/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/MIX50_s4/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/MIX50_s5/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/MIX50_s6/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/MIX50_s7/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/MIX50_s8/campaigns.jsonl`
- `paired_comparison_3/collected_2026-10-14/MIX50_s9/campaigns.jsonl`
- `paired_comparison_3/pair_sites.json`
- `paired_comparison_3b/collected_2026-10-14/CEDT_P4_s0/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/CEDT_P4_s1/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/CEDT_P4_s10/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/CEDT_P4_s11/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/CEDT_P4_s12/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/CEDT_P4_s13/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/CEDT_P4_s14/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/CEDT_P4_s2/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/CEDT_P4_s3/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/CEDT_P4_s4/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/CEDT_P4_s5/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/CEDT_P4_s6/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/CEDT_P4_s7/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/CEDT_P4_s8/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/CEDT_P4_s9/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/DTDE_s0/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/DTDE_s1/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/DTDE_s10/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/DTDE_s11/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/DTDE_s12/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/DTDE_s13/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/DTDE_s14/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/DTDE_s2/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/DTDE_s3/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/DTDE_s4/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/DTDE_s5/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/DTDE_s6/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/DTDE_s7/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/DTDE_s8/campaigns.jsonl`
- `paired_comparison_3b/collected_2026-10-14/DTDE_s9/campaigns.jsonl`
- `paired_comparison_3b/pair_sites.json`
- `scale_validation/results/scale_fit.json`
- `stacking_study/followup_2026-10-12/results/relax/results_16208957/outcomes.jsonl`
- `stacking_study/followup_2026-10-12/results/socph/results_16208958/soc_outcomes.jsonl`
- `stacking_study/results_2088961/outcomes.jsonl`

## Reproducing the results

`REPRODUCE.md` gives the software versions, the commands that regenerate
every summary and figure from the committed records, and the commands and
resources that re-run each computation (density functional theory jobs,
paired campaigns, simulations). `requirements.txt` lists the Python
packages; Quantum ESPRESSO 7.4 is needed only to re-run the DFT jobs.

## Licences and citation

Code (Python and shell scripts) is released under the MIT License
(`LICENSE`); data, result files and regenerated figures under the Creative
Commons Attribution 4.0 International licence (`LICENSE-CC-BY-4.0.txt`). The
candidate structures keep CC0-1.0. Cite with `CITATION.cff`.

`FILES.sha256` lists a SHA-256 hash for every file. Pseudopotentials are not
redistributed (GBRV 1.5 and PSlibrary 1.0.0 keep their own licences); file
names and hashes are given in `REPRODUCE.md` and in the experiment folders.
User identifiers were removed from the committed outputs.
