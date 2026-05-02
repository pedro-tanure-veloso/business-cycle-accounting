# BCA — Business Cycle Accounting

A Python implementation of the Business Cycle Accounting methodology of
Chari, Kehoe & McGrattan (2007) and Brinca, Chari, Kehoe & McGrattan
(2016, NBER WP 22663). Decomposes aggregate fluctuations into four
structural wedges — efficiency (A), labor (1−τ_l), investment (1+τ_x),
and government (g) — by Kalman-filter MLE on a VAR(1), with
counterfactual simulations and Brinca-Iskrev-Loria-style f-statistics.

The project goal is a usable BCA toolkit for **arbitrary US time
windows**. The BCKM 2016 US 1980Q1–2014Q4 replication is the locked-in
regression test that pins fidelity; new windows (e.g. COVID
2010Q1–2023Q4) are validated against narrative priors.
**Cross-country support is explicitly out of scope.**

## Status

- **Layer 1 — BCKM 2016 US 1980-2014 replication**: closed
  (2026-05-01). f-statistics match Table 11 to ≤ 0.01 in every channel
  at BCKM-θ on our dataset; full wrap-up in
  [`bckm_replication/REPORT.md`](bckm_replication/REPORT.md).
- **Layer 2 — generalizability validation (COVID 2010Q1–2023Q4)**:
  passed (2026-05-02). 6/6 narrative-prior rubric checks pass under
  both trend variants; labor wedge accounts for 49% / 77% of the
  output / hours f-stat weight over the 2019Q4–2022Q4 window —
  pipeline correctly identifies COVID as a primarily labor-driven
  contraction with a counter-cyclical investment wedge. Side-by-side
  COVID-vs-Great-Recession analysis in
  [`covid_analysis/REPORT.md`](covid_analysis/REPORT.md).
- **MLE result caching**: shipped (2026-05-02). Driver re-runs in
  ~3 seconds with cached pickles; cold runs ~10–15 minutes.
- **Layer 3 — cross-country**: explicitly out of scope.

## Quick start

```bash
# Python 3.10+; create a venv to keep deps isolated
python3.12 -m venv .venv
.venv/bin/pip install -e ".[dev]"

# FRED API key required only for the first data fetch; cached after.
# Get one at fred.stlouisfed.org/docs/api/api_key.html
echo "FRED_API_KEY=your_key" > .env

# Run the BCKM 1980-2014 regression scripts (cached parquet committed)
.venv/bin/python bckm_replication/scripts/eval_bckm_fstats.py

# Run the COVID 2010-2023 smoke test (~3s with cached MLE pickles,
# ~20 min for a fresh re-fit via --no-cache-mle)
set -a && source .env && set +a
.venv/bin/python covid_analysis/scripts/run_covid_analysis.py
```

Cached series live in `~/.bca_cache/fred/`; the BCKM-replication
parquets in `bckm_replication/data/`; COVID parquets and figures in
`covid_analysis/`.

## Repository layout

```
bca_core/                   pure-Python computational core (no web deps)
  params.py                 CalibrationParams (BCKM Table 1 defaults)
  constants.py              canonical BCKM Table 8/9/10 parameters (code convention)
  bckm_reference.py         loader for bckm_replication/matlab_reference/worktemp.mat
  model.py                  Prototype model: SS, log-linearization, decision rules
  klein.py                  Klein (2000) QZ solver for linear RE models
  bckm_lom.py               BCKM-faithful capital LOM (fixexpadj.m)
  wedges.py                 Static wedge extraction from data
  var_estimation.py         Kalman-filter MLE (BCKM mleqadj.m logic) + result caching
  counterfactuals.py        Wedge-alone simulations + f-statistics
  data/
    fred.py                 FRED fetcher with disk caching (incl. CPS / BLS labor)
    bea.py                  BEA NIPA + Fixed-Asset fetcher (diagnostic-only)
    adjustments.py          Durables PIM, sales-tax, per-capita, calgz detrending,
                            BLS-faithful labor construction (employment_cps × hours)
    pipeline.py             End-to-end fetch → adjust → detrend (build_us_dataset)

scripts/                    generic pipeline drivers (window-agnostic)
  run_var_counterfactuals.py  Main entry point — full pipeline + Figure 2B
  diagnose_counterfactuals.py PASS/FAIL diagnostics for the CF decomposition
  solow_residual_check.py     Solow residual sanity check

bckm_replication/           BCKM 2016 US 1980-2014 ground-truth artifacts
  REPORT.md                 final replication wrap-up (closed 2026-05-01)
  BCKM_RESULTS.md           paper Parts III & IV — regression-target tables
  BCKM_DIFF_GUIDE.md        element-wise diff guide vs worktemp.mat
  DATA_FORENSICS.md         BEA NIPA migration walkdown
  DIVERGENCE_ANALYSIS.md    earlier head-to-head analysis
  bca_paper.pdf             full paper PDF
  matlab_reference/         BCKM 2016 matlab code (paper ground truth)
  octave_output/            octave dumps used as fixtures
  data/                     pinned parquets + bootstrap/sensitivity dumps
  scripts/                  diag_*, eval_bckm_*, compare_*, plot_*, bootstrap_*, sensitivity_*
  figures/                  Figure 2A–E reproductions, Solow residual, etc.

covid_analysis/             Layer 2 — COVID 2010Q1–2023Q4 smoke test
  scripts/run_covid_analysis.py   one-shot driver (datasets → MLE → CFs → figures)
  data/                     cached parquets + .meta.json + (gitignored) .mle.pkl
  figures/                  figure_A/B/2C/2D/2E_covid.png + wedges_us_2010_2023.png
                            (BCKM Part III style, both trend variants)
  REPORT.md                 narrative results + COVID-vs-Great-Recession comparison
  Diary.md                  session-by-session journal

tests/                      pytest unit tests (91 tests)
  test_bckm_*               BCKM 1980-2014 regression (marked @pytest.mark.bckm)
  test_covid_analysis.py    COVID dataset shape + structural identities + smoke
  test_*                    window-agnostic structural tests
```

## Documentation map

- [CLAUDE.md](CLAUDE.md) — working rules and methodology constraints
  (read before changing anything in `bca_core/`); also carries the
  append-only "Findings" journal
- [BCA_info.md](BCA_info.md) — paper summary, **live methodology**
  (Parts I & II): prototype economy, accounting procedure, data and
  calibration. Applies to any US window the toolkit is run on
- [bckm_replication/REPORT.md](bckm_replication/REPORT.md) — final
  BCKM-replication wrap-up: bugs found and fixed, methodology
  decisions, residual issues
- [bckm_replication/BCKM_RESULTS.md](bckm_replication/BCKM_RESULTS.md)
  — paper Parts III & IV (cross-country findings + US 1980-2014
  empirical tables): the regression-test ground truth
- [covid_analysis/REPORT.md](covid_analysis/REPORT.md) — COVID
  smoke-test results and Part III-style comparison with the Great
  Recession
- [covid_analysis/Diary.md](covid_analysis/Diary.md) — Layer-2
  session journal (append-only, persists across context resets)
- [COVID_PLAN.md](COVID_PLAN.md) — original approved plan for the
  COVID smoke test (now closed; superseded by `covid_analysis/REPORT.md`)
- [Diary.md](Diary.md) — Layer-1 session log (BCKM replication era)
- [bca-web-app-instructions.md](bca-web-app-instructions.md) —
  eventual web-app spec (out of scope until Layer 2 lands; now
  unblocked)
- [bckm_replication/bca_paper.pdf](bckm_replication/bca_paper.pdf) —
  full paper

## Methodology in one screen

State `[k_hat, A_hat, taul_hat, taux_hat, g_hat]`; observables
`[y_hat, l_hat, x_hat, g_hat]` as log-deviations from the model
deterministic steady state. Model is log-linearized around its SS and
solved via Klein's QZ method.

Estimation jointly fits the VAR(1) transition matrix **P** (4×4), the
unconditional-mean reparametrization **Sbar** so that
`P_0 = (I − P) · Sbar`, and the shock Cholesky factor **Q** (10
lower-triangular elements) by maximizing the Kalman-filter log-
likelihood of the four observables. Latent wedge paths come from a
Rauch–Tung–Striebel smoother. Stationarity is enforced by a single
spectral-radius soft penalty (`5e5 · max(|eig P| − 0.995, 0)²`,
matching `mleqadj.m:134`). The Kalman filter uses the steady-state
DARE-derived constant gain, recomputed at every objective evaluation.
The full optimizer-result dict is content-addressed-cached as a
pickle under `*.mle.pkl` (gitignored) so cold runs of ~10–15 minutes
become re-runs of ~3 seconds when inputs are unchanged.

Counterfactuals keep the full 4-D VAR for rational expectations and
zero out *inactive* wedge columns in the structural equations. Capital
evolves endogenously via `k' = P_k @ state`. Per-wedge
counterfactuals are **incremental**: `(C_j − C_0)` is the right
operator (BCKM `gwedges2.m:90-115`), so the four single-wedge CFs
sum to ≈ all-active ≈ data.

## Calibration

BCKM Table 1 values: `α = 1/3, ψ = 2.5, δ = 0.05/yr, β = 0.975/yr,
σ = 1`. Growth rates are data-derived per call (calgz fsolve for `gz`,
working-age-pop slope for `gn`). `g_share` is set from data
(`mean(g) / mean(y)`); adjustment-cost coefficient `a` from BGG (1999)
elasticity of q wrt x/k = 0.25. For non-BCKM windows the labor
construction is anchored to the BCKM-empirical hours-per-capita mean
0.24279 (model `ss["l"]` ≈ 0.29 vs raw FRED hours/pop ≈ 23 — the
anchor avoids an 80× scale mismatch in the wedge extraction).

## Data

FRED defaults (with BEA branches as diagnostic-only opt-in — see
[bckm_replication/REPORT.md](bckm_replication/REPORT.md)
§"BEA NIPA migration"):

- National accounts: `GDP`, `PCE`, `PCDG`, `GPDI`, `GCE`, `NETEXP`,
  `GDPDEF`
- Sales tax: `ASLSTAX`
- Population: `LFWA64TTUSQ647N` (working-age 15–64)
- Hours / employment:
  - **BLS-faithful (preferred when available)**: `CE16OV`
    (BLS LNS12000000, CPS civilian employment, ages 16+) and `AWHAETP`
    (avg weekly hours, total private, all employees) — combined as
    `employment × avg_weekly_hours × 13 weeks/qtr`, mirroring BCKM
    `usdata.m:hours.dat`.
  - **Legacy fallback**: `PAYEMS × AWHNONAG`, then `HOANBS`, then
    `PRS85006023` index.

Adjustments follow BCKM Section 4: durables reclassified from C to X
via perpetual-inventory stock with a 4% imputed service flow added to
both C and Y; sales tax subtracted from Y and C; series deflated and
converted to per-working-age-capita; calgz-style trend removed from
`log(y)`.

## Tests

```bash
# Full suite (91 tests; some are slow — load parquet, run MLE)
.venv/bin/python -m pytest tests/ -v

# Fast suite — skip the long MLE-driven tests
.venv/bin/python -m pytest tests/ -m "not slow" -q

# Layer-2 work — also skip BCKM-1980-2014-specific regression tests
.venv/bin/python -m pytest tests/ -m "not bckm and not slow" -q
```

Marker conventions:
- `bckm` — tests pertinent only to the BCKM 1980Q1–2014Q4 regression
  (load `worktemp.mat` or pin BCKM-published values).
- `slow` — tests that load a parquet and run MLE (~10 min each).

## What this repo is not

- The `BCKM/` directory (gitignored, never modified) contains the
  original Matlab replication files. Reference-only — see
  [CLAUDE.md](CLAUDE.md).
- This is the computational core only. A FastAPI service and a React
  front end are sketched in
  [bca-web-app-instructions.md](bca-web-app-instructions.md) but are
  not implemented; the Layer-2 validation that gated them is now
  complete, so the web layer is the next major milestone.

## References

- Brinca, Chari, Kehoe, McGrattan (2016), "Accounting for Business
  Cycles," NBER WP 22663.
- Chari, Kehoe, McGrattan (2007), "Business Cycle Accounting,"
  *Econometrica* 75(3).
- Brinca, Costa Filho, Loria (2024), "Business cycle accounting:
  What have we learned so far?" *Journal of Economic Surveys* 38(4).
- Brinca, Iskrev, Loria (2022), "On Identification Issues in Business
  Cycle Accounting Models."
- Klein (2000), "Using the generalized Schur form to solve a
  multivariate linear rational expectations model."
- Bernanke, Gertler, Gilchrist (1999), "The financial accelerator in
  a quantitative business cycle framework."
