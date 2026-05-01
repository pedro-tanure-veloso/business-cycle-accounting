# BCA — Business Cycle Accounting

A Python implementation of the Business Cycle Accounting methodology of
Chari, Kehoe & McGrattan (2007) and Brinca, Chari, Kehoe & McGrattan
(2016, NBER WP 22663). Decomposes aggregate fluctuations into four
structural wedges — efficiency (A), labor (1−τ_l), investment (1+τ_x),
and government (g) — by Kalman-filter MLE on a VAR(1), with
counterfactual simulations and Brinca-Iskrev-Loria-style f-statistics.

The end deliverable is a toolkit for **arbitrary country/period
combinations**. The BCKM 2016 US 1980Q1–2014Q4 replication is the
locked-in regression test that pins fidelity; new windows (e.g. COVID
2010Q1–2023Q4) and new countries are validated against narrative priors.

## Status

- **BCKM 2016 US 1980-2014 replication**: closed. f-statistics match
  Table 11 to ≤ 0.01 in every channel at BCKM-θ on our dataset; full
  wrap-up in [`bckm_replication/REPORT.md`](bckm_replication/REPORT.md).
- **Layer 2 — generalizability smoke tests**: active. COVID
  2010Q1–2023Q4 is the current target (see
  [Diary.md](Diary.md) for the live plan).
- **Layer 3 — cross-country (OECD MEI Tables III/IV)**: future.

## Quick start

```bash
pip install -e ".[dev]"

# First run: fetch from FRED, save processed dataset
FRED_API_KEY=... python scripts/run_var_counterfactuals.py \
    --save-data bckm_replication/data/us_1980_2014.parquet

# Subsequent runs: use the saved parquet, no API key needed
python scripts/run_var_counterfactuals.py \
    --data bckm_replication/data/us_1980_2014.parquet
```

A FRED API key is only needed for the initial data fetch. Get one at
[fred.stlouisfed.org/docs/api/api_key.html](https://fred.stlouisfed.org/docs/api/api_key.html).
Cached series live in `~/.bca_cache/fred/`; the BCKM-replication
parquets live in `bckm_replication/data/`; new datasets land in
`data/`.

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
  var_estimation.py         Kalman-filter MLE (BCKM mleqadj.m logic)
  counterfactuals.py        Wedge-alone simulations + f-statistics
  data/
    fred.py                 FRED fetcher with disk caching
    bea.py                  BEA NIPA + Fixed-Asset fetcher (diagnostic-only)
    adjustments.py          Durables PIM, sales-tax, per-capita conversion, calgz detrending
    pipeline.py             End-to-end fetch → adjust → detrend (build_us_dataset)

scripts/                    generic pipeline drivers (window-agnostic)
  run_var_counterfactuals.py  Main entry point — full pipeline + Figure 2B
  diagnose_counterfactuals.py PASS/FAIL diagnostics for the CF decomposition
  solow_residual_check.py     Solow residual sanity check

bckm_replication/           BCKM 2016 US 1980-2014 ground-truth artifacts
  REPORT.md                 final replication wrap-up (closed 2026-05-01)
  BCKM_DIFF_GUIDE.md        element-wise diff guide vs worktemp.mat
  DATA_FORENSICS.md         BEA NIPA migration walkdown
  DIVERGENCE_ANALYSIS.md    earlier head-to-head analysis
  counterfactual-debugging-summary.md   CF-fix session record
  matlab_reference/         BCKM 2016 matlab code (paper ground truth)
  octave_output/            octave dumps used as fixtures
  data/                     pinned parquets + bootstrap/sensitivity dumps
  scripts/                  diag_*, eval_bckm_*, compare_*, plot_*, bootstrap_*, sensitivity_*
  figures/                  Figure 2A–E reproductions, Solow residual, etc.

tests/                      pytest unit tests (79 tests; mostly window-agnostic)
data/                       cache for new (Layer 2/3) datasets
```

## Documentation map

- [CLAUDE.md](CLAUDE.md) — working rules and methodology constraints
  (read before changing anything in `bca_core/`); also carries the
  append-only "Findings" journal
- [BCA_info.md](BCA_info.md) — paper summary with target tables
- [bckm_replication/REPORT.md](bckm_replication/REPORT.md) — final
  BCKM-replication wrap-up: bugs found and fixed, methodology
  decisions, residual issues
- [Diary.md](Diary.md) — session-by-session log; current next steps
  live here
- [bca-web-app-instructions.md](bca-web-app-instructions.md) —
  eventual web-app spec (out of scope until Layer 2 lands)
- [bca_paper.pdf](bca_paper.pdf) — full paper

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
elasticity of q wrt x/k = 0.25.

## Data

FRED defaults (with BEA branches as diagnostic-only opt-in — see
[bckm_replication/REPORT.md](bckm_replication/REPORT.md)
§"BEA NIPA migration"):

- National accounts: `GDP`, `PCE`, `PCDG`, `GPDI`, `GCE`, `NETEXP`,
  `GDPDEF`
- Sales tax: `ASLSTAX`
- Population: `LFWA64TTUSQ647N` (working-age 15–64)
- Hours / employment: `PAYEMS`, `AWHNONAG`, `PRS85006023`

Adjustments follow BCKM Section 4: durables reclassified from C to X
via perpetual-inventory stock with a 4% imputed service flow added to
both C and Y; sales tax subtracted from Y and C; series deflated and
converted to per-working-age-capita; calgz-style trend removed from
`log(y)`.

## Tests

```bash
pytest tests/ -v
```

79 tests covering the model SS and log-linearization, the Klein solver,
Kalman-filter behaviour, wedge extraction, counterfactual
signs/decompositions, the BCKM Table 12 peak-trough decomposition pin,
and the algebraic-identity invariants (T1–T7) on the cached
1980-2014 parquet. Run them before reporting any change as done.

## What this repo is not

- The `BCKM/` directory (gitignored, never modified) contains the
  original Matlab replication files. Reference-only — see
  [CLAUDE.md](CLAUDE.md).
- This is the computational core only. A FastAPI service and a React
  front end are sketched in
  [bca-web-app-instructions.md](bca-web-app-instructions.md) but are
  not implemented (premature until Layer 2 generalizability is
  validated).

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
