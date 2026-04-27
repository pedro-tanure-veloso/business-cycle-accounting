# BCA — Business Cycle Accounting for the US

A Python port of the Business Cycle Accounting methodology of Chari, Kehoe & McGrattan (2007) and Brinca, Chari, Kehoe & McGrattan (2016, NBER WP 22663). Decomposes US aggregate fluctuations into four structural wedges — efficiency (A), labor (1−τ_l), investment (1+τ_x), and government (g) — by Kalman-filter MLE on a VAR(1) and then runs counterfactual simulations.

The validation target is the BCKM Section 7 US MLE estimates (Tables 8–11 in [BCA_info.md](BCA_info.md)), in particular the Great Recession decomposition where the investment wedge dominates the fall in investment. Once those tables are matched, the model is treated as correctly implemented and used for ongoing US analysis.

## Status

Active development. The model solves correctly, the pipeline runs end-to-end, and the investment wedge correctly worsens during the Great Recession. The MLE results are qualitatively right but quantitative targets (P matrix diagonal, P₀ vector, φ-statistics) are not yet matched. Open items and the current "next steps" list live in [Diary.md](Diary.md); a head-to-head against the BCKM Matlab reference is in [DIVERGENCE_ANALYSIS.md](DIVERGENCE_ANALYSIS.md).

## Quick start

```bash
pip install -e ".[dev]"

# First run: fetch from FRED, save processed dataset
FRED_API_KEY=... python scripts/run_var_counterfactuals.py --save-data data/us_1980_2014.parquet

# Subsequent runs: use the saved parquet, no API key needed
python scripts/run_var_counterfactuals.py --data data/us_1980_2014.parquet
```

Output: VAR(1) estimates, smoothed wedge paths, φ-statistics, peak-to-trough Great Recession decomposition, and `figure_2B_mle.png`.

A FRED API key is only needed for the initial data fetch. Get one at [fred.stlouisfed.org/docs/api/api_key.html](https://fred.stlouisfed.org/docs/api/api_key.html). Cached series live in `~/.bca_cache/fred/`; processed datasets live in `data/`.

## Repository layout

```
bca_core/                   pure-Python computational core (no web deps)
  params.py                 CalibrationParams (BCKM Table 77 defaults)
  model.py                  Prototype model: SS, log-linearization, decision rules
  klein.py                  Klein (2000) QZ solver for linear RE models
  wedges.py                 Static wedge extraction from data
  var_estimation.py         Kalman-filter MLE (BCKM mleqadj.m logic)
  counterfactuals.py        Wedge-alone simulations + φ-statistics
  data/
    fred.py                 FRED fetcher with disk caching
    adjustments.py          Durables PIM, sales-tax, per-capita conversion
    pipeline.py             End-to-end fetch → adjust → detrend
scripts/
  run_var_counterfactuals.py    Main entry point — full pipeline + Figure 2B
  diagnose_counterfactuals.py   PASS/FAIL diagnostics for the CF decomposition
  plot_wedges.py                Static wedge plots (BCKM Figure 2 style)
tests/                      pytest unit tests (52 tests)
data/                       processed parquet datasets + metadata sidecars
```

## Documentation map

- [CLAUDE.md](CLAUDE.md) — working rules and methodology constraints (read before changing anything in `bca_core/`)
- [BCA_info.md](BCA_info.md) — paper summary with the target parameter tables (Sections 1–7)
- [REPORT.md](REPORT.md) — implementation write-up, results, "things not to do"
- [DIVERGENCE_ANALYSIS.md](DIVERGENCE_ANALYSIS.md) — head-to-head Python vs. BCKM Matlab
- [Diary.md](Diary.md) — session-by-session log; current next steps live here
- [counterfactual-debugging-summary.md](counterfactual-debugging-summary.md) — past CF-decomposition bug fix log
- [bca-web-app-instructions.md](bca-web-app-instructions.md) — eventual web-app spec (out of scope for now)
- [bca_paper.pdf](bca_paper.pdf) — full paper

## Methodology in one screen

State `[k_hat, A_hat, taul_hat, taux_hat, g_hat]`; observables `[y_hat, l_hat, x_hat, g_hat]` as log-deviations from the model deterministic steady state. The model is log-linearized around its SS and solved via Klein's QZ method.

Estimation jointly fits the VAR(1) transition matrix **P** (4×4), the unconditional-mean reparametrization **Sbar** so that `P_0 = (I−P)·Sbar`, and the shock Cholesky factor **Q** (10 lower-triangular elements) by maximizing the Kalman-filter log-likelihood of the four observables. Latent wedge paths come from a Rauch–Tung–Striebel smoother. Stationarity is enforced by spectral-radius and per-diagonal soft penalties. The Kalman filter is initialized via DARE (not Lyapunov, not diffuse) using BCKM Table 77 parameters — held fixed during optimization for speed, recomputed at final parameters for the smoother.

Counterfactuals keep the full 4-D VAR for rational expectations and zero out *inactive* wedge columns in the structural equations. Capital evolves endogenously via `k' = P_k @ state`.

## Calibration

Calibration uses BCKM Table 77 values (annualized): α=0.35, ψ=2.24, δ=0.0464, ρ=0.0286, γ=1.9%, n=0.98%. `g_share` is set from data (mean(g)/mean(y)). Adjustment-cost coefficient `a` is pinned by an elasticity of q with respect to x/k of 0.25 (BGG 1999).

## Data

US quarterly data, currently 1980Q1–2014Q4 (the BCKM Section 7 sample) — fetched from FRED:

- National accounts: `GDP`, `PCE`, `PCDG` (durables), `GPDI`, `GCE`, `NETEXP`, `GDPDEF`
- Sales tax: `ASLSTAX`
- Population: `LFWA64TTUSQ647S` (working-age 15–64)
- Hours / employment: `PAYEMS`, `AWHNONAG`, `PRS85006023`

Adjustments follow BCA_info Section 4: durables reclassified from C to X via perpetual-inventory stock with a 4% imputed service flow added to both C and Y; sales tax subtracted from Y and C; series deflated and converted to per-working-age-capita; linear trend removed from `log(y)`.

## Tests

```bash
pytest tests/ -v
```

52 tests covering the model SS and log-linearization, the Klein solver, Kalman-filter behaviour, wedge extraction, and counterfactual signs/decompositions. Run them before reporting any change as done.

## What this repo is not

- The `BCKM/` directory (gitignored, never modified) contains the original Matlab replication files. They are reference-only — see [CLAUDE.md](CLAUDE.md).
- This is the computational core only. A FastAPI service and a React front end are sketched in [bca-web-app-instructions.md](bca-web-app-instructions.md) but are not implemented.

## References

- Brinca, Chari, Kehoe, McGrattan (2016), "Accounting for Business Cycles," NBER WP 22663.
- Chari, Kehoe, McGrattan (2007), "Business Cycle Accounting," *Econometrica* 75(3).
- Brinca, Costa Filho, Loria (2024), "Business cycle accounting: What have we learned so far?" *Journal of Economic Surveys* 38(4).
- Brinca, Iskrev, Loria (2022), "On Identification Issues in Business Cycle Accounting Models."
- Klein (2000), "Using the generalized Schur form to solve a multivariate linear rational expectations model."
- Bernanke, Gertler, Gilchrist (1999), "The financial accelerator in a quantitative business cycle framework."
