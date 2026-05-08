---
title: "BCA — Business Cycle Accounting"
topic: "project-overview"
layer: "all"
status: "active"
last_updated: "2026-05-08"
---

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

## Quick start

```bash
# Python 3.10+; create a venv to keep deps isolated
python3.12 -m venv .venv
.venv/bin/pip install -e ".[dev]"

# Run BCA on any US window.
# A free FRED API key is required for the first data fetch;
# get one at https://fred.stlouisfed.org/docs/api/api_key.html
export FRED_API_KEY=your_key_here

python scripts/run_bca.py \
    --start 2000Q1 --end 2015Q4 --base 2007Q4 \
    --window 2007Q4 2009Q2 \
    --output-dir ./output
```

On the first run the data is fetched from FRED and cached automatically.
Every subsequent run with the same window takes **~3 seconds** (MLE
result is also cached). After the first run no API key is needed:

```bash
python scripts/run_bca.py \
    --start 2000Q1 --end 2015Q4 --base 2007Q4 \
    --data ./output/us_2000Q1_2015Q4_2007Q4.parquet \
    --output-dir ./output
```

The script produces 10 files in `--output-dir`:

| File | Description |
|---|---|
| `wedges.png` | 4-panel: efficiency, labor, investment, government (index base=100) |
| `figure_A.png` | Output, labor, investment data overlay |
| `figure_B.png` | Output + three wedges (BCKM Figure 2B style) |
| `figure_2C/2D/2E.png` | Per-wedge counterfactual decompositions |
| `wedges.csv` | Wedge level series |
| `fstats.csv` | φ/f-statistics table |
| `counterfactuals.csv` | Per-quarter CF paths for all wedges × variables |
| `peak_to_trough.csv` | Log-change decomposition from base to window end |

Run `python scripts/run_bca.py --help` for the full argument reference.

## Status

| Layer | Description | Status |
|---|---|---|
| 1 | BCKM 2016 US 1980Q1–2014Q4 replication | ✅ Closed (2026-05-01) |
| 2 | Generalizability — COVID 2010Q1–2023Q4 | ✅ Passed (2026-05-02) |
| 3 | Cross-country | 🚫 Out of scope |

**Layer 1**: f-statistics match Table 11 to ≤ 0.01 in every channel at
BCKM-θ on our dataset. Full wrap-up in
[`bckm_replication/REPORT.md`](bckm_replication/REPORT.md).

**Layer 2**: 6/6 narrative-prior rubric checks pass (both trend
variants). Labor wedge accounts for 49%/77% of the output/hours f-stat
weight over 2019Q4–2022Q4 — pipeline correctly identifies COVID as
primarily a labor-driven contraction. Details in
[`covid_analysis/REPORT.md`](covid_analysis/REPORT.md).

## Repository layout

```
scripts/                    User-facing entry points
  run_bca.py                ★ Main CLI — any US window, all outputs in one command
  build_quarterly_data.py   Stage 2 data pipeline (runs bca_core + Gemini LLM)
  generate_events.py        Per-quarter macro event log via Gemini + Google Search grounding
  run_var_counterfactuals.py  Legacy BCKM-1980-2014 driver (parquet + figures)
  diagnose_counterfactuals.py PASS/FAIL diagnostics for the CF decomposition
  solow_residual_check.py     Solow residual sanity check

data/                       Curated inputs to the dashboard pipeline
  events.md                 Per-quarter macro event log (auto-generated, frozen)

bca_core/                   Pure-Python computational core (no web deps)
  params.py                 CalibrationParams (BCKM Table 1 defaults)
  constants.py              Canonical BCKM Table 8/9/10 parameters (code convention)
  model.py                  Prototype model: SS, log-linearization, decision rules
  klein.py                  Klein (2000) QZ solver for linear RE models
  bckm_lom.py               BCKM-faithful capital LOM (fixexpadj.m port)
  wedges.py                 Static wedge extraction from data
  var_estimation.py         Kalman-filter MLE (BCKM mleqadj.m) + result caching
  counterfactuals.py        Wedge-alone simulations + φ/f-statistics
  data/
    fred.py                 FRED fetcher with disk caching (incl. CPS / BLS labor)
    bea.py                  BEA NIPA + Fixed-Asset fetcher (diagnostic-only)
    adjustments.py          Durables PIM, sales-tax, per-capita, calgz, BLS labor
    pipeline.py             End-to-end fetch → adjust → detrend (build_us_dataset)

bckm_replication/           BCKM 2016 US 1980-2014 ground-truth artifacts
  REPORT.md                 Final replication wrap-up (closed 2026-05-01)
  BCKM_RESULTS.md           Paper Parts III & IV — regression-target tables
  BCKM_DIFF_GUIDE.md        Element-wise diff guide vs worktemp.mat
  DATA_FORENSICS.md         BEA NIPA migration walkdown
  bca_paper.pdf             Full paper PDF
  matlab_reference/         BCKM 2016 Matlab code (paper ground truth, read-only)
  octave_output/            Octave dumps used as test fixtures
  data/                     Pinned parquets + bootstrap/sensitivity results
  scripts/                  diag_*, eval_bckm_*, compare_*, bootstrap_*, sensitivity_*
  figures/                  Figure 2A–E reproductions, Solow residual, etc.

covid_analysis/             Layer 2 — COVID 2010Q1–2023Q4 smoke test
  scripts/run_covid_analysis.py  One-shot driver (datasets → MLE → CFs → figures)
  data/                     Cached parquets + .meta.json + .mle.pkl (gitignored)
  figures/                  figure_A/B/2C/2D/2E_covid*.png + wedges_*.png
  REPORT.md                 Narrative results + COVID-vs-Great-Recession comparison
  Diary.md                  Session-by-session journal

bca_web/                    Stage 3 — React Dashboard (Static Web App)
  public/data/              Static JSON files consumed by the UI
  src/                      React components and CSS design system

tests/                      pytest unit tests (91 tests)
  test_bckm_*               BCKM 1980-2014 regression (marked @pytest.mark.bckm)
  test_covid_analysis.py    COVID dataset shape + structural identities + smoke
  test_*                    Window-agnostic structural tests
```

## CLI reference

```
python scripts/run_bca.py --start YYYYQN --end YYYYQN --base YYYYQN [options]

Required:
  --start YYYYQN          Sample start, e.g. 2000Q1
  --end   YYYYQN          Sample end,   e.g. 2015Q4
  --base  YYYYQN          Normalization period (index = 100), e.g. 2007Q4

Data:
  --data PATH             Load pre-built parquet (no API key needed)
  --fred-api-key KEY      FRED key (alternative to FRED_API_KEY env var)

Options:
  --window START END      Evaluation window for f-statistics (default: base→end)
  --mle-window START END  Restrict trend fitting to a sub-window (e.g. pre-COVID)
  --output-dir DIR        Output directory (default: cwd)
  --no-cache-mle          Force MLE re-run (bypass pickle cache)
  --labor-target FLOAT    Hours/pop anchor (default: 0.24279, BCKM-empirical)
  --quiet                 Suppress verbose MLE output
```

**No FRED key?** The script prints a step-by-step guide and the exact
`--data` path to reuse once the data is cached.

## Documentation map

- [`README_AI.md`](README_AI.md) — AI/agent context map: full glossary,
  pipeline logic flow, key data structures, critical conventions.
  Start here if you are an LLM working on this codebase.
- [`llms.txt`](llms.txt) — LLM project index (llmstxt.org standard):
  links every file with a one-line description.
- [`BCA_info.md`](BCA_info.md) — Paper summary and **live methodology**
  (Parts I & II): prototype economy, accounting procedure, data and
  calibration. Applies to any US window the toolkit runs on.
- [`vision-doc.md`](vision-doc.md) — Project roadmap: current milestone
  status and the next steps (FastAPI service + React front end).
- [`CLAUDE.md`](CLAUDE.md) — Working rules, methodology constraints, and
  the append-only "Findings" debugging journal. Read before changing
  anything in `bca_core/`.
- [`bckm_replication/REPORT.md`](bckm_replication/REPORT.md) — Final
  BCKM-replication wrap-up: bugs found and fixed, methodology decisions,
  residual issues.
- [`bckm_replication/BCKM_RESULTS.md`](bckm_replication/BCKM_RESULTS.md)
  — Paper Parts III & IV (regression-test ground truth tables).
- [`covid_analysis/REPORT.md`](covid_analysis/REPORT.md) — COVID
  smoke-test results and Part III-style comparison with the Great Recession.

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

FRED defaults (with BEA branches as diagnostic-only opt-in):

- National accounts: `GDP`, `PCE`, `PCDG`, `GPDI`, `GCE`, `NETEXP`,
  `GDPDEF`
- Sales tax: `ASLSTAX`
- Population: `LFWA64TTUSQ647N` (working-age 15–64)
- Hours / employment:
  - **BLS-faithful (preferred when available)**: `CE16OV`
    (BLS LNS12000000, CPS civilian employment, ages 16+) ×`AWHAETP`
    (avg weekly hours, total private) × 13 weeks/qtr, mirroring BCKM
    `usdata.m:hours.dat`. Active when both series cover the sample.
  - **Legacy fallback**: `PAYEMS × AWHNONAG`, then `HOANBS`.

Adjustments follow BCKM Section 4: durables reclassified from C to X
via perpetual-inventory stock with a 4% imputed service flow added to
both C and Y; sales tax subtracted from Y and C; series deflated and
converted to per-working-age-capita; calgz-style exponential trend
removed from `log(y)`.

## Tests

```bash
# Full suite (91 tests; some are slow — run MLE from scratch)
.venv/bin/python -m pytest tests/ -v

# Fast suite — skip long MLE-driven tests (~1 second)
.venv/bin/python -m pytest tests/ -m "not slow" -q

# Layer-2 work — also skip BCKM-1980-2014-specific regression tests
.venv/bin/python -m pytest tests/ -m "not bckm and not slow" -q
```

Marker conventions:
- `bckm` — tests tied to the BCKM 1980Q1–2014Q4 regression
  (load `worktemp.mat` or pin BCKM-published values).
- `slow` — tests that load a parquet and run MLE cold (~10 min each).

## What this repo is not

- The `BCKM/` directory (gitignored, never modified) contains the
  original Matlab replication files. Reference-only — see
  [CLAUDE.md](CLAUDE.md).
- This is the computational core only. The web-layer roadmap (FastAPI
  service + React front end) is in
  [`vision-doc.md`](vision-doc.md); the Layer-2 validation that
  gated it is complete and the web layer is the next major milestone.

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
