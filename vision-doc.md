---
title: "Business Cycle Accounting — Project Vision"
topic: "vision"
layer: "all"
status: "active"
last_updated: "2026-05-03"
---

# Business Cycle Accounting — Project Vision

## Goal

Build a **US-only web application** that performs Business Cycle Accounting (BCA) in the style of
Chari–Kehoe–McGrattan (2007) and Brinca–Chari–Kehoe–McGrattan (2016, NBER WP 22663). The user
picks a US time window, the app runs the full BCA procedure end-to-end, and returns the four
wedges (efficiency, labor, investment, government) plus counterfactual decompositions of output,
labor, and investment.

Cross-country support is explicitly **out of scope**.

## Primary references

- Brinca, Chari, Kehoe, McGrattan (2016), "Accounting for Business Cycles," NBER WP 22663. <https://www.nber.org/papers/w22663>
- Full technical appendix and replication files: <https://pedrobrinca.pt/2016-accounting-for-business-cycles/>
- Brinca, Costa Filho, Loria (2024), "Business cycle accounting: What have we learned so far?" *Journal of Economic Surveys* 38(4) — for the modern synthesis.
- Brinca, Iskrev, Loria (2022), "On Identification Issues in Business Cycle Accounting Models" — read before trusting MLE estimates.

## Architecture

Three layers, kept loosely coupled so each can be tested independently:

1. **`bca_core`** ✅ **Complete** — pure-Python package implementing the prototype model, data
   adjustments, wedge extraction, Kalman-filter VAR MLE, and counterfactuals. No web dependencies.
2. **`bca_api`** — FastAPI service that exposes `bca_core` over HTTP. Long-running MLE jobs run
   asynchronously; results are cached by `(sample_hash, calibration_hash)`.
3. **`bca_web`** — Single-page front end (React + Recharts, or Streamlit for a faster prototype).

---

## Stage 1 — `bca_core` ✅ Done

`bca_core` is a fully validated Python port of the BCKM MATLAB computational core. It provides:

- **Prototype model**: closed-economy one-sector growth model with four structural wedges
  (efficiency A, labor 1−τ_l, investment 1+τ_x, government g), log-linearized via Klein QZ,
  with BGG adjustment costs calibrated to elasticity 0.25.
- **Data pipeline**: FRED-backed US data with consumer durables reclassification, sales-tax
  adjustment, per-capita deflation, and calgz exponential detrending. Supports arbitrary US
  windows with window-specific γ/n calibration and an optional `mle_window` pre-anomaly anchor.
- **VAR(1) MLE**: steady-state Kalman filter (DARE-based constant gain), RTS smoother, 30-parameter
  θ = [P₀(4), P(16), Q_lower_tri(10)], spectral-radius penalty matching BCKM.
- **Counterfactuals + φ-statistics**: BCKM-style incremental single-wedge decompositions
  (C_j − C₀ per `gwedges2.m`); level-space SSR φ-statistic summing to one across wedges per variable.

**Validation status:**
- Layer 1 (BCKM US 1980Q1–2014Q4): f-stats match Table 11 to ≤0.01 in every channel. ✅
- Layer 2 (COVID-era 2010Q1–2023Q4): all six qualitative priors pass (labor wedge dominates
  trough; investment counter-cyclical; A spike in 2021 recovery). ✅

For implementation details, calibration constants, and known pitfalls see `CLAUDE.md` and
`bckm_replication/REPORT.md`.

---

## Stage 2 — `bca_api`: FastAPI service

Endpoints (minimum viable):

- `POST /runs` — body: `{start, end, calibration_overrides?}`. Returns `{run_id, status}`. Kicks
  off estimation asynchronously.
- `GET /runs/{run_id}` — returns status and, when complete, the full result: wedges,
  counterfactuals, φ-statistics, metadata.
- `GET /data` — returns the raw and adjusted US data series for transparency.

Persist results in SQLite or Postgres keyed by `(sample_hash, calibration_hash)` so repeated
runs are served from cache. MLE estimation takes ~2–12 minutes cold; a job queue (RQ or Celery)
is required for anything beyond a toy deployment.

`bca_core.data` already implements a FRED adapter with disk caching. The API layer calls
`build_us_dataset()` and passes the result to `estimate_var_mle()`.

---

## Stage 3 — `bca_web`: front end

Minimum viable UI:

1. **Selector panel.** US sample start/end pickers, "Run" button. Advanced section (collapsed by
   default) for calibration overrides (γ, n, adjustment cost, detrend method, mle_window).
2. **Results tab 1 — Wedges.** Four time-series plots, one per wedge, with recession shading.
3. **Results tab 2 — Counterfactuals.** For each of output, labor, investment: actual series plus
   the four single-wedge counterfactual paths overlaid.
4. **Results tab 3 — Summary.** φ-statistic table and peak-to-trough decomposition. Let the user
   specify the recession window interactively.
5. **Data tab.** Adjusted input series so the user can sanity-check the data construction.

Recharts or Plotly are both fine. For React: TanStack Query for server state.

---

## Known pitfalls

- **Do not** identify VAR innovations with wedges in the counterfactuals. Keep the distinction
  between the underlying state `s_t` and the wedge-alone counterfactual paths. This is the
  single biggest conceptual error in the literature.
- **Do not** skip the durables reclassification. It changes results materially and is required
  for consistency with BCKM.
- **Do not** use a perfect-foresight procedure and claim it's BCKM. BCKM Section 3.B.iv shows
  this inflates the US labor-wedge contribution from 46% to 92%.
- **Do not** assume the MLE has a unique global mode. Run from multiple starting points and
  report. Identification can be fragile — see Brinca–Iskrev–Loria (2022).
- The US result is sample- and adjustment-cost-sensitive (BCKM footnote 2). Extending the
  sample or changing `a` will move results — that is correct behavior, not a bug.

---

## Milestones

| # | Milestone | Status |
|---|-----------|--------|
| 1 | Prototype model + steady state + log-linear solution | ✅ Done |
| 2 | US data pipeline — four observables as a clean DataFrame | ✅ Done |
| 3 | Wedge extraction — plots match BCKM Figure 2 | ✅ Done |
| 4 | VAR MLE — parameters match BCKM Tables 8/9/10 | ✅ Done |
| 5 | Investment wedge extraction using solved VAR | ✅ Done |
| 6 | Counterfactuals + φ-statistics — f-stats match BCKM Table 11 to ≤0.01 | ✅ Done |
| 7 | Wrap `bca_core` in FastAPI (`bca_api`) | ⬜ Next |
| 8 | Build the front end (`bca_web`) | ⬜ |
| 9 | End-to-end validation on arbitrary US windows in the deployed app | ⬜ |
