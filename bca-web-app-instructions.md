# Business Cycle Accounting Web App — Build Instructions

## Project goal

Build a web application that performs Business Cycle Accounting (BCA) in the style of Chari–Kehoe–McGrattan (2007) and Brinca–Chari–Kehoe–McGrattan (2016, NBER WP 22663). The app should let a user select a country and sample period, run the BCA procedure end-to-end, and return the four wedges (efficiency, labor, investment, government) plus counterfactual decompositions of output, labor, and investment.

The canonical reference implementation is Pedro Brinca's MATLAB toolbox `BCAppIt!` and the replication files for BCKM (2016). We are porting the computational core to Python and wrapping it in a web app. We are **not** porting the MATLAB GUI.

## Primary references

- Brinca, Chari, Kehoe, McGrattan (2016), "Accounting for Business Cycles," NBER WP 22663. <https://www.nber.org/papers/w22663>
- Full technical appendix and replication files: <https://pedrobrinca.pt/2016-accounting-for-business-cycles/>
- `BCAppIt!` MATLAB toolbox (for cross-checking): <https://pedrobrinca.pt/software/bcappit-2/>
- Brinca, Costa Filho, Loria (2024), "Business cycle accounting: What have we learned so far?" *Journal of Economic Surveys* 38(4) — for the modern synthesis and extensions.
- Brinca, Iskrev, Loria (2022), "On Identification Issues in Business Cycle Accounting Models" — read before trusting your MLE estimates.

## Architecture

Three layers, kept loosely coupled so each can be tested independently:

1. **`bca_core`** — a pure-Python package that implements the prototype model, data adjustments, wedge extraction, VAR estimation, and counterfactuals. No web dependencies. This is the port of the MATLAB core.
2. **`bca_api`** — a FastAPI service that exposes `bca_core` over HTTP. Long-running jobs (MLE estimation) should run asynchronously; consider a job queue (Celery + Redis, or RQ for simplicity) if estimation takes more than a few seconds.
3. **`bca_web`** — a front end. A React + Recharts single-page app is the default; Streamlit or Dash are acceptable simpler alternatives if the user wants a prototype rather than a production app.

## Stage 1 — `bca_core`: the computational engine

### 1.1 Prototype model

Implement the closed-economy one-sector stochastic growth model from BCKM 2016, Section 1.A:

- Production: `y_t = A_t * k_t^α * ((1+γ)^t * l_t)^(1-α)`
- Utility: `U(c, l) = log(c) + ψ * log(1 - l)`
- Capital accumulation with adjustment cost: `(1+n) k_{t+1} = (1-δ) k_t + x_t - φ(x_t / k_t) * k_t`, with `φ(x/k) = (a/2)(x/k - b)^2` and `b = γ + δ + n`.
- Four wedges: `A_t` (efficiency), `1 - τ_{lt}` (labor), `1/(1 + τ_{xt})` (investment), `g_t` (government).
- Government wedge in the data = government consumption + net exports (per the closed-economy equivalence in CKM 2005).

### 1.2 Calibration (defaults; expose as parameters)

- Capital share `α = 1/3`
- Time allocation `ψ = 2.5`
- Annualized: depreciation `δ = 5%`, rate of time preference `ρ = 2.5%` → `β = 1/(1+ρ)^(1/4)` for quarterly
- Population growth `n` and technology growth `γ` — estimate from the data for the selected country
- Adjustment cost `a` calibrated so elasticity of `q` w.r.t. `x/k` is 0.25 (following BGG 1999), i.e., `a = 0.25 / b`

### 1.3 Data adjustments (get these right — they matter)

- **Consumer durables reclassification.** Strip durables expenditure out of consumption; add it to investment. Build a durables capital stock by perpetual inventory with the same `δ`. Add an imputed service flow (4% of the durables stock) to **both** consumption and output.
- **Sales taxes.** Subtract sales tax revenue from both consumption and output.
- **Per capita and real.** Divide all nominal aggregates by the GDP deflator, then by working-age population.
- **Trends.** Remove a country-specific linear trend from log output, log investment, and the log government wedge before estimation.

### 1.4 Log-linearized solution

Log-linearize the model around the deterministic steady state. Two acceptable implementations:

- **Preferred:** hand-code the state-space matrices. The model is small (one endogenous state — capital — plus the four exogenous wedges). Write out the log-linearized Euler equation, the labor FOC, the resource constraint, and the production function; solve using Klein (2000) or Sims' `gensys`. Python ports exist in QuantEcon-adjacent code.
- **Alternative:** use `dolo` or `econpizza` to parse a model file. Overkill for a model this small; only consider if you anticipate extending the model class.

Validate the solution by checking steady-state values against BCKM's reported calibration and by checking impulse responses to each wedge.

### 1.5 Wedge extraction

Two steps, in this order:

1. **Static wedges (no model solution needed).**
   - Efficiency wedge: invert the production function. `A_t = y_t / (k_t^α * l_t^(1-α))`.
   - Labor wedge: invert the labor FOC. `1 - τ_{lt} = (ψ * c_t / (1 - l_t)) / ((1-α) * y_t / l_t)`.
   - Government wedge: read directly from data.
2. **Investment wedge (requires the solved decision rules).** The investment Euler equation has an expectation, so you need the estimated VAR (Stage 1.6) to take that expectation. Solve for `τ_{xt}` recursively given the solution and the other three wedges.

For the capital stock series: use the capital accumulation law, data on investment, and an initial `k_0`. Follow BCKM's appendix for the initial condition convention.

### 1.6 VAR(1) estimation by MLE

Estimate the VAR on the latent wedge vector `s_t = (A_t, τ_{lt}, τ_{xt}, g_t)`:

```
s_{t+1} = P_0 + P * s_t + ε_{t+1},   ε ~ N(0, V),   V = Q Q'
```

Estimate `P_0`, `P`, and `Q` (lower triangular) jointly by maximum likelihood using the log-linearized decision rules as the measurement equation and the four observables (y, l, x, g).

Implementation: use `statsmodels.tsa.statespace.MLEModel`. Subclass it and define:

- `transition` matrix (from the VAR and the capital accumulation law)
- `design` matrix (from the log-linearized decision rules)
- `start_params` — start from BCKM's published estimates when replicating; move outward for new samples

Optimizer: `scipy.optimize` via `model.fit()`. Try BFGS first, fall back to Nelder-Mead or Powell if it fails. **Identification can be fragile** (see Brinca–Iskrev–Loria 2022). Run from multiple starting points and report whether they converge to the same mode.

### 1.7 Counterfactual simulations (the "wedge-alone" experiments)

This is the step people get wrong. Read BCKM 2016 Section 2.C carefully.

For each experiment (e.g., "efficiency wedge alone"):

- Keep the full 4-dimensional VAR governing agents' expectations.
- Re-solve the model under the restriction that only the wedge of interest fluctuates; the other three wedges are held at their means as **constant functions of the state** — not with their innovations zeroed.
- Feed the realized historical state sequence `s^d_t` through these restricted decision rules.
- Compute the resulting paths for output, labor, and investment.

**Do not** do the naive thing of zeroing innovations in the VAR before feeding it through the unrestricted decision rules — that is the Christiano–Davis (2006) mistake the BCKM paper explicitly warns against.

### 1.8 Summary statistics

Implement:

- **ϕ-statistic** per BCKM equation on p. 39: for each variable (output, labor, investment), `ϕ_i = (1/Σ(y - y_i)^2) / Σ_j (1/Σ(y - y_j)^2)`. Lies in [0,1], sums to one across wedges.
- **Peak-to-trough decomposition** for any user-specified recession window.

## Stage 2 — `bca_api`: FastAPI service

Endpoints (minimum):

- `POST /runs` — body: `{country, start, end, calibration_overrides?}`. Returns `{run_id, status}`. Kicks off estimation asynchronously.
- `GET /runs/{run_id}` — returns status and, when complete, the full result: wedges, counterfactuals, ϕ-statistics, metadata.
- `GET /data/{country}` — returns the raw and adjusted data series used for that country, for transparency.
- `GET /countries` — list of supported countries.

Persist results in SQLite or Postgres keyed by `(country, sample, calibration_hash)` so repeated runs are cached.

For data ingestion: build a small adapter layer so new data sources can be added without touching `bca_core`. Start with:
- **US:** FRED (`fredapi`) for NIPA aggregates, BLS for hours, Census for population.
- **OECD countries:** OECD SDMX API.

Cache raw data on disk; refresh quarterly.

## Stage 3 — `bca_web`: front end

Minimum viable UI:

1. **Selector panel.** Country dropdown, sample start/end pickers, "Run" button. Advanced section (collapsed by default) for calibration overrides.
2. **Results tab 1 — Wedges.** Four time-series plots, one per wedge, with recession shading.
3. **Results tab 2 — Counterfactuals.** For each of output, labor, investment: the actual series plus the four single-wedge counterfactual paths, overlaid.
4. **Results tab 3 — Summary.** ϕ-statistic table and peak-to-trough decomposition table. Let the user specify the recession window interactively.
5. **Data tab.** Show the adjusted input series so the user can sanity-check.

Recharts or Plotly are both fine. For React, TanStack Query for server state; Zustand or Redux Toolkit only if the app grows beyond this.

## Validation plan (do this before shipping)

Before claiming the port works, validate against Brinca's MATLAB output:

1. Download Brinca's replication files from his site.
2. Run the MATLAB code on US data for 1980:Q1–2014:Q4 (the BCKM sample).
3. Run your Python code on the same data with the same calibration.
4. Check that the following match to ~3–4 decimals:
   - Recovered wedges, period by period
   - VAR parameter estimates (`P_0`, `P`, `V`)
   - Log-likelihood at the converged parameters
   - ϕ-statistics for output, labor, investment
5. Verify you reproduce BCKM's headline finding for the US: labor wedge ≈ 46% of output movement, investment wedge secondary, efficiency wedge small.

Small discrepancies typically come from: optimizer tolerances, trend-removal conventions (level vs. log, linear vs. HP), initial capital stock choice, VAR initial-condition handling. All are fixable.

## Known pitfalls to avoid

- **Do not** identify VAR innovations with wedges in the counterfactuals. Keep the distinction between the underlying state `s_t` and the wedges as functions of that state. This is the single biggest conceptual error in the literature.
- **Do not** skip the durables reclassification. BCKM show it changes results materially for some countries; even when it doesn't, consistency matters.
- **Do not** use a perfect-foresight procedure (as in Ohanian–Raffo 2012) and claim it's BCKM. BCKM Section 3.B.iv shows this inflates the US labor wedge contribution from 46% to 92%. If you want to offer both procedures as options, label them clearly.
- **Do not** assume the MLE has a unique global mode. Run multiple starts and report.
- **Footnote 2 of BCKM:** the US result is sample- and adjustment-cost-sensitive. Extending the sample back to 1948 or lowering `a` raises the labor-wedge role and lowers the investment-wedge role. If you offer sample selection in the UI, the results will move — that's correct behavior, not a bug.

## Suggested milestones

1. Prototype model + steady state + log-linear solution. Unit tests against analytical steady state.
2. Data pipeline for the US. Produce the four observables as a clean `pandas.DataFrame`.
3. Static wedge extraction. Plot against BCKM Figure 2 for the US.
4. VAR MLE. Match BCKM parameter estimates.
5. Investment wedge extraction using the solved VAR.
6. Counterfactuals + ϕ-statistics. Match BCKM Figure 4 and Table 1 for the US.
7. Wrap in FastAPI.
8. Extend data pipeline to OECD countries.
9. Build the front end.
10. Validation against MATLAB replication files.

Ship each milestone behind tests. The model and wedge-extraction code are small enough that >90% test coverage is reasonable.
