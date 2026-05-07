---
title: "BCA — AI Context Map"
topic: "project-overview"
layer: "all"
status: "active"
last_updated: "2026-05-04"
---

# BCA — AI Context Map

This file is a machine-readable context map for AI agents and LLMs. It defines
every acronym, walks through the computation pipeline step by step, describes all
key data structures, and documents the conventions that are most likely to cause
silent errors if misunderstood.

Start here before reading any other file in the repository.

---

## 1. What this project does (no acronyms)

This project implements **Business Cycle Accounting** — a macroeconomic technique
that decomposes aggregate fluctuations (drops in output, hours worked, and investment)
into four structural distortions called "wedges." Each wedge corresponds to a type
of friction: how efficiently the economy produces (efficiency wedge), how much the
labor market is distorted (labor wedge), how expensive it is to invest (investment
wedge), and what the government takes out of private resources (government wedge).

The pipeline takes quarterly US macroeconomic time series, fits a statistical model
by maximum likelihood, and then runs counterfactual simulations to answer: "If only
the labor wedge had moved — and everything else stayed at its long-run average —
how much of the output decline would we explain?" The answer for the 2008–2009
recession is roughly 46% from the labor wedge.

The project is a Python port of the published Matlab replication code by Brinca,
Chari, Kehoe, and McGrattan (2016). It is validated against the paper's published
results (Layer 1) and against a COVID-era out-of-sample window (Layer 2).

A high-fidelity web application (**BCA Monitor**) is the primary interface for this project. 
It uses a static React frontend to visualize structural findings and an AI Hypothesis Layer 
(Gemini) to interpret them for human analysts. An automated quarterly pipeline (GitHub Actions) 
ensures the dashboard is always up to date with the latest BEA/FRED releases.

---

## 2. Glossary

Every acronym, shorthand variable name, and "inside" term used in the codebase,
fully expanded.

### 2.1 Author and paper abbreviations

| Symbol | Full expansion | Where defined |
|--------|---------------|---------------|
| **BCA** | Business Cycle Accounting | methodology name |
| **BCKM** | Brinca, Chari, Kehoe, McGrattan (2016 NBER paper) | `BCA_info.md` |
| **CKM** | Chari, Kehoe, McGrattan (2007 *Econometrica* paper) | `BCA_info.md` |
| **BGG** | Bernanke, Gertler, Gilchrist (1999) — source of the adjustment-cost elasticity formula | `bca_core/params.py` |
| **BIL** | Brinca, Iskrev, Loria (2022) — identification paper | `vision-doc.md` |

### 2.2 Algorithmic and statistical abbreviations

| Symbol | Full expansion | Where used |
|--------|---------------|-----------|
| **MLE** | Maximum Likelihood Estimation | `bca_core/var_estimation.py` |
| **VAR** | Vector Autoregression — the joint process governing the four wedges | `bca_core/var_estimation.py` |
| **KF** | Kalman Filter — computes the likelihood of the observables given parameters | `bca_core/var_estimation.py` |
| **RTS** | Rauch-Tung-Striebel smoother — backward pass that extracts smoothed wedge paths | `bca_core/var_estimation.py` |
| **DARE** | Discrete Algebraic Riccati Equation — solved once per parameter evaluation to get the steady-state Kalman gain | `bca_core/var_estimation.py` |
| **QZ** | Generalized Schur decomposition — used by Klein (2000) to solve the linear rational expectations model | `bca_core/klein.py` |
| **SS** | Steady State — the deterministic long-run equilibrium around which the model is log-linearized | `bca_core/model.py` |
| **LOM** | Law of Motion — the equation describing how the capital stock evolves each period | `bca_core/bckm_lom.py` |
| **CF** | Counterfactual — a simulated path where only one wedge moves and the rest stay at their long-run mean | `bca_core/counterfactuals.py` |
| **LL** | Log-Likelihood — the objective maximized by the MLE | `bca_core/var_estimation.py` |
| **PIM** | Perpetual Inventory Method — used to construct the capital stock from investment data | `bca_core/wedges.py` |
| **FOC** | First-Order Condition — the optimality condition used to extract labor and investment wedges | `bca_core/wedges.py` |
| **OLS** | Ordinary Least Squares — used as a warm-start for the MLE (not the final estimator) | `bca_core/var_estimation.py` |

### 2.3 Model and estimation parameters

| Symbol | Full expansion | Typical value |
|--------|---------------|---------------|
| **α (alpha)** | Capital share in Cobb-Douglas production | 1/3 |
| **ψ (psi)** | Leisure weight in the household utility function | 2.5 |
| **δ (delta)** | Annual capital depreciation rate | 0.05 |
| **β (beta)** | Annual household discount factor (`1 / (1 + ρ)`) | 0.975 |
| **ρ (rho)** | Annual rate of time preference | 0.025641 |
| **γ (gamma)** | Annual technology growth rate (data-derived via `calgz`) | ~0.019 for US 1980–2014 |
| **n** | Annual population growth rate (data-derived from working-age population slope) | ~0.0098 for US 1980–2014 |
| **a (adj_cost)** | Adjustment-cost curvature coefficient (`adj_cost_elasticity / b`) | ~12.88 |
| **b** | Steady-state investment-to-capital ratio: `(1+γ)(1+n) − (1−δ)` | ~0.0199 |
| **g_share** | Government consumption as a share of output (data mean `g/y`) | ~0.115 for US 1980–2014 |

### 2.4 Matrix and vector names

| Symbol | Expansion | Shape | Where |
|--------|-----------|-------|-------|
| **P** | VAR(1) transition matrix for the four wedges | 4×4 | `bca_core/constants.py`, `estimate_var_mle` |
| **P_0** | VAR intercept vector; equals `(I − P) · Sbar` | (4,) | `estimate_var_mle` return dict |
| **Q** | Lower-triangular Cholesky factor of the VAR shock covariance `V = Q·Qᵀ` | 4×4 | `bca_core/constants.py`, `estimate_var_mle` |
| **V** | VAR shock covariance matrix (`V = Q·Qᵀ`) | 4×4 | `bca_core/var_estimation.py` |
| **Sbar** | Unconditional mean of the VAR state in BCKM coordinates `[log z̄, τ̄_l, τ̄_x, log ḡ]` | (4,) | `bca_core/constants.py` |
| **F** | State transition matrix in state-space form (5×5) | 5×5 | `bca_core/bckm_lom.py` |
| **H** | Observation (design) matrix mapping state to observables | 4×5 | `bca_core/bckm_lom.py` |
| **Gamma (Γ)** | Capital law-of-motion coefficients `[γ_k, γ_z, γ_τl, γ_τx, γ_g]` | (5,) | `bca_core/bckm_lom.py` |
| **theta (θ)** | Full 30-parameter MLE vector `[Sbar(4), P.flatten(16), Q_lower_tri(10)]` | (30,) | `bca_core/var_estimation.py` |

### 2.5 Variable naming conventions

| Convention | Meaning | Example |
|-----------|---------|---------|
| `x_hat` | Log-deviation of variable x from its steady-state value: `log(x_t / x_ss)` | `k_hat`, `A_hat`, `g_hat` |
| `taul` | Labor wedge distortion τ_l (in levels, not logs, because it can be negative) | `bca_core/bckm_lom.py` |
| `taux` | Investment wedge distortion τ_x (in levels) | `bca_core/bckm_lom.py` |
| `log_z` | Log efficiency level `log(A_t)` in BCKM coordinates | `Sbar[0]` |
| `log_g` | Log government spending `log(g_t)` in BCKM coordinates | `Sbar[3]` |
| `obs_hat` | The observable matrix `(T × 4)` in log-deviation form `[y_hat, l_hat, x_hat, g_hat]` | `prepare_observables` output |
| `obs_offset` | Steady-state log-level offsets for each observable (Sbar-dependent) | `_build_ss` output |
| `ss_new` | Steady-state dict recomputed from current Sbar at each optimizer step | `_build_ss` output |
| `states` | RTS-smoothed latent state matrix `(T × 5)` `[k_hat, A_hat, taul_hat, taux_hat, g_hat]` | `estimate_var_mle` return |
| `bind` | "Base-period index" — the quarter used as the Y0 anchor for level-ratio statistics (2008Q1 in the BCKM window) | `counterfactuals.py` |

### 2.6 Data pipeline abbreviations

| Symbol | Expansion |
|--------|-----------|
| **FRED** | Federal Reserve Economic Data — the data source API used in production |
| **BEA** | Bureau of Economic Analysis — NIPA tables; used in diagnostic-only BEA branches |
| **BLS** | Bureau of Labor Statistics — source for hours and employment data |
| **NIPA** | National Income and Product Accounts — BEA's main accounting framework |
| **calgz** | "Calibrate growth rates" — BCKM's fsolve procedure that fits the exponential trend growth rate `gz` so that `log(y_pc)` has mean zero on the MLE window; mirrors `calgz.m` |
| **PCE** | Personal Consumption Expenditure |
| **GPDI** | Gross Private Domestic Investment |
| **rSTX** | Real sales-tax revenue (subtracted from output and investment per BCKM `usdata.m`) |
| **rKCD** | Real consumer-durable-goods capital stock (perpetual inventory) |
| **CE16OV** | FRED ticker for BLS civilian employment (CPS, ages 16+) |
| **AWHAETP** | FRED ticker for average weekly hours of all private-sector employees (BLS) |

### 2.7 File suffixes

| Suffix | Meaning |
|--------|---------|
| `.parquet` | Cached dataset produced by `build_us_dataset`; committed for the BCKM window |
| `.mle.pkl` | Content-addressed pickle of the full `estimate_var_mle` result dict; gitignored (large, regeneratable in ~12 min) |
| `.meta.json` | Metadata sidecar for the parquet — records `start`, `end`, `detrend_method` |

### 2.8 Validation layer terminology

| Term | Meaning |
|------|---------|
| **Layer 1** | BCKM 2016 US 1980Q1–2014Q4 regression test — the pipeline must reproduce Table 11 f-stats to ≤0.01 |
| **Layer 2** | Generalizability validation on a different US window (COVID 2010Q1–2023Q4) |
| **Layer 3** | Cross-country extension — explicitly out of scope |
| **BCKM-θ** | The published BCKM parameter triple `(Sbar, P, Q)` from Tables 8/9/10, used as a probe/warm-start |

---

## 3. Technical stack

| Component | Technology |
|-----------|-----------|
| Language | Python 3.10+ |
| Numerical | numpy ≥ 1.24, scipy ≥ 1.10 |
| Data | pandas ≥ 2.0, pyarrow ≥ 14 (parquet) |
| State-space model | statsmodels ≥ 0.14 (`MLEModel` subclass for Kalman machinery) |
| Optimizer | `scipy.optimize.minimize` with L-BFGS-B method |
| QZ decomposition | `scipy.linalg.qz` (Klein 2000 solver) |
| DARE | `scipy.linalg.solve_discrete_are` |
| Data API | `fredapi` ≥ 0.5 (FRED) |
| **Web UI** | **React 18, Vite, Recharts, Lucide Icons, Vanilla CSS** |
| **Hypothesis Layer**| **Gemini 2.5 Flash** (via `google-genai` Python SDK) |
| **Hosting** | **Vercel** (Automatic rebuilds on push) |
| **Automation** | **GitHub Actions** (Triggered 1st of every month or manually) |
| Caching | Content-addressed pickle keyed by SHA-256 of inputs; 90-day FRED disk cache under `~/.bca_cache/fred/` |

---

## 4. Core computation pipeline

The full pipeline from raw FRED data to φ-statistics in five numbered steps.

### Step 1 — Data construction (`bca_core/data/`)

**Entry point**: `build_us_dataset(start, end, ...)`

```
FRED API  →  raw quarterly NIPA series [GDP, PCE, GPDI, GCE, NETEXP, GDPDEF, ...]
          +  BLS labor [CE16OV, AWHAETP or fallback PAYEMS×AWHNONAG]
          +  population [LFWA64TTUSQ647N]
          ↓
adjustments.py:
  1. reclassify_durables()     — strip PCE durables from C, add to X; build
                                  durables capital stock by perpetual inventory;
                                  add 4% service flow to both C and Y
  2. subtract_sales_tax()      — subtract rSTX from Y and C
  3. to_real_per_capita()      — divide by GDPDEF deflator, then by working-age pop
  4. compute_labor_input()     — employment × avg_weekly_hours × 13/quarter,
                                  normalized to BCKM anchor 0.24279
  5. compute_government_wedge()— government consumption + net exports
  6. remove_trend() / calgz    — fit exponential trend gz by fsolve so log(y_pc)
                                  has mean zero on the MLE window; detrend all series
          ↓
Output: DataFrame with columns [y, l, x, g] — quarterly log-level series (T rows)
```

### Step 2 — Observable construction (`estimate_var_mle` → `prepare_observables`)

```
DataFrame[y, l, x, g]  →  obs_hat (T × 4) matrix
                            columns: [y_hat, l_hat, x_hat, g_hat]
                            each = log(series_t)  (raw log levels, not re-centered)
                            centering is done implicitly by obs_offset = log(ss_new)
                            inside the Kalman filter
```

### Step 3 — VAR MLE (`estimate_var_mle`)

```
obs_hat  →  [initialization]
              _fsolve_sbar_initmle():
                fsolve for Sbar such that model SS means match data means
                (mirrors BCKM initmle.m — gives optimizer a valid starting basin)
           →  [L-BFGS-B optimization over theta = (Sbar, P, Q_chol), 30 params]
              At each function evaluation:
                _build_ss(Sbar, P, Q):
                  1. Solve model SS from Sbar via nonlinear system
                  2. bckm_state_space() → (F, H) matrices  [mleqadj.m:167-232]
                  3. obs_offset = log(ss_new[var])          [Sbar-dependent intercept]
                _steady_state_kalman(F, H, Q_proc):
                  DARE → constant gain K, innovation covariance S
                _kf_ll(F, H, Q_proc, obs_offset):
                  Forward Kalman pass → log-likelihood
                penalty if spectral_radius(P) > 0.995
           →  [At convergence]
              _rts(x_filt, P_filt, ...):
                Rauch-Tung-Striebel backward smoother → smoothed states (T × 5)
           →  Result dict: {Sbar, P, Q, H, F, ss_new, states, log_likelihood, obs_offset}
           →  Cached to *.mle.pkl (SHA-256 content key)
```

### Step 4 — Counterfactuals (`bca_core/counterfactuals.py`)

```
states (T × 5) + (P, Sbar, ss_new, H)  →  run_all_counterfactuals()

For each wedge j in {efficiency, labor, investment, government}:
  solve_counterfactual(model, P_var, active_wedges=[j], ss=ss_new):
    bckm_state_space_cf(..., As=[0,0,0,0]):        → C_0 (no-wedge baseline policies)
    bckm_state_space_cf(..., As=[active]):          → C_j (single-wedge policies)
    Incremental policy: (C_j - C_0)                [BCKM gwedges2.m:90-115]
  run_counterfactual(states, cf_policy):
    Simulate k' = P_k @ state, y = P_y @ state
    using the incremental decision rules

Output: dict {wedge_name → {"y": (T,), "l": (T,), "x": (T,)}}
  — four single-wedge CF paths that sum to ≈ all-active ≈ data
```

### Step 5 — φ-statistics (`phi_statistics`) and f-statistics (`f_statistics_bckm`)

```
data_hat + counterfactual paths  →  phi_statistics()
  For each observable v in {y, l, x}:
    SSR_j = Σ_t (data_v_t - cf_j_v_t)^2   for each wedge j
    φ_j(v) = (1/SSR_j) / Σ_k (1/SSR_k)   [level-space SSR decomposition]
  Returns DataFrame (rows=wedges, cols=observables), each column sums to 1.0

                                         →  f_statistics_bckm()  [Table 11 formula]
  Level-ratio SSR over the Great Recession window (2008Q1–2011Q4):
    fY_j = (1/SSR_j_level) / Σ_k (1/SSR_k_level)

### Step 6 — Dashboard Generation (`scripts/build_quarterly_data.py`)
```
stats_payload  →  Gemini 2.5 Flash (LLM)  →  Qualitative structural hypotheses
               →  Latest FRED headline data  →  Real-time growth optics
               →  JSON Export (latest_quarter.json)
```

### Step 7 — Frontend Visualization (`bca_web/`)
```
latest_quarter.json  →  React State  →  Recharts Line/ComposedCharts
                                    →  Glassmorphism Metric Cards
                                    →  AI Hypothesis Layer Component
```
```

---

## 5. Key data structures

### 5.1 State vector (5-element, log-deviation form)

```
[k_hat, A_hat, taul_hat, taux_hat, g_hat]
  [0]     [1]     [2]       [3]      [4]

- k_hat  : log-deviation of capital stock from steady state
- A_hat  : log-deviation of total factor productivity (efficiency wedge)
- taul_hat: deviation of labor wedge distortion τ_l from its SS value
- taux_hat: deviation of investment wedge distortion τ_x from its SS value
- g_hat  : log-deviation of government consumption from steady state

Only wedges 1-4 receive VAR(1) shocks. Capital (index 0) evolves
endogenously through the capital law-of-motion.
```

### 5.2 Observable matrix (`obs_hat`, shape T×4)

```
Columns: [y_hat, l_hat, x_hat, g_hat]
  y_hat = log(y_t)      (log output, detrended quarterly per-capita real GDP)
  l_hat = log(l_t)      (log labor input, hours-per-working-age-pop, normalized)
  x_hat = log(x_t)      (log investment, detrended quarterly per-capita real)
  g_hat = log(g_t)      (log government wedge = gov consumption + net exports)
```

### 5.3 Parameter vector theta (30 elements)

```
theta = [Sbar(4), P.flatten()(16), Q_lower_tri(10)]

Sbar: unconditional VAR mean in BCKM coordinates [log z̄, τ̄_l, τ̄_x, log ḡ]
      Bounds: [-1,-1,-1,-5] to [1,1,1,1]  (mleqadj.m lines 109/119)

P:    4×4 VAR transition in CODE convention: state_{t+1} = P · state_t
      Row i = equation for wedge i; column j = coefficient on wedge j_{t}
      WARNING: BCKM Table 8 in the paper prints P TRANSPOSED.
               Always import from bca_core/constants.py, never re-transcribe.

Q:    4×4 lower-triangular Cholesky factor of shock covariance V = Q·Qᵀ
      10 lower-triangular elements packed row-by-row: Q[0,0], Q[1,0], Q[1,1], ...
```

### 5.4 `estimate_var_mle` return dict

```python
{
    "Sbar":          np.ndarray (4,)      # VAR unconditional mean
    "P":             np.ndarray (4, 4)    # VAR transition matrix (code convention)
    "Q":             np.ndarray (4, 4)    # Cholesky factor of shock covariance
    "H":             np.ndarray (4, 5)    # Observation matrix (policy rows [y, l, x, g])
    "F":             np.ndarray (5, 5)    # State transition matrix
    "ss_new":        dict                 # Steady-state dict recomputed from Sbar
    "states":        np.ndarray (T, 5)   # RTS-smoothed latent states
    "log_likelihood":float               # Log-likelihood at converged θ (our convention)
    "obs_offset":    np.ndarray (4,)     # SS log-level offsets for observables
    "P_0":           np.ndarray (4,)     # VAR intercept = (I-P) · Sbar
}
```

### 5.5 Counterfactual result dict

```python
{
    "efficiency":  {"y": np.ndarray(T,), "l": np.ndarray(T,), "x": np.ndarray(T,)},
    "labor":       {"y": ..., "l": ..., "x": ...},
    "investment":  {"y": ..., "l": ..., "x": ...},
    "government":  {"y": ..., "l": ..., "x": ...},
    "all_active":  {"y": ..., "l": ..., "x": ...},  # sum ≈ data
}
```

---

## 6. File-by-file roles

| File | Role | Key public API |
|------|------|----------------|
| `scripts/run_bca.py` | **Primary CLI** — any US window via `--start/--end/--base`; 10 output files, FRED key guide, ~3s on cache hit | `--help` |
| `bca_core/params.py` | Calibration constants with BCKM Table 1 defaults | `CalibrationParams` dataclass |
| `bca_core/constants.py` | Published BCKM parameter tables **in code convention** | `P_BCKM_TABLE8`, `SBAR_BCKM_TABLE8`, `QCHOL_BCKM_TABLE10` |
| `bca_core/model.py` | Steady-state computation, log-linearization, Klein QZ solution | `PrototypeModel`, `ModelSolution` |
| `bca_core/klein.py` | Generalized Schur QZ decomposition for linear rational expectations | `klein_solve`, `KleinSolution`, `BlancharKahnError` |
| `bca_core/bckm_lom.py` | BCKM-faithful capital law-of-motion and state-space matrices | `bckm_state_space`, `bckm_state_space_cf` |
| `bca_core/var_estimation.py` | Kalman-filter MLE, RTS smoother, content-addressed caching | `estimate_var_mle`, `prepare_observables`, `BCAStateSpace` |
| `bca_core/counterfactuals.py` | Incremental wedge-alone simulations + φ/f-statistics | `run_all_counterfactuals`, `phi_statistics`, `f_statistics_bckm` |
| `bca_core/wedges.py` | Static wedge extraction from time series | `extract_wedges_bckm_style`, `build_capital_stock` |
| `bca_core/bckm_reference.py` | Loader for BCKM `worktemp.mat` reference artifacts | `load_bckm_reference` |
| `bca_core/data/pipeline.py` | End-to-end data entry point | `build_us_dataset` |
| `bca_core/data/adjustments.py` | BCKM data adjustments (durables, sales tax, labor, calgz) | `reclassify_durables`, `compute_labor_input`, `remove_trend` |
| `bca_core/data/fred.py` | FRED series fetcher with 90-day disk cache | `FredDataFetcher` |
| `bca_core/data/bea.py` | BEA NIPA + Fixed Assets fetcher (diagnostic opt-in) | `BeaDataFetcher` |
| `scripts/build_quarterly_data.py` | **Automated Data Exporter** — Orchestrates Layer 2 + LLM + FRED for dashboard usage | `generate_hypotheses`, `main` |
| `bca_web/src/App.tsx` | **Main Dashboard Controller** — Responsible for data fetching, chart rendering, and AI logic | `App` Component |
| `.github/workflows/update_data.yml` | **CI/CD Workflow** — Schedules monthly runs of `build_quarterly_data.py` | — |

---

## 7. Validation status

| Layer | Window | Status | Gate |
|-------|--------|--------|------|
| Layer 1 | US 1980Q1–2014Q4 (BCKM replication) | ✅ Closed | f-statistics match BCKM Table 11 to ≤0.01 in every channel |
| Layer 2 | US 2010Q1–2023Q4 (COVID era) | ✅ Closed | 6/6 qualitative priors pass under both trend variants |
| Layer 3 | Cross-country | 🚫 Out of scope | — |

---

## 8. Critical conventions for AI agents

These are the conventions most likely to cause silent, hard-to-diagnose errors.

### 8.1 The P-matrix transpose trap

**BCKM Table 8 in the published paper prints P transposed** relative to the standard
textbook VAR convention and relative to what BCKM's own Matlab code uses.

- **Code convention** (used everywhere in this codebase): `state_{t+1} = P · state_t`
  — `P[i, j]` is the coefficient of `state_j_t` in the equation for `state_i_{t+1}`.
  Row `i` = "what determines wedge `i`'s update."
- **Paper convention** (Table 8 only): Row `i` = "what wedge `i` does to others" =
  the transpose of the code convention.

**Rule**: Always import P from `bca_core/constants.py`. Never re-transcribe Table 8.
If you need the paper orientation, take the transpose explicitly.

Consequence of getting this wrong: 501 nats of log-likelihood loss at published θ.

### 8.2 Counterfactuals are incremental, not single-wedge absolute

A "single-wedge counterfactual" simulates the **incremental contribution** of
activating wedge `j` relative to the **no-wedge baseline** (where all wedges are
held at their unconditional mean):

```
CF_j path = (C_j − C_0) · state_t    [BCKM gwedges2.m:90-115]
```

**Not**: `C_j · state_t` (which includes the no-wedge baseline and makes the four
CFs over-count the data, summing to ~−10.9% vs data −7.5% in the Great Recession).

The four incremental single-wedge CFs sum to approximately the all-active CF ≈ data.

### 8.3 The Y0 anchor is the base period, not the sample start

Level-ratio statistics (BCKM Table 11/12) are anchored at `bind` — the base period
index, which is **2008Q1** in the BCKM window. Not at 1980Q1 (sample start).

Using the sample start as anchor collapses the 1980–2008 cumulative drift into the
sum-of-squared-residuals and produces wrong f-statistics.

### 8.4 Sbar parameterizes the steady state, not just the VAR mean

`Sbar` is not just the unconditional mean of the wedge VAR. It parameterizes the
*physical steady state* of the economy through the per-call `_build_ss` → `bckm_state_space`
chain. Changing Sbar re-linearizes the model (different `H` matrix), which changes
both the policy functions and the observable-intercept `obs_offset`.

A fixed `obs_offset` while Sbar drifts causes the Kalman intercept to be
systematically wrong — this was the "phi0 bug" that kept the optimizer at a
data-independent attractor (`log_g ≈ −1.2` vs BCKM's `−1.94`).

### 8.5 Log-likelihood conventions

Two different LL conventions are used in the literature:

- **Our convention** (printed by `estimate_var_mle`): standard log-density, includes
  `−0.5 · T · n · log(2π)` ≈ +514 nats for T=140, n=4.
- **BCKM convention** (stored in `worktemp.mat` as `mle.Likelihood`): minimization
  objective `L = 0.5·(T·log|Ω| + tr(Ω⁻¹Σ_innov))`, no `log(2π)` term.
  `mle.Likelihood` ≈ −2403 (negative, because `T·log|Ω| ≈ −5398` dominates).

When comparing log-likelihoods between our output and BCKM's Matlab output,
always account for the 514-nat offset and the sign.

### 8.6 Labor input normalization

Raw FRED labor hours per population is ~23 (hours/week scale), while the model
steady-state `l_ss ≈ 0.29` (fraction-of-time scale) — an 80× scale mismatch.

The fix: `compute_labor_input` normalizes so the **sample mean matches 0.24279**
(the BCKM-empirical anchor from their `usdata.m`). The resulting series is fed
as raw log-levels (`log(l_t)`) to the Kalman filter; the steady-state offset
`obs_offset[1] = log(ss_new["l"])` handles the level on the prediction side.

Do **not** re-scale to `l_ss` before or after this normalization.
