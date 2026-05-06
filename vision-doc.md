---
title: "Business Cycle Accounting — Project Vision"
topic: "vision"
layer: "all"
status: "active"
last_updated: "2026-05-06"
---

# Business Cycle Accounting — Project Vision

## Goal

Build a **US-only web application** that performs Business Cycle Accounting (BCA) in the style of
Chari–Kehoe–McGrattan (2007) and Brinca–Chari–Kehoe–McGrattan (2016, NBER WP 22663). The user
picks a US time window, the app runs the full BCA procedure end-to-end, and returns the four
wedges (efficiency, labor, investment, government) plus counterfactual decompositions of output,
labor, and investment.

The flagship use case is a **quarterly-updated Business Cycle Monitor** that combines three things
no standard macro dashboard does together:

1. **GDP decomposition by demand component** — the standard macro picture.
2. **BCA wedge decomposition** — which structural margin is driving the business cycle.
3. **LLM-generated hypothesis layer** — given the wedge configuration, what does the literature
   say could be causing it?

The value proposition is the juxtaposition of all three. Showing that investment collapsed 30%
*and* the investment wedge barely moved tells you something a GDP dashboard alone never would:
the investment decline was not driven by an investment financing distortion, so look elsewhere.
That inference requires the BCA layer. The LLM layer then surfaces the relevant literature for
wherever the wedges actually point.

Cross-country support is explicitly **out of scope**.

## Primary references

- Brinca, Chari, Kehoe, McGrattan (2016), "Accounting for Business Cycles," NBER WP 22663. <https://www.nber.org/papers/w22663>
- Full technical appendix and replication files: <https://pedrobrinca.pt/2016-accounting-for-business-cycles/>
- Brinca, Costa Filho, Loria (2024), "Business cycle accounting: What have we learned so far?" *Journal of Economic Surveys* 38(4) — for the modern synthesis.
- Brinca, Iskrev, Loria (2022), "On Identification Issues in Business Cycle Accounting Models" — read before trusting MLE estimates.

## Architecture

Four layers, kept loosely coupled so each can be tested independently:

1. **`bca_core`** ✅ **Complete** — pure-Python computational core: data pipeline, wedge
   extraction, Kalman-filter VAR MLE, and counterfactuals. No web dependencies.
2. **`bca_api`** — FastAPI service that exposes `bca_core` over HTTP. Long-running MLE jobs run
   asynchronously; results are cached by `(sample_hash, calibration_hash)`. Also hosts the
   quarterly LLM call to the Anthropic API (hypothesis generation), cached until the next
   NIPA release.
3. **`bca_web`** — Single-page front end (React + Recharts). Three screens: Macro Overview,
   Wedge Decomposition, Hypothesis Layer.
4. **Nowcasting extension** *(planned)* — MIDAS regressions that extend the BCA to monthly
   frequency using theoretically motivated high-frequency indicators. Served as an additional
   panel in `bca_web` once implemented.

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
- **CLI entry point** (`scripts/run_bca.py`): any US window via `--start/--end/--base`;
  produces 10 output files (4 figures + 4 CSVs) with NBER recession shading and MLE caching.
  ~3 seconds on a warm cache.

**Validation status:**
- Layer 1 (BCKM US 1980Q1–2014Q4): f-stats match Table 11 to ≤0.01 in every channel. ✅
- Layer 2 (COVID-era 2010Q1–2023Q4): all six qualitative priors pass. ✅

For implementation details, calibration constants, and known pitfalls see `CLAUDE.md` and
`bckm_replication/REPORT.md`.

---

## Stage 2 — `bca_api`: FastAPI service

Endpoints (minimum viable):

- `POST /runs` — body: `{start, end, base, calibration_overrides?}`. Returns `{run_id, status}`.
  Kicks off estimation asynchronously.
- `GET /runs/{run_id}` — returns status and, when complete, the full result: wedges,
  counterfactuals, φ-statistics, metadata.
- `GET /data` — returns the raw and adjusted US data series for transparency.
- `GET /hypotheses/{quarter}` — returns the cached LLM hypothesis output for a given quarter.

Persist results in SQLite or Postgres keyed by `(sample_hash, calibration_hash)` so repeated
runs are served from cache. MLE estimation takes ~2–12 minutes cold; a job queue (RQ or Celery)
is required for anything beyond a toy deployment.

`bca_core.data` already implements a FRED adapter with disk caching. The API layer calls
`build_us_dataset()` and passes the result to `estimate_var_mle()`.

---

## Stage 3 — `bca_web`: Business Cycle Monitor

**Update cadence:** quarterly, triggered by the BEA NIPA release. No intra-quarter nowcasting
(until Stage 4 lands). All panels reflect the most recently completed quarter.

---

### Screen 1 — Macro Overview

**Purpose:** orient the user to current macroeconomic conditions before going into the structural
decomposition.

**Panel 1A — GDP and demand components (time series).**
Real GDP, consumption, investment, government consumption, and net exports indexed to 100 at a
user-selectable base period (default: the most recent NBER peak). Quarterly frequency. Recession
shading from NBER dates. Overlaid line chart with toggles per component; investment and net
exports on a secondary axis if scale differences are large.

Key juxtaposition to surface: investment vs. output. When investment falls faster than output,
the investment wedge is a candidate driver. When labor income falls alongside consumption, the
labor wedge is a candidate. Build this as an annotated callout, not just a raw chart.

**Panel 1B — Current quarter snapshot.**
A single-row summary table: most recent quarter's QoQ and YoY growth rates for GDP and each
demand component, plus their contributions to GDP growth in percentage points (standard BEA
decomposition).

**Panel 1C — Historical context.**
For each series, where the current level sits in its empirical distribution since 1980: percentile
rank and standard deviations from mean.

---

### Screen 2 — Wedge Decomposition

**Purpose:** the BCA layer. Which structural margin is most responsible for the current behavior
of output, labor, and investment.

**Panel 2A — Wedge time series.**
Four subplots, one per wedge (efficiency, labor, investment, government), in standard-deviation
units from their sample means. Recession shading. Annotate the current quarter value and its
historical percentile.

**Panel 2B — φ-statistics: current recession window.**
A stacked bar chart showing, for each of output, labor, and investment, the share of movement
explained by each wedge. The recession window is user-selectable: default to the most recent NBER
recession if ongoing, otherwise let the user select start and end quarters.

This is the central panel of the app.

**Panel 2C — Counterfactual paths.**
For each of output, labor, and investment: the actual data series plus the four single-wedge
counterfactual paths, as in BCKM Figure 2. Indexed to 100 at the start of the selected window.

**Panel 2D — Historical comparison.**
A table showing the φ-statistics for the current episode alongside the same statistics for major
historical US recessions. Lets the user answer: does the current wedge configuration look more
like 2008 or more like 2020?

| Episode | φ_y(eff) | φ_y(labor) | φ_y(inv) | φ_y(gov) |
|---------|----------|------------|----------|----------|
| 1981–82 | — | — | — | — |
| 1990–91 | — | — | — | — |
| 2001    | — | — | — | — |
| 2008–09 | ~0.26 | ~0.44 | ~0.23 | ~0.07 |
| 2020    | — | — | — | — |
| **Current** | **—** | **—** | **—** | **—** |

Values populated by the BCA engine at each quarterly update.

---

### Screen 3 — Hypothesis Layer

**Purpose:** given the wedge configuration, surface candidate micro-foundations from the
literature with an LLM that knows the current state of the world.

**Panel 3A — Wedge summary for LLM context (shown to user).**
A compact display of what is being passed to the LLM: the current φ-statistics, the direction
and percentile of each wedge, and the most analogous historical episode. Shown to the user so
they understand what the LLM is reasoning from.

Example:
> **Current wedge configuration (2024Q4):**
> Labor wedge accounts for 41% of output movement (φ_y = 0.41). Labor wedge is 1.6 SD below its
> mean — its weakest reading since 2010Q2. Most analogous historical episode: 2008Q4–2009Q2.

**Panel 3B — Micro-foundations reference table.**
A structured table mapping each wedge to its candidate micro-foundations from the literature.
This is the intellectual core of the prompt — it grounds the LLM in the actual BCA literature.

**Efficiency Wedge (A_t)**

| Mechanism | Key References | Signature |
|-----------|---------------|-----------|
| TFP / technology shocks | Kydland & Prescott (1982); Prescott (1986) | Persistent A_t decline across all sectors |
| Input-financing frictions | Chari, Kehoe, McGrattan (2007) Prop. 1; Mendoza (2010) | A_t falls when credit to working capital tightens |
| Misallocation across firms/sectors | Hsieh & Klenow (2009); Restuccia & Rogerson (2008) | A_t falls as dispersion in firm-level TFP rises |
| Trade disruptions / import costs | Kehoe & Ruhl (2008) | A_t falls with trade collapse; recovers with trade |
| Entry/exit distortions | Bilbiie, Ghironi & Melitz (2012) | A_t falls as firm entry drops in downturns |
| Energy price shocks | Hamilton (1983); Finn (2000) | A_t falls with oil price spikes |

**Labor Wedge (1 − τ_lt)**

| Mechanism | Key References | Signature |
|-----------|---------------|-----------|
| Search and matching frictions | Mortensen & Pissarides (1994); Shimer (2005) | Wedge falls as job-finding rate falls, vacancies collapse |
| Wage rigidity / sticky wages | Hall (2005); Shimer (2004) | Wedge opens as wages fail to clear the labor market |
| Household financial frictions | Guerrieri & Lorenzoni (2017) | Wedge falls as household deleveraging depresses labor supply |
| Labor income taxes | Prescott (2004); Ohanian, Raffo & Rogerson (2008) | Wedge tracks average marginal tax rate on labor income |
| Markup shocks / imperfect competition | Gali, Gertler & Lopez-Salido (2007) | Wedge opens as firm markups rise in downturns |
| Union bargaining / insider-outsider | Blanchard & Summers (1986) | Persistent wedge after shocks; hysteresis |
| Working capital / payroll financing | Neumeyer & Perri (2005); Christiano, Trabandt & Walentin (2011) | Wedge co-moves with credit spreads |

**Investment Wedge (1/(1 + τ_xt))**

| Mechanism | Key References | Signature |
|-----------|---------------|-----------|
| Financial frictions / collateral constraints | Bernanke, Gertler & Gilchrist (1999); Kiyotaki & Moore (1997) | Wedge co-moves with credit spreads, bank health |
| Investment-specific technical change | Greenwood, Hercowitz & Krusell (1997) | Wedge trends; not cyclical |
| Uncertainty shocks | Bloom (2009); Bloom et al. (2018) | Wedge spikes with VIX; recovers as uncertainty resolves |
| Bank balance sheet constraints | Gertler & Kiyotaki (2010); He & Krishnamurthy (2013) | Wedge worsens when bank equity/capital ratios fall |
| Corporate tax / depreciation policy | McGrattan & Prescott (2005) | Wedge tracks changes in tax treatment of investment |
| Credit crunch / lending standards | Lown & Morgan (2006); Gilchrist & Zakrajsek (2012) | Wedge co-moves with Senior Loan Officer tightening |

**Government Wedge (g_t)**

| Mechanism | Key References | Signature |
|-----------|---------------|-----------|
| Government spending shock | Ramey (2011); Barro (1981) | g_t rises in wars and fiscal expansions |
| Fiscal multiplier / crowding out | Woodford (2011); Christiano, Eichenbaum & Rebelo (2011) | g_t movement has differential effect at ZLB |
| Trade balance shock | Kehoe & Ruhl (2008) | Net exports swing driven by external demand or exchange rate |
| External demand collapse | Eaton et al. (2016) | g_t falls with global trade volume in synchronized recessions |

**Panel 3C — LLM-generated hypotheses.**
The output of the LLM call, displayed as a structured narrative:

1. **Pattern identification.** Which wedge(s) dominate and how the current configuration compares to historical episodes.
2. **Candidate mechanisms.** For each dominant wedge, 2–3 mechanisms from the reference table most consistent with the current configuration, with citations and current-conditions reasoning.
3. **What to watch.** For each candidate mechanism, one observable indicator that would confirm or disconfirm it over the next one to two quarters.

The output is clearly labeled: **"Hypotheses generated by AI — not structural findings. Treat as a starting point for analysis, not a conclusion."**

---

### LLM Prompt Specification

**System prompt**

```
You are an expert macroeconomist specializing in business cycle accounting (BCA)
in the tradition of Chari, Kehoe, and McGrattan (2007) and Brinca, Chari, Kehoe,
and McGrattan (2016).

You will be given:
1. The current BCA estimation results for the US economy: phi-statistics, wedge
   levels in standard deviation units, and a comparison to historical episodes.
2. A reference table of micro-foundations for each wedge from the academic literature.

Your task is to generate structured hypotheses about what economic mechanisms could
be driving the current wedge configuration. You have knowledge of the current state
of the US and global economy. Use this knowledge to assess which mechanisms from
the reference table are most plausible given observable conditions right now.

Format your response in three sections:
- Pattern Identification
- Candidate Mechanisms (for each dominant wedge: 2-3 mechanisms with citations
  and current-conditions reasoning)
- What to Watch (one falsifiable indicator per mechanism)

Be specific. Cite the papers from the reference table by author and year. Do not
cite references not in the table. If a mechanism from the table is implausible
given current conditions, say so and explain why. Acknowledge uncertainty.
```

**User prompt (assembled programmatically each quarter)**

```
## Current BCA Estimation Results — [QUARTER]

phi-statistics (share of output movement explained by each wedge):
  Efficiency wedge:  [phi_y_eff]
  Labor wedge:       [phi_y_lab]
  Investment wedge:  [phi_y_inv]
  Government wedge:  [phi_y_gov]

Current wedge levels (deviation from sample mean, in standard deviations):
  Efficiency:  [val] SD  ([percentile]th percentile since 1980)
  Labor:       [val] SD  ([percentile]th percentile since 1980)
  Investment:  [val] SD  ([percentile]th percentile since 1980)
  Government:  [val] SD  ([percentile]th percentile since 1980)

Direction of change (current quarter vs. prior quarter):
  Efficiency:  [improved / worsened / flat]
  Labor:       [improved / worsened / flat]
  Investment:  [improved / worsened / flat]
  Government:  [improved / worsened / flat]

Most analogous historical episode: [episode]
Similarity metric: [cosine similarity in phi-statistic space]

Peak-to-trough decomposition (recession window: [start] to [end]):
[Table]

## Micro-Foundations Reference Table
[Insert full table]

## Task
Generate structured hypotheses as specified.
```

**Design notes**

- **No world summary is included.** The LLM is assumed to have current knowledge. This keeps the prompt lean and avoids stale summaries of conditions the model already knows.
- **The reference table is always included in full** (~600 tokens). It is small enough that retrieval adds complexity without benefit.
- **The prompt asks for falsifiable indicators.** A hypothesis that cannot be disconfirmed is not useful.
- **One call per quarter, response cached.** Do not re-run the LLM on demand — cache the response until the next quarterly NIPA trigger.
- **Model:** use the most current available Anthropic Sonnet model, updated at each quarterly run so the world-knowledge assumption remains as valid as possible. Display the model name and approximate knowledge cutoff in the UI.

---

### Data Flow

```
NIPA Release (quarterly)
        │
        ▼
bca_api — data pipeline
  ├── Pull FRED: GDP, C, I, G, NX, hours, population, deflator
  ├── Apply BCKM adjustments (durables, sales tax, per-capita, calgz)
  └── Store adjusted series → database
        │
        ▼
bca_core — BCA engine
  ├── Extract efficiency wedge (production function inversion)
  ├── Extract labor wedge (labor FOC inversion)
  ├── Estimate VAR(1) by Kalman-filter MLE
  ├── Extract investment wedge (Euler recursion + estimated VAR)
  ├── Run 4 single-wedge counterfactual simulations
  ├── Compute φ-statistics and peak-to-trough decompositions
  ├── Compute historical episode comparison (cosine similarity in φ-space)
  └── Store results → database keyed by (country, quarter, calibration_hash)
        │
        ▼
bca_api — LLM call (Anthropic API)
  ├── Assemble prompt: φ-stats + wedge levels + reference table
  ├── Call current Anthropic Sonnet model
  ├── Parse and store structured response
  └── Cache until next NIPA release
        │
        ▼
bca_api (FastAPI) → serves all panels to bca_web
        │
        ▼
bca_web (React)
  ├── Screen 1: Macro Overview
  ├── Screen 2: Wedge Decomposition
  └── Screen 3: Hypothesis Layer
```

---

### Front-End Component Map

| Screen | Component | Data source |
|--------|-----------|-------------|
| 1A | `GDPChart` | FRED via bca_api |
| 1B | `QuarterlySnapshot` | FRED via bca_api |
| 1C | `HistoricalPercentiles` | bca_api computed |
| 2A | `WedgeTimeSeries` | bca_core output |
| 2B | `PhiStatBar` | bca_core output |
| 2C | `CounterfactualChart` | bca_core output |
| 2D | `HistoricalEpisodesTable` | bca_core output, pre-computed per episode |
| 3A | `WedgeSummaryBox` | bca_core output, formatted |
| 3B | `MicroFoundationsTable` | Static — hardcoded in front end |
| 3C | `HypothesisPanel` | LLM output, cached quarterly |

The `MicroFoundationsTable` is static and hardcoded in the front end. It does not change with new
data and should not be fetched from the API. Adding a new paper to the reference table is a code
change, not a data change.

---

## Stage 4 — Nowcasting: Monthly BCA Wedge Monitor *(planned)*

### Motivation

The BCA procedure is inherently quarterly — it requires NIPA aggregates released with a lag.
This creates a monitoring gap: if you want to know in real time whether a developing downturn
looks like a labor-wedge recession (US 2008) or an efficiency-wedge recession, you have to wait
one to two quarters.

The nowcasting extension exploits the structural interpretation of each wedge to motivate
theoretically grounded monthly predictors rather than black-box indicators. The goal is a working
monitoring tool, not a publishable econometric contribution.

---

### Wedge-by-Wedge Strategy

#### Efficiency Wedge (A_t)

**Structural basis.** A_t is a static inversion:
`A_t = y_t / (k_t^α · l_t^(1−α))`. If you can nowcast quarterly output and hours, A_t drops out
immediately. The capital stock moves slowly and can be held at its last observed quarterly value
without material error.

**Approach:** indirect MIDAS. Nowcast y_t and l_t at monthly frequency, then invert the
production function.

**Monthly indicators:**
- INDPRO (industrial production) — primary; tracks the production side at monthly frequency
- ISM Manufacturing PMI (composite + new orders sub-index) — forward-looking
- Chicago Fed CFNAI — pre-built factor model across 85 monthly series; usable as a single predictor

**Specification:** `A_{q(t)} = α + β(L^{1/3}) · IP_t + ε_t` with Almon or beta MIDAS weights.

---

#### Labor Wedge (τ_lt)

**Structural basis.** Also a static inversion from the labor FOC:
`1 − τ_lt = [ψ · c_t / (1 − l_t)] / [(1−α) · y_t / l_t]`. Given monthly nowcasts of output,
hours, and consumption, it drops out immediately.

**Approach:** indirect MIDAS on hours and claims.

**Monthly indicators:**
- BLS Employment Situation: nonfarm payrolls, average weekly hours, average hourly earnings
- JOLTS: job openings, hires, quits rates — forward-looking labor market tightness
- Initial and continued jobless claims (weekly, aggregated to monthly) — highest frequency signal
- Real retail sales (RSAFS, deflated) — proxy for consumption within-quarter

**Specification:** `τ_{l,q(t)} = α + β₁(L^{1/3}) · hours_t + β₂(L^{1/3}) · claims_t + ε_t`

---

#### Investment Wedge (τ_xt)

**Structural basis.** Unlike the other wedges, τ_xt is *not* a static inversion — it requires
the estimated VAR to handle expectations. It maps, via BCKM Propositions 1 and 2, to financial
frictions and collateral constraints.

**Approach:** direct MIDAS on the Gilchrist–Zakrajšek excess bond premium (GZ EBP). The GZ EBP
strips out the firm-level default-risk component of corporate spreads, leaving the financial
intermediary health component — the closest available market price to the shadow price of the
collateral constraint in the BCKM bank model.

**Monthly indicators:**
- GZ excess bond premium (Gilchrist & Zakrajšek 2012; available monthly from the Fed)
- BAA–AAA Moody's spread (FRED: BAA, AAA) — simpler baseline
- Chicago Fed NFCI — composite financial conditions (credit + leverage + risk sub-indices)
- VIX (FRED: VIXCLS) — options-implied volatility; proxy for uncertainty shocks

**Specification:** `τ_{x,q(t)} = α + β(L^{1/3}) · EBP_t + ε_t`

**Honest caveat.** The GZ-based nowcast is a partial nowcast of the financial-friction component
only. During episodes driven by investment-specific technical change, the spread will have less
predictive power. Estimate the model separately for NBER recession and non-recession subsamples
and report both R² values — if recession R² is high and non-recession R² is low, that is the
correct finding, not a failure.

---

#### Government Wedge (g_t)

**Structural basis.** In the closed-economy BCKM framework, `g_t = government consumption + net
exports`. This is directly observable, not a model residual.

**Approach:** direct data. Monthly Treasury Statement for federal outlays; BEA/Census Trade
release for net exports. Reasonable simplification: hold g_t at its last observed quarterly value
within-quarter and update only at the quarterly NIPA release. Flag this assumption in the output.

---

### Implementation Plan

1. **Quarterly baseline** — run BCA on the full quarterly sample; store extracted wedge series as MIDAS dependent variables.
2. **Efficiency and labor wedge nowcasts** — estimate MIDAS regressions; evaluate in pseudo-OOS from 2000 onward with a rolling window.
3. **Investment wedge nowcast** — estimate MIDAS on GZ EBP; run recession vs. non-recession subgroup check.
4. **Government wedge** — wire up MTS and trade data; hold constant within-quarter.
5. **Integration with `bca_api`** — expose nowcasted wedges alongside the quarterly BCA output; add "current quarter estimate" panel to `bca_web` with confidence bands and the date of the most recent monthly update.

**Stack:** `fredapi` for FRED data, `scipy.optimize` for NLS MIDAS estimation (~50 lines, more transparent than a library), `statsmodels.tsa.dynamic_factor_mq` for DFM extension if pursued.

**Evaluation per wedge:**
- Pseudo-OOS RMSE relative to AR(1) baseline
- Directional accuracy (% of quarters with correct sign of change)
- Recession vs. non-recession RMSE separately
- Within-quarter vintage analysis (how the nowcast evolves as month 1, 2, 3 arrive)

---

### Honest Caveats

- **Residuals of residuals.** The wedges already absorb everything the model cannot explain. Nowcasting them adds a second layer of approximation. Communicate uncertainty bands.
- **Data revisions.** BCA wedges are sensitive to precise values of output, hours, and consumption. First-vintage data is revised substantially. Consider a true real-time evaluation using the Philadelphia Fed Real-Time Data Set.
- **Investment wedge is partial by construction.** The GZ spread captures only the financial-friction component of τ_xt. This is correct behavior — explain it that way.
- **Use for directional monitoring, not point inference.** The primary value is identifying the *type* of recession in real time, not pinning exact wedge magnitudes.

---

## Known pitfalls

- **Do not** identify VAR innovations with wedges in the counterfactuals. Keep the distinction
  between the underlying state `s_t` and the wedge-alone counterfactual paths.
- **Do not** skip the durables reclassification. It changes results materially.
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
| 7 | `bca_api` — FastAPI service + async job queue | ⬜ Next |
| 8 | `bca_web` Screen 1 — Macro Overview (GDP chart, snapshot, percentiles) | ⬜ |
| 9 | `bca_web` Screen 2 — Wedge Decomposition (φ-stats, CFs, historical table) | ⬜ |
| 10 | `bca_web` Screen 3 — Hypothesis Layer (LLM integration, micro-foundations table) | ⬜ |
| 11 | Quarterly trigger — automated NIPA-release pipeline + LLM cache invalidation | ⬜ |
| 12 | Nowcasting extension — MIDAS regressions, monthly wedge monitor panel | ⬜ |
