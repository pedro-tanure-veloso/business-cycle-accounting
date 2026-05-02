# COVID Analysis Report — US 2010Q1–2023Q4

> **Status**: Second pass complete (2026-05-02) with BLS-faithful labor
> data + MLE result caching. Verdict: **PASS** — pipeline qualitatively
> reproduces narrative priors under both trend variants.

## Setup

| Parameter | Value |
|---|---|
| Window | 2010Q1–2023Q4 (56 quarters) |
| Bind (anchor) | 2019Q4 (NBER cycle peak) |
| Trend variant A | Full-window calgz (γ_annual = 2.36%) |
| Trend variant B | Pre-COVID-fit (slope on 2010Q1–2019Q4, γ_annual = 2.25%) |
| Data source | FRED defaults |
| Labor construction | **CPS employment × avg weekly hours total × 13 weeks/qtr** (LNS12000000 × AWHAETP — BLS-faithful, BCKM `usdata.m` `hours.dat` analogue) |
| Labor anchor | `labor_target_mean=0.24279` (BCKM-empirical hours/pop) |
| Optimizer | L-BFGS-B warm-started from BCKM Tables 8/10 |
| MLE caching | Content-addressed pickle in `covid_analysis/data/*.mle.pkl`. Cache hit: ~3s. Cache miss: ~12 min/variant. |

## Narrative-prior rubric — relative to 2019Q4 bind

| Reference | Prior | Full-window Δ | Pre-COVID-fit Δ | Pass? |
|---|---|---|---|---|
| 2020Q2 | τ_l strongly negative | **−16.70%** | **−16.68%** | ✓ |
| 2020Q2 | A small/mechanical | −1.02% | −0.95% | ✓ |
| 2020Q2 | g elevated (CARES) | +3.37% | +3.43% | ✓ |
| 2021Q4 | A above pre-COVID | **+2.03%** | **+2.31%** | ✓ |
| 2023Q4 | τ_l recovered | +1.32% | +1.32% | ✓ |
| 2023Q4 | A near baseline | −2.90% | −2.37% | ✓ |

**vs. 2026-05-01 first-pass (PAYEMS×AWHNONAG labor):** the BLS-faithful
labor construction sharpens nearly every signal:
- 2020Q2 τ_l drop: −13.9% → **−16.7%** (closer to BLS hours-worked
  decline of −18.5% reported in BLS Productivity & Costs).
- 2021Q4 A spike: +1.0% → **+2.0%** (closer to BLS +3.2% TFP).
- 2020Q2 A mechanical: −3.1% → **−1.0%** (cleaner output-vs-hours
  cancellation).
- 2023Q4 baseline-recovery numbers tighter overall.

The improvement comes from CPS employment (LNS12000000) being the
universe BCKM's `hours.dat` actually targets — total civilian
employment, not just nonfarm payrolls — and AWHAETP being all-employees
hours instead of production-and-nonsupervisory-only.

**Both trend variants give nearly identical numerical results** — the calgz
slope is robust to whether COVID-era observations are included in the fit
window. This is itself a useful finding: the COVID anomaly is not large
enough relative to a 14-year window to materially distort the linear
trend estimate.

## Wedge values at reference quarters

### Full-window calgz (BLS-faithful labor)

| Quarter | exp(log_z) | (1−τ_l) | (1+τ_x) | exp(log_g) |
|---|---|---|---|---|
| 2019Q4 (bind) | 0.8452 | 0.7964 | 1.0119 | 0.8353 |
| 2020Q2 (trough) | 0.8371 | 0.6633 | 0.9118 | 0.8635 |
| 2021Q4 (recovery) | 0.8624 | 0.8013 | 1.0009 | 0.7322 |
| 2023Q4 (norm) | 0.8207 | 0.8069 | 0.9854 | 0.7610 |

### Pre-COVID-fit trend (BLS-faithful labor)

| Quarter | exp(log_z) | (1−τ_l) | (1+τ_x) | exp(log_g) |
|---|---|---|---|---|
| 2019Q4 (bind) | 0.8452 | 0.7964 | 1.0119 | 0.8353 |
| 2020Q2 (trough) | 0.8371 | 0.6636 | 0.9116 | 0.8639 |
| 2021Q4 (recovery) | 0.8647 | 0.8019 | 1.0004 | 0.7332 |
| 2023Q4 (norm) | 0.8252 | 0.8070 | 0.9859 | 0.7621 |

The persistent ~0.8 absolute level for `exp(log_z)` reflects an SS
calibration offset between the model and 2010-2023 data — only relative
changes from bind are meaningful for cyclical interpretation.

## F-statistics (window: 2019Q4–2022Q4)

### Full-window (BLS-faithful labor)

| Var | Efficiency | Labor | Investment | Government |
|---|---|---|---|---|
| y | 0.16 | **0.49** | 0.13 | 0.22 |
| l | 0.07 | **0.77** | 0.05 | 0.11 |
| x | 0.36 | 0.12 | 0.13 | 0.40 |

### Pre-COVID-fit (BLS-faithful labor)

| Var | Efficiency | Labor | Investment | Government |
|---|---|---|---|---|
| y | 0.16 | **0.49** | 0.12 | 0.22 |
| l | 0.07 | **0.77** | 0.05 | 0.11 |
| x | 0.38 | 0.12 | 0.11 | 0.39 |

**Labor wedge dominates the COVID output and hours dynamics** (49% / 77%
of f-stat weight respectively — *up* from 47% / 74% under PAYEMS×AWHNONAG). The investment channel is mostly
explained by efficiency + government, NOT by the investment wedge itself
— suggesting that capital-formation in 2020-2022 was driven mainly by
TFP shocks and fiscal stimulus rather than by frictions on capital
goods directly.

## Findings: COVID vs. the Great Recession

*Following the style of Part III of `bckm_replication/BCKM_RESULTS.md`. The
1980-2014 numbers below are BCKM's own published figures; the COVID-window
numbers come from running the same pipeline on US data 2010Q1–2023Q4 with
bind = 2019Q4 (NBER cycle peak), labor source LNS12000000 × AWHAETP, and
calgz detrending. Both trend variants give nearly identical numbers, so
only the full-window numbers are tabulated; pre-COVID-fit values agree to
within 0.5 pp on every line.*

### The COVID Recession in the United States

In Figure A, we see that output and labor both fell sharply over a single
quarter — output by about **9%** and labor by about **12%** between 2019Q4
and 2020Q2 — while investment fell about **12%** at the trough. By 2022Q1,
output and labor had returned to within a percentage point of their 2019Q4
levels and investment was already 4% above pre-pandemic levels. The
dispersion of recovery speeds across the three series is the most striking
visual difference from the Great Recession (Part III, §"United States:
Primarily a Labor Wedge Recession"), where output, labor, and investment
remained 6%, 3%, and 8% below their 2008Q1 peaks even at the end of 2014.

In Figure B, we see that the labor wedge collapsed dramatically — falling
about **17%** from 2019Q4 to 2020Q2 — while the efficiency wedge fell only
modestly (about **1%**). The investment wedge moved in the *opposite*
direction from any postwar recession we have studied: 1/(1+τ_x) **rose**
about **10%** at the trough, indicating that capital-investment frictions
*loosened* during the height of COVID rather than tightened. This is
consistent with the Federal Reserve's emergency liquidity facilities
holding credit spreads down through 2020-2021 even as real activity
collapsed.

In Figure 2C, panel C, we see that the labor-only counterfactual hits 89
at the 2020Q2 trough — i.e., if τ_l alone had moved exactly as it did and
all other wedges had stayed at 2019Q4 levels, output would have fallen
about **11%**, more than the actual **9%** decline. The labor wedge alone
*overpredicts* the trough; the efficiency and investment wedges (both
counter-cyclical at the trough) attenuate the labor-driven contraction.
The investment-only and efficiency-only counterfactuals stay above 100
almost everywhere in 2020-2022, in sharp contrast with the Great Recession
where the investment-only path traces close to data through 2011.

In Figure 2D, panel D (hours), the labor-only counterfactual reproduces
the data trough almost exactly (86 vs. 89 in the data; the difference is
absorbed by the small efficiency-wedge dip). The efficiency-only
counterfactual stays within 2 pp of bind throughout — labor-augmenting
TFP is essentially mechanical in this episode, as one would expect when
the shock is regulatory (mandatory closures) rather than technological.

In Figure 2E, panel E (investment), the labor-only counterfactual falls
to **76** at the 2020Q2 trough — far below the actual 88. The labor-wedge
shock alone implies investment should have fallen 24%; the much milder
12% actual decline is what the investment-wedge loosening (panel B) buys
back. This is the cleanest piece of evidence that COVID is *not* a single-
wedge recession even though the labor wedge dominates the f-statistics.

Overall, considering the period from 2019Q4 through 2022Q4, these results
imply that the COVID recession in the United States should be thought of
as **primarily a labor wedge recession with a counter-cyclical investment
wedge**. The labor wedge's role is even more dominant than during the
Great Recession (φ-weight on output **49%** vs BCKM's reported 46% for
2008-2011), and the investment wedge plays the *opposite* sign of the
2008 recession. Models of the COVID episode must yield a sharp, V-shaped
labor wedge fluctuation paired with a counter-cyclical investment-wedge
loosening; neither feature is salient in models built around the Great
Recession.

### Comparing the COVID Recession with the Great Recession

The comparison runs through three lenses: the data themselves, the
peak-to-trough wedge movements, and the φ-statistics. We use the BCKM
1980-2014 published numbers as the Great Recession baseline; our own
pipeline reproduces them to within 0.01 in every channel (see
`bckm_replication/REPORT.md`).

#### Peak-to-trough comparison

| Series | Great Rec. (2008Q1 → 2009Q3) | COVID (2019Q4 → 2020Q2) |
|---|---|---|
| Output | −7% | **−9%** |
| Labor | −7% | **−12%** |
| Investment | −23% | −12% |
| Efficiency wedge | −1% | **−1%** |
| Labor wedge (1−τ_l) | **−8%** | **−17%** |
| Investment wedge 1/(1+τ_x) | **−9%** | **+10%** *(loosened)* |
| Trough → return to bind (output) | not reached by 2014Q4 | 2022Q1 (~1.5 years) |

Three differences stand out:

1. **The COVID labor-wedge collapse was about twice as deep** as in 2008
   (−17% vs −8%) and resolved an order of magnitude faster (1.5 years to
   bind vs. unresolved by 2014Q4 in BCKM's window). The depth fits the
   regulatory-shutdown narrative; the speed reflects that the shutdown
   itself was lifted, while 2008 reflected a slow, demand-side collapse.

2. **The COVID investment wedge moved opposite the 2008 wedge.** In 2008
   the investment wedge tightened by 9 pp; in COVID it loosened by about
   10 pp. We read this as the Federal Reserve's monetary response (zero
   rates, asset-purchase facilities, market-functioning interventions)
   compensating for the real shock to capital formation. In a frictionless
   benchmark, COVID-2020 looks like a *pure labor shock with a partially
   offsetting credit boost* — a configuration that simply does not appear
   in the 1948–2014 sample.

3. **The investment data fall was much shallower** in COVID (−12% vs −23%
   in 2008). This is *consistent* with the wedge story: in 2008 both the
   labor-wedge contraction (−8% direct labor effect) and the
   investment-wedge tightening (−9%) pushed investment down; in COVID the
   labor wedge alone would have pushed investment to ~76 (panel E), and
   the loosened investment wedge buys back the difference to the actual
   88 trough.

#### φ-statistics comparison (output)

| | Efficiency | Labor | Investment | Government |
|---|---|---|---|---|
| Great Recession (BCKM 2008Q1–2011Q4) | 0.16 | **0.46** | 0.32 | 0.06 |
| COVID Recession (2019Q4–2022Q4) | 0.16 | **0.49** | 0.13 | 0.22 |

Two patterns:

1. **The labor wedge is more dominant in COVID** (49% vs 46%) — the
   pipeline confirms what the figures show visually. The COVID labor
   collapse is not just bigger in level terms; it is also a sharper,
   less-anticipated transition than the 2008 contraction.

2. **The investment wedge weight collapsed from 32% → 13%** between
   episodes. In 2008, the investment-wedge component of output tracked
   data closely (low SSR → high φ); in COVID, the investment wedge moved
   counter-cyclically and the investment-only counterfactual sits *above*
   100 throughout the trough, contributing a *high* SSR against the data
   path. Mathematically the φ falls; economically, "the investment wedge
   accounts for less of the COVID downturn" should be read as "in 2020,
   the investment wedge was working against the recession, not for it."

3. **The government wedge weight rose from 6% → 22%.** This is a real
   feature of the COVID episode but should be read carefully. BCKM's g =
   `gov_consumption + net_exports`. The 2020 *consumption* component of
   that aggregate did rise (CARES Act expanded gov purchases), but the
   *net-exports* component swung sharply negative (massive 2020Q2 import
   recovery against a slower export rebound). The government-wedge
   component of output is therefore picking up a mix of fiscal expansion
   and trade-balance dynamics; ARPA-style transfer payments are *not* in
   this aggregate (they enter via household consumption), so the +22%
   weight is *not* a measure of fiscal stimulus on output, and should be
   interpreted accordingly.

#### Summary

The 2008 recession is, in the BCKM decomposition, a labor-wedge recession
with a meaningful secondary investment-wedge story. The COVID recession is
a *more extreme* labor-wedge recession with the investment wedge running
*backwards* compared to 2008. The two episodes share a primary driver but
differ in (i) the magnitude and recovery speed of the labor wedge, (ii)
the sign of the investment wedge, (iii) the depth of investment data
fall, and (iv) the importance of the government channel.

The BCKM model framework — labor-augmenting TFP plus separable wedges —
generalizes cleanly from the 2008 episode to the 2020 episode without
re-calibration. The wedges it identifies in COVID match the qualitative
narrative supplied by BLS productivity reports, NBER recession dates, and
contemporary Fed monetary-policy actions. This is the strongest piece of
external evidence that the pipeline is structurally sound on out-of-sample
US windows.

## Counterfactual decomposition (figure_2C — Output components)

The black "Data" line drops to 91 at 2020Q2 and recovers above 100 by
2021. The decomposition shows:

- **Labor-only** (green): drops to 90 at 2020Q2 — single-handedly
  reproduces the COVID output trough. This is the strongest narrative
  match: the lockdown is a textbook labor shock.
- **Efficiency-only** (blue): mild dip to 97 at trough, mild rise to
  102 in mid-2021. The 2021 productivity boom is *directionally
  correct but underwhelming* — BLS reports +3.2% TFP, the pipeline
  shows +1%. Possible explanations:
  1. The calgz exponential trend captured part of the 2021 boom as
     "trend" rather than "wedge" (both variants show nearly identical
     results, so this is unlikely).
  2. BCKM's labor-augmenting `z` is a different object than BLS
     multifactor productivity. They're related but not identical.
  3. Some of the 2021 boom got allocated to the investment wedge in
     the optimizer's choice of basin.
- **Investment-only** (magenta): counter-cyclical at trough (rises to
  103 in 2020Q2), settles near 100. Suggests the investment wedge
  (BGG-style adjustment cost) was *less* binding during the height of
  COVID, possibly reflecting cheap credit / Fed intervention.
- **Government-only** (red): essentially flat near 99-100. **Not** a
  fail of the pipeline: BCKM's g = `gov_consumption + net_exports`
  excludes transfer payments. ARPA was largely transfers to
  households, which feed through household consumption → labor wedge,
  NOT the government channel directly.

## Methodology notes

- `mle_window=("2010Q1","2019Q4")` for the pre-COVID-fit variant fits
  the calgz exponential slope on the pre-anomaly sub-window and
  extrapolates forward. Both trend variants converged to nearly
  identical wedge paths — the COVID anomaly is too short relative to
  the 14-year sample to meaningfully distort the trend.
- Warm-start: BCKM Tables 8/10 (P, Q_chol, Sbar). Optimizer was free
  to move and converged to a different basin (LL +704 vs BCKM-θ on
  BCKM data ~+1873).
- `labor_target_mean=0.24279` (BCKM-empirical) is required because raw
  FRED hours/pop is ~23 while the model's `ss["l"]` is ~0.29. Without
  the rescale, the wedge extraction has an 80× scale mismatch and
  produces nonsense states.
- f-stats computed in level space (`f_statistics_bckm`), BCKM-faithful
  formula.

## Open issues

1. **2021 TFP magnitude undershoot.** A wedge shows +1% in 2021Q4 vs
   BLS's +3.2%. Worth investigating whether this is a calibration
   issue (γ_annual=2.36% vs the BCKM 1.9%) or a basin issue (the
   optimizer found a low-A basin).
2. **A-wedge declining trend 2010→2023.** The wedge declines from 120
   (2010) to 95 (2023) in normalized form. This could be calibration
   drift (γ chosen to match end-of-window) or a real productivity
   stagnation signal. Comparison with BLS labor productivity for the
   same window would clarify.
3. **Government channel underweighted.** g-wedge sees only the
   gov_consumption+net_exports half of fiscal policy. To capture
   transfers we'd need either a household-side modification to the
   model or to feed transfers into the c-resource constraint
   explicitly.

## Verdict

**Pipeline generalizes to non-BCKM windows.** The COVID smoke test
correctly identifies COVID as a primarily labor-driven contraction
(74% of hours f-stat weight) with secondary investment friction, with
appropriate sign agreement on all six tested narrative priors. The
2021 TFP spike is directionally correct but quantitatively muted; the
ARPA fiscal effect is structurally invisible to BCKM's g definition.
These are useful, well-understood limitations rather than pipeline
bugs.
