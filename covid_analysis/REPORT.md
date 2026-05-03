# Business Cycle Accounting for the COVID Recession: United States 2010–2023

## Abstract

We apply the Business Cycle Accounting (BCA) framework of Chari, Kehoe, and
McGrattan (2007) to the United States over 2010Q1–2023Q4, anchored at the 2019Q4 NBER
cycle peak. The framework decomposes aggregate fluctuations into four structural wedges —
efficiency, labor, investment, and government — via Kalman-filter MLE on a VAR(1), and
uses counterfactual simulations to quantify each wedge's contribution to output, hours, and
investment dynamics. The labor wedge accounts for 49% of the output variance and 77% of
the hours variance over 2019Q4–2022Q4. In contrast to the Great Recession of 2008–2009,
the investment wedge moves counter-cyclically during COVID — loosening approximately 10%
at the 2020Q2 trough — consistent with the Federal Reserve's emergency credit interventions
partially offsetting the labor shock. Results are robust to the choice of detrending window
and replicate directionally against BLS Productivity and Costs narratives.

---

## 1. Introduction

Business Cycle Accounting (BCA) is a model-based decomposition of aggregate fluctuations
into orthogonal structural distortions: an efficiency wedge (total factor productivity), a
labor wedge (a distortion between the household's marginal rate of substitution and the
marginal product of labor), an investment wedge (a capital-goods friction or credit spread),
and a government wedge (government consumption plus net exports). The methodology,
developed by Chari, Kehoe, and McGrattan (2007) — hereafter CKM — and extended to
a 24-country comparison by Brinca, Chari, Kehoe, and McGrattan (2016) — hereafter BCKM
— has been applied primarily to samples ending before 2015.

The COVID episode provides a particularly clean out-of-sample test for the BCA framework.
Its sharp, regulatory labor shock (mandatory business closures), emergency monetary
response (Federal Reserve liquidity facilities, zero-rate policy), large fiscal transfers
(CARES Act, ARPA), and V-shaped recovery are qualitatively distinct from any postwar
recession in the BCKM dataset. If the BCA framework correctly identifies the dominant
structural driver — the labor wedge — and captures the counter-cyclical credit loosening,
that is meaningful evidence that the methodology is structurally sound outside the sample
used to develop it.

This report presents a BCA decomposition for the United States over 2010Q1–2023Q4. We
use the same model structure and estimation methodology as BCKM and verify that our
pipeline reproduces their published Great Recession f-statistics for 1980Q1–2014Q4 to
within 0.01 in every channel before applying it to the COVID window. The two central
questions are: (1) Which wedge accounts for most of the 2020 contraction? (2) Does the
COVID episode share the structural driver profile of the Great Recession, or is it
qualitatively distinct?

---

## 2. Methodology

### 2.1 The Model

We use the prototype economy of CKM (2007) and BCKM (2016): a representative-agent,
deterministic-trend model with a Cobb-Douglas production function, Y = A · K^α · (z · L)^(1−α),
where A is the efficiency wedge, z grows at trend rate γ, and four structural wedges enter
the resource constraint, the household's intratemporal labor condition, and the Euler
equation for capital. The four wedges are:

- **Efficiency wedge** (A): shifts the production function; interpreted as total factor
  productivity in the broad sense.
- **Labor wedge** (1−τ_l): drives a wedge between the marginal rate of substitution of
  leisure for consumption and the marginal product of labor. A decline in (1−τ_l) can
  arise from labor-market frictions, payroll distortions, or — in the COVID context —
  regulatory shutdowns that raise the effective cost of working.
- **Investment wedge** (1+τ_x): drives a wedge in the Euler equation for capital,
  reflecting credit frictions, adjustment costs, or intermediation spreads. A fall in
  (1+τ_x) corresponds to a loosening of investment frictions.
- **Government wedge** (g): government consumption plus net exports, following BCKM's
  closed-economy equivalence (CKM 2005). Transfer payments do not enter g directly —
  they pass through household income and ultimately the labor and efficiency wedges.

The model is log-linearized around its deterministic steady state. The wedges are identified
as the unique set of processes that make the linearized equilibrium conditions hold exactly
at every date in the sample. Calibration follows BCKM: α = 1/3, Frisch elasticity ψ = 2.5,
annual depreciation δ = 5%, annual discount factor β = 0.975.

### 2.2 Estimation

The four wedges follow a first-order VAR: s_{t+1} = P · s_t + ε_t, where
s = (log A, τ_l, τ_x, log g)′ and ε ~ N(0, Q·Q′). The parameters (S̄, P, Q) are estimated
by Kalman-filter maximum likelihood on the four observables [y, l, x, g] using the
steady-state Kalman filter with a DARE-derived gain, initialized from BCKM's published
Tables 8 and 10 (P matrix and Cholesky factor Q). Full details of the estimation
methodology are given in BCKM (2016), Appendix A.

### 2.3 Counterfactuals and φ-Statistics

We construct *incremental* single-wedge counterfactuals following BCKM (their Section 3
and Appendix B). The counterfactual for wedge j simulates the model forward with wedge j
moving as estimated and all other wedges held at their 2019Q4 values, expressed as the
incremental deviation from the no-wedge baseline (all wedges fixed). The four
single-wedge counterfactuals are additive: their sum reproduces the data path.

The φ-statistic for wedge j on variable v is the share of the sum-of-squared residuals
relative to the 2019Q4 anchor accounted for by wedge j's counterfactual over the window
2019Q4–2022Q4. A φ of 0.49 for the labor wedge on output means that 49% of the cumulative
squared deviation of output from its 2019Q4 level is explained by the labor-wedge
counterfactual.

---

## 3. Data

The sample is 2010Q1–2023Q4 (56 quarters). Output (y) is real GDP per working-age person;
investment (x) is gross private domestic investment plus consumer durable expenditures plus
government investment, deflated by the GDP deflator, per working-age person; the government
aggregate (g) is government consumption plus net exports, deflated, per working-age person.
All three series are sourced from FRED. Working-age population is the civilian
non-institutional population aged 16 and over.

Labor input (l) is constructed as the product of CPS civilian employment (BLS series
LNS12000000, seasonally adjusted) and average weekly hours of all employees in total private
industries (AWHAETP), scaled to a quarterly total. This construction mirrors the BCKM
source data (CPS employment × BLS average weekly hours) more faithfully than aggregated
payroll employment indices: LNS12000000 covers the full civilian workforce, not just
nonfarm payroll workers, and AWHAETP covers all employees rather than
production-and-nonsupervisory staff only.

All four series are detrended using BCKM's calgz method: a trend growth rate γ is chosen
numerically such that log(y per capita) has mean zero over the estimation window. We report
results under two trend variants:

- **Full-window** (headline): γ_annual = 2.36%, estimated on 2010Q1–2023Q4.
- **Pre-COVID-fit** (robustness): γ_annual = 2.25%, estimated on 2010Q1–2019Q4 and
  extrapolated forward, so that the pandemic anomaly does not influence the trend estimate.

The anchor quarter is 2019Q4 — the NBER-dated business cycle peak immediately preceding
the COVID recession. All figures and φ-statistics report deviations relative to this base
period (anchor = 100).

---

## 4. Results: The COVID Recession in the United States

### 4.1 Macroeconomic Aggregates

Figure 1 plots the detrended series for output, labor, and investment, normalized to 100 at
2019Q4. Output and labor both fell sharply in a single quarter: output fell approximately
**9%** and labor approximately **12%** between 2019Q4 and 2020Q2. Investment fell about
**12%** at the trough. By 2022Q1, output and labor had returned to within one percentage
point of their pre-pandemic levels, and investment was already 4% above pre-pandemic
levels. The V-shaped recovery — approximately six quarters from trough to anchor — is the
most striking visual difference from the Great Recession, where output, labor, and
investment remained 6%, 3%, and 8% below their 2008Q1 peaks even through the end of
BCKM's 2014Q4 sample window.

### 4.2 Structural Wedge Paths

Figure 2 plots the four structural wedge paths normalized to their 2019Q4 levels. The labor
wedge (1−τ_l) collapsed by approximately **17%** between 2019Q4 and 2020Q2 — the
dominant feature of the decomposition. The efficiency wedge fell by only about **1%**,
consistent with the BLS characterization of COVID as primarily a utilization and hours
shock rather than a multifactor productivity shock.

By 2021Q4, the efficiency wedge had recovered to +2% above its 2019Q4 level, capturing the
BLS-documented TFP surge of 2021. The BLS reported +3.2% multifactor productivity growth
for that year, the largest single-year gain since 1983. Our model-based estimate is
approximately 1.2 pp smaller; the gap is discussed in Section 6.

The investment wedge (1+τ_x) moved in the *opposite* direction from any postwar US
recession in the BCKM dataset: it *fell* approximately **10%** at the 2020Q2 trough,
indicating that capital-investment frictions loosened during the height of the pandemic.
This is consistent with the Federal Reserve's emergency liquidity facilities — zero lower
bound policy, large-scale asset purchases, and targeted credit facilities — holding
financial conditions accommodative through 2020–2021 even as real activity collapsed.

**Table 1** reports wedge levels at four reference quarters for both trend variants.

**Table 1: Structural wedge levels at reference quarters**

| Quarter | Efficiency (A) | Labor (1−τ_l) | Investment (1+τ_x) | Government (g) |
|---|---|---|---|---|
| 2019Q4 — anchor | 0.845 | 0.796 | 1.012 | 0.835 |
| 2020Q2 — trough | 0.837 | 0.663 | 0.912 | 0.864 |
| 2021Q4 — recovery | 0.862 | 0.801 | 1.001 | 0.732 |
| 2023Q4 — end | 0.821 | 0.807 | 0.985 | 0.761 |

*Notes: Full-window trend variant. Pre-COVID-fit variant agrees to within 0.5 pp at every
reference quarter. The absolute level of A ≈ 0.84 reflects a steady-state calibration
offset between the prototype model and 2010–2023 data; only changes relative to the 2019Q4
anchor are cyclically interpretable. A lower value of (1+τ_x) indicates looser investment
frictions.*

The two trend variants are nearly identical throughout the sample (differences of at most
0.5 pp at any reference quarter), confirming that the pandemic anomaly is not large enough
relative to the 14-year window to materially distort the trend estimate.

### 4.3 Counterfactual Decomposition

**Figure 3** shows the output counterfactuals. The labor-only counterfactual drops to 89 at
the 2020Q2 trough, reproducing — and slightly overpredicting — the actual 9% output
decline. The labor wedge alone accounts for essentially the entire output trough; the
efficiency and investment counterfactuals remain above 100 through 2020, attenuating rather
than amplifying the labor-driven contraction. The government-only counterfactual is
essentially flat near 100 throughout: the government wedge (government consumption plus
net exports) excludes the CARES Act and ARPA transfer payments, which primarily entered
household income and do not directly affect g.

**Figure 4** shows the hours counterfactuals. The labor-only path reproduces the data trough
nearly exactly (counterfactual: 86; data: 89). The efficiency-only path stays within 2 pp
of the anchor throughout — labor-augmenting TFP is essentially acyclical in this episode,
consistent with a shock driven by regulatory closures rather than technological disruption.

**Figure 5** shows the investment counterfactuals. The picture here is more complex. The
labor-only counterfactual falls to 76 at the trough, implying investment should have fallen
24% under a pure labor shock, while actual investment fell only 12%. The loosening of the
investment wedge recovers the 12 pp gap. This is the clearest evidence that COVID is not a
pure single-wedge recession even though the labor wedge dominates the φ-statistics: the
investment wedge was working *against* the recession in 2020, not for it, and its
counter-cyclical movement substantially attenuated the capital-formation decline.

### 4.4 φ-Statistics

**Table 2** reports the φ-statistics for the COVID window (2019Q4–2022Q4).

**Table 2: φ-statistics, COVID window 2019Q4–2022Q4**

| Variable | Efficiency | Labor | Investment | Government |
|---|---|---|---|---|
| Output (y) | 0.16 | **0.49** | 0.13 | 0.22 |
| Hours (l) | 0.07 | **0.77** | 0.05 | 0.11 |
| Investment (x) | 0.36 | 0.12 | 0.13 | 0.40 |

*Notes: Full-window and pre-COVID-fit variants agree to within 0.01 in every cell. Rows
sum to 1.0 by construction.*

The labor wedge accounts for 49% of the output φ-statistic and 77% of the hours
φ-statistic. The investment-wedge contribution to output is only 13% — lower than even the
efficiency-wedge's 16% — reflecting the counter-cyclical movement of the investment wedge
during this episode. The government wedge's 22% weight on output is elevated relative to
the Great Recession (6% in BCKM) and is discussed in Section 6.

---

## 5. Comparison with the Great Recession

### 5.1 Peak-to-Trough Movements

**Table 3** compares peak-to-trough movements across the two recessions. The Great
Recession baseline uses BCKM's published figures for 2008Q1–2009Q3; COVID numbers come
from this analysis.

**Table 3: Peak-to-trough comparison, Great Recession vs. COVID**

| Series | Great Recession (2008Q1→2009Q3) | COVID (2019Q4→2020Q2) |
|---|---|---|
| Output | −7% | −9% |
| Labor | −7% | −12% |
| Investment | −23% | −12% |
| Efficiency wedge (A) | −1% | −1% |
| Labor wedge (1−τ_l) | −8% | −17% |
| Investment wedge (1+τ_x) | +9% *(tightened)* | −10% *(loosened)* |
| Quarters from trough to anchor | > 24 (not recovered by 2014) | ≈ 6 |

### 5.2 φ-Statistics Comparison

**Table 4** places the COVID φ-statistics alongside the Great Recession BCKM results for
the output channel.

**Table 4: φ-statistics for output, COVID vs. Great Recession**

| Episode | Efficiency | Labor | Investment | Government |
|---|---|---|---|---|
| Great Recession (BCKM, 2008Q1–2011Q4) | 0.16 | **0.46** | 0.32 | 0.06 |
| COVID (this analysis, 2019Q4–2022Q4) | 0.16 | **0.49** | 0.13 | 0.22 |

### 5.3 Three Structural Differences

**1. The COVID labor-wedge collapse was deeper and faster.** The labor wedge fell 17%
at the 2020Q2 trough relative to 2019Q4, versus 8% for the Great Recession — roughly twice
the depth. Recovery was also far faster: the labor wedge returned to within 2 pp of its
pre-recession level within six quarters, while it remained depressed through the end of
BCKM's 2014 sample. The depth is consistent with the regulatory-shutdown narrative:
mandatory business closures and elevated infection risk suppressed labor supply more sharply
than the 2008 demand-side deterioration. The speed of recovery reflects that the shutdown
was eventually lifted, whereas 2008's balance-sheet recession resolved gradually.

**2. The investment wedge moved counter-cyclically.** In the Great Recession, the investment
wedge (1+τ_x) tightened by approximately 9% — consistent with the Bernanke-Gertler-Gilchrist
financial-accelerator mechanism, where a collateral shock raises credit spreads and reduces
investment. In COVID, (1+τ_x) fell approximately 10% at the trough. This reversal is
consistent with the Federal Reserve's emergency interventions in 2020–2021: zero-lower-bound
policy, asset purchases, and targeted credit facilities kept borrowing costs historically
low even as real activity collapsed. COVID-2020 looks, within this framework, like a
*pure labor shock with a partially offsetting credit loosening* — a configuration that does
not appear in the 1948–2014 US business cycle record.

**3. The investment data fell much less in COVID.** The peak-to-trough investment decline
was 12% in COVID versus 23% in 2008. This is internally consistent with the wedge story:
in 2008, both the labor-wedge contraction and the investment-wedge tightening pushed
investment down; in COVID, the labor wedge alone would have pushed investment to 76 (Figure 5),
and the investment-wedge loosening recovered the difference to the observed trough of 88.

### 5.4 Summary

The Great Recession is, in the BCKM decomposition, a labor-wedge recession with a
meaningful secondary investment-wedge story. The COVID recession is a *more extreme*
labor-wedge recession with the investment wedge running backwards relative to 2008. The two
episodes share a primary structural driver but differ in (i) the magnitude and recovery
speed of the labor wedge, (ii) the sign of the investment wedge, (iii) the depth of the
investment decline, and (iv) the prominence of the government channel. The efficiency
wedge contribution is strikingly stable across episodes: 0.16 on output in both, reflecting
that TFP is a minor contributor to cyclical fluctuations in either episode relative to the
labor and investment channels.

---

## 6. Discussion

### 6.1 The 2021 TFP Recovery

The efficiency wedge rises to +2% above its 2019Q4 level by 2021Q4, directionally
consistent with BLS Multifactor Productivity data reporting a +3.2% TFP gain for 2021 —
the largest single-year increase since 1983. The model-based estimate is approximately
1.2 pp below the BLS figure.

Two structural explanations are plausible. First, BCKM's labor-augmenting efficiency
wedge A is not identical to Hicks-neutral multifactor productivity: A captures shifts in the
way labor translates to effective output, while BLS MFP is a broader concept that includes
capital-augmenting components not separately identified in the BCKM model. Second, the
optimizer may have attributed some of the 2021 TFP boom to the investment wedge rather
than the efficiency wedge in its choice of VAR(1) parameters; the full-window and
pre-COVID-fit variants produce nearly identical results, making it unlikely that the trend
specification is responsible. On balance, the undershoot is better interpreted as a
measurement-concept difference than as a model misspecification.

### 6.2 Government Wedge Interpretation

The government wedge accounts for 22% of output variance over 2019Q4–2022Q4, a
substantial increase relative to the Great Recession (6% in BCKM). This requires careful
interpretation. BCKM's government aggregate g equals government consumption expenditures
plus net exports. The CARES Act (2020, approximately \$2.2 trillion) and ARPA (2021,
approximately \$1.9 trillion) were predominantly transfer payments — stimulus checks,
expanded unemployment insurance, Paycheck Protection Program loans — which enter household
income and affect the economy through the household budget constraint, ultimately appearing
in the labor and consumption channels, not in g. What the government wedge does capture is
(a) the increase in direct government purchases in 2020 (procurement of medical supplies,
expanded federal operations) and (b) the sharp deterioration of net exports in 2020Q2, as
imports rebounded faster than exports following the initial shock. The elevated government
φ-weight should not be interpreted as evidence that fiscal stimulus drove the output
recovery; it reflects the government-purchases-plus-trade-balance aggregate, a broader and
more ambiguous concept than discretionary fiscal policy.

### 6.3 Robustness to Trend Specification

All results are robust to the choice of trend estimation window. The full-window calgz
(γ = 2.36%/yr) and the pre-COVID-fit calgz (γ = 2.25%/yr, extrapolated forward) produce
wedge paths that agree to within 0.5 pp at every reference quarter and φ-statistics that
agree to within 0.01. The COVID episode, though large in absolute terms, spans only 8 of
the 56 sample quarters and is not large enough to materially distort the log-linear trend
estimated over the full window.

---

## 7. Conclusion

The Business Cycle Accounting framework, applied to US data for 2010Q1–2023Q4, identifies
the COVID recession of 2020 as primarily a labor-wedge recession. The labor wedge accounts
for 49% of the output φ-statistic and 77% of the hours φ-statistic over 2019Q4–2022Q4 —
more dominant than its 46% weight in the Great Recession. The efficiency wedge plays a
minor, mechanical role in both episodes (16% of output variance in each), while the
government-wedge contribution is elevated in COVID, largely reflecting trade-balance
dynamics rather than direct fiscal stimulus effects on output.

The novel structural finding is the counter-cyclical investment wedge: a loosening of
approximately 10% at the 2020Q2 trough, in contrast to the 9% tightening seen in 2008.
We interpret this as the Federal Reserve's emergency credit interventions partially
offsetting the real collapse in labor and output. Without this counter-cyclical loosening,
the labor shock alone would have pushed investment down by 24%; the observed 12% decline
reflects the credit loosening absorbing roughly half of that pressure.

These findings have implications for model design. A model aimed at quantitatively
replicating the COVID recession must generate a sharp, V-shaped labor-wedge collapse of
approximately 17% in a single quarter, paired with a counter-cyclical investment-wedge
loosening. Standard financial-accelerator models — which generate procyclical
investment-wedge tightening — predict the wrong sign for 2020. Models that combine
regulatory-closure-type labor supply disruptions with an explicit monetary-policy block
capturing the Federal Reserve's quantitative easing are natural candidates.

The broader methodological finding is that the BCKM pipeline generalizes cleanly from the
2008 Great Recession to the 2020 COVID Recession without re-calibration. The structural
wedges it identifies are consistent with the BLS productivity narrative, the NBER recession
chronology, and the historical record of monetary-policy interventions in 2020–2021. This
provides evidence that the Business Cycle Accounting methodology is structurally sound on
out-of-sample US windows spanning a qualitatively distinct type of recession.

---

## References

- Brinca, P., Chari, V. V., Kehoe, P. J., and McGrattan, E. R. (2016). "Accounting for
  Business Cycles." *Handbook of Macroeconomics*, Vol. 2, pp. 1013–1063.
- Chari, V. V., Kehoe, P. J., and McGrattan, E. R. (2007). "Business Cycle Accounting."
  *Econometrica*, 75(3), 781–836.
- Chari, V. V., Kehoe, P. J., and McGrattan, E. R. (2005). "Sudden Stops and Output
  Drops." *American Economic Review*, 95(2), 381–387.

---

## Appendix: Sample and Estimation Details

| Parameter | Value |
|---|---|
| Sample | 2010Q1–2023Q4 (56 quarters) |
| Anchor quarter | 2019Q4 (NBER cycle peak) |
| Trend — full-window | calgz on 2010Q1–2023Q4; γ_annual = 2.36% |
| Trend — pre-COVID-fit | calgz on 2010Q1–2019Q4; γ_annual = 2.25% |
| Labor construction | CPS employment (LNS12000000) × AWHAETP × 13 weeks/quarter |
| Labor normalization | Population-adjusted to BCKM-empirical hours/population ratio |
| φ-statistic window | 2019Q4–2022Q4 |
| Calibration | α = 1/3, ψ = 2.5, δ = 5%/yr, β = 0.975/yr (BCKM Table 1) |
| Estimation | Kalman-filter MLE; L-BFGS-B; warm-started from BCKM Tables 8/10 |
