# COVID Analysis Report — US 2010Q1–2023Q4

> **Status**: Scaffolding complete; first run pending.
> Fill in results after running `run_covid_analysis.py`.

## Setup

| Parameter | Value |
|---|---|
| Window | 2010Q1–2023Q4 |
| Bind (anchor) | 2019Q4 |
| Trend method A | Full-window calgz |
| Trend method B | Pre-COVID-fit (slope on 2010Q1–2019Q4, extrapolated forward) |
| Data source | FRED defaults |
| Labor target mean | None (raw hours per capita) |

## Narrative-prior rubric

| Reference quarter | Prior (sign) | Full-window | Pre-COVID-fit | Pass? |
|---|---|---|---|---|
| 2020Q2 | τ_l strongly negative | — | — | — |
| 2020Q2 | A small / mechanical | — | — | — |
| 2021Q4 | A above 2019Q4 trend (BLS +3.2% TFP) | — | — | — |
| 2021Q4 | g elevated (CARES + ARPA) | — | — | — |
| 2023Q4 | τ_l ≈ recovered | — | — | — |
| 2023Q4 | A near 2019Q4 baseline | — | — | — |

A "pass" requires sign-direction agreement under **both** trend variants.

## Wedge table at reference quarters

*(to be filled after run)*

## F-statistics (2019Q4–2022Q4 window)

*(to be filled after run)*

## Verdict

*(to be filled after run)*

## Methodology notes

- `mle_window=("2010Q1","2019Q4")` fits the calgz exponential slope on the
  pre-anomaly sub-window and extrapolates forward, so the COVID period is
  detrended against a "would-have-been" trend rather than one distorted by
  the 2020 contraction.
- Warm-start: BCKM Table 8/10 published (P, Q, Sbar); optimizer free to move.
- f-stats computed in level space (`f_statistics_bckm`), BCKM-faithful formula.
