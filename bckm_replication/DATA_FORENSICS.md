# Data Forensics: BCKM 2016 US Replication

_Last updated: 2026-05-01_

This document records the full data-layer forensics performed while
trying to close the residual log-likelihood and f-statistic gap between
our Python BCA pipeline and BCKM 2016's published US results
(Tables 8–11). It exists because the lessons here are easy to lose to
context compaction and because they matter for **future periods and
countries**, not just for the 1980Q1–2014Q4 US replication.

The TL;DR is at the top. Mechanisms, decisions, and rejected
hypotheses are below it.

---

## TL;DR

1. **Every BEA-NIPA toggle makes the LL gap to BCKM `worktemp.mat`
   monotonically *worse*.** Our FRED-default pipeline is already the
   best operating point we have on this dataset and time window.
2. **The y-channel BEA migration passes the Stage-1 gate**
   (mean|diff| = 0.0074 vs FRED 0.0068; both well under the 0.025
   gate). y is essentially a wash.
3. **The x-channel BEA migration fails the Stage-1 gate by 23×**
   (BEA 0.1535 vs FRED 0.0065). The cause is **BEA NIPA vintage drift**
   accumulated through the 2018 and 2023 comprehensive revisions —
   not a construction bug in our BEA branch.
4. **The g-channel BEA migration fails the Stage-1 gate by 6×**
   (BEA 0.2100 vs FRED 0.0345). The cause is the **bind-ratio
   coupling** in BCKM's `maketrend.m`: every observable is anchored
   at `y_pc(2008Q1)`, so a per-channel level gate is order-coupled
   with the y construction.
5. **FRED's single-deflator approach is "accidentally robust"** to BEA
   vintage drift in the x and g channels because nominal numerator and
   GDPDEF denominator are revised in tandem.
6. **All three BEA branches (`y_source`, `x_source`, `g_source`) are
   preserved as opt-in diagnostic infrastructure**, not deleted. They
   are well-documented and gate-deferred so a future country/period
   exercise can re-enable them and re-run the joint walkdown.
7. **Productionization decision (user-validated 2026-05-01)**: keep
   `{y, x, g}_source = "fred"` as the default. Confidence in
   generalizability matters more than chasing the last few nats of
   numerical match to BCKM's specific 2014-vintage data pull.

---

## Goal reframing

The original framing of this work was "match Tables 8–11 to ~3–4
decimals." That framing is wrong for a generalizable code base.

The realised framing is:
**we want confidence that our code produces credible BCA decompositions
for arbitrary country/period combinations.**

BCKM's 2014 numerical results are the validation target for the *one*
period for which we have published ground truth (`worktemp.mat`). They
are not the goal. If our pipeline reproduces BCKM's *qualitative*
decomposition (which channel drives which downturn) and matches their
quantitative results to f-stat tolerance ~0.01–0.05 across the four
wedges, the pipeline is trusted and can be applied to other periods.

This reframing is consistent with the bigger project goal in
[`bca-web-app-instructions.md`](bca-web-app-instructions.md): a
three-layer web app (`bca_core` Python computational engine + `bca_api`
FastAPI service + `bca_web` React/Recharts front end) that performs BCA
for any country/period combination. That goal does not require an
exact match to BCKM 2014 numerics; it requires that the pipeline be
**right in mechanism** so a 2025 user can analyze, e.g., the 2020
COVID recession or a non-US country with a different vintage of data.

---

## The replication targets

For US 1980Q1–2014Q4 (BCKM Tables 8–11):

| Target | Source in `worktemp.mat` | Status as of 2026-05-01 |
|---|---|---|
| Sbar (4-vec) | `mle.sbar` | matched at BCKM-θ; our optimizer drifts to a nearby attractor |
| P (4×4) | `mle.P` | matches Table 8 to ≤4.3e-5 (canonical in `bca_core/constants.py`) |
| Q (4×4) | `mle.Q_chol` | matches Table 10 to ≤4.3e-5 (canonical) |
| LL at BCKM-θ | `mle.likelihood = −2402.88` | our LL (BCKM-form) at BCKM-θ on **BCKM `Y_raw`** = −2386.28 — 16.6 nats off |
| LL gap on **our df** | data-construction-driven | +1719.66 (FRED defaults) — 683 nats short of BCKM-form |
| Table 11 f-stats (A, τ_l, τ_x, g) | `tables[...]` = (0.16, 0.46, 0.32, 0.06) | at BCKM-θ on our df: (0.151, 0.485, 0.306, 0.058) — **matches to ~0.01** |
| Table 11 f-stats at our converged θ | — | (0.128, 0.630, 0.181, 0.061) — fY[τ_l] over by 0.17, fY[τ_x] under by 0.14 |

**Key observation**: at BCKM-θ on our df, we already match Table 11 to
~0.01 in every cell. The remaining problem is that our optimizer walks
to a different (Sbar, P, Q) basin where the f-stat distribution
re-weights from τ_x toward τ_l. That's a basin-escape / identification
problem, *not* a structural pipeline bug. (See `Diary.md` and the
`Findings (most recent first)` section of `CLAUDE.md` for the basin
diagnostics.)

---

## The data layer: FRED vs BEA

BCKM's `usdata.m` constructs five real per-capita observables from BEA
NIPA tables (T10105, T10106, T10109, T30904, T30905, FAAt101, FAAt103)
and BLS series (hours, civilian population, armed forces). It writes
the result to `worktemp.mat` as `Y_raw`.

Our pipeline initially used FRED tickers as a more convenient proxy
(GDPC1, GDPDEF, etc.). The FRED approach is "BCKM-adjacent": it
re-aggregates the same underlying NIPA components but does so by
dividing the **nominal aggregate** by a **single deflator (GDPDEF)**,
where BCKM divides each component by its own chain-real deflator
before summing.

For most channels these two approaches give numerically similar
results because the nominal numerator and GDPDEF denominator move
together. But they can diverge in specific cases:

- **Terms-of-trade movements**: when pIM ≠ pEX (e.g. 2008 oil shock),
  BEA's chain-real `g = rEX − rIM` can differ from FRED's
  `(EX − IM)/GDPDEF` by a few percent of GDP.
- **Composition shifts**: when a chain-real sub-aggregate's mix
  changes over time (e.g. consumer durables shifting from cars to
  electronics), BEA's chain-real series picks up the index drift
  while a single-deflator series does not.
- **Vintage revisions**: see the "BEA vintage drift" deep dive below.

We tested all three relevant channels (y, x, g) by adding
`y_source`/`x_source`/`g_source` flags to `bca_core/data/pipeline.py`,
each switching between a FRED legacy branch and a BCKM-faithful BEA
branch.

---

## Per-channel migration outcomes

### y-channel — PASSES Stage-1 gate

**BCKM construction** (`usdata.m:51`):
```
Y = rGDP − rSTX + 0.04·rKCD + rDCD
```
where `rKCD` is consumer-durables real net stock (FAAt101 line 15,
quarterly via log-linear interpolation) and `rDCD` is durables real
depreciation (FAAt103 line 15).

**Result** (`scripts/diag_gate_y_channel.py`):

| Branch | mean\|diff\| vs `bckm.Y_raw[:,0]` |
|---|---|
| FRED  | **0.0068** |
| BEA   | **0.0074** |
| gate  | 0.025 |

Both pass. BEA is +0.0006 worse than FRED — well within
floating-point and quarterization noise. **y-channel data
construction is not a leverage point for the LL gap.**

### x-channel — FAILS Stage-1 gate by 23×

**BCKM construction** (`usdata.m:53`):
```
X = rCD + rGPDI + rGI − (rCD/(rCND+rCS+rCD))·rSTX
```

**Result** (`scripts/diag_gate_x_channel.py`, bind-anchored
log-deviation form to cancel the `maketrend.m:15` y_pc(by) anchor
constant):

| Branch | mean\|diff\| vs `bckm.Y_raw[:,1]` |
|---|---|
| FRED  | **0.0065** |
| BEA   | **0.1535** |
| gate  | 0.025 |

BEA fails the gate by 23×. Component decomposition
(`scripts/diag_x_components.py`):

| Component | log-growth 1980→2008 (our BEA pull) | implied %/yr |
|---|---|---|
| rCD       | +1.55 | ~5.7%/yr |
| rGPDI     | +1.02 | ~3.6%/yr |
| rGI       | +0.77 | ~2.7%/yr |

Total BEA-X log-growth: **+1.06 nats**. BCKM-implied total: **+0.679
nats**. Our pull shows ~14pp of cumulative excess real growth versus
BCKM's 2014-vintage pull of the same NIPA series.

This is **not a construction bug** — all five components match
`usdata.m:30-38,53` formulas exactly (verified by
`scripts/diag_bea_fa_lines.py` against BCKM's pre-quarterized .dat
snapshots). The gap is BEA NIPA vintage drift (next section).

### g-channel — FAILS Stage-1 gate by 6×

**BCKM construction** (`usdata.m:34-35`, `usdata.m:56`):
```
g = (rEX − rIM) + rGC + worktemp.adjg·rGI   [chain-real]
```

**Result**:

| Branch | mean\|diff\| vs `bckm.Y_raw[:,3]` |
|---|---|
| FRED  | **0.0345** |
| BEA   | **0.2100** |
| gate  | 0.025 |

BEA fails the gate by 6×. Root cause: `maketrend.m:15` anchors
**every** observable at `y_pc(2008Q1)`, so the g-observable's level
depends on **both** g(2008Q1) **and** y(2008Q1). At 2008Q1:

| Pull | g/y |
|---|---|
| BEA chain-real (rEX−rIM, rGC, rGI all chain-real) | **0.118** |
| FRED single-deflator ((NX_nominal+GC_nominal+GI_nominal)/GDPDEF) | **0.097** |
| BCKM 2014 vintage | **0.102** |

The ~$300B gap between FRED's `NX_nominal/GDPDEF` and BEA's chain-real
`rEX − rIM` comes from terms-of-trade movement (2008 oil shock pushed
pIM > pEX, breaking the chain ≡ nominal÷GDPDEF identity). BEA is
faithful to BCKM's *formula*; FRED is closer to BCKM's *numerical
result* on this dataset.

**The gate is order-coupled.** Per-channel gating cannot work for g
alone while y stays on FRED — the level mismatch propagates through
the y_pc(2008Q1) anchor. We documented this in
`bca_core/data/pipeline.py` and **deferred** the g-channel gate to
after a joint y+x+g BEA migration test (next section).

---

## The joint walkdown

We added `--y-source`, `--x-source`, `--g-source` flags to
`scripts/diag_worktemp_compare.py` so the full Cartesian over BEA
toggles could be evaluated. The LL at BCKM-θ vs
`worktemp.mle.likelihood`:

| config (y, x, g) | LL ours | gap (ours − \|bckm\|) | x bias | g bias |
|---|---|---|---|---|
| (fred, fred, fred) — **default** | **+1719.66** | **−683**  | −0.023 | −0.035 |
| (bea,  fred, fred)               | +1683.68    | −719  (−36)   | +0.046 | +0.034 |
| (bea,  fred, bea)                | +1173.42    | −1229 (−547)  | +0.046 | +0.279 |
| (bea,  bea,  bea)                | −882.64     | −3286 (−2603) | −0.212 | +0.279 |

**Every BEA toggle makes the LL gap monotonically worse.** The fully
BEA-faithful configuration is **2603 nats further from BCKM** than the
FRED default.

Three mechanisms drive the monotone-worse pattern:

1. **`calgz` trend-slope coupling.** BCKM's `calgz.m` fsolves the
   trend slope `gz` against the MLE-window mean of detrended log y_pc.
   When y switches to BEA, y_pc changes, gz refits, and the new gz is
   then applied to x and g detrending — flipping their bias signs vs
   the BCKM target. The g and x channels carry no useful "BEA-faithful"
   signal once the trend is anchored on a different y.

2. **BEA NIPA vintage drift.** BEA chain-real series have been
   back-revised across the 2018+ comprehensive revisions, shifting
   historical levels 1–3% per series and cumulating to ~14pp over
   1980–2014 (see deep dive below). FRED's single-deflator approach
   happens to cancel most of this drift because numerator and
   denominator are revised in tandem.

3. **Bind-ratio dependency.** `maketrend.m:15` anchors every real
   series at `y_pc(2008Q1)`, so the level-gate against `bckm.Y_raw`
   depends on the y_pc(2008Q1) value — not just on the channel under
   test.

---

## BEA vintage drift: deep dive

The single most consequential finding from this forensics exercise.

### What is BEA vintage drift

BEA NIPA tables undergo two kinds of revision:

1. **Annual revisions** every July/August for the prior 3–4 years,
   incorporating new source data (tax returns, Census data, BEA
   surveys).
2. **Comprehensive revisions** roughly every 5 years (2013, **2018**,
   **2023**, ~2028), which:
   - Re-base the chain-real index (e.g., 2017=100 → 2022=100).
   - Re-classify component categories (e.g., consumer durables splits,
     fixed-asset reclassifications).
   - Update the underlying input-output tables, source-data
     methodology, and price-index methodology.
   - Apply silent retrospective revisions to historical levels going
     back 30+ years.

A single comprehensive revision can shift historical real levels by
1–3% per series. **Two comprehensive revisions (2018 + 2023) compound
to 14pp over a 28-year window** for sub-aggregates with high
methodology sensitivity (consumer durables, NX, gov investment).

### Why chain-real components drift more than headline aggregates

Chain-real series are non-additive by construction. The BEA chain
Fisher index for a sub-aggregate (e.g., rCD) is computed as a Fisher
average of Laspeyres and Paasche indices using **adjacent-period
prices**. When a comprehensive revision:

- Re-classifies what falls into "consumer durables" (e.g., software
  reclassified from durables to investment),
- Updates the price-index source for a component (e.g.,
  housing-services rent imputation methodology),
- Re-bases the chain index from year_old to year_new,

the historical Fisher chain-real series is **recomputed from scratch**
with the new classifications/prices/base. The result is that
`rCD[1980Q1]` in a 2014 BEA pull and `rCD[1980Q1]` in a 2026 BEA pull
are **different numbers** — even though both are labeled "real,
chained 2017 dollars" or similar.

For x-channel components, this means:

- **rCD** (consumer durables) had its definition revised in 2018 to
  reclassify software from goods to services and to update the
  hedonic adjustment for electronics.
- **rGPDI** (gross private domestic investment) had its
  intellectual-property-products component re-scoped multiple times.
- **rGI** (government investment) had its R&D capitalization treatment
  updated in 2013 and again in 2023.

Each of these revisions cumulates into ~5–6%/yr nominal-real growth in
our 2026 pull versus the ~3.5%/yr that BCKM's 2014 pull recorded.

### Why FRED's single-deflator approach is "accidentally robust"

FRED's approach to constructing real x is roughly:
```
real_x ≈ (PCE_durables_nominal + GPDI_nominal + Gov_invest_nominal) / GDPDEF
```

When BEA does a comprehensive revision:

- The **nominal** components (numerator) are revised — but nominal
  revisions are smaller because they don't depend on chain-Fisher
  recomputation; they're just $-value updates.
- **GDPDEF** (denominator) is revised in tandem — every comprehensive
  revision updates GDPDEF using the *same* methodology updates.
- The **ratio** is therefore *more stable* than either component
  alone, because methodology shifts in the numerator and denominator
  partially cancel.

This is why FRED-x matches BCKM's 2014 pull to mean|diff| = 0.0065
even though the *components* of FRED-x have all been revised since
2014. The cancellation isn't perfect, but it's much better than
recomputing each chain-real component from scratch.

### Why y was clean but x and g weren't

y is dominated by **rGDP**, the headline real-GDP series. rGDP is the
single series BEA puts the most effort into stabilizing across
revisions because it's the headline number. Sub-aggregate revisions
that average to zero at the rGDP level can still have large effects on
x (which is *not* dominated by any single stable headline) and g
(which is a sum of small chain-real components where rounding and
re-classification effects compound).

**Empirical check**: our BEA-y pull's log-growth 1980→2008 matches
BCKM-y to within 0.7pp; our BEA-x pull's log-growth 1980→2008 exceeds
BCKM-x by ~14pp.

### What we could do about it (and didn't)

Three options exist for matching BCKM's 2014 vintage exactly:

1. **ALFRED archive snapshot**: pull the as-of-2014 BEA NIPA tables
   from FRED's ALFRED (Archived FRED). This would give us the exact
   data BCKM saw. Rejected because (a) ALFRED only has FRED tickers,
   not BEA NIPA Table line items; (b) it would lock us to a specific
   vintage forever, defeating the goal of generalizability.
2. **Level-correction multipliers**: estimate a series-by-series
   multiplicative correction from a 2014-vintage anchor and apply it
   to current pulls. Rejected because it's brittle, hides the issue,
   and would produce wrong results for non-2014 periods.
3. **NIPA Underlying Detail Tables**: pull the more granular UDT
   series and reconstruct BCKM's component definitions from scratch.
   Rejected as not worth the engineering effort for a 14pp gap that
   doesn't change f-stat structure.

**The decision** (user-validated 2026-05-01): accept the residual
data-layer gap. Our FRED-default pipeline matches BCKM's f-stat
*structure* at BCKM-θ to ~0.01; the LL gap is a level offset, not a
structural disagreement.

---

## Three failure mechanisms (recap)

| # | Mechanism | Where it lives |
|---|---|---|
| 1 | **calgz trend coupling** — gz refits on y, then is applied to x and g | `matlab_reference/calgz.m`, `bca_core/data/pipeline.py` (calgz path) |
| 2 | **BEA chain-real vintage drift** — components revised 1–3%/series/comprehensive revision | BEA NIPA tables since 2018 |
| 3 | **Bind-ratio dependency** — y_pc(2008Q1) is the universal anchor | `matlab_reference/maketrend.m:15` |

These three are independent. Fixing any one of them in isolation
does not close the gap; you'd need an ALFRED archive pull (mechanism
2) **and** a joint y+x+g BEA migration with re-fit gz (mechanisms 1
and 3) to do it cleanly.

---

## Where the residual LL gap likely lives

After ruling out y-channel (clean), x-channel (vintage drift), and
g-channel (bind-ratio coupling), the most likely remaining
data-construction leverage point is:

- **labor / working-age-population**: variance ratio of `df["l"]` to
  `bckm.Y_raw[:,2]` is **1.31** — the only channel still above 1.2.
  We use OECD MEI 15–64 working-age population; BCKM uses BLS Civilian
  Population + DoD Armed Forces. The series have similar levels but
  different short-run noise structure, which could shift the
  Kalman-innovation distribution by tens of nats.

We have **not** migrated labor to a BLS-faithful pull because BLS
post-2014 dropped the BCKM-canonical universe and the available
substitutes (CPS-LF, CES-Total) carry their own measurement-system
discontinuities. This is logged as the highest-leverage remaining
data-layer item.

The non-data-layer remainder of the gap — the basin-escape problem
where our optimizer drifts away from BCKM's (Sbar, P, Q) once it
starts walking — is documented in `Diary.md` and is not a data
forensics issue.

---

## Implications for the web-app project

For the broader project (a web app that runs BCA on arbitrary
country/period data), the lessons from this forensics exercise are:

1. **Vintage discipline is a first-class concern.** Any reproducibility
   guarantee a user might want from the web app should explicitly
   record the data vintage (year + month) of each NIPA / OECD pull, so
   the app can either re-pull from the same vintage or warn the user
   that current data may differ from a reference computation.

2. **Single-deflator vs chain-real is a documented choice.** The web
   app's data layer should expose both options per channel, with
   default = single-deflator (FRED-style) for robustness to vintage
   drift, and BCKM-faithful chain-real as an opt-in for replication
   exercises against published papers.

3. **The Stage-1 gate is methodology, not just a one-off check.**
   `mean|diff| ≤ 0.025` against a published reference is a sensible
   gate for any country/period combination where we have a reference.
   For combinations without a reference (most of them), per-channel
   sanity checks (variance ratios, log-growth rates, terms-of-trade
   bias) replace the gate.

4. **Calibration-vs-data coupling matters.** BCKM's `calgz.m` shows
   that the trend slope is fit jointly with the data, not separately.
   Any web-app implementation must do this fsolve, not hardcode a
   per-country γ from a literature table.

5. **The bind-period anchor is methodology, not data.** It must be
   chosen per-period (peak of the cycle of interest) and propagated
   consistently through detrending, normalization, and counterfactual
   construction. A web-app UI should let the user pick `bind` and warn
   when the choice materially shifts results.

---

## What we did not try

For completeness, options considered and explicitly rejected:

- **ALFRED archive pull** to recover 2014-vintage BEA tables. Rejected
  because it locks us to a specific vintage and doesn't help future
  periods.
- **Level-correction multipliers** to anchor against BCKM's published
  numbers. Rejected as brittle and methodology-hiding.
- **Single-deflator approach for x using the GPDI deflator instead of
  GDPDEF.** Not tested; would shift the GPDI sub-component but leave
  rCD/rGI vintage drift untouched. Marginal benefit unclear.
- **Reconstructing labor from the BLS Productivity & Costs release.**
  Not tested; would require a separate data adapter and is gated on
  the user's higher-priority items (basin escape, web-app
  productionization).
- **Switching to a perfect-foresight (Ohanian-Raffo) approach** for
  wedge extraction. Explicitly out of scope per `bca-web-app-instructions.md`
  ("known pitfalls: don't use perfect-foresight").

---

## Diagnostic scripts (durable artifacts)

These scripts are kept in `scripts/` as regression and re-run
infrastructure for any future θ change or country/period extension:

| Script | Purpose |
|---|---|
| `diag_gate_y_channel.py` | Stage-1 gate for y |
| `diag_gate_x_channel.py` | Stage-1 gate for x (bind-anchored log-deviation form) |
| `diag_x_components.py` | Component decomposition of BEA-x (localizes vintage drift) |
| `diag_bea_fa_lines.py` | Recon utility: dumps BEA Fixed Assets line numbers |
| `diag_worktemp_compare.py` | End-to-end walkdown vs `worktemp.mat` (with `--y-source/--x-source/--g-source` flags) |

Each is independently runnable and self-documenting. None require
ALFRED, all use the BEA NIPA API + FRED for current pulls.

---

## Cross-references

- Method finding journal: [`CLAUDE.md`](CLAUDE.md), section
  "Findings (most recent first)".
- Per-session debugging journal: [`Diary.md`](Diary.md).
- Original project framing: [`bca-web-app-instructions.md`](bca-web-app-instructions.md).
- BCKM paper summary: [`BCA_info.md`](BCA_info.md), Section 7
  (US MLE estimates).
- Pipeline source flags: [`bca_core/data/pipeline.py`](bca_core/data/pipeline.py),
  `y_source`, `x_source`, `g_source` parameters.
