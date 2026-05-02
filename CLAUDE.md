# BCA Project — Claude Working Instructions

## Goal

Build a **usable Business Cycle Accounting toolkit** that produces credible
four-wedge decompositions — efficiency (A), labor (1−τ_l), investment
(1+τ_x), government (g) — for **arbitrary country/period combinations**.
The end deliverable is a three-layer web app (`bca_core` model +
`bca_api` service + `bca_web` UI) where users pick a country and time
window and get back a wedge decomposition + counterfactuals.

### What "credible" means and how we validate it

**Layer 1 — BCKM 2016 US 1980Q1–2014Q4 regression test.** The pipeline
must reproduce the paper's headline results to a documented tolerance:

- **Table 8** (P matrix): 4×4 VAR transition matrix, τ_l diagonal ≈ 1.001
- **Table 9** (P₀ vector): [0.0140, 0.0008, 0.0129, −0.0137]
- **Table 10** (Q matrix): lower-triangular Cholesky factor of V (V = Q·Q′)
- **Table 11** (f-statistics): fYA=0.16, fYτL=0.46, fYτx=0.32

**Status (2026-05-01):** at BCKM-θ on our df, f-stats match Table 11
to ≤0.01 in every channel. The pipeline structure is correct (verified
end-to-end, see Findings). The ~150-nat residual LL gap is concentrated
in BEA NIPA vintage drift (chain-real series back-revised post-2018 to
levels that differ ~14pp cumulatively from BCKM's 2014 vintage); every
BEA-vs-FRED toggle on our window makes the gap **worse** (joint walkdown,
2026-05-01 finding). FRED defaults are the production setting.

This is the regression target — kept passing forever — not the live
target. The 1980-2014 cached parquet
(`bckm_replication/data/us_1980_2014_calgz.parquet`) and the
BCKM-window tests (`tests/test_bckm_table12.py`, parts of
`tests/test_bckm_reference.py`) pin this fidelity. Final replication
wrap-up: `bckm_replication/REPORT.md`.

**Layer 2 — narrative-prior smoke tests on other US windows.** Validate
generalizability by running the pipeline on windows where we have strong
qualitative priors (e.g. COVID 2010Q1–2023Q4 with bind=2019Q4: τ_l
collapses 2020Q2; A spikes 2021 with record TFP; g elevated through
ARPA). A "pass" is sign-direction agreement on the wedge decomposition
under both full-window and pre-anomaly-fit detrending. This is the
*current* live target (post-2026-05-01 pivot).

**Layer 3 — cross-country (future).** OECD MEI data layer for the
24-country Tables III/IV of the paper. Out of scope until US
generalizability lands.

### Scope of BCKM-replication-specific rules

Several rules in this document — most prominently the
"Follow BCKM exactly" methodology rule, the BCKM-anchored constants
(`labor_target_mean=0.24279`, the BEA-branch
`start_year=1980, end_year=2014` hardcodes in
`bca_core/data/bea.py:195-273`), and the data-source priority
("BEA/BLS first, FRED only as fallback") — apply to the **Layer 1
regression** only. They are deliberate fidelity choices that make our
1980-2014 output BCKM-faithful.

For Layer 2 and Layer 3, the relevant defaults are:

- **`labor_target_mean=None`** (raw hours-per-capita on the sample
  window; BCKM-anchored 0.24279 is opt-in for the regression path).
- **FRED defaults** for `y_source/x_source/g_source` (BEA branches are
  diagnostic-only opt-in; their 1980-2014 hardcodes are out-of-scope
  to remove until a separate migration sprint).
- **Data-derived γ and n** (calgz fsolve and working-age-pop slope —
  generic, no change needed).
- **`mle_window`** parameter to fit the calgz slope on a sub-window
  (e.g. pre-COVID) when the full sample contains a structural anomaly
  that would distort the trend.

When a Layer-2 or Layer-3 deliverable would benefit from violating a
BCKM-replication rule (e.g. using an OECD-MEI population series instead
of FRED `LFWA64TTUSQ647N`), the rule yields. The Layer-1 regression
test is the gate: as long as the cached BCKM parquet still produces the
pinned f-stats, the rule's purpose has been served.

## Primary Reference Materials

1. **`BCA/BCA_info.md`** — paper summary with target parameter tables (Sections 1–7); primary spec
2. **`BCA/bckm_replication/matlab_reference/`** — original Matlab replication files from the paper; ground truth for methodology
3. **`BCA/bca_paper.pdf`** — full paper PDF
4. **`BCA/bckm_replication/REPORT.md`** — final BCKM 1980-2014 replication wrap-up

**When implementing or changing anything, check the Matlab files first.** Key files in `bckm_replication/matlab_reference/`:
- `mleqadj.m` — Kalman-filter MLE estimation
- `kfilter.m` — steady-state Kalman filter (DARE-based constant gain)
- `datamine.m` — calibration parameters and data construction
- `runall.m` / `fstats3.m` — counterfactuals and f-statistics

## BCKM Folder — Never Commit, Never Modify

**NEVER read, modify, commit, or push anything inside `BCA/BCKM/`.** This folder contains the original Matlab replication files from the paper. It is reference-only and is excluded from version control via `.gitignore`. If you need to consult the Matlab methodology, read the files but do not stage or alter them in any way.

## Loops — stop them when you report back

When the user invokes `/loop` and you self-pace via `ScheduleWakeup`,
**stop the loop the moment you report back to the user**, unless something
specific is still running in the background that genuinely needs the next
tick (a long-running fit, a CI build, a file watch). Concretely: at the end
of any turn that delivers a substantive update to the user, **omit the
`ScheduleWakeup` call** so no further wakeups fire.

Why: the user is not always at the keyboard. A wakeup that fires while
they're afk runs autonomously, makes decisions, and surprises them when
they come back. That's bad. Loops are for "watch this background task";
they are not a license for autonomous iteration during idle gaps.

Operational rule: if you have nothing materially new to do without user
input — e.g. you're awaiting approval on a fix, or no observable event has
happened — **the loop ends now** and the user re-arms it later if they want
another tick.

If you scheduled a wakeup before realizing the user was afk and the wakeup
fires while still pending: when it fires, recognize the conversation has
returned to interactive use and DO NOT take further actions or
reschedule — just confirm the loop is dead.

## Debugging Findings — append-only journal

**Working rule for Claude (you):** every time a debugging session establishes
a non-trivial fact about this codebase or about BCKM — something that
falsified a hypothesis, ruled out a suspect, fixed a sign/anchor/convention
bug, or proved two quantities are (or are NOT) the same — **append it to this
section before the session ends**. One bullet per finding. Format: short
verdict, where it was verified (script + line / table-mat field), and what it
rules out. These accumulate; do not delete or rewrite earlier entries unless
a later finding directly invalidates them (in which case mark the old one
`SUPERSEDED` and add the replacement below).

The point: every iteration we lose facts to context compaction. This section
is the durable memory.

### Findings (most recent first)

- **2026-05-01 — Joint BEA migration walkdown: every BEA toggle makes
  the LL gap WORSE. FRED defaults are the best operating point.
  The whole "BEA over FRED" priority is the wrong leverage point for
  this dataset and time window.** `scripts/diag_worktemp_compare.py`
  now takes `--y-source / --x-source / --g-source` flags. Walking the
  full Cartesian gives the LL gap at BCKM-θ vs `worktemp.mle.likelihood`:

  | config (y, x, g) | LL ours | gap (ours − \|bckm\|) | x bias | g bias |
  |---|---|---|---|---|
  | (fred, fred, fred) baseline | **+1719.66** | **−683** | −0.023 | −0.035 |
  | (bea,  fred, fred)          | +1683.68    | −719  (−36)  | +0.046 | +0.034 |
  | (bea,  fred, bea)           | +1173.42    | −1229 (−547) | +0.046 | +0.279 |
  | (bea,  bea,  bea)           | −882.64     | −3286 (−2603)| −0.212 | +0.279 |

  Three mechanisms explain the monotone-worse pattern:
  1. **calgz trend coupling**: `gz` is fsolved against y_pc on the MLE
     window. Switching y to BEA changes y_pc, which re-fits gz, which
     is then applied to x and g detrending — flipping their bias signs
     vs the BCKM target. The g and x channels carry no useful "BEA-
     faithful" signal once the trend is anchored on a different y.
  2. **BEA vintage drift**: chain-real BEA series have been back-revised
     across the 2018+ comprehensive revisions, shifting historical
     levels 1-3% per series and cumulating to ~14pp over 1980-2014
     (see x-channel finding below). FRED's single-deflator approach
     happens to cancel most of this drift.
  3. **Bind-ratio dependency**: BCKM's `maketrend.m:15` anchors every
     real series at `ypc(by)`, so the level-gate against `bckm.Y_raw`
     depends on the y_pc(2008Q1) value — not just on the channel under
     test. (Per the g-channel order-coupling finding below.)

  **Decision (user-validated 2026-05-01)**: keep `{y,x,g}_source="fred"`
  as the production default; preserve all three BEA branches as opt-in
  ablation infrastructure (well-documented, gate-deferred). The BEA
  migration is **complete** in the sense that BCKM-faithful constructions
  are available for diagnostic A/B testing, but **not used** because
  every toggle reduces fidelity to BCKM's `worktemp.mat` ground truth.
  Implication for closing the residual data-layer LL gap: the leverage
  is NOT in BEA-vs-FRED data sources. Most likely candidates remain:
  labor/working-age-pop (var ratio 1.31, only channel still above 1.2),
  detrending sub-conventions, and the off-by-one obs differencing
  reflected in the 16.6-nat objective-only gap from 2026-04-30. The
  basin-escape problem (Sbar drift away from BCKM's once optimization
  starts) is the real bottleneck for f-stat replication.

- **2026-05-01 — x-channel BEA migration FAILS Stage-1 gate due to BEA
  vintage drift; FRED-x is already an exceptionally tight match to BCKM
  growth and is the right default. Don't migrate x to BEA.** Step 3 of
  the migration plan: `x_source="bea"` flag plumbed in
  (`bca_core/data/pipeline.py`); BEA branch reuses
  `fetch_real_components` per `usdata.m:53`:
  ``X = rCD + rGPDI + rGI − (rCD/(rCND+rCS+rCD))·rSTX``. Stage-1 gate
  (`scripts/diag_gate_x_channel.py`, bind-anchored log-deviation form
  to cancel detrending-anchor constants):
  FRED mean|diff|=**0.0065**, BEA mean|diff|=**0.1535** — BEA is
  **23× worse than FRED**. Component decomposition
  (`scripts/diag_x_components.py`): our BEA rCD log-grew +1.55 nats
  1980→2008 (≈+5.7%/yr nominal-real), rGPDI +1.02 nats (+3.6%/yr),
  rGI +0.77 nats (+2.7%/yr). BCKM-implied total X growth over the
  same window is +0.679 nats; ours is +1.06 nats — a ~14pp
  cumulative excess. **Not a construction bug.** All five components
  match BCKM `usdata.m:30-38,53` formulas exactly (recon
  `scripts/diag_bea_fa_lines.py` + `bca_core/data/bea.py:354-357`
  comments verify the modern T10105/T10106/T10109/T30904/T30905 line
  mappings against BCKM's pre-quarterized .dat snapshots). The gap is
  **BEA NIPA vintage drift**: BCKM ran 2014; the chain-real series
  have been rebased and back-revised through 2018, 2023+
  comprehensive revisions, shifting historical levels by 1-3% per
  series and cumulating to ~14pp over 28 years. FRED-x's
  single-deflator approach (`(GPDI + pce_durables + gov_inv)/GDPDEF`)
  happens to cancel this drift cleanly because GDPDEF rebases and the
  nominal numerator's historical revisions partially offset.
  **Resolution**: keep `x_source="fred"` as default — already the
  best gate result of any channel. Commit `x_source="bea"` as
  diagnostic-only opt-in (parallel to g). The original Stage-1
  walkdown's −0.023 x bias (`scripts/diag_worktemp_compare.py`) is
  almost entirely a level offset, NOT a growth disagreement —
  bind-anchored growth gap is just 0.0065 (FRED) and falls inside
  numerical precision of BCKM's `maketrend.m:15` ypc(by) anchor.
  79/79 fast tests still pass with FRED default. Implication for the
  data-layer LL gap: **x is not a leverage point** — the residual gap
  must live in g (still on FRED) and l (variance-ratio 1.31, also
  FRED).

- **2026-05-01 — y-channel BEA migration PASSES Stage-1 gate; y is not
  where the data-construction gap lives. x and g are.** Step 4 of the
  BEA NIPA migration is complete: `BeaDataFetcher.fetch_durables_components`
  pulls FAAt101 line 15 (consumer-durables current-cost net stock,
  year-end annual) and FAAt103 line 15 (current-cost depreciation, SAAR),
  quarterizes log-linearly for the stock and constant-within-year for
  the flow, and divides by the quarterly pCD deflator. Pipeline exposes
  `y_source="fred"` (default, legacy) and `y_source="bea"` (BCKM-faithful
  per `usdata.m:51`: `Y = rGDP − rSTX + 0.04·rKCD + rDCD`). Stage-1 gate
  (`scripts/diag_gate_y_channel.py`) against `bckm.Y_raw[:,0]`:
  FRED mean|diff|=**0.0068**, BEA mean|diff|=**0.0074** — both well
  under the 0.025 gate. BEA is +0.0006 worse than FRED — basically a
  wash; the y-channel data-construction differences between FRED and
  BEA are second-order at the level we're measuring against BCKM.
  79/79 fast tests still pass with FRED default. **Implication for
  the data-layer gap**: the worktemp.mat walkdown
  (`scripts/diag_worktemp_compare.py`, 2026-05-01) shows Stage-1
  per-channel mean(ours − bckm) = (y: +0.000, l: −0.000, **x: −0.023,
  g: −0.035**). y carries essentially zero bias — the ~666-nat data-layer
  LL gap (per the 2026-04-30 LL-decomposition finding) is concentrated
  in x and g. **The next leverage point is the x-channel migration**
  (`usdata.m:53`: `X = rCD + rGPDI + rGI − (rCD/rCNDS)·rSTX`), then a
  joint y+x+g gate.

- **2026-05-01 — BEA NIPA migration is order-coupled at the bind ratio:
  per-channel level gate against `bckm.Y_raw` cannot pass for g alone
  while y stays on FRED.** Stage-1 element-wise comparison after Step 2
  (g→BEA): BEA mean|diff|=0.2100 vs FRED 0.0345 vs gate 0.025 → **FAIL**.
  Root cause traced to `matlab_reference/maketrend.m`: BCKM defines
  `mled[:,k] = var_pc/ypc(bind)*(1+gzt)^by`, so `Y_raw[:,3] = log(g/y at
  bind)` after detrending — the **level** of the g-observable depends on
  both g and y constructions at 2008Q1. At 2008Q1: BEA chain-real g/y =
  0.118 vs FRED nominal-deflated g/y = 0.097 vs BCKM 0.102. The ~$300B
  gap between FRED's `NX_nominal/GDPDEF` and BEA's chain-real `rEX−rIM`
  comes from terms-of-trade movement (2008 oil shock pushed pIM>pEX,
  breaking the chain ≡ nominal÷GDPDEF identity). **The BEA construction
  is BCKM-faithful per `usdata.m:34-35,56`; the gate failure is a
  consistency artifact, not a bug.**
  **Resolution (option d, user-approved 2026-05-01)**: keep
  `g_source="fred"` as default so the f-stat baseline doesn't move;
  expose `g_source="bea"` for the post-y joint-gate test. The
  per-channel gate is **DEFERRED** to after y migrates (Step 4); only
  then is the y-trend normalizer BCKM-faithful and the level-comparison
  meaningful. Documented in `bca_core/data/pipeline.py` `g_source`
  docstring. Pipeline behavior unchanged; 79/79 fast tests still pass.

- **2026-04-30 — Optimizer investigation: BCKM-θ is essentially a critical
  point of our LL on BCKM data; the +11-nat residual gap decomposes
  cleanly into 6 nats of routine L-BFGS-B convergence + 5 nats of
  shrinkage-induced basin drift. No structural bug.**
  `scripts/diag_optimizer_basin_v2.py` runs four checks at BCKM-θ on
  BCKM `Y_raw`:

  | Stage | LL | max\|eig(P)\| | ‖Sbar−Sbar_BCKM‖∞ | ‖P−P_BCKM‖∞ | ‖Q−Q_BCKM‖∞ |
  |---|---|---|---|---|---|
  | BCKM-θ (ref)              | +1887.73 | 0.9960 | — | — | — |
  | After penalized L-BFGS-B  | +1893.75 | 0.9952 | 7.6e-3 | 1.3e-3 | 1.6e-3 |
  | Full optimizer (cached)   | +1899.03 | 0.9950 | 5.6e-2 | 2.7e-2 | 3.0e-2 |

  • **Q1 (penalized walk)**: 23 L-BFGS-B iters from BCKM-θ → +6.02 nats.
    All four θ components drift by ≤ 1.6e-3. This is gradient noise / BCKM
    stopping early — not a structural disagreement.
  • **Q2 (shrinkage gain)**: extra +5.28 nats from multi-restart + 50-iter
    multiplicative-shrink loop. This walks to a meaningfully different
    basin (Sbar drift jumps 7x, P drift 20x, Q drift 19x).
  • **Q3 (Q sign flips)**: all 16 sign basins give bit-identical LL
    (+1893.7474) and bit-identical penalty. Our walked Q is in the
    `++++` basin — same as BCKM (gap 1.58e-3 in our basin, vs ≥ 1.5e-2
    for any sign-flipped basin). **Sign flip is not the issue.**
  • **Q4 (per-quarter innovations)**: per-channel RMSE differences are
    ≤ 3e-4 across all 4 obs channels. Largest single-quarter improvement
    (1980Q1) is +0.037 nats — driven by our t=0 boundary (innov[0] =
    obs[0] − log(ss) is non-zero by construction; BCKM uses Y(2)−Y0
    after differencing, dropping this point). Remaining 139 quarters
    contribute < 0.005 each.

  **Implications:**
  1. The optimizer is structurally correct. BCKM-θ is ~6 nats below a
     nearby penalized local max — typical optimizer-stopping noise.
  2. The shrinkage loop (BCKM-faithful per `runmleadj.m`) finds a
     different ~5-nat-better basin where weight shifts from τ_l to τ_x
     in f-stats. This is at the limit of what's distinguishable from
     numerical fidelity, given finite-difference Jacobians and the off-
     by-one obs differencing convention.
  3. **Hunting structural bugs in the optimizer is closed.** No deeper
     fix lives there. The remaining ~150-nat data-construction gap on
     our df (FRED-vs-BEA-NIPA) is the higher-leverage open issue.

  Open methodological question for the user: should we disable the
  multiplicative-shrink loop (`pb=0.99, nps=50` in `var_estimation.py:
  815-828`)? It's BCKM-faithful but pushes us 5 nats / 13pp of f-stat
  weight away from BCKM's published θ on BCKM's own data. Default:
  keep it (BCKM uses it; the LL gain is real).

- **2026-04-30 — phi0 fix landed: `obs_offset_kf = obs_offset_wedge`
  reproduces BCKM Table 11 to ~0.01 with our Sbar warm-start.**
  Patch in `var_estimation.py:476-484` replaces the static
  `obs_offset_kf = obs_hat[0, :]` with the Sbar-dependent
  `obs_offset_wedge = log(ss_new[var])`. All 79 fast tests still pass.
  Pre-vs-post pivot table on our df:

  | scenario | LL_pre | LL_post | fstats_post (A,τ_l,τ_x,g) |
  |---|---|---|---|
  | BCKM-θ on our df | +1645.97 | **+1719.66** | 0.151, 0.485, 0.306, 0.058 |
  | BCKM(P,Q) + OUR Sbar_init | +1248.62 | **+1716.54** | **0.154, 0.469, 0.316, 0.060** |
  | Our converged MLE | +1826.47 | +1825.64 | 0.128, 0.630, 0.181, 0.061 |

  Table 11 target: (0.16, 0.46, 0.32, 0.06). The middle row matches to
  **0.01** in every channel — i.e., when our Sbar-fsolve seed is paired
  with BCKM's published (P, Q), our pipeline reproduces Table 11 almost
  exactly. **The phi0 hypothesis is confirmed and resolved.**

  On BCKM `Y_raw`:
  - Sbar `log_g` moved from −1.218 (pre) to **−1.991 (post)**, matching
    BCKM's −1.935 to 0.06. The data-independent attractor at
    log_g ≈ −1.2 is **gone**.
  - ‖ΔSbar‖∞ collapsed from 0.72 to **0.056**.
  - LL gap (ours converged − BCKM-θ-on-BCKM-data) dropped from +27 to
    **+11 nats**.

  Residual problems after the fix:
  1. Our optimizer still finds a (P, Q) basin ≈ +106 nats above BCKM-θ
     on our df (was +180); the landscape walk now shows a barrier
     (LL crashes to −1.7M around α=0.55–0.65, likely DARE failure in
     interpolated θ — the path through θ-space is pathological even
     though both endpoints are well-conditioned).
  2. f-stats at our MLE-converged θ on our df are now (0.128, 0.630,
     0.181, 0.061) — substantial improvement from (0.050, 0.877,
     0.053, 0.021), but still not Table 11. Reshuffles in fY[τ_l] vs
     fY[τ_x] are now driven by P/Q drift, not Sbar.
  3. ~150-nat data-construction gap to BCKM `Y_raw` remains (BEA NIPA
     todo, still BLOCKED).

  Implication: the dominant Track A bug is fixed. Remaining work splits
  into (a) understanding why our optimizer escapes BCKM's (P, Q) basin
  on top of our (P, Q)-correct warm-start, and (b) the BEA NIPA data
  switch.

- **2026-04-30 — Basin gap is ~83% data-construction, ~17% objective.**
  Two diagnostics together pin the cause:
  1. `scripts/diag_ll_landscape.py` walks θ_bckm → θ_ours along the
     convex line in θ-space, scoring our LL at 21 grid points.
     **LL is monotone-increasing across all 20 steps**: from +1645.97
     (BCKM-θ) to +1826.47 (our θ), Δ=+180 nats, no barrier. f-stats
     degrade smoothly along the path: (0.151, 0.485, 0.306) → (0.050,
     0.877, 0.053). Verdict: BCKM-θ is *not* a local max of our LL
     surface — the optimizer correctly walks uphill on a smooth slope
     to a different argmax. So we are **NOT in a "basin escape"
     situation**; we are on a different objective function.
  2. `scripts/diag_mle_on_bckm_data.py` runs our optimizer on BCKM's
     `Y_raw` (permuted [0,2,1,3]) instead of our `df`. Result: LL
     converges to +1898.64, beating BCKM's stored θ by **+27 nats on
     BCKM's own data** (vs +180 nats on our df). f-stats jump from our
     own-data converged (0.05, 0.88, 0.05, 0.02) to (0.115, 0.614,
     0.235, 0.037) — clear directional move toward Table 11.
  Decomposition: 180 − 27 = **~150 nats of the LL gap is the data
  source** (FRED-vs-BEA-NIPA construction); **~27 nats is residual
  objective-formulation difference** independent of data. Track B
  (data) is the higher-leverage fix.

- **2026-04-30 — Basin disagreement is concentrated in Sbar, especially
  `log_g`. P and Q converge to BCKM.** When fed BCKM's `Y_raw`, our
  optimizer lands at ‖ΔP‖∞ = 0.025 and ‖ΔQ‖∞ = 0.007 vs BCKM Tables
  8/10 — close. But ‖ΔSbar‖∞ = **0.72**, dominated by log_g: ours
  −1.22 (g/y ≈ 0.30) vs BCKM −1.94 (g/y ≈ 0.14, matching data). On our
  `df` the converged log_g is also ≈ −1.21 — **so log_g ≈ −1.2 is a
  data-independent attractor of our objective**, not a data-fit
  artifact. Concrete suspect: `var_estimation.py:476` sets
  `obs_offset_kf = obs_hat[0, :]` (constant in θ) where BCKM
  `mleqadj.m:232` sets `phi0 = Y(:,1) − C·X0(1:5)` with X0 = SS state
  in absolute coords — making BCKM's intercept Sbar-dependent through
  C·X_ss. In Path A our state is mean-zero by construction so we set
  C·X_ss = 0 and drop that term. If BCKM's X0 is non-zero in their
  coords, our LL loses an Sbar-dependent term that pulls Sbar toward
  data fit. **Highest-leverage Track A candidate** — verify next by
  reading `mleqadj.m:225-235` and checking whether `Y0 = X0(1:5)` in
  BCKM is zero or not at the linearization.

- **2026-04-30 — At full BCKM-θ, our f-stats match Table 11 to ~0.01.
  This SUPERSEDES the "structure is broken" framing.** Running
  `scripts/eval_bckm_fstats.py` with `(SBAR_BCKM_TABLE8, P_BCKM_TABLE8,
  QCHOL_BCKM_TABLE10)` from `bca_core/constants.py` produces
  fY[A]=**0.1513**, fY[τ_l]=**0.4853**, fY[τ_x]=**0.3057**, fY[g]=0.0577
  vs Table 11 target (0.16, 0.46, 0.32, —). Max cell gap ≈ 0.025.
  LL at BCKM-θ = +1645.97 (our convention; BCKM-form ≈ +1131). Our MLE
  converges to LL=+1826 with very different f-stats (fY[A]=0.05, fY[τ_l]=0.88,
  fY[τ_x]=0.05) — i.e. the optimizer finds a higher-LL basin that disagrees
  with BCKM's f-stat decomposition. **Diagnosis: pipeline structure is
  correct at published θ; remaining problem is OPTIMIZER BASIN choice.**
  The earlier diagnostic framing ("LL=+1233 vs BCKM +2706 ≈ 1500-nat gap")
  conflated our LL convention with BCKM's `mle.likelihood` sign convention
  (see the `mleqadj.m:257` finding below). Cached at `/tmp/fstats_clean.txt`
  (2026-04-30 16:26). Implication: **stop hunting structural bugs at BCKM-θ;
  the structure is fine. Hunt the basin instead** — multi-start strategy,
  Q parameterization, or stationarity-penalty interaction with the optimizer.

- **2026-04-30 — Open paradox in `bckm.components["mlx"]` regression:
  +0.082 on k is structurally impossible from our LOM, but R²=1.0.**
  OLS β = [+0.082, +0.286, −1.718, +0.348, +0.109] (rank=5, cond=14.9,
  R²=1.000000). Our LOM gives [0.000, +0.368, −1.550, +0.457, +0.136].
  Our and BCKM's algorithms both produce Δgammak=0 by construction, so
  β[k]=+0.082 cannot come from either pipeline's `bckm_capital_lom`.
  Best remaining hypotheses: (a) BCKM stored `mlx` was generated with a
  different Y0 anchor than `gwedges2.m:21` claims, (b) BCKM applied an
  additional row operation we haven't found, or (c) a numerical artifact
  in the matlab CF loop that doesn't survive analytical scrutiny. Not
  blocking — given the f-stat finding above, the labor→x cell mismatch
  contributes <0.025 to the f-stat gap, well within the noise budget.

- **2026-04-30 — `params.py:adj_cost_elasticity` has a wrong comment, but the
  value is not the labor→x bug.** The comment claims "Calibrated to match
  BCKM adja=12.88 (mleqadj.m, datamine.m)". BCKM's actual formula
  (`bca_steady2.m:58`) is `adja = 0.25/(n + γ + δ)`, which gives ≈ 12.56 at
  US growth rates, not 12.88. Setting `adj_cost_elasticity=0.25` would
  reproduce BCKM's `a` exactly. Sweep over `{0.250, 0.256448, 0.260, 0.270,
  0.280}` (a ∈ {12.56, 12.88, 13.06, 13.56, 14.06}) shows max|gap| in
  `(C_lab − C_0)[x, :]` is monotone in `a` but stays ≥ 0.15 at all values
  and never produces P_x[k] ≠ 0. The hypothesis "the labor→x gap is
  miscalibrated `a`" is **falsified**. Comment is still wrong and should be
  fixed. (`scripts/diag_capital_lom_labor.py` sweep block.)

- **2026-04-30 — `worktemp.adjc` is a SWITCH, not a coefficient value.**
  `datamine.m:54: worktemp.adjc = 2; %1 for no, 2 for BGG, 3 for 4*BGG`.
  Briefly mistook `adjc=2` as `a=2` (vs our 12.88) — that was wrong; the
  switch selects which formula `adja` is computed from. Always check the
  comment in `datamine.m` before treating an `adjc/adjk/...` integer as a
  parameter value.

- **2026-04-30 — Our `bckm_capital_lom` matches `fixexpadj.m` line-for-line.**
  Verified: `res_adjust` (bca_core/bckm_lom.py) ≡ `res_adjust2.m`
  (Z indexing, As-pinning at lines 46-53, Newton iter on labor FOC, residual
  formula); `bckm_capital_lom` (lines 198-266) ≡ `fixexpadj.m:50-74` (Gamma
  quadratic + linear solve); `bckm_C_matrix:339` `C[1,:] =
  [phixk,0,0,0,0] + phixkp*Gamma` ≡ `fixexpadj.m:106`. `a0, a1, a2`
  partials are bit-identical across As ∈ {[0,1,0,0], [0,0,0,0], [1,1,1,1]}
  → Δgammak = 0 by construction. Our LOM is structurally fine.

- **2026-04-30 — Linearized state path matches BCKM's linearized state path
  to machine precision.** `scripts/diag_state_path_compare.py`: log_k
  matches exactly; log_g eps; log_z 3.65e-3; τ_l 1.11e-2; τ_x 2.18e-3.
  Cross-check: applying BCKM's nonlinear formula `(1−Tault)/(1−Tault(Y0))`
  to our extracted state reproduces `bckm.wedges["tault"]` to **4.44e-16**,
  and `(Zt/Zt(Y0))^(1−θ)` reproduces `bckm.wedges["zt"]` to 1.84e-4. The
  τ_l "1.1pp gap" is **pure linearization curvature**, not a state bug.
  Conclusion: **CF math should use the linearized state, not `bckm.wedges`**;
  state extraction is not the bug.

- **2026-04-30 — `bckm.wedges` (from `worktemp.mat`) are NONLINEAR;
  `Xt0` in `gwedges2.m:80` is LINEARIZED. Don't mix them.** Specifically
  `gwedges2.m:70-77`: lowercase `lzt, tault, lkt, lgt` (LIN, used in
  `Xt0 = [lkt, lzt, tault, tauxt, lgt, ones]` line 80, fed to `(C2−C0)*Xt0`)
  vs capital `Zt, Tault` (NL, used only in published `w.zt`,
  `w.tault` at lines 196-197). Direct subtraction of `bckm.wedges` from
  our linearized states is apples-to-oranges. Use `Y_raw` + Cobb-Douglas
  inversion to reconstruct the linearized series for any state-level
  comparison.

- **2026-04-30 — OLS β recovered in `diag_labor_x_cell.py` Layer 3 IS BCKM's
  true `(C_lab − C_0)[x, :]`.** OLS of `bckm.components["mlx"]` (=
  `exp(YMl(:,2) − YMl(Y0,2))*100`, x-component of bind-normalized labor CF;
  per `gwedges2.m:171,205`) on bind-centered linearized states is rank-5
  full-rank (cond=14.9, R²=1.000000, residual std 8e-6). β =
  [+0.082, +0.286, −1.718, +0.348, +0.109] in (k, z, τ_l, τ_x, g) order.
  Our derived `phixkp · (G_lab − G_zero)` =
  [0.000, +0.368, −1.550, +0.457, +0.136]. The +0.082 entry on k is
  **structurally impossible** for any LOM with Δgammak=0 — but our
  algorithm and the `fixexpadj.m` algorithm both yield Δgammak=0. Open
  paradox; suspect anchor/Y0 convention or a stored-vs-recomputed BCKM
  components mismatch.

- **2026-04-29 — Counterfactual decomposition is INCREMENTAL `(C_j − C_0)`,
  not single-wedge `C_j`.** Already in main "Methodology Rules" section
  above; logged here too because it was the dominant CF bug for several
  sessions. Without the `−C_0` subtraction, the four single-wedge CFs
  over-count the data drop (−10.92% sum vs data −7.48%) and the inverse-SSR
  f-stat decomposition collapses onto whichever wedge tracks the no-wedge
  baseline. Verified end-to-end in `tests/test_bckm_table12.py`.

- **2026-04-29 — Y0 anchor is `bind` (2008Q1), NOT sample start.** Already
  in main section; logged here because it caused multiple silent bugs
  (`f_statistics_bckm` with anchor=0 collapsed pre-2008 drift into SSR;
  Table 12 "peak" was being read as local argmax 2007Q4 instead of 2008Q1).

- **2026-04-29 — Single canonical source for BCKM Table 8 / 9 / 10.**
  `bca_core/constants.py` exports `P_BCKM_TABLE8`, `SBAR_BCKM_TABLE8`,
  `QCHOL_BCKM_TABLE10` in code convention (verified element-wise to
  octave_output to ≤4.3e-5). The paper prints P transposed; using paper
  orientation as code orientation costs **501 nats of LL** at published θ.
  Never re-transcribe Table 8 — always import from `constants.py`.

- **2026-04-29 — Labor normalization: feed raw labor, do NOT rescale to
  `l_ss`.** Already in main section; logged here because it was the
  dominant LL-gap source at BCKM θ before the 2026-04-29 fix.

- **2026-04-29 — `mle.likelihood` in `worktemp.mat` is `L`, not `−L`.**
  BCKM `mleqadj.m:257` minimizes `L = 0.5(T·log|Ω| + tr(Ω⁻¹Σ_innov)) +
  penalty`, no `n_obs·log(2π)`. Stored value is signed; at converged θ it's
  negative (≈ −2402.88) because `T·log|Ω|` dominates. Our printed `ll_ours`
  adds `0.5·T·n·log(2π) ≈ 514.6`. Always state the convention when quoting
  an LL gap. Verified end-to-end in `scripts/diag_bckm_data_isolation.py`:
  same-data, same-θ gap is **16.6 nats**, not 757.

## Methodology Rules

### Follow BCKM exactly (Layer 1 regression scope)

This rule applies to the **BCKM 1980Q1–2014Q4 US regression test**
(see "Goal → Scope of BCKM-replication-specific rules"). Outside that
window the rule yields to whatever is correct for the active country/
period — model structure stays BCKM-faithful, but data-construction
and labor-anchor choices are window-appropriate.

- If our approach differs from `mleqadj.m` on the regression path, ask
  the user before proceeding
- Do not introduce approximations or shortcuts that are not in the paper
- Document any unavoidable deviations clearly

### Data sources — BEA / BLS first, FRED only as fallback (Layer 1 scope)

**McGrattan's working rule (passed down from advisor): always pull data from
the original source — BEA NIPA, BLS Productivity & Costs, BLS CPS — and use
FRED only when the original source has no equivalent.** BCKM's `usdata.m` does
exactly this: every series in their construction is sourced from BEA NIPA
tables or BLS series-level files, not aggregated FRED tickers.

**Apparent OECD-vs-NIPA contradiction with the paper text — resolved.**
`BCA_info.md` §4 says "data come mainly from OECD" — true, but specifically
about the 24-country **cross-country** exercise (paper's Tables III/IV).
For the **US-specific MLE block** (Tables 8/9/10/11, our validation target),
McGrattan used NIPA + BLS instead. Three pieces of evidence pin this:
1. `matlab_reference/usdata.m` loads `nipa115/116/119/32/33/394/395.dat`
   (BEA NIPA Tables 1.1.5/1.1.6/1.1.9/3.2/3.3/3.9.4/3.9.5), `atab10d.dat`
   and `btab100d.dat` (BEA Fixed Asset tables), `hours.dat` (BLS), and
   `civpop.dat`/`armed.dat` (BLS Current Population Survey + DoD). No
   OECD codes appear anywhere.
2. `usdata.m:101` writes `worktemp.mat` directly. There is no
   `oecddata.m` or any other US-data construction script in
   `matlab_reference/`.
3. `worktemp.mat`'s `(Sbar, P, Q)` plugged through our pipeline produces
   f-stats matching Table 11 to ~0.01 in every channel — so `usdata.m`
   IS what produced the paper's US headline tables.

`BCA_info.md` §4 line 175 acknowledges this implicitly: "while the U.S.
NIPA accounts have quarterly data on consumer durable expenditures for
the 1980:1–2014:4 sample we use, the OECD has more limited data."

Implication: do **not** redo the US data layer in OECD — it would
diverge from `worktemp.mat:Y_raw` element-wise. If we ever extend BCA
to other countries (Tables III/IV), we'll need a separate OECD pipeline,
but that's out of scope until the US Table 11 replication is closed.

Why this matters for replication:
- FRED republishes BEA/BLS series but sometimes truncates them when BEA
  rebases the chain index (e.g. `PCDGCC96`, `PCNDGC96`, `PCESVC96`,
  `A782RX1Q020SBEA` only have data from 2007Q1+ on FRED after the 2017
  chain rebase, even though BEA's NIPA API serves the full 1947+ history).
- FRED aggregates can paper over construction details that matter at the
  margin (e.g. PCE chain-type price index vs component-specific deflators;
  state-only vs federal+state+local sales tax).
- "Faithful to BCKM" means the same NIPA line item, not a FRED ticker that
  happens to have a similar name.

Concrete order of preference for each US series we need:
1. **BEA NIPA API** for any nominal/real/deflator NIPA quantity (national
   accounts, gov consumption/investment, PCE breakdown, sales tax aggregates).
   `bca_core/data/bea.py` already implements `BeaDataFetcher` with disk cache.
2. **BLS series-level files** for hours, employment, productivity (BLS
   `Productivity & Costs` releases, `CES` employment, `CPS` labor force).
3. **OECD MEI** for working-age population (the BCKM-canonical universe is
   gone from BLS post-2014; OECD MEI 15-64 is the cleanest current proxy).
4. **FRED** only if 1–3 don't carry the series, AND only after verifying the
   FRED ticker has full history (not truncated at a chain-rebase boundary).

When a FRED-default series is being used as a stand-in for a BEA/BLS source,
note the original source explicitly in the FRED dictionary docstring so the
fallback is auditable.

### Key methodological choices (already verified against paper)

**Observables**: `[y_hat, l_hat, x_hat, g_hat]` are time-detrended log levels. `prepare_observables` centers them against `obs_offset` so the prediction side sees zero at SS. BCKM `mleqadj.m:237` keeps observables uncentered and folds the SS into the obs-equation intercept `phi0 = Y0 − C(:,1:5)·X0(1:5)` — mathematically equivalent. (P₀ in BCKM is the VAR drift `(I−P)·Sbar`, not an SS-offset absorber; the SS lives in `phi0` on the model side.)

**Labor normalization**: Feed raw labor (`df["l"]`) directly — **do not rescale** to model `l_ss`. BCKM `mleqadj.m:237-238` uses `Y_l = log(hpc)` with no centering, and the SS-vs-data level gap is absorbed by `Sbar[2]` (τ_l) at MLE time. The earlier "rescale `df["l"]` so sample mean = `l_ss`" rule (and the matching `l_hat = log(l/l_ss)` centering in `prepare_observables`) was a bug: it injected a ~+0.20 phantom innovation per quarter at any Sbar where `l_new(Sbar) ≠ l_ss_calib`, which dominated the LL gap to BCKM at fixed θ. Fixed 2026-04-29 (Option A in the BCKM-θ diagnostic): `prepare_observables` now produces `l_hat = log(df["l"])` raw, and `obs_offset[1] = log(ss_new["l"])` cancels it on the prediction side.

**Government share**: Set `g_share = mean(g_dt) / mean(y_dt)` from data (≈ 0.115 for US 1980–2014). The paper does *not* publish a target g/y; this is data-derived under BCKM's convention `g = gov_consumption + net_exports` (BCA_info.md §3, §4 — closed-economy equivalence per CKM 2005). The earlier 0.166 figure in this doc was incorrect — it corresponded to gov_consumption alone, without the net-exports subtraction.

**Growth rates**: Use BCKM Table 2 (US row) values: `γ = 1.9%/yr`, `n = 0.98%/yr`. These are the published averages that BCKM's `calgz.m` fsolve (numerical `gz` such that detrended log-y has mean zero on the MLE window) and population-growth formula `gn = (iP_end/iP_1)^(1/(T-1))−1` produce on US 1980Q1–2014Q4. Hardcoding them is a shortcut around BCKM's per-call calibration step. Do not substitute HP-filtered, VAR-derived, or other ad-hoc growth-rate estimates.

**Cross-country calibration constants** (BCKM `BCA_info.md` Table 1 + `matlab_reference/datamine.m`): `α = 1/3`, `ψ = 2.5`, `δ_annual = 0.05`, `β_annual = 0.975` (so `ρ_annual = 1/0.975 − 1 ≈ 0.02564`), `σ ≈ 1` (log utility). These are the values currently in `bca_core/params.py`. **Do not change them to other "BCKM" values** — a previous session moved them to `(0.35, 2.24, 0.0464, 0.9722)` thinking those were Table 77 / Table 7-7 values; they are not, and the resulting basin failed Phase B/D against `worktemp.mat`.

**Estimation**: Kalman-filter MLE (`mleqadj.m` logic), 30-parameter theta = [P₀(4), P(16), Q_lower_tri(10)].

**P matrix convention — the paper transposes**. `BCA_info.md` §7 Table 8 is printed in a "rows = drivers, columns = receivers" *narrative* convention: row 0 reads "what z does" (its self-persistence in column 0, its outgoing spillover to τ_l in column 1, …). This is the **transpose** of the textbook convention `state_{t+1} = P · state_t` that BCKM's matlab code (`mleqadj.m:222`) and our codebase actually use. From 2026-04 through 2026-04-30 the same Table-8 matrix was hardcoded in **nine** independent places (var_estimation, eval scripts, tests) all in paper convention but used as code convention — a silent transpose at every BCKM-θ probe, the warm-start, and the counterfactual decomposition. Quantified cost at published θ: LL = +1195 (wrong) vs +1697 (correct), **a 501-nat gap from a single transposed matrix**. The fix: a single canonical module **`bca_core/constants.py`** exports `P_BCKM_TABLE8`, `SBAR_BCKM_TABLE8`, `QCHOL_BCKM_TABLE10` in CODE convention (verified element-wise against `octave_output/P_bckm.csv` to ≤4.3e-5). **Always import from there. Never re-transcribe Table 8.** If you need the paper's printed orientation, take the transpose.

**Kalman initialization / steady-state filter**: Use the DARE-derived steady-state covariance, evaluated at the *current* VAR parameters on every objective call (steady-state Kalman, BCKM `mleqadj.m` style). On a 5×5 system the per-call DARE is cheap (~ms) and eliminates the optimized-vs-final-smoother LL mismatch that arises when Σ₀ is frozen at BCKM Table 77 params while the optimizer drifts. The `_steady_state_kalman` helper returns a constant gain K, innovation cov S, and Σ_filt; the RTS smoother uses the same constants. Frozen-Σ₀ is no longer used.

**LL formula and the BCKM `mle.likelihood` sign convention**. BCKM `mleqadj.m:257` is

```
L = 0.5·(T·log|Ω| + tr(Ω⁻¹·Σ_innov)) + penalty
```

where `Ω` is the steady-state innovation covariance and `Σ_innov = Σ_t innov_t innov_t'`. **There is no `n_obs·log(2π)` constant** — `mleqadj.m` is a minimizer of `L`, not a probability density. Our `_kf_ll` in `var_estimation.py:546` computes

```
ll_ours = -0.5·(T·n_obs·log(2π) + T·log|Ω| + quad)
        = -L − 0.5·T·n_obs·log(2π)
```

so we add a `0.5·T·n_obs·log(2π) ≈ 514.6` nat normalization (T=140, n_obs=4) that BCKM omits. **This is harmless for parameter estimation — the constant doesn't change argmax — but it makes printed LL values incomparable.** When stating an LL gap, always specify the convention: `"BCKM form"` (positive-minimization L) vs `"our form"` (log-density LL).

**Critical sign trap**: `worktemp.mat`'s `mle.Likelihood` field stores `L` directly — it is **not** `−L`. The stored value `−2402.88` is *negative* because at the converged θ the `T·log|Ω| ≈ −5398` term dominates the positive `quad ≈ 538` term; `L = 0.5·(−5398 + 538) ≈ −2430`. So the conversions are:

```
L_bckm                   = -2402.88   (stored mle.Likelihood, may be negative)
L_ours_at_BCKM_θ_BCKMdata = -2386.28   (our mleqadj.m:257-form L on BCKM Y_raw)
ll_ours_at_BCKM_θ_BCKMdata = -L − 514.6 = +1871.67   (printed by ``estimate_var_mle``)
```

Verified by `scripts/diag_bckm_data_isolation.py`: at identical (Sbar, P, Q, data) the gap between our `L` and BCKM's stored `mle.Likelihood` is **16.6 nats**, not 757. Earlier confusion was from comparing `ll_ours` (with 2π) directly to `|mle.Likelihood|` (BCKM form), conflating two different L definitions and the sign convention.

**Stationarity constraints**: Spectral radius penalty only — `5e5 * max(|eig(P)| − 0.995, 0)²`, matching `mleqadj.m:134`. **Do not** add per-diagonal bounds. BCKM allows individual diagonals to exceed 0.995 (Table 8 has τ_l = 1.001); only the eigenvalue constraint is BCKM-faithful. If L-BFGS-B deadlocks at the penalty boundary, switch optimizer or clip the warm-start diagonal — do **not** reintroduce a per-diagonal penalty.

**Counterfactuals**: Keep full 4D VAR for rational expectations; zero out inactive wedge columns in structural equations. Capital evolves endogenously through `k' = P_k @ state`.

**Counterfactual decomposition is INCREMENTAL, not single-wedge** (BCKM `gwedges2.m:90-115`). A "single-wedge CF" is the **incremental contribution** of activating that wedge relative to the no-wedge baseline `C0` (As=[0,0,0,0]). Concretely, BCKM constructs `YMz = (Xt0 − Xt0[Y0])(C1 − C0)' + YM0[Y0]`, **not** `YMz = (Xt0 − Xt0[Y0])C1' + YM0[Y0]`. The `−C0` subtraction is what makes the per-wedge CFs **additive** (sum of four single-wedge CFs ≈ all-active CF ≈ data) and what makes Table 12 work: peak-trough drop = −7.0% with components −1.9, −3.4, −4.5, +2.7 summing to −7.1. Without `−C0`, each single-wedge CF carries the no-wedge baseline plus its wedge effect, so the four CFs over-count the data drop (we hit −10.92% sum vs data −7.48% on 2026-04-29) and the inverse-SSR f-stat decomposition collapses onto whichever wedge happens to track the baseline closest. **The BCKM additive identity is the ground truth that pins this — verified end-to-end against `worktemp.mat` per-quarter CF paths in `tests/test_bckm_table12.py`.**

**The Y0 anchor is `bind` (= 2008Q1 in our 1980-2014 dataset), NOT the sample start (1980Q1)**. BCKM `gwedges2.m:21` sets `Y0 = worktemp.bind` and every level-ratio (`w.yt`, `w.mzy`, etc.) is anchored at this base period. `fstats3.m` then computes f-stats on slices of these already-Y0-anchored series. Using anchor=0 (sample start) in `f_statistics_bckm` collapses the cumulative pre-2008 drift into the SSR and is wrong. The `peak` in BCKM Table 12 is also `bind` (2008Q1 = 1.0 by construction), not the actual local argmax 2007Q4.

## Testing Standards

- Write unit tests for every non-trivial function in `bca_core/`
- Test files live in `tests/`
- Run `pytest tests/ -v` before reporting any result
- Critical tests to maintain:
  - Model solves correctly (Klein QZ), decision rules have right signs
  - `solve_counterfactual` with all-wedges active reproduces `_solve_with_var`
  - IRF for positive τ_x shock → investment and output fall (not rise)
  - Kalman filter log-likelihood is finite and well-behaved at BCKM parameters
  - `prepare_observables` returns zero-centered obs when fed SS-constant data; `phi0` carries raw SS offsets (`phi0[1] = log(l_ss)` per Option A — no labor rescaling, BCKM `mleqadj.m:237` `Y_l = log(hpc)` convention)

## Do Not

- Use OLS VAR as the final estimator (warm-start only)
- Use Lyapunov or diffuse covariance for Kalman initialization
- Add per-diagonal penalties on P (BCKM imposes only the spectral-radius bound)
- Run a *transient* time-varying Kalman recursion with a frozen Σ₀ during optimization — that produces a 100+-unit LL gap vs the final smoother and biases counterfactuals. Use the steady-state Kalman with DARE-per-call instead.
- Use `std > 0.01` for random BCKM perturbations, or leave diagonal unclipped
- Skip the BCKM Table 77 warm-start
- Rescale `df["l"]` to model `l_ss` before `prepare_observables` (see Labor normalization above — that rescale was the dominant LL-gap source at BCKM θ; raw labor is what BCKM feeds)

Details on why these are wrong are in `bckm_replication/REPORT.md` → "Things not to do".

## Asking Before Deviating

If you are considering any of the following, **stop and ask the user first**:
- Changing the state-space formulation (F, H matrices)
- Changing the penalty structure or stationarity bounds
- Switching from DARE to any other initialization method
- Adding or removing parameters from the theta vector
- Changing the data series used (e.g., investment definition)
- Changing calibration constants (α, β, δ, g_share formula)

## Progress Tracking

When you reach ~90% of your context window, write your progress to `Diary.md` with:
- Date/session identifier
- What was completed since the last entry
- What is in progress
- What is blocked and why
- Exact next step (file, function, what to change)

This ensures continuity across context resets.

### Layer 2 diary: `covid_analysis/Diary.md`

**`covid_analysis/Diary.md`** is the session journal for all Layer 2 work
(COVID-era smoke test and any subsequent non-BCKM window work). Append a
dated entry at the end of every session that makes meaningful progress on
this work. The entry should record: what was done, what is in progress,
what is blocked, and the exact next step.

**Commit `covid_analysis/Diary.md` whenever you append to it.** Treat it
as durable memory that persists across context resets — the same way
`CLAUDE.md`'s Findings section works for the BCKM replication. Do not let
a session end without a Diary entry and a commit if anything material
happened.
