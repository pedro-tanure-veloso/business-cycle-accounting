# BCA Project — Claude Working Instructions

## Goal

Build a **usable Business Cycle Accounting model for the US economy** that lets us analyze the US business cycle through the lens of four structural wedges: efficiency (A), labor (1−τ_l), investment (1+τ_x), and government (g).

The concrete validation target is correctly replicating the **Great Recession decomposition** — specifically the results in **`BCA/BCA_info.md` Section 7 (United States — MLE Estimates)**:

- **Table 8** (P matrix): 4×4 VAR transition matrix, τ_l diagonal ≈ 1.001
- **Table 9** (P₀ vector): [0.0140, 0.0008, 0.0129, −0.0137]
- **Table 10** (Q matrix): lower-triangular Cholesky factor of V (V = Q·Q′); BCKM displays it with mirrored upper triangle but only the 10 lower-tri elements are free parameters
- **Table 11** (f-statistics): fYA=0.16, fYτL=0.46, fYτx=0.32

Once these are matched, the model is trusted as correctly implemented and can be used for ongoing US economic analysis.

## Primary Reference Materials

1. **`BCA/BCA_info.md`** — paper summary with target parameter tables (Sections 1–7); primary spec
2. **`BCA/matlab_reference/`** — original Matlab replication files from the paper; ground truth for methodology
3. **`BCA/bca_paper.pdf`** — full paper PDF

**When implementing or changing anything, check the Matlab files first.** Key files in `matlab_reference/`:
- `mleqadj.m` — Kalman-filter MLE estimation
- `kfilter.m` — steady-state Kalman filter (DARE-based constant gain)
- `datamine.m` — calibration parameters and data construction
- `runall.m` / `fstats3.m` — counterfactuals and f-statistics

## BCKM Folder — Never Commit, Never Modify

**NEVER read, modify, commit, or push anything inside `BCA/BCKM/`.** This folder contains the original Matlab replication files from the paper. It is reference-only and is excluded from version control via `.gitignore`. If you need to consult the Matlab methodology, read the files but do not stage or alter them in any way.

## Methodology Rules

### Follow BCKM exactly
- If our approach differs from `mleqadj.m`, ask the user before proceeding
- Do not introduce approximations or shortcuts that are not in the paper
- Document any unavoidable deviations clearly

### Key methodological choices (already verified against paper)

**Observables**: `[y_hat, l_hat, x_hat, g_hat]` are time-detrended log levels. `prepare_observables` centers them against `obs_offset` so the prediction side sees zero at SS. BCKM `mleqadj.m:237` keeps observables uncentered and folds the SS into the obs-equation intercept `phi0 = Y0 − C(:,1:5)·X0(1:5)` — mathematically equivalent. (P₀ in BCKM is the VAR drift `(I−P)·Sbar`, not an SS-offset absorber; the SS lives in `phi0` on the model side.)

**Labor normalization**: Feed raw labor (`df["l"]`) directly — **do not rescale** to model `l_ss`. BCKM `mleqadj.m:237-238` uses `Y_l = log(hpc)` with no centering, and the SS-vs-data level gap is absorbed by `Sbar[2]` (τ_l) at MLE time. The earlier "rescale `df["l"]` so sample mean = `l_ss`" rule (and the matching `l_hat = log(l/l_ss)` centering in `prepare_observables`) was a bug: it injected a ~+0.20 phantom innovation per quarter at any Sbar where `l_new(Sbar) ≠ l_ss_calib`, which dominated the LL gap to BCKM at fixed θ. Fixed 2026-04-29 (Option A in the BCKM-θ diagnostic): `prepare_observables` now produces `l_hat = log(df["l"])` raw, and `obs_offset[1] = log(ss_new["l"])` cancels it on the prediction side.

**Government share**: Set `g_share = mean(g_dt) / mean(y_dt)` from data (≈ 0.115 for US 1980–2014). The paper does *not* publish a target g/y; this is data-derived under BCKM's convention `g = gov_consumption + net_exports` (BCA_info.md §3, §4 — closed-economy equivalence per CKM 2005). The earlier 0.166 figure in this doc was incorrect — it corresponded to gov_consumption alone, without the net-exports subtraction.

**Growth rates**: Use BCKM Table 2 (US row) values: `γ = 1.9%/yr`, `n = 0.98%/yr`. These are the published averages that BCKM's `calgz.m` fsolve (numerical `gz` such that detrended log-y has mean zero on the MLE window) and population-growth formula `gn = (iP_end/iP_1)^(1/(T-1))−1` produce on US 1980Q1–2014Q4. Hardcoding them is a shortcut around BCKM's per-call calibration step. Do not substitute HP-filtered, VAR-derived, or other ad-hoc growth-rate estimates.

**Cross-country calibration constants** (BCKM `BCA_info.md` Table 1 + `matlab_reference/datamine.m`): `α = 1/3`, `ψ = 2.5`, `δ_annual = 0.05`, `β_annual = 0.975` (so `ρ_annual = 1/0.975 − 1 ≈ 0.02564`), `σ ≈ 1` (log utility). These are the values currently in `bca_core/params.py`. **Do not change them to other "BCKM" values** — a previous session moved them to `(0.35, 2.24, 0.0464, 0.9722)` thinking those were Table 77 / Table 7-7 values; they are not, and the resulting basin failed Phase B/D against `worktemp.mat`.

**Estimation**: Kalman-filter MLE (`mleqadj.m` logic), 30-parameter theta = [P₀(4), P(16), Q_lower_tri(10)].

**P matrix convention — the paper transposes**. `BCA_info.md` §7 Table 8 is printed in a "rows = drivers, columns = receivers" *narrative* convention: row 0 reads "what z does" (its self-persistence in column 0, its outgoing spillover to τ_l in column 1, …). This is the **transpose** of the textbook convention `state_{t+1} = P · state_t` that BCKM's matlab code (`mleqadj.m:222`) and our codebase actually use. From 2026-04 through 2026-04-30 the same Table-8 matrix was hardcoded in **nine** independent places (var_estimation, eval scripts, tests) all in paper convention but used as code convention — a silent transpose at every BCKM-θ probe, the warm-start, and the counterfactual decomposition. Quantified cost at published θ: LL = +1195 (wrong) vs +1697 (correct), **a 501-nat gap from a single transposed matrix**. The fix: a single canonical module **`bca_core/constants.py`** exports `P_BCKM_TABLE8`, `SBAR_BCKM_TABLE8`, `QCHOL_BCKM_TABLE10` in CODE convention (verified element-wise against `octave_output/P_bckm.csv` to ≤4.3e-5). **Always import from there. Never re-transcribe Table 8.** If you need the paper's printed orientation, take the transpose.

**Kalman initialization / steady-state filter**: Use the DARE-derived steady-state covariance, evaluated at the *current* VAR parameters on every objective call (steady-state Kalman, BCKM `mleqadj.m` style). On a 5×5 system the per-call DARE is cheap (~ms) and eliminates the optimized-vs-final-smoother LL mismatch that arises when Σ₀ is frozen at BCKM Table 77 params while the optimizer drifts. The `_steady_state_kalman` helper returns a constant gain K, innovation cov S, and Σ_filt; the RTS smoother uses the same constants. Frozen-Σ₀ is no longer used.

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

Details on why these are wrong are in `REPORT.md` → "Things Not To Do".

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
