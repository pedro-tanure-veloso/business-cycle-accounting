# BCA Project — Claude Working Instructions

## Goal

Build a **usable Business Cycle Accounting model for the US economy** that lets us analyze the US business cycle through the lens of four structural wedges: efficiency (A), labor (1−τ_l), investment (1+τ_x), and government (g).

The concrete validation target is correctly replicating the **Great Recession decomposition** — specifically the results in **`BCA/BCA_info.md` Section 7 (United States — MLE Estimates)**:

- **Table 8** (P matrix): 4×4 VAR transition matrix, τ_l diagonal ≈ 1.001
- **Table 9** (P₀ vector): [0.0140, 0.0008, 0.0129, −0.0137]
- **Table 10** (V matrix): symmetric shock covariance
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

**Observables**: `[y_hat, l_hat, x_hat, g_hat]` as log-deviations from model SS (not sample means). This is how BCKM constructs observables — P₀ absorbs the offset between model SS and data mean.

**Labor normalization**: After detrending, rescale `df["l"]` so its sample mean equals the model's `l_ss`. The pipeline normalizes over the full fetch window (1947–2024), so post-hoc rescaling is required to make `mean(l_hat) ≈ 0`.

**Government share**: Set `g_share = mean(g_dt) / mean(y_dt)` from data (≈ 0.166 for US 1980–2014), not the model default of 0.20.

**Growth rates**: Use BCKM Table 2 (US row) values: `γ = 1.9%/yr`, `n = 0.98%/yr`. Do NOT use data-estimated growth rates for the structural calibration.

**Cross-country calibration constants** (BCKM `BCA_info.md` Table 1 + `matlab_reference/datamine.m`): `α = 1/3`, `ψ = 2.5`, `δ_annual = 0.05`, `β_annual = 0.975` (so `ρ_annual = 1/0.975 − 1 ≈ 0.02564`), `σ ≈ 1` (log utility). These are the values currently in `bca_core/params.py`. **Do not change them to other "BCKM" values** — a previous session moved them to `(0.35, 2.24, 0.0464, 0.9722)` thinking those were Table 77 / Table 7-7 values; they are not, and the resulting basin failed Phase B/D against `worktemp.mat`.

**Estimation**: Kalman-filter MLE (`mleqadj.m` logic), 30-parameter theta = [P₀(4), P(16), Q_lower_tri(10)].

**Kalman initialization / steady-state filter**: Use the DARE-derived steady-state covariance, evaluated at the *current* VAR parameters on every objective call (steady-state Kalman, BCKM `mleqadj.m` style). On a 5×5 system the per-call DARE is cheap (~ms) and eliminates the optimized-vs-final-smoother LL mismatch that arises when Σ₀ is frozen at BCKM Table 77 params while the optimizer drifts. The `_steady_state_kalman` helper returns a constant gain K, innovation cov S, and Σ_filt; the RTS smoother uses the same constants. Frozen-Σ₀ is no longer used.

**Stationarity constraints**: Spectral radius penalty only — `5e5 * max(|eig(P)| − 0.995, 0)²`, matching `mleqadj.m:134`. **Do not** add per-diagonal bounds. BCKM allows individual diagonals to exceed 0.995 (Table 8 has τ_l = 1.001); only the eigenvalue constraint is BCKM-faithful. If L-BFGS-B deadlocks at the penalty boundary, switch optimizer or clip the warm-start diagonal — do **not** reintroduce a per-diagonal penalty.

**Counterfactuals**: Keep full 4D VAR for rational expectations; zero out inactive wedge columns in structural equations. Capital evolves endogenously through `k' = P_k @ state`.

## Testing Standards

- Write unit tests for every non-trivial function in `bca_core/`
- Test files live in `tests/`
- Run `pytest tests/ -v` before reporting any result
- Critical tests to maintain:
  - Model solves correctly (Klein QZ), decision rules have right signs
  - `solve_counterfactual` with all-wedges active reproduces `_solve_with_var`
  - IRF for positive τ_x shock → investment and output fall (not rise)
  - Kalman filter log-likelihood is finite and well-behaved at BCKM parameters
  - `prepare_observables` gives `mean(l_hat) ≈ 0` after labor rescaling

## Do Not

- Use OLS VAR as the final estimator (warm-start only)
- Use Lyapunov or diffuse covariance for Kalman initialization
- Add per-diagonal penalties on P (BCKM imposes only the spectral-radius bound)
- Run a *transient* time-varying Kalman recursion with a frozen Σ₀ during optimization — that produces a 100+-unit LL gap vs the final smoother and biases counterfactuals. Use the steady-state Kalman with DARE-per-call instead.
- Use `std > 0.01` for random BCKM perturbations, or leave diagonal unclipped
- Skip the BCKM Table 77 warm-start
- Normalize observables by sample means (use model SS instead)

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
