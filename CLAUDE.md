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
2. **`BCA/BCKM/Multicountry - End/USAN2/`** — original Matlab replication files from the paper; ground truth for methodology
3. **`BCA/bca_paper.pdf`** — full paper PDF

**When implementing or changing anything, check the Matlab files first.** Key files:
- `mleqadj.m` — Kalman-filter MLE estimation
- `protmod.m` / `solve*.m` — model and linearization
- `counterfactual*.m` — wedge counterfactuals

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

**Growth rates**: Use BCKM Table 77 values: `γ = 1.9%/yr`, `n = 0.98%/yr`. Do NOT use data-estimated growth rates for the structural calibration.

**Estimation**: Kalman-filter MLE (`mleqadj.m` logic), 30-parameter theta = [P₀(4), P(16), Q_lower_tri(10)].

**Kalman initialization**: Use DARE (not Lyapunov, not diffuse) with BCKM Table 77 parameters as prior; freeze during optimization for speed, recompute at final parameters for smoother.

**Stationarity constraints**: Spectral radius penalty + per-wedge diagonal bounds `[0.995, 1.005, 0.995, 0.995]` for [A, τ_l, τ_x, g]. The τ_l bound is relaxed to 1.005 because BCKM's estimated τ_l diagonal is ≈ 1.001.

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
- Use spectral-radius penalty without per-diagonal penalty
- Solve DARE inside the optimization loop (only at final parameters)
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
