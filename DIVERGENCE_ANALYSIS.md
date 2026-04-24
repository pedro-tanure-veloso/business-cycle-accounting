# BCKM Replication — Sources of Divergence from the Paper

## Context

Our current MLE run (`scripts/run_var_counterfactuals.py` on FRED data, 1980Q1–2014Q4) produces:

| Quantity    | Ours                                        | BCKM Section 7 target         |
|-------------|---------------------------------------------|-------------------------------|
| ll          | 1806                                        | —                             |
| P diag      | [0.89, 0.99, 0.96, 0.99]                    | [0.989, 1.001, 0.967, 0.995]  |
| P₀          | [−0.006, 0.020, −0.005, −0.003]             | [0.014, 0.001, 0.013, −0.014] |
| φᵧ          | A 0.23, L 0.28, X 0.25, G 0.23              | A 0.16, L 0.46, X 0.32, G ~0  |
| PtT x (GR)  | CF-X −1.05 (actual −0.35)                   | CF-X similar range            |

The τₓ wedge correctly worsens during the Great Recession, and the peak-to-trough investment response has the right sign and dominant contributor. So the model is qualitatively correct. But the VAR itself, the intercept vector P₀, and φ-statistics diverge materially.

Below is a head-to-head comparison of our Python implementation vs. the Matlab files in `BCA/BCKM/Multicountry - End/USAN2/`, ranked by expected impact on results.

---

## Tier 1 — Likely high impact

### 1.1 Data construction (usdata.m / datamine.m)

BCKM constructs the observables as follows:

```matlab
% usdata.m — adjustments to NIPA
Y = rGDP - rSTX + 0.04*rKCD + rDCD              % + durables service flow + durables cons
C = rCND + rCS - share_cnd*rSTX + 0.04*rKCD + rDCD
X = rCD + rGPDI + rGI - (rCD/rCNDS)*rSTX        % durables + GPDI + gov investment
G = rGC + rEX - rIM                              % gov CONSUMPTION + net exports
```

```matlab
% datamine.m — trend removal
gz = 1.018^(1/4) - 1                             % CALIBRATED
mled = [t, ypc, xpc, hpc, gpc]                   % per-capita levels
Y = log(mled(:,2:5)) - log([(1+gz)^t, (1+gz)^t, 1, (1+gz)^t])
% observables are log-levels minus log((1+gz)^t) growth trend
% hours (col 3) is NOT detrended
```

Our implementation ([bca_core/data/pipeline.py](bca_core/data/pipeline.py), [bca_core/data/adjustments.py](bca_core/data/adjustments.py)):

- Uses raw FRED `GDP` without subtracting sales tax
- Does not add imputed durables service flow (`0.04·rKCD` — 4% annual return on the stock of consumer durables)
- Government wedge = `GCE + NETEXP` (close to BCKM), but BCKM splits `GC` from `GI` and puts `GI` in investment
- Investment = `GPDI + durables` but omits government investment `GI`
- Growth rate `γ` is **estimated** by OLS on log(y); BCKM uses **calibrated** γ = 1.01934^(1/4) = 0.00478/qtr
- Labor normalizer: we use `target_mean=0.25` (later rescaled to `l_ss`); BCKM uses `hpc/1300` (hours per person per quarter / 1300), a fixed divisor tied to the definition of "fraction of time endowed worked"

**Impact on estimates**: MLE fits the VAR to whatever observables you feed it. Different Y/X/G definitions produce different shock magnitudes and different P₀. This is almost certainly the largest single source of the φ-statistic gap.

**What to do**: match the NIPA adjustments in `adjustments.py`. The imputed durables service flow and the GC/GI split are the two we are missing. `rKCD` (consumer-durables stock) is not on FRED; it comes from Flow of Funds `btab101d`. For a first pass we can approximate using `PCDG` cumulated with BEA-style depreciation (ask before implementing).

### 1.2 Sample period

- BCKM MLE: 1948Q1–2014Q3 (~267 quarters) in `datamine.m`; 1969Q1–2014Q2 (~182 quarters) in `runall.m` for the wedge decomposition
- Ours: 1980Q1–2014Q4 (140 quarters)

A sample shorter by 40–50% attenuates evidence on rare events (e.g., 1980s disinflation, early-80s recession) that helps pin down the efficiency and labor wedge covariance. **Extending to 1948Q1 start should be straightforward** in our FRED pipeline — most series go back that far. Ask before changing since it changes all figures.

### 1.3 Observables are log-levels with a fixed gz trend, not log-deviations from model SS

In `mleqadj.m`:

```matlab
Y = log(ZVAR) - log([(1+gz)^t, (1+gz)^t, 1, (1+gz)^t])
Ybar = Y(2:T,:) - Y(1:T-1,:)*D'                  % = Y(2:T,:) since D=0
```

The observables are log-levels (of per-capita, base-period-normalized, sales-tax-adjusted series) with the deterministic growth trend subtracted. The filter then matches them to the model via a **constant offset**:

```matlab
phi0 = Y0 - C*X0(1:5)                            % Y0 = [log(ys); log(xs); log(ls); log(gs)]
C    = [C, phi0]                                 % C is 4x6
```

`phi0` is a deterministic constant in the observation equation that absorbs the gap between the model's steady state and the data mean.

Our implementation ([bca_core/var_estimation.py](bca_core/var_estimation.py)):

- `prepare_observables` returns log-deviations from model SS (`log(l/l_ss)`, `log(x / x_ss_over_y_ss)` etc.)
- We post-hoc rescale `df["l"]` so `mean(l_hat) ≈ 0`
- No `phi0` column in our observation matrix

**Consequence**: BCKM's P₀ absorbs the SS-vs-sample-mean drift inside the filter; ours has that drift pre-removed by rescaling, so our estimated P₀ is measuring something different (and smaller in magnitude). Our P₀ ≈ [−0.006, 0.020, −0.005, −0.003] vs BCKM's [0.014, 0.001, 0.013, −0.014] is not surprising under this setup.

### 1.4 Parametrization: Sbar vs. P₀ directly

BCKM optimizes over `Sbar` (the unconditional mean of log wedges) and derives `P₀ = (I − P)·Sbar` ([mleqadj.m:29](BCKM/Multicountry%20-%20End/USAN2/mleqadj.m), [runmleadj.m:14](BCKM/Multicountry%20-%20End/USAN2/runmleadj.m)):

```matlab
Sbar    = [log(1); .05; .0; log(.07)]            % defaults
P0      = (eye(4)-P)*Sbar
% Theta order: Sbar(1:4), P(:), Q(lower tri)
```

And critically, `initmle.m` uses `fsolve` to pick Sbar such that model SS matches data sample means:

```matlab
x0 = [ys-Ym(1); xs/ys-Ym(2); ls-Ym(3); gs/ys-Ym(4)]
% fsolve zeros these residuals → initial Sbar
```

So the **model steady state depends on Sbar** and is re-solved each MLE iterate. Our implementation keeps the model SS fixed at calibration values and optimizes P₀ directly. This means:

- BCKM's optimization landscape has a smooth, well-conditioned "match sample means" direction via Sbar.
- Ours has P₀ multiplied by `(I − P)⁻¹` to get Sbar, which for near-unit-root P blows up small changes in P₀ into huge changes in Sbar. That can produce degenerate local minima where P₀ is small but unconditional means are implausible.

**What to do**: reparametrize `theta = [Sbar(4), vec(P)(16), vec(Q_tri)(10)]` and recover `P₀ = (I − P)·Sbar` inside `_unpack`. Keep everything else the same. This is a modest refactor of [var_estimation.py](bca_core/var_estimation.py).

### 1.5 Calibration constants differ

| Parameter | Ours (params.py)         | BCKM (datamine.m)                |
|-----------|--------------------------|----------------------------------|
| θ / α (capital share) | 1/3 ≈ 0.3333 | 0.35                             |
| ψ (leisure weight)    | 2.5          | 2.24                             |
| δ annual              | 0.05         | 1 − (1 − 0.0464)^(1/4) = 0.0464 |
| β (quarterly)         | (1/1.025)^0.25 ≈ 0.9938 | 0.9722^(1/4) ≈ 0.9930     |
| σ (risk aversion)     | (not exposed) | 1.000001 (≈ log utility)        |

Small differences individually, but they compound in the steady state. For example ψ = 2.24 vs 2.5 changes `l_ss` by several percent, which changes the scale of `l_hat` and feeds into P₀.

**What to do**: add `theta` (α), `psi`, `delta_annual`, `rho_annual`, `sigma` as constants in `CalibrationParams` with BCKM defaults, and use them as the default when doing the replication exercise.

---

## Tier 2 — Likely medium impact

### 2.1 Kalman filter form: steady-state gain vs. time-varying

BCKM uses the **asymptotic / steady-state** Kalman likelihood ([mleqadj.m:244–258](BCKM/Multicountry%20-%20End/USAN2/mleqadj.m), [kfilter.m](BCKM/Multicountry%20-%20End/USAN2/kfilter.m)):

```matlab
[K, Sigma]  = kfilter(A, Cbar, B*B', Rbar, B*B'*C')      % DARE-derived gain, covariance
Omega       = Rbar + Cbar*Sigma*Cbar'                     % steady-state innovation cov
for i=2:T
  Xt(i,:)   = Xt(i-1,:)*A' + innov(i-1,:)*K'              % constant gain K throughout
  innov(i,:) = Ybar(i,:) - Xt(i,:)*Cbar'
end
sum1 = innov'*innov/T
L = 0.5 * (T*(log(det(Omega)) + trace(Omegai*sum1)) + penalty)
```

This is the asymptotic form:

```
L = (T/2) · [log|Ω| + tr(Ω⁻¹ · S_innov)] + penalty
```

assuming the filter has converged to the steady-state gain immediately. `Cbar = C·A` and `Rbar = R + C·B·B'·C'` shift by one period so that `Ybar[t] = Cbar·X[t] + noise` uses the predicted state.

Ours ([var_estimation.py](bca_core/var_estimation.py), `_kf_ll`) uses the standard time-varying Kalman filter with the exact log-likelihood `−(T/2)·k·log(2π) − (1/2)·Σₜ [log|Ωₜ| + innovₜ'·Ωₜ⁻¹·innovₜ]`.

Both are correct. In the limit T → ∞ and stationary process they agree. Differences:

- Ours is more sensitive to the `Sigma_0` prior. BCKM effectively ignores transient filter gain — the steady-state gain is reached at t=1 by construction.
- The gradient landscape is slightly different; optimizers may find different local minima.

**What to do**: probably not worth changing right now. The asymptotic form is a performance optimization. If we're seeing convergence issues it's worth trying.

### 2.2 Initial state X₀: SS vs. DARE prior

BCKM sets `X₀ = [log(k_s); log(z_s); τ_l,s; τ_x,s; log(g_s); 1]` — the exact model SS ([mleqadj.m:161](BCKM/Multicountry%20-%20End/USAN2/mleqadj.m)). No uncertainty.

Ours uses `x = zeros(5)` as the initial mean with `Sigma_0 = dare(F, H, Q_proc + εI, R_obs)` as the covariance.

For T=140 this probably matters little for the MLE but could matter for the smoother (early-period smoothed values). Worth checking.

### 2.3 Penalty structure

BCKM ([mleqadj.m:135](BCKM/Multicountry%20-%20End/USAN2/mleqadj.m)):

```matlab
penalty = 500000 * max(max(abs(eig(P))) - .995, 0)^2    % spectral radius only
```

Ours ([var_estimation.py](bca_core/var_estimation.py)):

```python
eig_max = max(abs(eig(P)))
penalty = 5e5 * max(eig_max - 1.005, 0)^2               # spectral radius
DIAG_BOUNDS = [0.995, 1.005, 0.995, 0.995]              # per-wedge
diag_excess = max(|diag(P)| - DIAG_BOUNDS, 0)
penalty += 5e5 * sum(diag_excess^2)                     # per-diagonal
```

We added the per-wedge diagonal penalty after observing τₓ,τₓ = 1.016 in a run without it. But BCKM tolerates τ_l,τ_l = 1.001 (Table 8), which our bounds explicitly allow at 1.005, so that's fine. The question is whether our tighter 0.995 cap on A, τ_x, g prevents BCKM's solution from being found.

Looking at BCKM's Table 8 (target): max diag = 1.001 (τ_l). Our per-diag at 1.005 is above that. So our penalty should not block BCKM's solution. But the per-diag penalty makes the landscape more restrictive than BCKM's.

**What to do**: worth trying a run with BCKM's exact penalty (spectral only, bound 0.995) to see if we find a different optimum. Ask before changing since it's a deliberate choice to avoid a specific failure mode (see [REPORT.md](REPORT.md) "Things Not To Do").

### 2.4 Adjustment cost parameter

BCKM uses `a = 12.88` (quarterly, from `runmleadj.m:113`, `adjc=2`, "BGG"). We calibrate `a = adj_cost_elasticity / b` where `b = (1+γ)(1+n) − (1−δ)` and `adj_cost_elasticity = 0.25`. For our calibration:

```
b ≈ 0.00478 + 0.00244 + 0.01220 ≈ 0.01942                         (investment/capital SS)
a = 0.25 / 0.01942 ≈ 12.87                                         (matches BCKM 12.88 ✓)
```

**OK, this matches.** Good.

---

## Tier 3 — Minor / confirmed not-the-issue

### 3.1 Counterfactual construction (`fixexpadj.m` vs `solve_counterfactual`)

BCKM ([fixexpadj.m](BCKM/Multicountry%20-%20End/USAN2/fixexpadj.m)) keeps the full P in the transition matrix (agents form expectations over ALL wedges), and zeros the STATIC-observation coefficients for inactive wedges via `As` switches:

```matlab
C = [ [phiyk, phiyz*As(1), phiyl*As(2), 0, phiyg*As(4)] + phiykp*Gamma(1:5)'; ... ]
```

The decomposition is then linear:
```matlab
YMz  = (Xt - Xt_base) * (C1 - C0)' + YM0_base           % isolates z's static contribution
YMn  = (Xt - Xt_base) * C0' + YM0_base                  % "no-wedge" (pure dynamic) baseline
```

Our `solve_counterfactual` zeros both the static C_wedge columns AND the expectations-through-D_wedge columns for inactive wedges, then re-solves for Phi_k, Phi_c. This is a different decomposition: ours is the "re-solve with restricted decision rules" whereas BCKM's is a "marginal static contribution, holding dynamics common."

Both are defensible BCA decompositions. BCKM's appendix (Section 2.C of the 2016 paper) describes the "fixed-expectations" approach we mostly implement. **Likely not a major source of the phi-stat discrepancy**, since if our VAR is different, our counterfactuals will be different regardless.

### 3.2 Phi-statistic formula

BCKM ([runall.m:68–78](BCKM/Multicountry%20-%20End/USAN2/runall.m)):

```matlab
mzye = Y(t)/Y(base) - mzy(t)/mzy(base)                  % level ratios
temp1 = mean(mzye^2, 1)                                 % MSE
phi_z = (1/temp1(1)) / sum(1/temp1)
```

Ours ([counterfactuals.py](bca_core/counterfactuals.py)):

```python
residual = data_hat - cf_hat                            # log-deviations
ssr = sum(residual^2)
phi = (1/ssr) / sum(1/ssr_i)
```

For small deviations `log(x/x₀) ≈ x/x₀ − 1`, so the formulas are approximately equivalent. For the Great Recession (investment −35%) the approximation breaks down. This is at most a second-order effect.

### 3.3 Counterfactuals base period

BCKM normalizes all counterfactual paths to a base period (1969Q1 or 1948Q1) before computing errors. We don't — we compute errors over the entire sample. This matters for the phi-statistics but only for the φ**value** (same MSE ordering).

---

## Tier 4 — Already OK

- **30-parameter θ layout**: we match BCKM (4 Sbar/P₀ + 16 P + 10 Q_tri). ✓
- **Klein QZ model solution**: we use the closed-form undetermined-coefficients solution (equivalent to BCKM's `roots` of the quadratic in `gammak`). ✓
- **Model equations** (resource constraint, labor FOC, Euler with BGG adj costs): checked vs. `res_adjust.m`, ours matches. ✓
- **Adjustment cost parameter `a = 12.88`** once calibrated from elasticity = 0.25: matches. ✓
- **DARE for initial covariance**: we use DARE; BCKM's steady-state Kalman gain is also DARE-derived. ✓

---

## Recommended order of fixes

1. **Extend sample back to 1948Q1 (or at least 1969Q1)** — single biggest leverage, trivial change. Ask first.
2. **Add imputed durables service flow to Y and C** — next biggest data fix. Consumer-durables stock is not on FRED; discuss with user whether to approximate via cumulated PCDG or pull the FOF tables from the BCKM/BEA folder.
3. **Use calibrated γ = 1.01934^(1/4) instead of estimating** — one-line change in run script.
4. **Update `CalibrationParams` defaults to BCKM values** (θ=0.35, ψ=2.24, δ=0.0464, β=0.9722 annual, σ=1) — one commit.
5. **Reparametrize θ from P₀ to Sbar, and add `phi0` column to the observation matrix** — the two cleanly-separable methodological fixes that would bring our setup in line with `mleqadj.m`. Medium-sized refactor.
6. **Revisit per-diagonal penalty** after items 1–5; see if we can drop it back to spectral-only now that other issues are resolved.
7. Consider switching to asymptotic Kalman likelihood last, only if convergence is still problematic.

Items 1–3 can likely be done without changing the algorithm — we'd expect measurable movement toward BCKM's numbers just from the data/calibration side. Items 4–5 complete the alignment.

---

## Files to examine more closely

| Matlab | What it does | Our counterpart |
|--------|--------------|-----------------|
| [mleqadj.m](BCKM/Multicountry%20-%20End/USAN2/mleqadj.m) | Likelihood, state-space build, SS | [var_estimation.py](bca_core/var_estimation.py) |
| [initmle.m](BCKM/Multicountry%20-%20End/USAN2/initmle.m) | Solves for Sbar matching sample means | — (we skip this) |
| [runmleadj.m](BCKM/Multicountry%20-%20End/USAN2/runmleadj.m) | Multi-restart L-BFGS driver | estimate_var_mle in var_estimation.py |
| [res_adjust.m](BCKM/Multicountry%20-%20End/USAN2/res_adjust.m) | Euler residual with BGG adj cost | PrototypeModel in [model.py](bca_core/model.py) |
| [fixexpadj.m](BCKM/Multicountry%20-%20End/USAN2/fixexpadj.m) | Counterfactual policy construction | solve_counterfactual in [counterfactuals.py](bca_core/counterfactuals.py) |
| [usdata.m](BCKM/Multicountry%20-%20End/USAN2/usdata.m) | Raw BEA→observables | [pipeline.py](bca_core/data/pipeline.py) + [adjustments.py](bca_core/data/adjustments.py) |
| [datamine.m](BCKM/Multicountry%20-%20End/USAN2/datamine.m) | Calibration + worktemp build | scripts/run_var_counterfactuals.py |
| [runall.m](BCKM/Multicountry%20-%20End/USAN2/runall.m) | Full pipeline + phi stats | scripts/run_var_counterfactuals.py |
| [kfilter.m](BCKM/Multicountry%20-%20End/USAN2/kfilter.m) / [doubalg.m](BCKM/Multicountry%20-%20End/USAN2/doubalg.m) | DARE via doubling algorithm | scipy.linalg.solve_discrete_are |
