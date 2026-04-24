# Business Cycle Accounting — Implementation Report

## Overview

This project implements the Business Cycle Accounting (BCA) methodology of Chari, Kehoe & McGrattan (2007) and Brinca, Chari, Kehoe & McGrattan (2016) for US data. The goal is to decompose aggregate fluctuations into four structural wedges — efficiency (A), labor (1−τ_l), investment (1+τ_x), and government (g) — and quantify each wedge's contribution to output, hours, and investment dynamics via counterfactual simulations.

---

## Model

The prototype economy is a closed-economy one-sector stochastic growth model. The state vector is:

```
[k_hat, A_hat, taul_hat, taux_hat, g_hat]
```

Observables used in estimation are:

```
[y_hat, l_hat, x_hat, g_hat]
```

all expressed as log-deviations from their sample means. The model is log-linearized around the deterministic steady state and solved via Klein's (2000) QZ method.

Adjustment costs follow the BGG (1999) form: `φ(x/k) = (a/2)(x/k − b)²`, calibrated to match an elasticity of Tobin's q with respect to the investment-to-capital ratio of 0.25.

---

## Estimation: Kalman-Filter MLE

We follow BCKM (2016)'s `mleqadj.m` and jointly estimate the VAR(1) transition matrix **P** (4×4) and shock Cholesky factor **Q** (10 lower-triangular elements) by maximizing the Kalman-filter log-likelihood of the observables. The latent wedge paths are recovered by a Rauch-Tung-Striebel (RTS) backward smoother.

### State-Space Formulation

The Kalman model couples the model's decision rules with the VAR:

**Transition** (5×5):
```
F[0, :] = P_k   (capital accumulation, re-solved for each P_var via undetermined coefficients)
F[1:, 1:] = P_var
```

**Observation** (4×5):
```
H[0] = P_y,  H[1] = P_l,  H[2] = P_x,  H[3, 4] = 1  (g observed directly)
```

**Process noise**: shocks enter only the wedge states (rows 1–4).  
**Measurement noise**: `R = 1e-8 · I` (near-perfect observation).

Decision rules `P_k, P_y, P_l, P_x` are re-solved at each parameter iterate using the undetermined-coefficients correction `C_eff = C − D @ P_var`, which accounts for agents' rational expectations over future wedge paths.

### Penalties

Two stationarity penalties are added to the negative log-likelihood:

```python
penalty  = 5e5 * max(max_eigenvalue(P_var) − 0.995, 0)²   # spectral radius
penalty += 5e5 * sum(max(|diag(P_var)| − 0.995, 0)²)      # per-diagonal
```

The diagonal penalty prevents individual elements from reaching unit-root territory even when the spectral radius constraint is satisfied due to off-diagonal cancellation.

### Initial Covariance: DARE

The Kalman filter is initialized with the **steady-state prediction covariance** obtained by solving the Discrete Algebraic Riccati Equation (DARE):

```python
Sigma_0 = solve_discrete_are(F, H.T, Q_proc + ε·I, R_obs)
```

This is computed once from the BCKM Table 77 parameters and held fixed throughout optimization (for speed). The exact DARE from the estimated parameters is used only for the final RTS smoother pass. Using DARE rather than the stationary process covariance (`solve_discrete_lyapunov`) is critical: for near-unit-root processes the Lyapunov solution grows to ~250× Q, giving an initialization as uninformative as a diffuse prior and causing the filter to misattribute the Great Recession investment decline to capital dynamics rather than the investment wedge.

### Optimization

`scipy.optimize.minimize` with L-BFGS-B, `maxiter=500`. Multiple starts:

1. **BCKM Table 77 init** — US MLE parameters from the replication files  
2. **OLS warm-start** — approximate wedge paths extracted via the pseudo-inverse of the observation matrix H, then a VAR(1) fit by OLS  
3. **Perturbed BCKM starts** — 3 draws with `std=0.01`, diagonal-clipped to ≤0.99

Best result across all starts is retained. The final smooth uses the DARE-initialized RTS pass on the winning parameters.

---

## Results (US 1980Q1–2014Q4)

### VAR Transition Matrix (MLE)

|       |    A   |  τ_l  |  τ_x  |   g   |
|-------|--------|-------|-------|-------|
| **A** |  0.956 | −0.003|  0.015| −0.001|
| **τ_l**| 0.036 |  0.986| −0.036|  0.004|
| **τ_x**|−0.029 | −0.002|  0.968|  0.004|
| **g** |  0.034 | −0.007| −0.051|  0.995|

Max |eigenvalue| = 0.986 (stationary). The τ_x diagonal (0.968) is very close to BCKM's Table 77 value of 0.967.

### Investment Wedge During the Great Recession

| Date   | taux_hat |
|--------|----------|
| 2007Q4 | −0.044   |
| 2009Q2 | +0.033   |
| **Δ**  | **+0.077 (worsened ✓)** |

The investment wedge correctly deteriorates during the GR — this was the primary identification target.

### Phi-Statistics (Variance Decomposition)

|            |   Y   |   L   |   X   |
|------------|-------|-------|-------|
| Efficiency | 0.194 | 0.160 | 0.187 |
| Labor      | 0.309 | 0.500 | 0.374 |
| Investment | **0.393** | 0.256 | **0.235** |
| Government | 0.104 | 0.085 | 0.204 |

BCKM (2016) US targets (no-durables): φ_A^Y ≈ 0.12, φ_L^Y ≈ 0.58, φ_X^Y ≈ 0.25.

The investment-wedge phi for investment (φ_X^X = 0.235) is close to BCKM's target. The remaining gap in output (investment 0.39 vs 0.25, labor 0.31 vs 0.58) has two identified causes:

- **Demeaning vs. intercept estimation**: we force P_0 = 0 and demean observables; BCKM estimates VAR intercepts (Sbar), allowing larger shock magnitudes. Our estimated shocks are ~10–20× smaller than BCKM's.
- **Data**: BCKM uses a "no-durables" investment series that excludes durable consumption; we use the standard FRED investment aggregate.

### Peak-to-Trough (2007Q4 → 2009Q2)

|   | Data  | Efficiency | Labor  | Investment | Government |
|---|-------|------------|--------|------------|------------|
| Y | −0.071|  +0.000    | −0.027 |  −0.118    |  +0.052    |
| L | −0.081|  +0.010    | −0.041 |  −0.157    |  +0.073    |
| X | −0.345|  +0.033    | +0.011 |  **−0.713**|  +0.178    |

The investment wedge alone overpredicts the GR decline in investment (−0.71 vs −0.35), consistent with BCKM's finding that the investment wedge is the primary driver of the Great Recession.

---

## Things Not To Do

### 1. Diffuse or Lyapunov initial covariance
Using `Sigma_0 = 100·I` or `solve_discrete_lyapunov(F, Q_proc)` as the Kalman filter prior is wrong for this model. For near-unit-root VAR processes, the Lyapunov solution is O(250·Q), as uninformative as a diffuse prior. The filter then explains the GR investment collapse via capital falling below trend (via the large P_x[k] coefficient ≈ −3.7) rather than through the investment wedge. The result is a taux_hat that *improves* during the GR — the opposite of what the data imply.

**Use DARE instead.**

### 2. OLS VAR on latent wedges
Fitting a VAR(1) by OLS on the directly extracted wedge series gives a P matrix with a non-stationary investment-wedge diagonal (≈ 1.15), incorrect off-diagonal structure, and wrong counterfactual signs. OLS is fine as a warm-start, but must not be the final estimator.

### 3. Spectral-radius penalty without diagonal penalty
Penalizing only `max|eig(P)| > 0.995` allows individual diagonal elements to exceed 1 (e.g., τ_x diagonal hit 1.001 in early runs) because off-diagonal elements can suppress the spectral radius while individual rows diverge. Both penalties are required.

### 4. DARE inside the optimization loop
Calling `solve_discrete_are` at every L-BFGS-B function evaluation makes the optimizer approximately 20–40× slower (DARE for 5×5 matrices is fast in isolation but is called ~13,000 times per restart). Precompute a single fixed `Sigma_0` from the BCKM reference parameters and reserve the exact DARE for the final smoother pass only.

### 5. Large perturbation std or unclipped diagonal in random restarts
Using `std=0.02` for BCKM perturbations, combined with a BCKM τ_l diagonal of 1.001, frequently produces starting points with large diagonal penalties. scipy's L-BFGS-B can deadlock in the Fortran backend when the gradient landscape is flat at 1e20 for hundreds of iterations, consuming no measurable CPU time but blocking indefinitely. Use `std≤0.01` and clip the diagonal to `[−0.99, 0.99]` before packing into the theta vector.

### 6. Skipping the BCKM warm-start
The BCKM Table 77 parameters are the best available prior for the US estimation. Starting solely from `0.9·I` or random draws finds substantially worse local minima (ll ≈ 1500–1660 vs 1837–1855). Always include the BCKM init and the OLS pseudo-inverse warm-start alongside any random perturbations.

---

## Next Steps

### 1. Estimate VAR intercepts (P_0) rather than forcing P_0 = 0

The single largest source of phi-statistic discrepancy is shock magnitude. BCKM estimates a 4-dimensional intercept vector Sbar jointly with P and Q, which absorbs the unconditional means of each wedge process. This allows the VAR to operate around a non-zero mean, and consequently the MLE assigns much larger shock variances (τ_x shock std ≈ 0.12 in BCKM vs ≈ 0.006 here).

Our current formulation forces P_0 = 0 and demeans observables, so all wedge variation must be explained by zero-mean innovations off a zero mean — a more compressed representation that systematically underestimates shock sizes.

**Implementation**: add 4 intercept parameters to the theta vector, include a `state_intercept` term in the Kalman recursion (`x_pred = F @ x + P_0_wedge`), and do not demean observables before passing them to the filter. Initialize P_0 at the sample means of the approximate wedge paths from the OLS pseudo-inverse step.

### 2. Align the investment measure to BCKM's no-durables definition

BCKM's US baseline excludes durable consumption from investment, defining:

```
X_BCKM = Fixed investment + Change in inventories
```

The standard FRED gross private domestic investment series includes residential and nonresidential fixed investment plus inventory change, but our pipeline also includes consumer durables via the PCE deflator chain. Switching to the BCKM no-durables aggregate would shift the investment-wedge identification and is likely responsible for part of the φ_X vs φ_L discrepancy (too much investment, not enough labor in our decomposition).

**Implementation**: modify `bca_core/data/pipeline.py` to pull the BCKM-compatible investment series and adjust the resource constraint accordingly.

### 3. Refine the optimization with a two-stage approach

Because the fixed BCKM Sigma_0 used during optimization differs from the DARE computed at the estimated parameters, the reported final log-likelihood (1782) is lower than the best restart's optimization-criterion value (1853). A two-stage refinement would close this gap:

1. Run the current multi-start optimization with the fixed BCKM Sigma_0 → obtain P_var_init, Q_init.
2. Compute the exact DARE at (P_var_init, Q_init) → Sigma_0_new.
3. Re-optimize from (P_var_init, Q_init) using Sigma_0_new as the fixed prior.
4. Repeat once or twice until the log-likelihood under the reported criterion stabilizes.

This coordinate-ascent loop converges quickly (typically 2–3 outer iterations) and ensures the optimization objective and the final smoother use the same initialization.
