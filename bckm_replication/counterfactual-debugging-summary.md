---
title: "Counterfactual Decomposition — Debugging Summary"
topic: "debugging"
layer: "bckm-replication"
status: "archived"
last_updated: "2026-05-03"
---

# Counterfactual Decomposition — Debugging Summary

## Problem Statement

The counterfactual decomposition in `bca_core/counterfactuals.py` produced wrong results:

- **Investment-only counterfactual** went to 280 (index, 2008Q1=100) for investment — should be 90-100, tracking data downward.
- **ALL single-wedge counterfactuals** showed the wrong sign during the Great Recession: output and investment went UP instead of DOWN.
- **Peak-to-trough** for investment-only: x went UP by +1.04 while actual data went DOWN by -0.35.

## Diagnostic Approach

We wrote `scripts/diagnose_counterfactuals.py` — a 10-step diagnostic that runs in order, printing PASS/FAIL for each step. The steps test IRF signs, all-wedges reproduction of data, smoothed wedge paths, policy vector signs, internal consistency, D\_wedge construction, cross-wedge effects, one-period verification, and sign conventions.

## Bug Found and Fixed

### `data_hat` constructed from endogenous capital evolution (FIXED)

**File:** `scripts/run_var_counterfactuals.py`, line 85

**Before (buggy):**
```python
cf_all_pol = solve_counterfactual(proto, P_var, active_wedges=[0, 1, 2, 3])
data_hat = run_counterfactual(smoothed_states, cf_all_pol)
```

**After (fixed):**
```python
obs = var_result["obs_hat"]
data_hat = {
    "y": obs[:, 0],
    "l": obs[:, 1],
    "x": obs[:, 2],
}
```

**Why this was wrong:** `run_counterfactual` evolves capital endogenously via `k' = P_k @ state`, but the Kalman smoother's capital includes gain corrections (backward smoothing adjustments) that the deterministic simulation doesn't have. Over 140 quarters, this causes cumulative divergence:

| Variable | Max error (hat-space) |
|----------|----------------------|
| y        | 0.007                |
| l        | 0.022                |
| x        | **0.083**            |
| k        | 0.022                |

Using the actual observables (which are exactly `D @ smoothed_state` by construction of the state-space model with `obs_cov = 0`) eliminates this divergence entirely.

**Impact:** The `data_hat` errors propagated into phi-statistics and peak-to-trough decompositions, making all downstream comparisons unreliable.

## Full Diagnostic Results

| Step | Test | Result | Detail |
|------|------|--------|--------|
| 1 | Existing tests | **PASS** | All 52 tests pass |
| 2 | IRF sign check | **PASS** | +taux shock: x=-0.109 (correct <0), y=-0.014 (correct <0), c=+0.017 (correct >0) |
| 3 | All-wedges reproduces data | **PASS** (after fix) | Policy vectors match obs via smoothed states: max\_err < 0.003 |
| 4 | Smoothed taux\_hat during recession | **WARNING** | taux\_hat fell from +0.088 to -0.034 (investment wedge *improved*) |
| 5 | Policy vector signs | **PASS** | P\_x\[3\]=-9.63, P\_y\[3\]=-1.46, P\_k\[3\]=-0.19 (all negative, correct) |
| 6 | CF matches \_solve\_with\_var | **PASS** | All policy vectors match exactly |
| 7 | D\_wedge sign check | **INFO** | C\_eff\[1,2\]=0.043 — investment wedge has tiny net effect in Euler due to (1+tau\_x) appearing on both sides |
| 8 | Diagonal P\_var | **INFO** | Cross-wedge interactions significant (P\_x differs by up to 0.23) |
| 9 | One-period verification | **PASS** | Exact match at mid-sample |
| 10 | Wedge sign conventions | **WARNING** | (1+tau\_x) decreased from 1.091 to 0.966 during recession |

## Why the "Wrong Sign" Is Not a Code Bug

### The investment wedge improved during the Great Recession

The Kalman smoother estimates that `taux_hat` (log-deviation of 1+tau\_x) **decreased** from +0.088 at 2007Q4 to -0.034 at 2009Q2. This means the investment wedge *improved* — investment became less distorted.

Since the policy coefficient P\_x\[3\] = -9.63 is correctly negative (higher tax -> less investment), the investment-only counterfactual correctly shows investment **increasing** when the wedge improves. The "wrong sign" is the model's answer, not a code error.

### Wedge movements during the Great Recession (smoothed)

| Wedge | 2007Q4 | 2009Q2 | Change | Interpretation |
|-------|--------|--------|--------|----------------|
| A (efficiency) | 0.987 | 0.975 | -0.011 | Worsened slightly |
| 1-tau\_l (labor) | 1.095 | 1.002 | -0.093 | **Worsened significantly** |
| 1+tau\_x (investment) | 1.091 | 0.966 | -0.122 | Improved (counter-intuitive) |
| g (government) | 0.873 | 1.054 | +0.188 | Increased (stimulus) |

The labor wedge absorbed most of the recession, consistent with findings in the BCA literature (CKM 2007).

### Phi-statistics (variance decomposition)

After the fix, the phi-statistics are:

|            | y    | l    | x    |
|------------|------|------|------|
| Efficiency | 0.52 | 0.49 | 0.61 |
| Labor      | 0.23 | 0.25 | 0.09 |
| Investment | 0.05 | 0.07 | 0.03 |
| Government | 0.20 | 0.20 | 0.27 |

The efficiency wedge dominates (52% of output variation), and the investment wedge is small (5%). This is consistent with the standard BCA literature.

## Why Investment-Only CF Goes to 280

Three compounding factors:

1. **taux\_hat dropped by -0.122 during the recession.** With P\_x\[3\] = -9.63, the direct effect is -9.63 x (-0.122) = +1.18 in x\_hat — a massive increase in investment.

2. **Capital divergence over 140 quarters.** The investment-only CF evolves capital endogenously with only taux active. By 2007Q4 (t=111), CF capital has diverged by 0.18 from data capital. This shifts the CF to a completely different part of state space.

3. **Resource constraint amplification.** Investment is a residual: x = (y - c - g) / x\_ss. Small changes in consumption (driven by the Euler equation) get amplified into large investment swings. The coefficient x\_c (consumption's effect on investment) is approximately -7, creating a large multiplier.

## Why C\_eff Is Small for the Investment Wedge

The investment wedge `(1+tau_x)` appears on **both sides** of the Euler equation (CKM 2007, eq. 8):

- **LHS (today):** `(1+tau_x_t)` — cost of investing today
- **RHS (tomorrow):** `(1-delta) * (1+tau_x_{t+1})` — resale value includes future tax

When the wedge is persistent (rho\_taux = 0.976), these nearly cancel:

```
C_eff[Euler, taux] = 1.0 - c_delt * rho_taux
                   = 1.0 - 0.974 * 0.976
                   = 0.043
```

This means a persistent investment wedge has almost no net effect on the Euler equation — the cost of investing today is offset by the higher resale value of capital tomorrow. The investment-only counterfactual is therefore dominated by capital dynamics rather than the wedge itself.

## Files Changed

| File | Change |
|------|--------|
| `scripts/run_var_counterfactuals.py` | Fixed `data_hat` to use actual observables; removed unused imports |
| `scripts/diagnose_counterfactuals.py` | **New** — 10-step diagnostic script |

## Potential Next Steps

1. **Investigate why taux\_hat improves during the recession.** This could be a modeling/identification issue. The 4D VAR with cross-wedge interactions may be misattributing the investment decline to other wedges.

2. **Consider alternative counterfactual presentations.** Starting CFs from the recession peak (rather than t=0) would reduce capital divergence and may produce more interpretable recession-window comparisons.

3. **Check sensitivity to P\_var structure.** The cross-wedge interactions in the estimated P\_var are significant (Step 8). A restricted VAR (e.g., block-diagonal) might produce more economically sensible wedge paths.
