# BCA Project Diary

---

## Session: 2026-04-23

### What was completed since the last entry

**Phase 1 — CalibrationParams BCKM defaults** ✅  
Updated `bca_core/params.py` to match BCKM (2016) Table 77 calibration:
- `alpha`: 1/3 → 0.35 (BCKM θ=0.35)
- `psi`: 2.5 → 2.24 (BCKM ψ=2.24)
- `delta_annual`: 0.05 → 0.0464
- `rho_annual`: 0.025 → 0.02860 (so β_annual = 0.9722)

Updated `tests/test_params.py` and `tests/test_var_estimation.py` to match new defaults. The `test_at_steady_state` test had a wrong premise — `prepare_observables` computes `y_hat = log(y_dt)` without subtracting `log(y_ss)`, so at SS `y_hat = x_hat = g_hat = log(y_ss) ≠ 0`. Fixed the test to assert this correct behavior. All 52 tests now pass.

**Phase 2 — Sample period extension** ✅  
Changed `start="1980Q1"` → `start="1969Q1"` in `scripts/run_var_counterfactuals.py`. Effective sample is 1976Q1–2014Q4 (T=152) because FRED `LFWA64TTUSQ647S` (OECD working-age population) only starts ~1976. Also added `Sbar` diagnostic print to the script.

**Phase 3 — Sbar reparametrization** ✅  
Most impactful change. BCKM optimizes Sbar (unconditional wedge mean) and derives `P_0 = (I−P)·Sbar`. Near the unit root, this gives a much better-conditioned optimization landscape than optimizing P₀ directly.

Changes to `bca_core/var_estimation.py`:
- `_unpack`: extracts `Sbar = theta[:4]`, computes `P_0 = (I−P) @ Sbar`
- `_pack`: first 4 elements are now Sbar, not P_0
- BCKM warm-start: converts `_P_0_bckm` to `_Sbar_bckm` via `solve(I−P, P_0)`
- OLS warm-start: converts `P_0_ols` to `Sbar_ols`
- Perturbed restarts: perturb Sbar instead of P_0
- Return dict includes both `"P_0"` and `"Sbar"` keys

**Phases 4 and 5 — attempted and reverted**

Phase 4 (GI reclassification) was attempted twice and reverted both times:
- First attempt: incorrectly identified `A955RC1Q027SBEA` as government investment; it is actually government CONSUMPTION (≈77% of GCE). This gave g ≈ GI+NETEXP ≈ 0.7% of GDP — obviously wrong.
- Second attempt: correctly derived GI = GCE − gov_consumption, added it to x_adj. Pipeline ran but x_hat mean = +0.147, phi-stats dominated by investment (71%). Root cause: model SS x/y ≈ 25.45% but data x/y ≈ 30% with GI included — the Kalman smoother absorbs the 5% gap through taux_hat, making the investment wedge nearly flat (std=0.03) and producing degenerate phi-stats. This requires BCKM's approach of re-solving the model SS at each iterate (using Sbar/phi0), which our implementation does not support.

Phase 5 (level-ratio phi-stats) was attempted and reverted: the `(Y(t)/Y(0) − CF(t)/CF(0))²` formula over a 38-year sample amplifies secular drift, making the government wedge dominate (phi_G=80%) regardless of business-cycle variation.

Both phases reverted to original code.

### Final pipeline results (Phases 1–3 only)

```
T = 152 (1976Q1–2014Q4)
g_share = 0.1678
SS: y=1.3528, l=0.3343, x/y=0.2545, g/y=0.1678

P diagonal: [0.9954, 0.9909, 0.9889, 0.9874]   (A improved from 0.886 → 0.9954!)
Max |eigval|: 1.0051  (barely non-stationary)
Sbar: [1.3668, -4.2908, -1.9985, -2.8325]
P_0:  [-0.0009, 0.0033, -0.0051, -0.0049]        (target: [0.014, 0.001, 0.013, -0.014])

Log-likelihood: optimized=1990.58, final=1846.11  (144-unit gap)

Smoothed wedge std: A=0.015, τ_l=0.093, τ_x=0.113, g=0.122

phi-stats (log-deviation):
  efficiency:  [0.015, 0.059, 0.015]
  labor:       [0.034, 0.132, 0.032]
  investment:  [0.140, 0.120, 0.049]
  government:  [0.811, 0.689, 0.904]   ← dominating (target: ~0)

Peak-to-trough 2007Q4→2009Q2:
  actual y: -0.0708
  CF-investment: -0.1463
  CF-government: +0.0410
```

### What is in progress

Nothing — all planned phases either completed or cleanly reverted. Tests pass. Pipeline runs.

### What is blocked and why

**1. 144-unit log-likelihood gap (highest priority)**  
`optimized ll = 1990` vs `final ll = 1846`. The optimizer uses a frozen Σ₀ computed once at startup from BCKM parameters; the true Σ₀ at the optimal parameters is very different. BCKM's `kfilter.m` uses an asymptotic (DARE-derived) constant Kalman gain K, which is Σ₀-independent. We use a time-varying Kalman filter with a frozen Σ₀.  
Fix: Switch to steady-state Kalman with constant K throughout. K comes from DARE at each candidate parameter set (run once per objective evaluation, not inside the filter loop). This is feasible — scipy has `solve_discrete_are`.

**2. Government wedge dominates phi-stats (phi_G ≈ 81%)**  
Root cause: `compute_government_wedge` = GCE + NETEXP, where GCE includes government investment (GI). BCKM uses only government consumption (GC = A955RC1Q027SBEA) in G, and puts GI into X. This would fix the composition of G and X.  
Fix is blocked by: we need to re-solve the model SS at each Sbar iterate (because moving GI from G to X changes g_share which changes SS which changes the observation matrix H). BCKM handles this via a phi0 constant in the state-space — not currently implemented.

**3. Max eigval = 1.0051 (barely non-stationary)**  
Soft penalty allows optimizer to find solutions just outside the spectral radius bound. The penalty `5e5 * max(eig_max - 1.005, 0)^2` = 0.005 at eigval=1.0051, which is negligible compared to ll≈1990.  
Fix: Tighten threshold from 1.005 to ≤0.999, or use a hard constraint projection.

**4. T=152 not ~182 (sample too short)**  
FRED `LFWA64TTUSQ647S` starts ~1976, not 1969. BCKM's US sample is 1969Q1–2014Q4 (182 quarters).  
Fix: Replace with a series that goes back to 1969 — e.g., construct from `LNS11300060` (BLS civilian labor force 16-64) or interpolate annual Census data. Or use total population / scaling factor as a proxy.

### Exact next step

**File:** `bca_core/var_estimation.py`  
**Function:** `_kf_ll` (the Kalman filter loop, ~line 350–430)  
**What to change:** Replace the time-varying Kalman filter with a steady-state filter using a constant gain K:

```python
# Before optimization loop, solve DARE at current theta to get K:
from scipy.linalg import solve_discrete_are
Sigma_inf = solve_discrete_are(F.T, H.T, Q_state, R_obs)
K = Sigma_inf @ H.T @ np.linalg.inv(H @ Sigma_inf @ H.T + R_obs)
# Use constant K throughout the filter; no Sigma update inside loop
```

The key point: run DARE once per objective evaluation (not inside the T-step loop), use the resulting K as the fixed gain for all T steps. This is O(n³) once per objective call vs. O(T·n³) for the time-varying filter — same cost asymptotically but eliminates the Σ₀ sensitivity entirely.

Reference: `BCA/BCKM/Multicountry - End/USAN2/kfilter.m` lines 1–40 for the exact implementation pattern BCKM uses.

---

## Session: 2026-04-27

### Pickup state

Re-read `BCA_info.md` (especially Section 4 on data construction and Section 7 on the US MLE target tables) and traced the data pipeline end-to-end. No code changes this session — this entry is a scoping note.

### Findings (data side, vs BCA_info Section 4 / BCKM Matlab)

Already implemented correctly:
- Durables PIM stock + 4% service flow added to both C and Y (`bca_core/data/adjustments.py:14`).
- Sales tax subtracted from Y and C (`bca_core/data/adjustments.py:58`).
- Per-capita / GDP-deflator conversion.

Missing or wrong:
- **GI is in G, not in X.** `compute_government_wedge` (`adjustments.py:140`) sets `g = GCE + NETEXP`, but FRED `GCE` lumps gov consumption and gov investment. Matlab definition: `G = rGC + rEX − rIM` and `X = rCD + rGPDI + rGI − …`. Need to fetch a gov-consumption-only series (e.g. NIPA `A955RC1Q027SBEA`) or compute `gov_investment = GCE − gov_consumption`, subtract from G, add to X. This is the leading suspect for φ_G dominating at ~81% in the last run.
- **Durables expenditure (`rDCD`) not added to Y.** Matlab is `Y = rGDP − rSTX + 0.04·rKCD + rDCD`; current `adjustments.py:52` adds the service flow but not durables expenditure itself.
- **Sales tax not applied to X.** Matlab subtracts a `(rCD/rCNDS)·rSTX` term from X to remove the tax embedded in durables consumption that is being moved to investment. Current code does not.
- **Detrending uses OLS-estimated γ.** `pipeline.py:114` fits a linear trend to `log(y)` and uses that slope as γ for detrending. BCKM uses calibrated `gz = 1.018^(1/4) − 1` (`datamine.m`). The `params.gamma_annual` set in `run_var_counterfactuals.py:120` is only used by the model SS, not by the detrender — so we currently have *two* growth rates floating around: estimated for detrending, calibrated for the model.
- **Sample period.** `run_var_counterfactuals.py:110` uses `start="1969Q1"` (effective T=152 from 1976Q1 due to OECD pop series start). BCKM Section 7 target is **1980Q1–2014Q4** (T=140). The extension makes apples-to-apples comparison against Tables 8–11 harder.

### Architectural finding — Phase 4 blocker was overstated

The previous diary entry said GI reclassification was blocked because "BCKM re-solves the model SS at each iterate (using Sbar/phi0)." Re-tracing the code, this is wrong.

`var_estimation.py:317` calls `proto.steady_state()` once on entry to `estimate_var_mle` and then uses the resulting `ss` and `log_linearize(ss)` output throughout the optimization. The MLE loop only varies `P_var` and `Q_chol`. As long as `params.g_share` is set from corrected data in `main()` *before* constructing `params` (which it is — `run_var_counterfactuals.py:116`), the SS rebuilds cleanly and downstream nothing else depends on per-iterate SS.

BCKM's `phi0` (in `mleqadj.m`) is a constant offset added to the observation matrix to absorb the gap between model SS and data mean. It is computed *once* from the SS, not at each iterate. The Phase 4 attempts likely failed for a different reason (probably mismatch between corrected data ratios and the SS that was rebuilt from them, e.g. the model's `b = γ + δ + n` not being adjusted), not because of architectural impossibility.

**Implication:** The GI split is mechanical, not architectural — same difficulty class as Phases 1–3.

### Next steps (in priority order)

**Step 1 — Sample alignment + calibrated γ** *(15 min, low risk)*
- Change `start="1969Q1"` → `start="1980Q1"` in `scripts/run_var_counterfactuals.py:110`.
- Plumb a `gamma_annual` argument through `build_us_dataset` → `pipeline.py` so detrending uses the calibrated value (1.9%/yr) rather than OLS slope on log(y).
- Goal: clean baseline whose results map onto Tables 8–11.

**Step 2 — GI split** *(half day, biggest expected impact on φ_G)*
- Add a FRED series for gov investment in `bca_core/data/fred.py`. Candidates: `A782RC1Q027SBEA` (gov gross investment), or derive `gov_investment = GCE − A955RC1Q027SBEA` (gov consumption).
- In `compute_government_wedge`, replace `GCE + NETEXP` with `gov_consumption + NETEXP`.
- In `reclassify_durables` (or a new step), add `gov_investment` to `x_adj` alongside `gpdi + pce_durables`.
- `g_share_data` in `main()` will auto-recompute to BCKM-like ~0.10 (down from 0.166).
- Confirm SS rebuilds correctly at the new g_share; re-run tests; re-run pipeline; report φ-stats.

**Step 3 — NIPA fidelity tweaks** *(small, do while in `adjustments.py`)*
- Add `rDCD` (durables expenditure) to `y_adj` in addition to the service flow.
- Apply sales tax to X via the `(rCD/rCNDS)·rSTX` term.

**Step 4 — Steady-state Kalman gain** *(previously queued; still relevant)*
- The `_kf_ll` constant-gain rewrite from the 2026-04-23 entry — closes the 144-unit LL gap. Independent of Steps 1–3; can be done before, after, or in parallel.

### What is blocked and why

Nothing hard-blocked. Step 2's only subtlety is making sure `params.b = γ + δ + n` and the resulting model `x/y` ratio stays consistent with the corrected data `x/y` after GI is added — a sanity check, not a blocker.

### Exact next step

Start with Step 1 (sample lock + calibrated γ), then run the existing pipeline once to capture a clean pre-GI baseline before touching `adjustments.py`.
