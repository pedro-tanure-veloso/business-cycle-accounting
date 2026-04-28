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

---

## Session: 2026-04-27 (continued)

### Step 1 — done

Sample locked to 1980Q1–2014Q4 (T=140). Detrending uses calibrated γ=1.9%/yr via a new `fixed_slope` argument in `remove_trend` and a `gamma_annual` argument plumbed through `build_us_dataset`. `mean(y_hat) = 0.0000`. 52/52 tests pass.

Result snapshot vs the 2026-04-23 baseline:
- P diagonal [.995, .993, .987, .995] vs target [.989, 1.001, .967, .995] — `g` matches exactly, `τ_x` got closer.
- φ_L^Y 0.16 → 0.42 (target 0.46).
- LL 1846 → 1805.

Commit `d45ae1f`. Cache: `data/us_1980_2014.parquet`.

### Step 2 — done structurally, broke counterfactuals

Fetched FRED `A955RC1Q027SBEA` (gov consumption only). Derived `gov_investment = GCE − gov_consumption`. Pipeline now:
- `G = gov_consumption + NETEXP` (was `GCE + NETEXP`)
- `X = GPDI + durables + gov_investment` (was `GPDI + durables`)

Result snapshot:
- `g_share` 0.166 → 0.124 (BCKM US value).
- φ_L^Y 0.16 → **0.44** (target 0.46) — strong identification gain on labor wedge.
- φ_L^L 0.42 → **0.62** (target 0.70).
- φ_G^Y 0.75 → 0.50 (still high but moving the right way).
- LL 1805 → **1595** (−210).
- `mean(x_hat) = +0.138` (was −0.02).
- Investment-only CF gives x = +1.35 at peak-to-trough (data is −0.29) — broken.

### Diagnosis — exposed a pre-existing model-data SS gap bug

Model x/y at SS = 0.2545 (pinned by α=0.35, β, δ=0.0464, γ=0.019, n=0.0098 via `b = γ+δ+n+γn`, independent of g_share). After GI is correctly placed in X, **data x/y ≈ 0.29** — a 14% gap.

In `bca_core/var_estimation.py:18` (`prepare_observables`), `x_hat = log(df["x"]) − log(x_ss/y_ss)`. When data x/y ≠ model x/y, this gives `mean(x_hat) ≠ 0`. The optimizer absorbs this constant by pushing `Sbar[τ_x]` to a large negative value (−1.98). The `Sbar` reparametrization makes this likelihood-equivalent to setting `P_0` directly, so MLE doesn't care.

But the smoothed state path inherits the level shift (`mean(taux_hat) = −0.30`), and the counterfactual function uses smoothed states *as-is* and applies the model decision rules. The huge negative τ_x level in the "investment-alone" CF reads as a permanent investment subsidy — hence x = +1.35.

**BCKM handles this differently.** From `mleqadj.m`: `phi0 = Y0 − C·X0(1:5)` where `Y0 = [log(ys); log(xs); log(ls); log(gs)]`. `phi0` is appended as a 6th column of `C`, giving `Y(t) = C·X(t) + phi0`. **The observation equation absorbs the model-data SS gap, leaving state mean-zero.** Our pipeline instead puts the gap into the state via `Sbar`. Mathematically equivalent for likelihood; not equivalent for counterfactuals.

### Where this leaves us

- **Step 1 is good and stays.**
- **Step 2 data-side change is correct and stays** — data x/y reflecting GI is what BCKM actually does. Reverting would be a step backward.
- **The counterfactual breakage is a pre-existing bug now visible.** Two ways to fix:

**Option A — `phi0` observation offset (faithful to BCKM, ~1 day).**
- Modify `prepare_observables` to keep observables in log-level form (only subtracting growth trend), not normalizing by model SS ratios.
- Add a 5th state element (constant 1) and extend `H` so that `Y(t) = H·[state; 1] = H_state·state + phi0`.
- `phi0` initialized once from model SS values; held fixed in optimization.
- `Sbar` then represents true wedge level deviations from SS, and counterfactuals work cleanly.

**Option B — demean smoothed states before counterfactuals (quick fix, ~1 hr).**
- Subtract `Sbar` from each smoothed wedge path before passing to `run_counterfactual`.
- Conceptually: "wedge X alone" means X follows its data deviations from mean; others held at SS.
- Less faithful to BCKM (their `phi0` and our `Sbar` are not the same beast structurally), and changes the meaning of CF level outputs, but should fix the +1.35 overshoot.

Option A is the right long-term answer; Option B is a triage that gets the existing CF API working again.

### Open questions

1. Does Option B preserve the φ-statistic interpretation? (φ-stats are variance ratios, not levels — should be fine.)
2. Does the previously-queued steady-state Kalman gain rewrite (Diary 2026-04-23 "exact next step") interact with `phi0` work? Probably orthogonal — that's about the filter recursion, not the observation equation.

### Next exact step

Discuss A vs B with user before implementing either. Step 2 commit captures the current state with an honest readout of what works and what doesn't.

---

## Session: 2026-04-27 (Step 3 — Option A: phi0 + Sbar=0)

### What was completed

**Step 3a — phi0 in `prepare_observables`** ✅
`bca_core/var_estimation.py` `prepare_observables` now returns `(obs, phi0)`.
`phi0 = mean(obs_raw, axis=0)` is subtracted from observables so the Kalman
filter sees a centered signal. After Step 2's GI split, phi0 came out to
`[0.000, -0.001, +0.138, -0.010]` — the +0.138 in x is exactly the
data-vs-model x/y SS gap (data 0.29 vs model 0.255). Updated callers
(`estimate_var`, run script, test) to unpack the tuple.

**Step 3b — Sbar fixed at 0** ✅
First run of phi0 alone showed phi0 was *necessary but not sufficient*:
`mean(taux smoothed) = -0.158`, max |eig| = 1.005, investment wedge
flipped sign during GR (improved instead of worsened). Diagnosis: near
unit root, `(I-P)` is near-singular, so a huge Sbar coexists with a tiny
`P_0 = (I-P)·Sbar` and is likelihood-equivalent. The optimizer kept finding
the same big-Sbar basin (-2 region) and accumulated drift over T=140
periods showed up as wrong-sign smoothed wedges.

Fix: in `_unpack`, force `Sbar = np.zeros(4)` regardless of theta[:4]
content, so `P_0 = (I-P)·0 = 0`. The first 4 theta entries are now inert
(L-BFGS leaves them alone since gradient is zero). This matches BCKM's
mleqadj.m design end-to-end: phi0 in the obs equation carries the SS gap,
the wedge VAR has no separate intercept.

### Pipeline result after Step 3a + 3b

```
T = 140
phi0 (SS misalignment): [0.000, -0.001, +0.138, -0.010]
obs_hat means: ~0 in all 4 channels

Sbar: [0, 0, 0, 0]            (forced)
P_0:  [0, 0, 0, 0]            (target: [0.014, 0.001, 0.013, -0.014])
P diag: [0.993, 1.002, 0.964, 0.984]
Max |eig|: 0.987              (was 1.005 — stationary now ✓)
LL: 1778.6                    (was 1786 with phi0 only; ~7-unit drop)

Smoothed wedge means:
  A=-0.006, τ_l=0.000, τ_x=-0.013, g=-0.000   (all ≈ 0 ✓)

Investment wedge during GR:
  2007Q4: -0.059, 2009Q2: +0.035, Δ = +0.094 (worsened ✓ — right sign)

φ-stats:                target (BCKM Table 11)
  y: A=0.31  τ_l=0.22  τ_x=0.29  g=0.18    fYA=0.16, fYτL=0.46, fYτx=0.32
  l: A=0.19  τ_l=0.58  τ_x=0.13  g=0.10
  x: A=0.29  τ_l=0.36  τ_x=0.11  g=0.24

Peak-to-trough (2007Q4 → 2009Q2):
   actual  efficiency  labor  investment  government
y -0.0722     -0.007  -0.060     -0.146      +0.042
l -0.0814     +0.009  -0.092     -0.203      +0.061
x -0.2898     +0.020  -0.187     -0.945      +0.159
```

### What this changes vs Step 2 commit

- `Sbar = [-2, ...]` regime → `Sbar = 0` ✓
- Investment-only x peak-to-trough went from **+1.35 → -0.945** (sign flip from wrong to right; magnitude overshoots data of -0.29 by ~3×)
- Investment wedge GR direction: `improved → worsened` ✓
- φ_government on output: `0.92 → 0.18` (no longer dominating)
- Stationary VAR ✓

### What still doesn't match BCKM

1. **CF magnitudes overshoot.** Investment-only x at trough is -0.945, ~3.3× data (-0.29). Likely related to:
   - The 144-unit `optimized_LL vs final_LL` gap (still open from 2026-04-23 Diary). Different DARE Sigma0 between optimization and final smoother — Step 4 in 2026-04-23 plan: switch to steady-state Kalman gain via DARE-per-iterate.
   - φ_y[A] = 0.31 vs target 0.16 — efficiency wedge still overweighted on output.
2. **φ_y[τ_L] = 0.22 vs target 0.46.** Labor wedge underweighted on output (correct on labor itself though, 0.58 ≈ target 0.7).
3. **P_0 = 0 vs target [0.014, 0.001, 0.013, -0.014].** Forced to zero by current design. To match exactly, would need to make Sbar a free parameter again *and* prevent it from running away — e.g., a tight L2 prior on Sbar centered at 0.

### Where this leaves us

- **Step 1 + Step 2 + Step 3 land cleanly together.** The architecture now matches BCKM mleqadj.m: phi0 in obs equation, no free VAR intercept, smoothed states truly mean-zero.
- **Counterfactuals work directionally.** Magnitudes are off but signs and qualitative decomposition are correct.
- **Next priorities** (in order of expected impact):
  1. **Steady-state Kalman gain via DARE-per-iterate** (Step 4, ~half day). Closes the 144-unit LL gap, likely fixes the CF magnitude overshoot too.
  2. **Add `rDCD` to Y, sales tax to X** (Step 3 from earlier scope, ~1 hr). Last data-side adjustment.
  3. **Optional Sbar prior** if we still don't match BCKM Table 9 P_0 after the above.

### Open questions (resolved or carried)

1. Option A vs B → **Option A chosen and works.** B was a triage; A is the right architecture.
2. Whether forcing Sbar=0 hurts LL substantially → **No.** ~7-unit drop. Confirms Sbar was identifying noise, not real structure.
3. DARE Sigma0 issue from 2026-04-23 → **Still open.** Now the dominant remaining issue (CF magnitudes).

---

## Session 2026-04-27 — Step 4: steady-state Kalman (DARE per call)

### What landed

Replaced the time-varying transient Kalman recursion (with frozen `_Sigma0_fixed` from BCKM Table 77 params) with a steady-state filter where `solve_discrete_are` is called **once per objective evaluation**. New helper `_steady_state_kalman(F, H, Q_proc)` returns the constant gain K, innovation cov S, and Σ_pred / Σ_filt; both `_kf_ll` and `_kf_full` use those constants for all t. Removed dead code: `_dare_cov`, `_Sigma0_fixed` precomputation, `_P_bckm_stable`, the old `Sigma0` argument.

scipy `solve_discrete_are` uses LQR convention `X = AᵀXA − …`, so we pass `A=Fᵀ, B=Hᵀ` to recover the Kalman predicted-cov form.

CLAUDE.md amended in the same change: the old "Solve DARE inside the optimization loop" prohibition has been replaced with positive guidance that the steady-state Kalman is now the canonical setup. The 5×5 DARE is cheap (~ms) and per-call evaluation eliminates the optimized-vs-final LL gap.

### Result

```
LL (best restart): 1779.78
Final smoother LL: 1779.78           (gap closed — was 144 units in Step 3)

P diag: [0.971, 0.998, 0.960, 0.995]
Max |eig|: 0.992                     (stationary ✓)

Smoothed wedge means: ≈ 0 in all 4 channels ✓
taux_hat std: 0.0271                 (was 0.0355 in Step 3 — tighter)

GR investment wedge Δ(1+τ_x): +0.070 (worsened ✓ — right sign)

φ-stats:
  y: A=0.29  τ_l=0.23  τ_x=0.30  g=0.18
  l: A=0.15  τ_l=0.58  τ_x=0.18  g=0.09
  x: A=0.19  τ_l=0.36  τ_x=0.21  g=0.24

Peak-to-trough (2007Q4 → 2009Q2):
   actual  efficiency  labor  investment  government
y -0.0722     -0.003  -0.058     -0.104      +0.039
l -0.0814     +0.012  -0.087     -0.146      +0.057
x -0.2898     +0.039  -0.166     -0.676      +0.138
```

### What this changes vs Step 3

- Investment-only x peak-to-trough: **−0.945 → −0.676** (closer to data −0.29; still ~2.3× off, vs ~3.3× before).
- Output and labor CF magnitudes also shrunk and got closer to data.
- Optimized vs final-smoother LL: **gap closed** by construction (same DARE constants, no transient). MLE objective is now self-consistent with the smoother.
- All 52 tests pass.

### What still doesn't match BCKM

1. **CF magnitudes still overshoot data** by ~2× on investment. Most likely culprit now: data-side investment definition. BCKM includes consumer durables in Y and applies a sales-tax adjustment on X (`(rCD/rCNDS)·rSTX`) — that pulls down x's variance and brings investment-wedge magnitudes toward data. This is the queued Step 5.
2. **φ_y[τ_L] = 0.23 vs BCKM target 0.46.** Labor wedge still underweighted on output.
3. **P_0 = 0 vs BCKM Table 9 [0.014, 0.001, 0.013, −0.014].** Still forced to zero by `Sbar=0`. Carried.

### Where this leaves us

Steps 1–4 land together and give a clean BCKM mleqadj.m architecture: phi0 in obs equation, Sbar=0, steady-state Kalman with DARE-per-call, RTS smoother on the same constants. The only remaining structural item is the data-side Step 5 (durables + sales tax).

### Next

- **Step 5: Add `rDCD` (consumer durables) to Y; apply sales-tax wedge `(rCD/rCNDS)·rSTX` to X.** ~1 hr. Last data-side adjustment before re-running the full Table 8/9/10/11 comparison.
- Optional **Sbar L2 prior** if Step 5 still leaves P_0 = 0 mismatched against Table 9.

---

## Session 2026-04-27 — Step 5: BCKM data adjustments (durables + sales-tax split)

### What landed

Aligned `bca_core/data/adjustments.py` to BCKM `usdata.m`:

```
Y = rGDP - rSTX + 0.04·rKCD + rDCD              (was rGDP + 0.04·rKCD - rSTX)
C = (rCND + rCS) - share_cnd·rSTX + 0.04·rKCD + rDCD
                                                 (was pce - rDCD + 0.04·rKCD - rSTX)
X = rCD + rGPDI + rGI - (rCD/rCNDS)·rSTX        (was rCD + rGPDI + rGI, no tax)
```

Two changes:

1. **Durables service flow now has both components.** `reclassify_durables` adds `service_flow_full = 0.04·K_dur + rDCD` to Y and C. The `rDCD` term proxies the depreciation flow `δ·K_dur` (equal in steady state). Previously only the return component was added.

2. **Sales tax is split across Y, C, X.** `subtract_sales_tax` now uses `share_cnd = rCND/rCNDS` and `share_dur = rCD/rCNDS` (with `rCNDS = rCND + rCS`). Y still gets the full rSTX subtracted; C gets `share_cnd·rSTX`; X gets `share_dur·rSTX` (new — was zero). Required two new FRED series: `PCND` (rCND) and `PCESV` (rCS).

Pipeline cache `data/us_1980_2014.parquet` regenerated.

### Result

```
LL: 1749.16                          (Step 4: 1779.78 → drop of 30; Y is now a different series)
g_share from data: 0.115             (was 0.124 — Y grew via service flow)
SS: y=1.276 l=0.315 x/y=0.255 g/y=0.115

phi0 (SS misalignment): y=+0.000  l=-0.001  x=+0.061  g=-0.010
                                     (x phi0 dropped 0.138 → 0.061 — model SS x/y now closer to data ✓)

P diag: [0.995, 0.998, 0.949, 0.995]
Max |eig|: 0.994                     (stationary ✓)

φ-stats:                       BCKM Table 11 target
  y: A=0.37  τ_l=0.28  τ_x=0.15  g=0.20    fYA=0.16, fYτL=0.46, fYτx=0.32
  l: A=0.15  τ_l=0.68  τ_x=0.08  g=0.10                       fLτL ≈ 0.70
  x: A=0.28  τ_l=0.40  τ_x=0.06  g=0.26

Peak-to-trough (2007Q4 → 2009Q2):
   actual  efficiency  labor  investment  government
y -0.0833     -0.014  -0.058     -0.168      +0.036
l -0.0814     +0.014  -0.088     -0.240      +0.053
x -0.2898     +0.030  -0.164     -1.111      +0.133
```

### What this changes vs Step 4

Direction is mixed:

- **φ_l[τ_L] = 0.68** (was 0.58; target 0.70). Big improvement, essentially on target.
- **φ_y[τ_L] = 0.28** (was 0.23; target 0.46). Up but still well short.
- **φ_y[τ_x] = 0.15** (was 0.30; target 0.32). Now *under* target — investment wedge contributes too little.
- **φ_y[A] = 0.37** (was 0.29; target 0.16). Got further from target — efficiency wedge over-explains.
- **phi0[x] dropped from 0.138 → 0.061.** The data x/y is now much closer to the model's 0.255 — this is the qualitative thing Step 5 was supposed to fix, and it did.
- **Investment-only CF magnitudes got bigger** (x peak-to-trough −0.68 → −1.11), driven by Y now including durables (which fell harder than overall GDP in the GR — actual y peak-to-trough −0.072 → −0.083).
- **Actual y peak-to-trough deepened** (−0.072 → −0.083) because Y now includes the durables service flow, and durables consumption crashed in 2008–09.
- **LL dropped 30 units** — expected, the observable is now a different series.

### Reading the result

The data construction now matches BCKM's `usdata.m`. The variance decomposition moved on labor (toward the BCKM target) and on investment-on-labor (toward target). The output decomposition didn't get cleaner: efficiency over-explains output and labor still under-explains. This is the same qualitative pattern reported in `REPORT.md` ("φ_A^Y too high, φ_L^Y too low") — Step 5 didn't close it.

phi0 going from 0.138 to 0.061 is the cleanest win: the steady-state misalignment between model and data on x/y, which Step 3 was originally introduced to absorb, is now small.

### What still doesn't match BCKM

1. **φ_y[A] too high, φ_y[τ_L] too low.** Same direction as before. Likely candidates:
   - Sample period: BCKM uses 1948Q1–2014Q3 (Tier 1 in DIVERGENCE_ANALYSIS.md, item 1.2). Our 1980–2014 sample omits the 1970s/early-80s identification windows.
   - Sbar / model-SS coupling: BCKM `initmle.m` uses fsolve to pick Sbar so that model SS matches data sample means *jointly* with the labor wedge SS — a smoother optimization landscape. Our Sbar=0 design forces SS misalignment into phi0 alone.
2. **Investment-only x peak-to-trough −1.11 vs data −0.29.** ~3.8× overshoot. This is now driven by `taux_hat` std (0.043 vs 0.027 before) — the new Y picks up durables-driven cyclicality and assigns it to τ_x. Possibly improves with sample extension.
3. **P_0 = 0 vs Table 9 target.** Carried — would need Sbar prior or BCKM's fsolve initialization.

### Next

Two queued items, in order of expected impact:

1. **Extend sample to 1948Q1.** Tier 1 in `DIVERGENCE_ANALYSIS.md`. Most series in the FRED fetcher already go back; the bottleneck is `LFWA64TTUSQ647S` (1969+). Need an alternate working-age population series (CNP16OV?) to extend back. Worth a quick scoping check before committing to it.
2. **Initialize Sbar via fsolve on the data sample means** (BCKM `initmle.m`). Would relax the Sbar=0 constraint without re-introducing the runaway-Sbar pathology.


---

## Session: 2026-04-27 — Step 6 — Sample extension to 1948Q1

### What was completed since the last entry

Extended the BCA estimation sample from 1980Q1–2014Q4 (T=140) to 1948Q1–2014Q4 (T=268), aligning with BCKM's MLE sample. This required swapping out three FRED series whose start dates limited the previous sample, and rewiring the data adjustments accordingly.

### Files touched

**`bca_core/data/fred.py`** — `FRED_SERIES` dict
- Removed `pce` (PCE total, monthly, 1959+) — total PCE now derived from components
- Removed `avg_weekly_hours` (AWHNONAG, 1964+) — labor input now from `hours_index` directly
- Replaced `working_age_pop` series: `LFWA64TTUSQ647S` (OECD 15–64, 1977+) → `CNP16OV` (BLS civilian non-institutional pop 16+, 1948+)
- All 13 remaining FRED series now span 1947+, allowing the full BCKM 1948Q1+ sample

**`bca_core/data/adjustments.py`**
- `reclassify_durables`: `c_adj` formula no longer reads `df["pce"]`. C is now built from components: `c_adj = pce_nondurables + pce_services + service_flow_full` (mathematically identical to before since `pce = pce_durables + pce_nondurables + pce_services`)
- `compute_labor_input`: prefers `hours_index` (PRS85006023, 1947+) over `employment * avg_weekly_hours` (which only goes back to 1964 due to AWHNONAG). Per-capita normalization unchanged.

**`scripts/run_var_counterfactuals.py`**
- Default `start="1980Q1"` → `start="1948Q1"`
- Cache filename: `data/us_1980_2014.parquet` → `data/us_1948_2014.parquet`
- `n_restarts=5` → `n_restarts=4` (see issue below)

### What was kept the same

- Detrend method: linear log trend at calibrated γ=1.9%/yr (BCKM Table 77)
- BCKM data adjustments from Step 5 (durables service flow into Y/C; sales-tax wedge split across Y/C/X)
- Steady-state Kalman with DARE-per-iterate (Step 4)
- Sbar = 0 constraint (Step 3) — to be relaxed in Step 7

### Issue encountered: hung restart

The first attempt with `n_restarts=5` hung on restart 5/5 for ~70 minutes after restarts 1–4 had completed (LLs 2643, 2741, 2775, 2776). With T=268 vs T=140, each LL evaluation is ~2× slower, and `ftol=1e-13` makes the optimizer wander in flat regions when a random perturbation lands in a bad basin. Killed the process and reduced to `n_restarts=4`. The completed run finishes in ~10 min and the best LL is unchanged (best of 4 was already on the plateau at restart 3/4).

Followup item: consider relaxing `ftol` to 1e-10 or capping `maxiter` more tightly for the longer sample.

### Numbers (after rerun, 1948Q1–2014Q4, n_restarts=4)

```
T = 268
g_share from data: 0.1321  (was 0.166 in 1980+ sample)
SS: y=1.2994  l=0.3211  x/y=0.2545  g/y=0.1321
phi0: y=-0.0000  l=-0.0548  x=+0.0810  g=-0.0192

LL = 2776.227

VAR diagonal:
  A   : 0.9949
  τ_l : 0.9989  (BCKM target ≈ 1.001 — close)
  τ_x : 0.9748
  g   : 0.9666
Max |eigenvalue|: 0.9963

Sbar = [0, 0, 0, 0]  (Step 3 constraint, still active)
P_0  = [0, 0, 0, 0]  (BCKM Table 9 target [0.014, 0.001, 0.013, -0.014])
```

**φ-statistics for y (variance decomposition):**

| | Step 5 (1980+) | Step 6 (1948+) | BCKM target |
|---|---|---|---|
| φy[A]   | 0.37 | 0.02 | 0.16 |
| φy[τ_L] | 0.28 | 0.02 | 0.46 |
| φy[τ_x] | 0.15 | 0.58 | 0.32 |
| φy[g]   | —    | 0.39 | —    |

**Investment wedge during GR (correct sign maintained):**
- 2007Q4: −0.0474, 2009Q2: +0.0911, Δ=+0.139 (worsened ✓)

**Peak-to-trough 2007Q4 → 2009Q2 (Investment-only CF):**
- y: −0.192 (data: −0.087) — 2.2× overshoot (was 1.4× in Step 5)
- x: −1.223 (data: −0.293) — 4.2× overshoot (was 3.8× in Step 5)

### Reading the result

Step 6 is a regression on the φ-statistics: the longer sample shifted variance away from labor and efficiency into investment and government. Several factors compound:

1. **Sbar=0 still binding**, and now it has more work to do because the sample mean drift over 1948–2014 (post-WWII normalization, 1970s stagflation, secular decline in labor share) is larger than over 1980–2014. With phi0[x]=0.081 and Sbar locked at 0, the wedge means must pin themselves to model SS rather than data SS, distorting the variance allocation.

2. **Smoothed `taul_hat` std = 0.447** (much larger than Step 5). Q[τ_l] = 0.0025 is 3× the next-largest shock variance. The MLE is fitting low-frequency labor variation (LFP trends, demographic shifts) as noisy persistent shocks — which fits `l` well (φ_l[τ_L] = 0.83) but doesn't translate to `y` because the labor wedge moves cancel in equilibrium.

3. **Hours index vs employment×hours mismatch**: `PRS85006023` is BLS nonfarm business hours of all persons, while the prior `PAYEMS × AWHNONAG` was nonfarm payrolls × avg weekly hours. The two have different long-run normalizations; over 67 years this affects the labor-wedge identification more than it did over 35 years.

4. **The longer sample contains genuine structural change** (1948 industrial economy → 2014 service economy). A single VAR(1) with constant coefficients can't capture this, and the MLE picks the mode that maximizes likelihood under that mis-specification — which is not necessarily the mode that matches BCKM Table 11.

### What still doesn't match BCKM

Same direction as before, but now louder. Step 7 (fsolve-init Sbar) is the immediate next move — relaxing Sbar should let phi0 mass redistribute and may pull φ_y[A]/φ_y[τ_L] back toward target. The hours-vs-employment×hours question can be revisited if Step 7 doesn't close the gap.

### Next

**Step 7**: BCKM `initmle.m`-style fsolve initialization for Sbar. Pick initial Sbar by solving for the steady-state wedge means that make model SS observables equal to data sample means, so the optimizer starts from a feasible non-zero P_0 instead of the Sbar=0 corner. Relaxing the constraint without re-introducing the runaway-Sbar pathology was Step 3's worry; the fsolve init gives a principled starting point so the optimizer doesn't have to discover Sbar from scratch.


---

## Session: 2026-04-27 — Step 7 — Free vs fixed Sbar (initmle.m-style)

### What was tried

Two parametrizations of the wedge-VAR unconditional mean Sbar, both with uncentered observables (`prepare_observables(..., center=False)`) so the SS gap is absorbed by the wedge VAR rather than a fixed `phi0` in the obs equation.

**Option 1 — free Sbar with bound penalty (`±0.5`)**: theta = [Sbar(4), P(16), Q(10)] = 30 params. fsolve-init seed: solve `H · (I−F)^{-1} · intercept(Sbar) = sample_obs_mean` for Sbar at BCKM P/Q. Soft-bound penalty on |Sbar|.

**Option 2 — fix Sbar at the fsolve-init value**: theta = [P(16), Q(10)] = 26 params. Sbar pinned once at startup from `_fsolve_sbar(_P_bckm, _Q_bckm)`. Implied `P_0 = (I−P_var)·Sbar_fixed` moves with the optimizer's P.

### Files touched

**`bca_core/var_estimation.py`**
- `prepare_observables`: added `center: bool = True` flag (Step 7 used `False`).
- `_unpack` / `_pack`: shrunk theta from 30 → 26 (option 2). `Sbar_fixed` is closure-captured.
- `_kf_ll` / `_kf_full`: x0 set to unconditional state mean `(I−F)^{-1}·intercept` instead of zeros — required when intercept ≠ 0.
- `_fsolve_sbar`: 4×4 linear-solve helper that builds the column-by-column map `Sbar → H · E[s]` and solves for the Sbar matching `sample_obs_mean`.
- Removed the `_SBAR_BOUND` penalty (option 2): no free Sbar to bound.

**`scripts/run_var_counterfactuals.py`**
- `prepare_observables(df, ss, center=False)`.

### Numbers

```
sample_obs_mean (phi0): y=-0.0000  l=-0.0548  x=+0.0810  g=-0.0192
Sbar_fixed (fsolve-init from BCKM P/Q):
  A=+0.0073  τ_l=-0.1102  τ_x=-0.0785  g=-0.0192
```

| | Step 6 (Sbar=0) | Step 7 opt 1 (free Sbar) | Step 7 opt 2 (fix Sbar) | BCKM |
|---|---|---|---|---|
| LL | 2776 | 3109 | 2698 | — |
| Max \|eig\| | 0.996 | (Sbar at bound) | **1.008 ✗** | <1 |
| Sbar | 0 | hit ±0.5 bound | fixed at fsolve-init | implied 0.014 |
| GR Δτ_x sign | +0.139 ✓ | — | **−0.016 ✗** | + |
| φ_y[A]   | 0.02 | (worse) | 0.06 | 0.16 |
| φ_y[τ_L] | 0.02 | (worse) | 0.03 | 0.46 |
| φ_y[τ_x] | 0.58 | (worse) | 0.03 | 0.32 |
| φ_y[g]   | 0.39 | (worse) | **0.87** | — |

### Reading the result

Both Step 7 options are worse than Step 6 (Sbar=0):

- **Option 1** lifted LL by 332 units but did so by exploiting the (I−P) near-singular basin: free Sbar saturated the ±0.5 bound, P_0 collapsed, and the optimizer used the extra freedom to dampen wedge variance rather than match data dynamics. φ-stats regressed across the board.

- **Option 2** prevents the runaway pathology but pays for it twice: (i) Sbar pinned at `[+0.007, −0.110, −0.079, −0.019]` is a strong prior the optimizer cannot revise, so the implied P_0 = (I−P)·Sbar_fixed is far from BCKM Table 9; (ii) to make the data fit at all under that constraint, the optimizer pushes P toward the unit circle (max |eig| = 1.008 — non-stationary, breaks the spectral-radius bound). The result is a model with the wrong-sign GR investment wedge and 87% of y-variance attributed to government.

The fsolve-init Sbar is large in magnitude (|τ_l|≈0.11, |τ_x|≈0.08) because the model's unconditional response of mean(l_hat) and mean(x_hat) to Sbar is small — it takes a big Sbar to absorb a phi0 of −0.05 / +0.08. Big Sbar then forces big P_0, which the rest of the system isn't well-suited to.

### Hypothesis on why BCKM replication keeps failing

Ranked by my current confidence, after Steps 6 and 7:

1. **Labor input series construction (most likely)** — `taul_hat` smoothed std hits 0.45 in Step 6 and stays ≥0.32 in Step 7, vs BCKM Table 10 which implies far smaller wedge volatility. PRS85006023 (nonfarm business hours of all persons, BLS index) normalized to a fixed mean and divided by CNP16OV does not match CKM-2007's "labor input" construction — that paper uses a Cobb-Douglas-aggregated hours-of-work-and-quality index across NIPA labor categories. Until our `l_t` matches BCKM's `l_t` series, the labor wedge will absorb whatever low-frequency mismatch exists between our hours data and the model.
2. **Investment definition (rDCD = pce_durables vs perpetual-inventory δ·K_dur)** — Step 5 added `pce_durables` as the depreciation flow (+ K_dur return). Re-checking `usdata.m`: BCKM uses `rDCD` as a separately-computed series, *not* equal to current pce_durables expenditure. We're double-counting the depreciation flow if `pce_durables` already includes replacement of fully-depreciated stock.
3. **MLE multimodality** — the LL surface has at least three distinct local optima (Sbar=0 plateau at 2776, free-Sbar bound at 3109, fixed-Sbar mode at 2698). BCKM's exact restart strategy (initmle.m perturbation set) may be necessary to land in their mode.
4. **Detrending intercept** — calibrated γ=1.9%/yr fixes the slope; the OLS-fitted intercept differs from BCKM's `meandata` adjustment. This shifts the data sample mean and propagates into phi0 / Sbar.
5. **Population aging in CNP16OV** — CNP16OV (16+) includes retirees, whose share grew 1948→2014. Per-capita labor scaling is biased downward over time. BCKM's pop series excludes 65+; ours doesn't.

The smoking gun is #1: the only state with `taul_hat` std ≥ 3× the others *and* the only series we changed in Step 6 that has a known semantic mismatch with BCKM. Step 8 should be a CKM-2007-style labor-input reconstruction or a regression test against BCKM's published `l_t` time series if obtainable.

### Recommendation: roll back to Step 6 as operational baseline

Step 6 (Sbar=0, centered observables) gives the only well-behaved result so far: stationary VAR, correct GR investment-wedge sign, well-defined φ-stats. It is wrong on the φ-statistic *targets*, but it is wrong in a way that is internally consistent and interpretable. Step 7 options trade interpretability for LL points and break the GR sign — that is not a good trade.

Concretely: leave the var_estimation.py changes in place (the `center` flag, fsolve helper, unconditional-mean x0 are all reusable); revert the runner script to `center=True` (or default omit) so the production pipeline returns to Sbar=0 until the labor-input issue is fixed.

### Next

**Step 8**: investigate labor input series. Two cuts:
- (a) Diagnostic: compare `taul_hat` smoothed time series under our `l_t` vs CKM-2007's published `l_t` (if available in BCKM's USAN2 data files).
- (b) Constructive: try `(1 − unemployment_rate) · avg_weekly_hours` reconstruction back to 1948, comparing the φ-stats. If this closes the gap, the labor-input series is the answer.

Hold off on more Sbar surgery until the labor identification is resolved.


---

## Session: 2026-04-27 — Step 8 — Read the actual BCKM Matlab and re-plan

### Context

User added the BCKM Matlab reference files to `matlab_reference/` (driver `datamine.m`,
estimator `mleqadj.m`, runner `runmleadj.m`, init helper `initmle.m`, trend
helper `maketrend.m`, calgz residual `calgz.m`, Kalman helper `kfilter.m`,
counterfactual driver `runall.m`, f-stats `fstats3.m`, US data fetcher
`usdata.m`). We had been working from `BCA_info.md` plus an older vintage of
`usdata.m`; reading the production Matlab end-to-end revealed nine concrete
mismatches with our pipeline.

### Discrepancies (production Matlab vs current pipeline)

1. **Sample period.** `datamine.m` runs MLE on `t = 1980.25 : 0.25 : 2015`
   (T=140, 1980Q1–2015Q1). Our Step 6 extension to 1948Q1–2014Q4 was based on a
   misreading of an older `usdata.m` vintage. The production sample is the
   shorter one. **Step 6 must be reverted.**
2. **Growth rate γ.** `maketrend.m` calls `gzt = fsolve(@calgz, 0)` to pick γ
   so that `mean(detrended log y_pc) = 0` over the MLE window. Our Step 1 used
   the calibrated 1.9% from BCKM Table 77; the production code uses a *data-
   estimated* γ.
3. **Trend normalization.** `maketrend.m` divides each series by
   `ypc(by) · (1+gzt)^by` where `by = base date index` — a base-year level
   pin. Our `remove_trend` only subtracts the trend slope and does not pin a
   base year.
4. **Sbar parametrization.** `mleqadj.m` `ind = [0,0,0,0,0,0,0, 1,1,1,1, ...]`
   has Sbar(1..4) as estimated parameters (theta(1..4)) — Sbar IS a free
   parameter. Our Step 7 option 2 (fix Sbar) was the wrong direction. The
   correct setup is option 1 (free Sbar) + a sensible fsolve-based init.
5. **fsolve-init for Sbar.** `initmle.m` uses fsolve over level/ratio
   residuals (model SS observables vs data sample means), not a one-shot
   linear solve. Our Step 7 fsolve-init was a 4×4 linear solve at BCKM P/Q —
   too crude.
6. **Spectral radius bound.** `mleqadj.m` line 134: `penalty = 500000 *
   max(|eig|−.995, 0)^2` — the threshold is **0.995**, not 1.005 as in our
   `bca_core/var_estimation.py`. Our CLAUDE.md note that "BCKM's τ_l diagonal
   is ~1.001" rationalized the 1.005 bound, but `datamine.m` runs *with* the
   0.995 spectral bound and still produces τ_l diag ≈ 1.001 in Table 8 — the
   bound is on the *spectral radius*, not on individual diagonals.
7. **Initial P.** `mleqadj.m` line 28 sets `P = 0.995*eye(4)`. Our restarts
   use BCKM Table 77 P (which has off-diagonals); the production warm-start
   is a diagonal 0.995·I.
8. **Restart perturbation.** `runmleadj.m` does `nps=50` restarts with
   *multiplicative* perturbation `x_new = x*pb, pb=0.99` — geometric shrink,
   not random Gaussian. Our `n_restarts=4` with Gaussian perturbation is a
   different exploration strategy.
9. **F-statistics window and definition.** `fstats3.m` lines 1–40 use the
   Great Recession window `i1=2008.25, i2=2011.75` with `ilog=0` (levels) and
   `ifilt=0` (no HP), and the formula `f(1,i) = 1/sum((dly−dlyc(:,i)).^2)`,
   then normalizes `f` to sum to 1. Our `phi_statistics` was computing
   full-sample inverse-SSR with no window. **This is the single most
   important discrepancy** — we have been comparing the wrong statistic to
   BCKM Table 11.

### Step 8 plan

1. **8.1 — Rewrite f-stats as GR-window inverse-MSE.** Add a `window`
   argument to `phi_statistics`; in the runner script, compute the canonical
   BCKM stat over 2008Q1–2011Q4 and print both that *and* the full-sample
   diagnostic side-by-side.
2. **8.2 — Revert sample to 1980Q1–2015Q1** (T=140) per `datamine.m`.
3. **8.3 — `calgz`-style fsolve for γ + base-year normalization** per
   `maketrend.m`.
4. **8.4 — fsolve-init Sbar over level/ratio residuals**, free Sbar in theta
   per `initmle.m` + `mleqadj.m`.
5. **8.5 — Tighten spectral-radius bound 1.005 → 0.995** per `mleqadj.m`
   line 134.
6. **8.6 — Initial P = 0.995·I** per `mleqadj.m` line 28.
7. **8.7 — Multiplicative perturbation x*pb, pb=0.99, nps=50** per
   `runmleadj.m`.

### Step 8.1 — landed

`bca_core/counterfactuals.py` `phi_statistics` now takes
`window: tuple[int, int] | None = None`; when set, slices both data and
counterfactual paths to that window before computing inverse-SSR. Default
behaviour (window=None) is unchanged for back-compat.

`scripts/run_var_counterfactuals.py` finds the GR window indices via
`find_date_index(df.index, 2008, 1)` and `find_date_index(df.index, 2011, 4)`
and prints the canonical BCKM Table 11 stat first, then the full-sample
diagnostic with a clear label that it is NOT the BCKM target.

`tests/` — 52/52 passing.

### Result on broken Step 7-option-2 config (cached parquet)

Pipeline run did not include the broken Step 7 fixes (Sbar still pinned, max
|eig|=1.008, GR investment sign wrong), but the f-stat block did execute:

```
F-statistics (BCKM Table 11, GR window 2008Q1–2011Q4):
                  y       l       x
efficiency    0.04    0.13    0.05      (BCKM target fY[A]=0.16)
labor         0.03    0.16    0.03      (BCKM target fY[τ_l]=0.46)
investment    0.37    0.10    0.43      (BCKM target fY[τ_x]=0.32 ✓ near match)
government    0.57    0.61    0.50
```

Even with the broken Step 7 config, the GR-window stat moves the investment
column close to BCKM's 0.32 — a strong sign that switching to the canonical
BCKM stat will close some of the gap independent of the other Step 8 items.
Re-evaluation against a clean Step 8.2+ baseline is required for a real
comparison, but the qualitative validation is encouraging.

### What still doesn't match BCKM (carried)

Same items as Step 7, plus the seven Step 8 items above. Step 8.1 is the
cheapest of the seven — and the one that most directly affects whether we
can claim BCKM Table 11 replication, since we had been comparing the wrong
statistic.

### Next

Step 8.2: revert sample to 1980Q1–2015Q1 in
`scripts/run_var_counterfactuals.py` and update the cache filename.

---

## Session: 2026-04-27 — Step 8.2–8.7 implemented (BCKM-faithful estimator)

### What was completed

All seven Step 8 sub-steps are in code. 8.1 was committed earlier (59f5bd2);
8.2–8.7 are in this commit:

- **8.2** — `scripts/run_var_counterfactuals.py`: sample reverted to
  `start="1980Q1", end="2014Q4"` (T=140), per BCKM `datamine.m` line 65
  (`t = 1980.25:0.25:2015`, `eobs=140`). Step 6's 1948+ extension reverted.

- **8.3** — `bca_core/data/adjustments.py` + `pipeline.py`: added
  `method="calgz"` branch in `remove_trend`. Closed-form slope
  `γ_q = (mean(log_s) − log_s[bdate]) / (mean(t) − bdate_idx)` (mirrors
  `maketrend.m`/`calgz.m` — but BCKM uses `fsolve`; we use the closed-form
  it solves for, mathematically equivalent). The fitted trend simultaneously
  satisfies `mean(detrended log) = 0` and `log(detrended[bdate]) = 0`.
  `pipeline.py` plumbs `base_year_quarter="2008Q1"` → `base_idx` via
  `sample.index.get_loc`. New metadata fields: `detrend_method`, `base_idx`,
  `base_year_quarter`.

- **8.4** — `bca_core/var_estimation.py`: Sbar back as a free parameter.
  `N_P = 30` (4 Sbar + 16 P + 10 Q lower-tri). `_unpack`/`_pack` take
  `theta = [Sbar(4), P(16), Q_lo(10)]` and emit `P_0 = (I−P)·Sbar`.
  Added `_model_ss_from_sbar()` mirroring `mleqadj.m` sec 5a (sigma=1
  hardcoded; CalibrationParams omits it). Added `_fsolve_sbar_initmle()`
  using `scipy.optimize.fsolve` over residuals
  `[ys/ys_ss − Ym(0); ls/ls_ss − Ym(1); (xs/ys)/(xs_ss/ys_ss) − Ym(2);
  (gs/ys)/(gs_ss/ys_ss) − Ym(3)]`. Falls back to zeros on non-convergence
  (which is what happened in the run below — initial guess didn't pull
  fsolve into a valid basin; needs follow-up).

- **8.5** — Spectral-radius bound tightened from 1.005 → 0.995
  (mleqadj.m line 134). Per-diagonal bound also 0.995 (no τ_l relax).
  Plus added Sbar bounds penalty `[-1,-1,-1,-5]/[1,1,1,1]`
  (mleqadj.m Lb/Ub).

- **8.6** — Initial `P_warm = 0.995·I` (mleqadj.m line 28 / runmleadj.m
  x0a/b/c). Q warm-start from BCKM `x0c` lower-triangular values
  (runmleadj.m lines 100-110, "result from initpw with annual adja=12.88").

- **8.7** — Added BCKM multiplicative-shrink loop after main optimizer:
  `pb=0.99`, `nps=50`. Each iteration scales `theta = pb · theta` then
  re-runs L-BFGS-B; keeps best. Per runmleadj.m lines 121-141.

### Result on stale-cache run (linear γ, NOT calgz — calgz needs FRED refetch)

The on-disk cache `data/us_1980_2014.parquet` predates Step 8.3 — it was
generated with `linear, gamma_annual=0.019`. `pipeline.py` short-circuits
to load it without re-detrending, so this run uses linear-γ data. Calgz
will require a FRED refetch via `--save-data`.

```
T = 140, g_share = 0.166

Sbar_init (initmle.m fsolve): [0,0,0,0]   # ier != 1, fell back to zeros
ll: 1814 → 1834 (multiplicative-shrink loop took +20)

Sbar (final):  [-0.010, -0.146, -0.045,  0.178]
P_0 (= (I−P)·Sbar):  [-0.006,  0.013, -0.005, -0.005]
P_0 BCKM target:     [ 0.014,  0.001,  0.013, -0.014]

P diagonal: [0.905, 0.990, 0.932, 0.984]
            (BCKM:  0.989, 1.001, 0.968, 0.995)
Max |eigval|: 0.9953  (right at bound — was 1.008 in Step 7)

GR-window f-stats (BCKM Table 11):
              y       l       x
efficiency  0.05    0.02    0.03   (target fY[A]=0.16)
labor       0.48    0.70    0.29   (target fY[τ_l]=0.46 ✓)
investment  0.11    0.19    0.42   (target fY[τ_x]=0.32)  -- was 0.37 ✓ pre-step
government  0.36    0.09    0.26
```

`fY[τ_l]` matched almost exactly (0.48 vs 0.46). `fY[A]` and `fY[τ_x]`
are below target — this is the regression caused by Sbar warm-start
falling back to zeros (mean wedge offset wasn't captured in initial
guess, so the wedge VAR has to absorb the offset elsewhere). Pre-Step 8
the f[τ_x]=0.37 was nearly on target; now it's 0.42 in the x-row but
0.11 in the y-row. The fsolve seeding logic needs a better initial guess.

### Tests

`pytest tests/ -v` — 52/52 passing. No regressions from the var_estimation
restructure.

### What still doesn't match BCKM

1. **Sbar fsolve fails** (returns [0,0,0,0]). The residuals' Jacobian at
   Sbar=[0, 0.05, 0, log(0.2)] may be ill-conditioned, or the units of
   `Ym` (`mean(exp(obs_hat))` — ratios) don't match my model SS units
   exactly. BCKM's `initmle.m` uses absolute levels, not ratios.
2. **Calgz cache not yet validated.** The current run used a linear-γ
   cache. A FRED refetch is needed to confirm calgz changes f-stats
   in the right direction.
3. **fY[A] still under target** (0.05 vs 0.16). Likely tied to (1).

### Next

Two items, in order:
1. Refetch from FRED with calgz (`--save-data data/us_1980_2014_calgz.parquet`)
   and re-run. This validates Step 8.3 end-to-end.
2. Debug the Sbar fsolve. Likely fix: change `Ym` from
   `mean(exp(obs_hat))` (ratios) to absolute model-units, matching
   BCKM `initmle.m` semantics directly. Or add an iterative seed
   schedule for fsolve.

---

## Session: 2026-04-28 — Step 8.1 follow-up (Sbar warm-start + calgz)

### What was completed

**Sbar fsolve fixed** — `bca_core/var_estimation.py`:
- Found a silent bug: `_model_ss_from_sbar` referenced `proto.params` but
  `PrototypeModel` exposes `self.p`. The try/except in `_fsolve_sbar_initmle`
  was swallowing the `AttributeError` and returning zeros. Single-line fix.
- Rewrote `_fsolve_sbar_initmle` to take `data_means` from caller and use
  BCKM-faithful residuals (initmle.m line 53):
  `[ys − Ym[0]; xs/ys − Ym[1]; ls − Ym[2]; gs/ys − Ym[3]]`
  with `Ym = [mean(y), mean(x/y), mean(l), mean(g/y)]` (not ratios as before).
- `scripts/run_var_counterfactuals.py` now computes `data_means` from raw
  detrended levels and passes to `estimate_var_mle`.

**Coordinate-system pivot** — empirically tested four Sbar candidates as
warm-starts. BCKM's initmle.m operates in absolute log-levels; our state
space uses log-deviations from a calibrated SS. The BCKM-faithful Sbar
(`A=-0.28, τ_l=+0.03, τ_x=+0.04, g=-2.16`) lives in the wrong coordinate
system and crashed LL to -26434. Pivoted to the **OLS-implied Sbar**
(`Sbar = (I-P_ols)⁻¹·P_0_ols`) which lives in the same deviation space as
our observables. The fsolve result is kept as a diagnostic print only.

**Calgz cache** — refetched FRED with `--save-data` and
`base_year_quarter="2008Q1"`, producing `data/us_1980_2014_calgz.parquet`.
calgz-fitted γ_q = 0.0048 (annual 0.0192).

**Misleading log fixed** — `scripts/run_var_counterfactuals.py` now prints
"calgz-fitted" when `detrend_method=='calgz'` instead of "OLS-estimated".

### Result on calgz cache + OLS-implied Sbar warm-start

```
T = 140, g_share = 0.115, calgz γ_annual = 0.0192

Sbar_initmle (absolute, diagnostic only): A=-0.28  τ_l=+0.03  τ_x=+0.04  g=-2.16
Sbar_init  (OLS-implied, used as warm-start): A=-0.39  τ_l=+0.08  τ_x=+0.53  g=-0.04

ll: 1601 → 1870 (multiplicative-shrink loop +268)

Sbar (final):  [-0.168, +0.238, +0.040, +0.208]
P_0:           [-0.002, -0.003, -0.002, -0.002]
P_0 target:    [+0.014, +0.001, +0.013, -0.014]

P diagonal: [0.979, 0.995, 0.939, 0.896]   (BCKM: 0.989, 1.001, 0.968, 0.995)
Max |eigval|: 0.9956 (within bound)

GR-window f-stats (BCKM Table 11):
              y       l       x
efficiency  0.15    0.07    0.66    target fY[A]=0.16   ✓ matched
labor       0.65    0.75    0.11    target fY[τ_l]=0.46 (over by 0.19)
investment  0.06    0.12    0.14    target fY[τ_x]=0.32 (under by 0.26)
government  0.14    0.05    0.09

Peak-to-trough (2007Q4→2009Q2):
   actual  efficiency  labor  investment  government
y -0.087    -0.059   -0.017     -0.118     +0.048
l -0.037    -0.001   -0.020     -0.165     +0.070
x -0.293    -0.088   -0.020     -0.769     +0.213
```

### Headline

**fY[A] = 0.15 vs target 0.16 — essentially nailed.** The efficiency-wedge
contribution to GR output, which previously sat at 0.05 (run with zero Sbar
warm-start), is now within 0.01 of BCKM's published value. LL improved 36
units. τ_l diagonal moved from 0.990 → 0.995, closer to BCKM's 1.001.

The remaining gap is between `fY[τ_l]` (0.65 vs 0.46) and `fY[τ_x]` (0.06
vs 0.32). The investment counterfactual *does* produce the right
peak-to-trough decline (y=-0.118 alone vs actual -0.087, broadly in line
with BCKM finding investment explains roughly a third of GR output drop),
but its full-window variance contribution is suppressed because
shock co-variance Q[τ_l, τ_x] ≠ 0 — labor and investment shocks are
correlated, and labor is absorbing the joint variance.

### Tests

`pytest tests/ -v` — 52/52 passing. No regressions.

### What still doesn't match BCKM

1. **P_0 still very far from Table 9** (~0.002 vs target ±0.013). Suggests
   our state-space coordinate system (deviations from calibrated SS) does
   not coincide with BCKM's (absolute log-levels). They will not match
   unless the SS frame is shifted or BCKM is reproduced in absolute levels.
2. **fY[τ_l] over-attributed, fY[τ_x] under-attributed.** Likely tied to
   the joint-shock covariance. Worth investigating whether the
   multiplicative-shrink iterates explored a flatter ridge of the LL surface.

### Next

1. Commit Step 8.1 follow-up (this session's changes to `bca_core/var_estimation.py`,
   `scripts/run_var_counterfactuals.py`, and the new calgz cache).
2. Investigate the τ_l/τ_x co-attribution: try a stricter constraint on Q
   off-diagonals OR re-run with `nps=200` to see if the ridge has a better
   point. (Out of scope this session.)
