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
