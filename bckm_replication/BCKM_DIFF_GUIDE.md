---
title: "BCKM ↔ Our Pipeline — Data Manipulation Reference"
topic: "reference"
layer: "bckm-replication"
status: "reference"
last_updated: "2026-05-03"
---

# BCKM ↔ Our Pipeline — Data Manipulation Reference

When you cross-machine compare our pipeline against BCKM's `worktemp.mat`,
the two stacks use different conventions in places where neither is wrong
— they're just different. This guide catalogs every conversion you need
so the diff is apples-to-apples.

If you add a new diff target, update this file. If you skip a conversion
and get nonsense, the bug is almost always in this list.

## Source of truth

- **BCKM**: `octave_output/bckm_*.csv` and `octave_output/bckm_dump.mat`,
  produced by [matlab_reference/dump_worktemp_for_diff.m](matlab_reference/dump_worktemp_for_diff.m)
  on the Octave-equipped machine. Dumps the raw obs, the detrended `mled`
  matrix, and the full converged MLE solution under `worktemp.mle.{...}`.
- **Ours**: `data/us_1980_2014_calgz.parquet` for observables;
  pipeline outputs from `scripts/eval_bckm_fstats.py` for likelihoods,
  innovations, and counterfactual paths.

## Constants you must reuse

```
gz         = 0.004725592524     # quarterly trend, worktemp.params(2)
gn         = 0.0024511032       # quarterly population, worktemp.params(1)
base_year  = 2008.25            # 2008Q1 — pin point for calgz detrending
T          = 140                # quarters, 1980Q2..2015Q1
```

The base year ties to `maketrend.m`'s `by` index and to our calgz
detrender's anchor; both agree at 2008Q1 in the current pipeline. If
either side changes its base, every level number in the diff slides.

## Manipulations required to compare observables

### 1. BCKM's `mled` is rescaled levels, not detrended levels

`maketrend.m` builds each non-labor column of `mled` as
`series_pc / ypc(by) * (1 + gz)^by`. To recover detrended levels, divide
by `(1+gz)^t` where `t` is the calendar quarter index. This is the same
operation that `mleqadj.m:237-238` does inside the Kalman feed:
`Y = log(mled[:,2:5]) - log((1+gz)^t)`.

```python
T = len(bckm_raw)
trend = (1 + gz) ** np.arange(T)
for c in ["y", "x", "g", "c_real"]:
    bckm[c] = bckm_raw[c].values / trend
```

### 2. BCKM normalizes per series; we preserve the `c+x+g=y` identity

After detrending, BCKM rebases **each series independently** so its
value at the base year equals 1. This breaks the closed-economy identity
`c + x + g = y` by ~2.2 in level (each series gets its own scale).

We instead divide *all* series by the same `y` trend factor (calgz
convention), which preserves the identity. Both are correct internally;
they just cannot be compared without a common rescale.

To match BCKM's convention, do the same per-series rescale on our side:

```python
ours = ours_raw.copy()
for c in ["y", "x", "g", "c", "l"]:
    v = ours_raw[c].values
    ours[c] = v / v[base_idx]   # base_idx = quarter index of 2008Q1
```

After this step, the `(c+x+g)−y` residual will be **the same nonzero
value on both sides** (~−2.2). That is the diagnostic that the
normalization mismatch is fully accounted for. If ours sums to ~0 and
theirs sums to ~−2.2, you forgot this step.

### 3. Labor (`l` / `hpc`) carries no trend — leave it alone

`mled` column 4 is raw `hpc` (hours per capita ≈ 0.243), with no trend
or rebase applied. Our parquet's `l` column is rescaled to model SS so
its sample mean matches BCKM's published `l_ss`. **Do not** apply the
base-year normalization from step 2 to labor — it produces garbage
because the two sides start in different magnitudes (0.24 vs 1.0).

For labor, the meaningful diagnostic is correlation, not level
agreement. Mean(l) on both sides should already be ~0.243 by
construction (we rescale to that target on import).

## Manipulations required to compare the MLE solution

### 4. BCKM `Likelihood` sign convention

`worktemp.mle.Likelihood` stores the **value of the function being
minimized**, which is the standard Gaussian negative log-likelihood
**without** the `0.5·T·p·log(2π)` constant:

```
L_BCKM = 0.5·T·log|Ω| + 0.5·Σ_t innov_t' Ω⁻¹ innov_t  +  0.5·penalty
```

(See `mleqadj.m:257`.) Our pipeline reports a positive log-likelihood
that *includes* the `2π` constant. To compare:

```
LL_ours_full = -L_BCKM + 0.5·T·p·log(2π)
            ≈ -L_BCKM + 0.5·139·4·log(2π)
            ≈ -L_BCKM + 510.95
```

So `mle.Likelihood = -2402.876` ⇒ expected `LL_ours ≈ +2913.83` if our
pipeline matches BCKM exactly. Always strip or add the constant before
comparing.

### 5. BCKM uses 4 observables, not 5

`mleqadj.m:13` reads `mled(:, 2:5)` which is `[y, x, h, g]`. The 5-column
`worktemp.obs` (and the 5-column `mle.obs`) is just the raw input
preserved in the struct — it is **not** the Kalman input. Don't compare
against the 5-col version. Our pipeline's `H = np.zeros((4, 5))` already
matches the 4-observable layout.

### 6. Observable order: BCKM is `[y, x, h, g]`, ours is `[y, l, x, g]`

Different permutation of the same series. To compare innovations,
covariances, or any per-row quantity, permute one side:

```python
# BCKM rows in our order:   bckm[:, [0, 2, 1, 3]]   # y,x,h,g → y,h,x,g (=l)
```

(`l` and `h` are the same series — labor — different letter.)

### 7. State dimensionality: BCKM 6D, ours 5D

BCKM augments the state with a constant `1` so the SS intercept lives in
the last column of `C` (`mleqadj.m:231-232`). Ours uses a 5D state and
folds the SS gap into a separate `obs_offset` vector. Mathematically
identical; only matters when you read `worktemp.mle.P` (4×4 wedges) or
`mle.P0` (4×1) — those are wedge-block quantities, comparable directly
to ours, the augmented "1" never appears in P or P0.

### 8. Likelihood bookkeeping: BCKM drops the first observation

`mleqadj.m:239-240` slices `Y(2:T)` and decrements `T = T-1`, so the
likelihood is a sum over **139** innovations, not 140. Our pipeline uses
all 140. Per-quarter contribution ≈ 17 nats, so this accounts for ~17
nats of LL gap, not for thousands.

## Reference scripts

- [scripts/diff_bckm_data.py](scripts/diff_bckm_data.py) — observables diff (steps 1–3)
- [scripts/eval_bckm_fstats.py](scripts/eval_bckm_fstats.py) — likelihood + f-stat diff at BCKM θ (steps 4–8)
- [matlab_reference/dump_worktemp_for_diff.m](matlab_reference/dump_worktemp_for_diff.m) — dumper that produces all `bckm_*` files

## Adding a new diff target

When you build a new diagnostic, walk this checklist:

1. Trended? Apply step 1.
2. Comparing levels across series? Apply step 2.
3. Touching labor? Don't apply step 2 to it (step 3).
4. Comparing likelihoods? Reconcile sign + 2π constant (step 4).
5. Comparing innovations / Σ / Ω rows? Permute observable order (step 6).
6. Comparing state vectors / P matrices? Confirm 4-block vs 6-block alignment (step 7).
7. Counting observations? Match T or correct for the off-by-one (step 8).
