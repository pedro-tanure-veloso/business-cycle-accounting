---
title: "BCKM Computational Strategy"
topic: "reference"
layer: "bckm-replication"
status: "reference"
last_updated: "2026-05-03"
---

# BCKM Computational Strategy

This document describes the **exact** computational strategy implemented in the
Matlab code in this folder. It is meant to be the authoritative reference: any
re-implementation should match what is described here, file by file. Where the
code is silent, the doc says so explicitly.

The pipeline is a single chain:

```
usdata.m       → data.mat        (raw NIPA assembly, one-shot)
datamine.m     → worktemp.mat    (calibration + detrend + harness)
gmle.m         → worktemp.mle    (MLE of 30-θ state-space)
   runmleadj.m → uncmin × N      (10-restart minimization wrapper)
   mleqadj.m   → L(θ)            (likelihood objective)
gwedges2.m     → worktemp.w      (wedge time series + counterfactuals)
runall.m       → bca???.mat      (multi-country driver + f-stats collection)
```

---

## 0. Conventions and naming

**Quarterly frequency.** All time-series are quarterly; the calendar grid runs
on `t = 1980.25:0.25:2015` (i.e. 1980Q1, 1980Q2, ..., 2014Q4 — note `2015` here
denotes 2015Q0 = 2014Q4 in BCKM's `(Y + Q/4)` convention).

**Base date.** `bdate = 2008.25` (2008Q1) is hard-coded in `datamine.m`. All
"index = 100" series in `gwedges2.m` and `runall.m` use 2008Q1 as the base.
The Great Recession window is 2008Q1 → 2011Q4.

**Sample.** MLE sample is `1980.25:2015` (≈140 quarters). The base date is
inside the sample (`bind` between `iobs` and `eobs`).

**Variables (lower-case = data, no decoration = SS):**
- `ypc, xpc, hpc, gpc, cpc`: per-capita GDP, investment, hours, government,
  consumption
- `Y, X, H, G, C`: aggregate (un-divided) levels
- `y, x, h, l, g, c, k`: model SS values (also `ks, hs, ys` etc.)
- `lyt, lxt, lht, lgt, lct, lkt`: log-detrended time-series
- `lz, ltaul, ltaux, lg`: log of model-SS wedges (note `taul`, `taux` are
  *not* logged — they are levels in `[-1,1]`)
- `Sbar, P, Q`: stationary mean, transition matrix, lower-Cholesky shock
  covariance of the wedge VAR (4-dimensional)

**State-space dimensions inside the MLE:**
- State `X` is 6-dimensional augmented: `[log k, log z, taul, taux, log g, 1]`
- Observable `Y` is 4-dimensional: `[log y, log x, log h, log g]`
- Wedge VAR is 4-dimensional `S = [log z, taul, taux, log g]`

---

## 1. Raw data: `usdata.m`

**One-shot** script — produces `data.mat` from NIPA tables. Sample 1969:Q1
to 2014:Q2 (`T = 1:182`). Series:

```
rGDP, rPCE                                 from nipa116
rCD, rCND, rCS, rGPDI, rEX, rIM, rG        from nipa115/116
pCD                                        from nipa119
rGC, rGI                                   from nipa394/395
rSTX  = (nipa32(7,T) + sum(nipa33([9,11],T)))' / nipa119(4,T)' * 100
nKCD  = btab100d(T,9) / 1000               (consumer-durable stock, nominal)
nDCD  = atab10d(T,27) / 1000               (consumer-durable depreciation)
rKCD  = nKCD ./ pCD                        (real CD stock)
rDCD  = nDCD ./ pCD                        (real CD depreciation)
Pop   = 10^3·(civpop−civpop_inst) + 10^6·armed
H     = hours·10^9 / 4                     (annual→quarterly)
```

**Aggregates (BCKM convention):**

```
Y = rGDP - rSTX + 0.04·rKCD + rDCD                       sales-tax-adjusted GDP
                                                          + 4% CD service flow
                                                          + CD depreciation
C = rCND + rCS - (rCND+rCS)/(rCND+rCS+rCD)·rSTX
                + 0.04·rKCD + rDCD                        sales-tax-adj services
                                                          + CD service flow
X = rCD + rGPDI + rGI - rCD/(rCND+rCS+rCD)·rSTX           durables + GPDI + Gov-I
                                                          - CD share of sales tax
G = rGC + rEX - rIM                                       gov consumption + NX
                                                          (closed-economy CKM 2005)
```

`Y, X, G, C` are then converted to per-capita:
```
ypc  = Y./Pop·10^9, xpc = X./Pop·10^9, gpc = G./Pop·10^9, cpc0 = C./Pop·10^9
hpc  = H./Pop                                       (no further scaling)
```

The script truncates to start at `t == 1969.25` and HP-filters
`log[ypc, hpc, xpc, gpc]` with smoothing parameter 1600 to populate
`worktemp.lhpo` for the moment statistics in Table II.

**Output:** raw `data.mat` and a partially-built `worktemp.mat` (only the HP
moments — the rest is filled by `datamine.m`).

---

## 2. Calibration and detrending: `datamine.m` + `maketrend.m` + `calgz.m`

### 2.1 Calibration (`datamine.m`)

Hard-coded constants (independent of country/sample):
```
psi    = 2.5
sigma  = 1.000001                    (≈ log utility, with regularization)
theta  = 1/3                         (capital share α)
beta   = 0.975^(1/4)                 (annual β=0.975 → quarterly)
delta  = 1 - (1 - 0.05)^(1/4)        (annual δ=0.05 → quarterly)
```

Population growth from data:
```
gn = (iP_end / iP_1)^(1 / (T-1)) - 1
```
where `iP` is column 6 of `data.mat` (population index).

Tech growth `gz` is **solved numerically** (see §2.3).

MLE settings:
```
mlestart = 1980.25, mleend = 2015     (the MLE sample bracket)
nps      = 50                          (number of uncmin restarts in runmleadj)
pb       = 0.99                        (perturbation factor for restarts)
adjc     = 2                           (1=no, 2=BGG, 3=4×BGG capital adjustment cost)
```

`adja` (the actual adjustment-cost magnitude) is computed downstream in
`gwedges2.m` as
```
qadj = 0.25 / sum(params([gn, gz, delta]))
adja = adjcs(adjc),    adjcs = [0, qadj, 4·qadj]
```

### 2.2 Detrending (`maketrend.m`)

The detrending is **base-year normalization with a numerically calibrated
trend rate `gz`**. The function returns three things:

```
mled (T × 7):  [t, ypc·norm(by), xpc·norm(by), hpc, gpc·norm(by), cpc·norm(by), cpci·norm(by)]
                where norm(by) := (1/ypc(by)) · (1+gz)^by
                                                  ↑
                                base-year normalization scaled by trend at by

Y    (T × 6):  log(mled(:,2:7)) − [(1+gz).^t for cols 2,3,5,6,7;
                                   ones for col 4 (hours)]
gzt  (scalar): the calibrated quarterly tech growth rate
```

Conventions:
- All flow series (y, x, g, c, cpci) carry the trend `(1+gz)^t`.
- **Hours `hpc` is treated as stationary** — column 4 of `mled` is the raw
  per-capita hours series, and the corresponding column of `Y` subtracts
  `log(ones)` (i.e. nothing).
- The base-year scaling `(1+gz)^by · 1/ypc(by)` ensures that the level of
  each series at `t = bdate` is exactly the trend level
  `1·(1+gz)^by` for output, and `xpc(by)/ypc(by)·(1+gz)^by` for investment,
  etc.
- `cpci = ypc - xpc - gpc` is the *implied* consumption (resource constraint),
  separate from the data-measured `cpc`. Both are kept.

### 2.3 Trend rate calibration (`calgz.m`)

`gzt = fsolve(@calgz, 0)` finds the unique `gz` such that the
sample-period-mean of detrended log-output is zero:
```
mY = mean(  log(ypc·(1+gz)^(by-mles+1))
          - log(ypc(by-mles+1))
          - log((1+gz)^(0:T-1)') ),    over the MLE sample
```
This is BCKM's way of jointly choosing the trend rate AND the level so that
detrended log-`y` is centered at zero on the MLE window. As a consequence,
**`mean(Y(:,1)) = 0` exactly within the sample**.

For the US 1980Q1–2014Q4 sample, `gz ≈ 0.0047` per quarter (≈ 1.9%/year).

---

## 3. Steady state and harness: `bca_steady2.m` + `bca_params2.m`

`bca_steady2.m` reads `worktemp.mat`, computes the steady state from the
sample means of detrended observables, and **writes** `bca_params2.m` (a
data-file containing literal numeric assignments).

```
GDPs    = mean(exp(Y(:,1)));     Xs = mean(exp(Y(:,2)));     Hs = mean(exp(Y(:,3)));
gwedges = mean(exp(Y(:,4)));     Cs = mean(exp(Y(:,5)));

ks      = xs / (delta + grate - 1)             (capital from investment & δ, grate)
ewedges = ys / (ks^θ · hs^(1-θ))               (residual-defined efficiency)
adja    = 0.25 / (gn + gz + delta)             (capital adjustment cost coeff.)
```

The other constants (`beta, delta, eta = gn, gamma = gz, sigma, theta,
psi, lwedges = xwedges = 1`) are passed through unchanged from
`worktemp.params`. `bca_params2.m` is then `eval`-ed by `bca_simul2.m` and
`bca_wedges2.m` to seed the perfect-foresight simulator.

**Key fact:** `bca_steady2.m` is *for the perfect-foresight (PF) variant*
(`bca_simul2.m`), which is a sanity check — not the BCA core. The MLE itself
does NOT use `bca_steady2`/`bca_params2`. It computes its own steady state
inside `mleqadj.m` from `(I-P)\P0`.

---

## 4. MLE: `gmle.m` → `runmleadj.m` → `mleqadj.m`

### 4.1 Top-level harness (`gmle.m`)

`gmle.m` calls `runmleadj()` to get the 30-vector `Theta`, then re-shapes it
into named matrices and stores everything in `worktemp.mle`. **`Q` is stored
as a *symmetric* 4×4 matrix here** (Q(2,1) = Q(1,2), etc.) — but inside
`mleqadj.m`, only the *lower-triangular* 10 elements `Theta(21:30)` are
treated as free; the symmetric upper triangle is mirrored for the storage
representation in `gmle.m` only.

Wait — re-reading `gmle.m`: it actually does set `Q(1,2)=Theta(22)` etc., and
then mirrors `Q(2,1)=Q(1,2)`. **This is a storage convention bug
(BCKM-internal): `Theta(22)` is `Q(2,1)` inside `mleqadj.m` (lower-tri
Cholesky factor), but `gmle.m` stores it at `Q(1,2)` and mirrors. The
*objective* uses lower-tri, so the inconsistency is silent — `Theta(22)`
maps to `Q(2,1)` regardless.**

The mapped-but-mismatched stored `Q` (used as a Cholesky factor) is
`Q_lower`, with `Q·Q'` being the actual VAR shock covariance.

### 4.2 Restart loop (`runmleadj.m`)

10 starting points, BFGS at each. Initial Sbar is solved by
`fsolve(@initmle, [0; 0.05; 0; log(0.2)])` — this matches model SS levels to
sample means of the raw observables (see §4.4).

Three pre-coded P/Q starting matrices `x0a, x0b, x0c` differ in their Q
columns and correspond to `adjc = 1, 2, 3`. The selected `x0` uses the
adjc-indexed column.

```
x1 = uncmin(x0,        'mleqadj', adja)        % run 1
x1 = uncmin(x1,        'mleqadj', adja)        % run 2 (refine)
for i = 3:nps:                                  % runs 3..10
  x2 = uncmin(0.99·x2, 'mleqadj', adja)
end
[L, idx] = min(F);   x = X(:, idx)             % best of 10
```

The `0.99·x2` perturbation is intentional — it nudges the optimizer off the
local optimum to look for better basins.

The bootstrap-standard-error block (lines 161–225) is **commented out** in
production. Only the point estimate is computed.

### 4.3 Likelihood (`mleqadj.m`)

This is the workhorse. Reads `worktemp.mat` for data and parameters; receives
`Theta` (30×1) and `adja` (scalar) from the optimizer.

#### 4.3.1 Parameter unpacking

```
Sbar = Theta(1:4)
P    = reshape(Theta(5:20),  4, 4)   (column-major: Theta(5)=P(1,1), Theta(6)=P(2,1),...)
Q    = lower-tri 4×4 with Theta(21:30) filling [Q(1,1), Q(2,1), Q(3,1), Q(4,1),
                                                Q(2,2), Q(3,2), Q(4,2),
                                                Q(3,3), Q(4,3),
                                                Q(4,4)]
P0   = (I − P) · Sbar                 (so the implied unconditional mean is Sbar)
D, R = 0                              (no measurement-error AR or covariance)
```

Bounds: `P` entries clipped to `[-0.999, 0.999]`; Sbar to `[-1, -1, -1, -5]`
to `[1, 1, 1, 1]`. **Out-of-bounds returns L = 1e+20.**

Penalty: `500000 · max(|max(eig(P))| − 0.995, 0)^2`.
**No per-element bound on the eigenvalues — only the spectral radius.**

#### 4.3.2 Steady state (model)

```
[lz, taul, taux, lg] = (I − P) \ P0           (the implied unconditional mean of S)
zs = exp(lz);  gs = exp(lg)
beth = beta·(1+gz)^(-sigma)
kls  = ((1+taux)·(1 − beth·(1−delta)) / (beth·theta))^(1/(theta−1)) · zs
A    = (zs/kls)^(1−theta) − (1+gz)·(1+gn) + 1 − delta            (positive shift)
B    = (1−taul)·(1−theta)·kls^theta·zs^(1−theta) / psi
ks   = (B + gs) / (A + B/kls)
cs   = A·ks − gs
ls   = ks / kls
ys   = ks^theta · (zs·ls)^(1−theta)
xs   = ys − cs − gs
X0   = [log ks; log zs; taul; taux; log gs; 1]    (6-vector)
Y0   = [log ys; log xs; log ls; log gs]           (4-vector)
```

#### 4.3.3 Decision rule (linearization)

The Euler equation for capital (with adjustment costs `adja`) and labor
choice are linearized via numerical differentiation of `res_adjust(Z, param)`
at the SS values of `Z = [k_{t+2}, k_{t+1}, k_t, z_{t+1}, z_t, taul_{t+1},
taul_t, taux_{t+1}, taux_t, g_{t+1}, g_t]`:

```
del   = max(|Z|·1e-5, 1e-8)              (per-component step)
for i = 1:11:
  dR(i) = (res_adjust(Z+δ_i) − res_adjust(Z−δ_i)) / (2·del(i))
end
```

Then the saddle-path coefficients:
```
a0 = dR(1), a1 = dR(2), a2 = dR(3)                   (the k_{t+2}, k_{t+1}, k_t coefs)
b0 = dR(4:2:11)' = [dR(4), dR(6), dR(8), dR(10)]     (forward-state coefs)
b1 = dR(5:2:11)' = [dR(5), dR(7), dR(9), dR(11)]     (current-state coefs)

gammak = root of [a0, a1, a2] inside the unit circle
gamma  = -((a0·gammak + a1)·I + a0·P') \ (b0·P + b1)'      (4×1)
gamma0 = (1 − gammak)·log(ks) − gamma' · [log(zs); taul; taux; log(gs)]
```

So the linearized capital LoM is:
```
log k_{t+1} = gamma0 + gammak·log k_t + gamma' · [log z_t, taul_t, taux_t, log g_t]
```

**`res_adjust.m`** is the per-period Euler residual with quadratic capital
adjustment costs `phi(x/k) = (a/2)·(x/k − b)²`, where `b = (1+gz)·(1+gn) −
1+δ` is the SS investment-capital ratio. Hours `l` is solved
**internally** via 5 Newton steps on the labor FOC inside
`res_adjust.m` — there's no explicit closed-form. (Note: the residual code
hard-codes a starting `l = 0.99` for the iteration, which is a hold-over and
serves merely as initial condition for the Newton iteration.)

#### 4.3.4 Observation-equation coefficients (4 observables × 5 state cols)

Linearized reduced-form coefficients (no time-subscript = SS values) — for
output, investment, hours:

```
philh  = -(psi·ys·(1−θ) + (1−θ)·(1−taul)·ys·(1−ls)/ls·θ + (1−θ)·(1−taul)·ys)
philk  = ( psi·ys·θ + psi·(1−δ)·ks − (1−θ)·(1−taul)·ys·(1−ls)/ls·θ ) / philh
philz  = ( psi·ys·(1−θ) − (1−θ)²·(1−taul)·ys·(1−ls)/ls ) / philh
phill  = ( (1−θ)·(1−taul)·ys·(1−ls)/ls · 1/(1−taul) ) / philh
philg  = ( -psi·gs ) / philh
philkp = ( -psi·(1+gz)·(1+gn)·ks ) / philh

phiyk  = θ + (1−θ)·philk          (output ∂/∂ log k)
phiyz  = (1−θ)·(1+philz)          (output ∂/∂ log z)
phiyl  = (1−θ)·phill              (output ∂/∂ τ_l)
phiyg  = (1−θ)·philg              (output ∂/∂ log g)
phiykp = (1−θ)·philkp             (output ∂/∂ log k_{t+1})
phixk  = -ks/xs · (1-δ)           (investment ∂/∂ log k)
phixkp =  ks/xs · (1+gz)·(1+gn)   (investment ∂/∂ log k_{t+1})
```

The "kp" suffix denotes the coefficient on `log k_{t+1}`, which is
substituted out via `log k_{t+1} = Gamma · X` (so it shows up in the
final C as `+ phiykp · gamma(1:5)'`, etc.).

There is **no** `philx` or `phiyx` row in `mleqadj.m` — the production-side
investment wedge `taux` enters only through `philkp/phixkp` after the capital
LoM substitution. (Compare `mleqadj_check.m`, which is the *non-collapsed*
variant with explicit `philx, phiyx, ...` columns and is more readable.)

#### 4.3.5 State-space matrices (6D state, 4D obs)

```
A     = [ gammak,                            gamma',                          gamma0;
          [0; 0; 0; 0],                      P,                                  P0;
          [0,    0,      0,         0],      0,                                   1 ];      6×6

B     = [ 0, 0, 0, 0;
          Q;
          0, 0, 0, 0 ];                                                                      6×4

C     = [ [phiyk,  phiyz,  phiyl,    0,   phiyg]   + phiykp·Gamma(1:5)';                     output
          [phixk,    0,      0,      0,     0  ]   + phixkp·Gamma(1:5)';                     invest
          [philk,  philz,  phill,    0,   philg]   + philkp·Gamma(1:5)';                     hours
          [   0,      0,      0,     0,     1  ] ];                                          gov
                                                                                       (4×5)

phi0  = Y0 − C(:, 1:5) · X0(1:5)            (the "intercept" of the obs eq.)
C     = [C, phi0]                                                                            (4×6)
```

Two key observations:
1. The 4th observable (government) is identity in `log g` — it reads off the
   state directly.
2. The intercept `phi0` is folded into the 6th column of `C`, paired with
   the 6th state component which is the constant `1`. So the augmented state
   is **driverless** in its 6th coordinate — `[1] → A·X → ...; 1` always.

#### 4.3.6 Likelihood evaluation

```
T          = length(ZVAR)                                  (= eobs − iobs + 1)
Y          = log(ZVAR) − log[(1+gz)^t for cols 1,2,4; ones for col 3]
                                                          (re-applies the trend
                                                          subtraction inside MLE
                                                          since ZVAR is base-year
                                                          normalized but NOT
                                                          time-detrended)
Ybar       = Y(2:T, :) − Y(1:T-1, :) · D'      = Y(2:T, :)   (since D=0)
T          = T − 1                                          (drop first obs)
Cbar       = C·A − D·C                          = C·A         (since D=0)
Rbar       = R + C·B·B'·C'                      = C·B·B'·C'   (since R=0)

[K, Sigma] = kfilter(A, Cbar, B·B', Rbar, B·B'·C')           steady-state Kalman
                                                              with cross-cov

Omega      = Rbar + Cbar · Sigma · Cbar'                     innovation cov
Omegai     = inv(Omega)

Xt(1, :)   = X0'                                              initialization
innov(1,:) = Ybar(1,:) − X0'·Cbar'

for i = 2:T:
  Xt(i, :)    = Xt(i-1, :)·A' + innov(i-1, :)·K'              Kalman update
  innov(i, :) = Ybar(i, :)    − Xt(i, :)·Cbar'
end

sum1 = innov(1:T, :)' · innov(1:T, :) / T                     sample MSE

L    = 0.5·( T·(log|Omega| + tr(Omega⁻¹·sum1)) + penalty )    minimized!
```

**Key facts about this LL:**
- It is a **negative half log-likelihood without the `(2π)` constant** — i.e.,
  it is a (positive) cost function that BFGS minimizes.
- After dropping the first observation: `T_effective = T_data − 1` (e.g., 139
  for a 140-quarter sample).
- `Ω` is **time-invariant**: BCKM uses the steady-state Kalman filter with
  a constant gain `K`, never the transient time-varying recursion.
- The cross-covariance `V12 = B·B'·C'` is non-zero and is correctly accounted
  for in `kfilter.m` (it inverts `V2` into `A` and `V1` if `rank(V2) ≥
  size(C, 1)`, and otherwise iterates).

#### 4.3.7 Steady-state Kalman filter (`kfilter.m` + `doubalg.m`)

`kfilter(A, C, V1, V2, V12)`:
- **If `rank(V2) ≥ rows(C)`** (the typical case here, since `V2 = Rbar`
  inherits rank from `C·B·B'·C'`):
  ```
  A_new  = A − (V12/V2)·C
  V1_new = V1 − V12·(V2 \ V12')
  [k, s] = doubalg(A_new, C, V1_new, V2)            (DARE via doubling)
  k      = k + V12 / V2
  ```
- **Else**: iterative Riccati update with tolerance `1e-16`.

The doubling algorithm `doubalg` solves the discrete-time algebraic
Riccati equation `S = A·S·A' + V1 − A·S·C'·(V2 + C·S·C')⁻¹·C·S·A'` with
quadratic convergence — typically converges in 30–50 iterations.

### 4.4 Initial Sbar (`initmle.m`)

`initmle(Sbar)` is called via `fsolve` and returns the residual:
```
x0 = [ ys − Ym(1);
       xs/ys − Ym(2);
       ls − Ym(3);
       gs/ys − Ym(4) ]
```
where `Ym = mean(exp(ZVAR))` is the sample mean of the *level* observables
(so `Ym(1)` is mean output, `Ym(2)` is mean investment, etc.) and `ys, xs,
ls, gs` are the model SS values *given* the candidate `Sbar`. The fsolve
finds `Sbar` such that:
1. Model SS output equals sample mean output
2. Model SS investment-output ratio equals sample mean of `exp(ZVAR(:,2))`
3. Model SS hours equals sample mean hours
4. Model SS gov-output ratio equals sample mean of `exp(ZVAR(:,4))`

**The mismatched scaling between channels (level vs. ratio) is BCKM-original
— do not "fix" it.** It produces a reasonable warm-start that BFGS then
refines.

### 4.5 The `mleqadj_check.m` variant

This is a debug/verification version of `mleqadj.m` with:
- Hard-coded calibration constants for "the original CKM 2007 paper":
  `gn = 1.015^(1/4)−1`, `gz = 1.016^(1/4)−1`, `beta = 0.9722^(1/4)`,
  `delta = 1−(1−0.0464)^(1/4)`, `psi = 2.24`, `theta = 0.35` —
  **NOT the BCKM 2016 calibration used in `datamine.m`**. Useful for
  sanity-checking against CKM 2007 numbers but NOT to be confused with the
  production calibration.
- Explicit `philz`, `philx`, `philg` columns (not collapsed via `gamma·Gamma`)
  — a pedagogical version.
- Different decision-rule formulation: solves a quadratic in `gammak` whose
  coefficients are explicit functions of `coef1, coef2, coef3, coef4` (the
  Euler-equation linearization with adjustment costs).

The two should agree numerically when both use the same calibration; the
production code is the *less* readable but *more compact* `mleqadj.m`.

---

## 5. Optimizer (`uncmin.m` + `umlnmin.m` + `umstop.m` + `enorm.m`)

Custom BFGS implementation, no MATLAB optimization toolbox dependency.

**Tolerances (all hard-coded):**
```
epsa    = 1e-16     (absolute gradient-norm tolerance)
xtol    = sqrt(1e-10) ≈ 1e-5    (change in x)
ftol    = 1e-10     (change in f)
gtol    = 1e-7      (relative gradient norm: ||g|| / (1 + |f|))
steptl  = 1e-16     (step length minimum — was 1e-3·sqrt(epsm), changed)
maxits  = 10000
```

**Gradient:** finite-difference, forward by default. If a line-search fails
(`iretcd == 1`), fall back to central differences (`iagflg = -1`).
- Forward step: `fddev = sqrt(1e-16) · diag(max(|x|, xsize))` ≈ `1e-8·|x|`.
- Central step: `cddev = (1e-16)^(1/3) · diag(...) ≈ 5e-6·|x|`.

**Hessian:** maintained via Cholesky factor `L`, updated by the standard
BFGS formula:
```
L · S · L' (after update) where S = the new Hessian
```
The update is the rank-2 BFGS in factor form (Gill-Murray-Wright). The
positive-definiteness check is explicit: if `L · S · L'` has negative or
imaginary eigenvalues, **the update is skipped** (the Hessian stays
unchanged for that iteration).

**Line search (`umlnmin.m`):** cubic interpolation with safeguarded
back-tracking. Accepts when `f(x + λ·p) < f(x) + 1e-12·g'·p·λ` (very loose
sufficient-decrease).

**Stopping (`umstop.m`):** terminates if 3 of {gradient, x-change, f-change}
all hit tolerance simultaneously, OR if gradient drops below `epsa`, OR if
the line-search failed twice (`iretcd == 1`), OR if `maxits` reached, OR
if 5 consecutive `mxtake` (max-step-length) events occurred.

**`enorm.m`** is a vectorized "norm-each-element" helper — used in some
auxiliary scripts but not in the main pipeline.

---

## 6. Wedge extraction: `gwedges2.m`

After MLE finishes, `gwedges2.m` re-solves the state-space at the converged
`Theta` and back-computes the four wedge time series.

### 6.1 Re-solve at converged θ

```
[L, Sbar, P0, P, Q, A, B, C, param] = mleqadj(x0, adja)        % gives A, C, etc.
T   = length(ZVAR)
Y   = log(ZVAR) − [(1+gz)^t for cols 1,2,4,5; ones for col 3]    % time-detrend
                                                                  % (Y now has 5 columns
                                                                  % since ZVAR has 5 cols
                                                                  % including consumption)
lyt, lxt, llt, lgt, lgct = Y(:, 1:5)
```

### 6.2 Capital path (deterministic, linear)

```
lkt(1) = lk = log(ks)            % SS log-capital
Kt(1)  = exp(lk) = ks            % SS level capital
for i = 1:T:
  lktp(i) = lk + ((1-δ)·(lkt(i)-lk) + (xs/ks)·(lxt(i)-lx)) / ((1+gz)·(1+gn))
  lkt(i+1) = lktp(i)
  Ktp(i)   = ((1-δ)·Kt(i) + exp(lxt(i))) / ((1+gz)·(1+gn))
  Kt(i+1)  = Ktp(i)
end
```

The log-version is **linearized around the SS** (so it's NOT the same
trajectory as the level version `Kt` for large deviations); both are kept,
but the Solow-residual computation later uses the **level `Kt`**.

**Initial condition:** capital starts at the model SS `ks`. This means there
is a transient adjustment in the first ~20–30 quarters. `gwedges2.m` treats
this as fine because the wedge analysis is centered at `bind = 2008Q1`
(2008.25), well past the transient.

### 6.3 Wedge formulas

**Efficiency wedge (Solow residual, levels-form):**
```
Zt    = (Y_t / (Kt^θ · H_t^(1-θ)))^(1/(1-θ))       % level
w.zt  = (Zt / Zt(Y0))^(1-θ)                        % base 1 at Y0=2008.25
```

The `^(1-θ)` exponent at the end is crucial: BCKM defines the *wedge* as
`A_t^(1-θ)` rather than `A_t` itself, so it appears multiplicatively in the
production function as `Y = (A·H)^(1-θ)·K^θ` (Hicks-augmenting).

A *log-linear* version is also computed for the state-space projection:
```
lzt = lz + (lyt − ly − θ·(lkt − lk))/(1−θ) − llt + ll
```

These two should agree to first order at the SS but differ by curvature for
larger deviations.

**Labor wedge (levels-form, intratemporal FOC):**
```
Tault   = 1 - psi/(1-θ)·(Ct/Y_t)·(H_t/(1-H_t))      % τ_l in level [0, 1)
w.tault = (1 - Tault) / (1 - Tault(Y0))             % base 1 at Y0
```

A log-linear version for state-space:
```
tault = taul + (1 − taul)·(lyt − ly − lct + lc − 1/(1-l)·(llt − ll))
```

**Investment wedge (linear, solved from C-matrix row 2):**
There is **no closed-form** for `taux` — the investment Euler is forward-
looking and contains capital adjustment costs. BCKM uses the **linearized
investment equation** (row 2 of the obs-equation matrix `C`) and inverts it:

```
tauxt = (lxt − C(2,1)·lkt − C(2,2)·lzt − C(2,3)·tault − C(2,5)·lgt − C(2,6))
        / C(2,4)
```

`C(2, j)` are the SS-evaluated coefficients on (log k, log z, taul, taux,
log g, 1) in the linearized investment equation. The taux column `C(2,4)`
contains `phixkp · gamma(4)` (the capital LoM substitution), which is what
identifies it. **A consistency check `tauxchk` uses row 1 (output) instead
of row 2 — the two should agree at convergence.**

**Government wedge:** `w.gt = exp(lgt − lgt(Y0))`. Trivial — it is
directly observable (the 4th observable).

### 6.4 Counterfactual decomposition: `fixexpadj.m`

For each subset of "active" wedges `As ∈ {0, 1}^4`, solve a **separate**
state-space and produce projected observables.

The construction inside `fixexpadj.m`:
- Same SS computation as `mleqadj.m`.
- Same `gammak, gamma, gamma0` computation, using `res_adjust2.m` instead of
  `res_adjust.m`. **`res_adjust2.m` differs from `res_adjust.m` by adding
  the wedge-mask multipliers `Az, Al, Ax, Ag` and `s0` (initial state) to the
  argument list** — the residual is computed at a state where wedges
  `S_j = As(j)·s_j + (1-As(j))·Sbar(j)`. This is BCKM's "fixed expectations"
  device: the *decision rule* (`gammak, gamma, gamma0`) is computed at the
  full S-process P, but the *observable equation* C zeros out the columns
  corresponding to inactive wedges.
- The C matrix is constructed with the wedge multipliers `As` applied to the
  output and hours rows:
```
C = [ [phiyk, phiyz·As(1), phiyl·As(2), 0, phiyg·As(4)] + phiykp·Gamma(1:5)';   % output
      [phixk,           0,           0, 0,           0] + phixkp·Gamma(1:5)';   % invest
      [philk, philz·As(1), phill·As(2), 0, philg·As(4)] + philkp·Gamma(1:5)';   % hours
      [    0,           0,           0, 0,       As(4)] ];                       % gov
```

So when `As = [1, 0, 0, 0]` (only A active), `C(1, 3) = phiyz · 1 = phiyz`
but `C(1, 4) = phiyl · 0 = 0` etc. — the labor and gov wedges contribute
nothing to output or hours. The investment row is unchanged across As (no
direct wedge entry).

### 6.5 Wedge contributions to observables

Given the per-wedge state-space matrices `C_k` for `k ∈ {none, A, L, X, G,
all, ...}`, the contribution of each wedge bundle to observables is computed
via:
```
YM0   = Xt0 · C'                                            % full (all wedges)
YMz   = (Xt0 - Xt0(Y0,:))·(C1 - C0)' + YM0(Y0,:)            % only A active
YMl   = (Xt0 - Xt0(Y0,:))·(C2 - C0)' + YM0(Y0,:)            % only τ_l
YMx   = (Xt0 - Xt0(Y0,:))·(C3 - C0)' + YM0(Y0,:)            % only τ_x
YMg   = (Xt0 - Xt0(Y0,:))·(C4 - C0)' + YM0(Y0,:)            % only g
YMnox = (Xt0 - Xt0(Y0,:))·(C1+C2+C4 - 2·C0)' + YM0(Y0,:)    % all but τ_x
... etc.
```

The "subtract 2·C0" term in 3-wedge counterfactuals is because `C0` is the
"all-wedges-off" decision rule — adding three (Cj - C0)'s gives `C1 + C2 +
C4 - 3·C0`, and we add `+C0` once to recover the additive baseline; net `-
2·C0`.

This is the **linear superposition** of wedge contributions in BCKM. It
relies on the linearity of the obs-equation in the active-mask `As` — which
is exact for the linearized model.

`Xt0 = [lkt, lzt, tault, tauxt, lgt, ones]` is the 6-D state from §6.2-6.3.
The observable contributions are then exponentiated and base-normalized:
```
mz_y = exp(YMz(:,1) - YMz(Y0, 1)) · 100      % output, base-100 at Y0
mz_h = exp(YMz(:,3) - YMz(Y0, 3)) · 100      % hours
mz_x = exp(YMz(:,2) - YMz(Y0, 2)) · 100      % investment
mz_c = (Ymz_C / Ymz_C(Y0)) · 100             % consumption (from resource constraint)
```

### 6.6 Moment statistics (Tables II/III)

`gwedges2.m` lines 224+ HP-filter (smoothing `400·freq = 1600` for
quarterly) the wedges, observables, and per-wedge model components, then
compute:
- Table IIA: relative std of each wedge / std of output
- Table IIA1o: relative std of each observable / std of output
- Table IIB: cross-correlations among wedges
- Table IIIA: relative std of model-implied series / std of output (lines
  282–331)
- Table IIIB: cross-correlations among model components

These populate `worktemp.tableIIA1, IIA2, IIB, IIIA1, ...IIIB` etc.

---

## 7. F-statistics (the BCKM Table 11 measure)

The f-statistic is computed in **two places** with two different
implementations:

### 7.1 In-pipeline (`gwedges2.m` lines 365+ and `runall.m` lines 67+)

Identical code. For each variable v ∈ {Y, H, X, C} and each wedge w ∈
{z, l, x, g}, build the wedge-removed error series:
```
m{w}{v}e[t] = Y_norm[t] - mw_v_norm[t]
            = w.Y[t, v_col] / w.Y[ind, v_col]
              - w.{m{w}{v}}[t] / w.{m{w}{v}}[ind]
```
where `ind = bind` (start of GR window), `wsize = ind+wsize-1` is the GR
window end (= `wend`).

Then the f-statistics are the **inverse-MSE weights**:
```
temp = mean(m_we[:, 1:4].^2, 1)             % size-4 row of MSEs by wedge
fz   = (1/temp(1)) / sum(temp.^-1)
fl   = (1/temp(2)) / sum(temp.^-1)
fx   = (1/temp(3)) / sum(temp.^-1)
fg   = (1/temp(4)) / sum(temp.^-1)
```

This is the **level-ratio** formulation — series are indexed at base 1 at
the start of the window, errors are ratio-deviations, MSE is over the GR
window. **The 4 wedge contributions sum to 1.**

### 7.2 Standalone (`fstats3.m`)

Same idea, different implementation. Uses `worktemp.w.mzy`, `mly`, `mxy`,
`mgy` etc. (already 100-normalized at Y0), but with options:
- `ilog` ∈ {0, 1}: use levels (default) or logs.
- `ifilt` ∈ {0, 1}: HP-filter the series first (default no).

For levels, no HP filter:
```
f(1, i) = 1 / sum((dly - dlyc(:, i)).^2)         % output, wedge i
f(2, i) = 1 / sum((dll - dllc(:, i)).^2)         % hours
f(3, i) = 1 / sum((dlx - dlxc(:, i)).^2)         % investment
fstat   = f ./ (sum(f')' · ones(1, 4))           % normalize each row to 1
```

The window is `i1 = 2008.25 → i2 = 2011.75` (Great Recession). This
matches BCKM Table 11 exactly.

**Why two versions?** `fstats3.m` is the published Table 11 implementation;
`gwedges2.m` / `runall.m` is the same idea but uses the in-memory
`worktemp.w.*` variables for cross-country batch processing. They produce
the same numbers up to floating-point.

---

## 8. Perfect-foresight cross-check: `bca_simul2.m`

Independent, **non-linearized** simulation of the model under perfect
foresight given the wedges from `bca_wedges2.m`. This is the "PF" track —
not the BCA-MLE track. Used in `runnit.m` to cross-validate that the
linear MLE wedges, when fed back into a non-linear PF simulator, reproduce
the observed time series.

The PF simulator solves the full Euler/labor system as a *forward* boundary-
value problem via the secant method (`secant.m`), starting from the
`vinitial.dat` initial guess. Each one-wedge counterfactual is a separate
PF solve.

The output is `pftemp.mat` with simulated `yt, ht, xt` for each
counterfactual.

`bca_resid2.m` is the PF residual that `secant.m` zero-finds:
```
foc = [ Euler eqn for each t in 1..T-1 (capital);    %  T-1 equations
        terminal cond: kn(T) = deriv·kn(T-1);        %  1 equation
        labor FOC for each t (psi·c·h/y - (1-θ)·lwedge·(1-h) = 0) ]
                                                      %  T equations
                                                      %  → 2T equations, 2T unknowns
                                                      %  unknowns are [kn(1..T), h(1..T)]
```

This is a 2T-dimensional nonlinear system solved by Newton's method (with
finite-difference Jacobian) inside `secant.m`. **It uses no Kalman / no
likelihood — it is just deterministic forward simulation given the
wedges.**

---

## 9. Counter-check: `gacc.m` / `gacc_a.m`

Quick growth-accounting plot — applies the Solow-residual formula
`dA = dy − α·dk − (1-α)·dh` directly to detrended observables (with `K`
built by perpetual inventory from investment). Used to sanity-check the
efficiency wedge against a naive direct calculation. Not part of the
production pipeline.

---

## 10. What the production result looks like (USA, 1980-2014)

After the pipeline runs:

**Calibration (BCKM Table 1):**
- α = 1/3, β = 0.975 (annual), δ = 0.05 (annual), ψ = 2.5, σ ≈ 1
- γ ≈ 1.9% / year, n ≈ 1.0% / year (data-derived)

**MLE estimates (BCKM Table 8/9/10):**
- Sbar = (logz, taul, taux, logg) ≈ (0.13, 0.37, -0.05, -1.94)
- P diagonal ≈ (0.989, 1.001, 0.968, 0.995), max|eig(P)| ≈ 0.996
- Q lower-tri 10 elements (Cholesky of shock covariance) — numerical

**F-stats on Great-Recession window 2008Q1-2011Q4 (BCKM Table 11):**
- fY[A]   ≈ 0.16   (efficiency wedge)
- fY[τ_l] ≈ 0.46   (labor wedge — dominant)
- fY[τ_x] ≈ 0.32   (investment wedge)
- fY[g]   ≈ 0.06   (government wedge — minor)

These four sum to 1 by construction.

---

## 11. Inventory of every .m file in this folder

| File                  | Role                                         | Production path?          |
|-----------------------|----------------------------------------------|---------------------------|
| `usdata.m`            | Build `data.mat` from raw NIPA               | yes (one-shot)            |
| `datamine.m`          | Calibrate + detrend + harness `worktemp.mat` | yes (entry point)         |
| `datamine_l.m`        | Variant for `data_l.mat`                     | alt data path             |
| `maketrend.m`         | Detrend by base-year norm + `gz`             | yes                       |
| `calgz.m`             | Trend-rate calibration residual              | yes (subroutine of above) |
| `bca_steady2.m`       | SS + write `bca_params2.m`                   | only PF check             |
| `bca_params2.m`       | Numeric SS values (auto-generated)           | only PF check             |
| `bca_wedges2.m`       | Direct wedge calc for PF check               | only PF check             |
| `bca_simul2.m`        | PF simulation (per-wedge)                    | only PF check             |
| `bca_resid2.m`        | PF Euler residual                            | only PF check             |
| `gmle.m`              | Top-level MLE harness                        | yes                       |
| `runmleadj.m`         | 10-restart MLE driver                        | yes                       |
| `mleqadj.m`           | MLE objective                                | yes                       |
| `mleqadj_check.m`     | Pedagogical / debug variant                  | reference                 |
| `initmle.m`           | Sbar warm-start (fsolve target)              | yes                       |
| `kfilter.m`           | Steady-state Kalman gain                     | yes                       |
| `doubalg.m`           | DARE doubling algorithm                      | yes                       |
| `res_adjust.m`        | Per-period Euler residual                    | yes (mleqadj subroutine)  |
| `res_adjust2.m`       | Same w/ wedge masks                          | yes (fixexpadj subroutine)|
| `fixexpadj.m`         | Per-mask state-space construction            | yes (gwedges2 subroutine) |
| `gwedges2.m`          | Wedge extraction + counterfactuals + tables  | yes                       |
| `runall.m`            | Multi-country driver + f-stats collection    | yes                       |
| `fstats3.m`           | Standalone f-stat (Table 11)                 | yes (post-hoc)            |
| `runnit.m`            | PF cross-check driver (multi-country)        | only PF check             |
| `uncmin.m`            | Custom BFGS                                  | yes (optimizer)           |
| `umlnmin.m`           | Cubic line search                            | yes                       |
| `umstop.m`            | Convergence check                            | yes                       |
| `secant.m`            | Newton/secant fixed-point                    | only PF check             |
| `enorm.m`             | Element-wise norm helper                     | aux                       |
| `printit.m`           | Plot script (graphics)                       | post-hoc                  |
| `ptime.m`             | Time-axis builder                            | aux                       |
| `gacc.m`, `gacc_a.m`  | Direct growth-accounting check               | post-hoc                  |
| `Untitled.m`,         | Throwaway plot scripts                       | manual                    |
| `Untitled2.m`         |                                              |                           |
| `test.m`              | Manual debugging                             | manual                    |
| `dump_worktemp_for_diff.m` | (User-added) export worktemp for diff   | added by user             |
| `octave_fresh_run.m`,      | (User-added) Octave runners                | added by user             |
| `octave_multistart.m`      |                                            |                           |

---

## 12. Critical implementation details that are *easy to get wrong*

These are the points where re-implementations have historically diverged
from BCKM. Each one is a hard requirement.

1. **Hours is uncentered.** `mled(:, 4) = hpc` raw — no division by sample
   mean, no division by `ls`, no log-shift. The corresponding column of
   `Y` subtracts `log(ones)` (i.e., nothing). The sample mean of `log(hpc)`
   is *not* zero — it is whatever the data say (≈ log(0.243) for US).
   Consequently `Sbar(2)` (taul) absorbs the SS-vs-data hours level gap at
   MLE time.

2. **First observation is dropped.** `T_effective = T_data − 1`. The
   `Ybar` matrix has `T-1` rows. `sum1 = innov'·innov / T_effective`. The
   `T` factor in front of `log|Ω|` is the same `T_effective`.

3. **Likelihood is L = -loglike + const, minimized.** The (2π) constant is
   omitted. Going from BCKM `L` to a standard Gaussian log-likelihood `LL`:
   `LL = -L - 0.5·T_eff·k·log(2π)` where `k = 4` (number of observables) and
   `T_eff` is the effective sample size after the first-obs drop.

4. **Steady-state Kalman, per call.** The DARE is solved at every call to
   `mleqadj.m` (`kfilter` is invoked from inside the objective at the
   *current* P, Q values). The covariance `Σ` is NOT frozen at any reference
   parameter set.

5. **No per-element bound on `eig(P)`.** Only the spectral radius is bounded.
   Individual diagonals can exceed 0.995 (Table 8 has `taul` diagonal
   ≈ 1.001 — by design).

6. **Government is gov-cons + net exports.** `G = rGC + rEX - rIM` (CKM
   2005 closed-economy convention). NOT just government consumption.

7. **Investment includes durables.** `X = rCD + rGPDI + rGI - rCD/(...)·rSTX`.
   Excluding durables makes `xpc` ~3 percentage points smaller and shifts
   `Sbar(3)` (taux).

8. **The Cholesky `Q` is lower-triangular with 10 free parameters.**
   `Theta(21:30)` fills `Q(1,1), Q(2,1), Q(3,1), Q(4,1), Q(2,2), Q(3,2),
   Q(4,2), Q(3,3), Q(4,3), Q(4,4)`. The shock covariance is `Q·Q'`. Symmetric
   storage in `gmle.m` (mirroring upper-tri = lower-tri) is for display only.

9. **Decision rule uses *full* P, even in counterfactuals.** `fixexpadj`
   solves for `gammak, gamma, gamma0` at the *full* 4-D VAR P, not at a
   reduced one-wedge VAR. The counterfactual switch only zeros out the
   *observation-equation* columns (`As(j)·phi*` in C), never the state
   transition.

10. **Capital path in wedge extraction starts at SS.** `Kt(1) = ks` and
    `lkt(1) = lk`. The first ~20-30 quarters are a transient that the
    user is expected to discard before doing GR-window analysis.

11. **The base-year normalization happens inside `maketrend.m`, but the
    *time*-detrending happens twice — once in `mled(:,j+1) =
    series·(1+gz)^by/ypc(by)` and again in `mleqadj.m` line 237 as
    `Y = log(ZVAR) − log((1+gz).^[0:T-1])`. The two together give a
    series that is base-year-1 at `bdate` AND time-detrended, which is
    what the state-space wants.** This is not a bug — BCKM is using the
    base-year scaling to set the *level* of the trend and the in-objective
    subtraction to remove the time-component.

---

## 13. Closing note

This document is meant to be exhaustive. If your re-implementation produces
different `f`-statistics than the BCKM Table 11 numbers (0.16/0.46/0.32/0.06
for output), the divergence is **almost certainly** in:

(a) Data construction: the NIPA aggregates `Y, X, G, C` and the per-capita
    deflators are heavily convention-dependent; a mismatch of a few percent
    here propagates linearly through `Sbar` and into the wedge time series.
(b) The `Y_l = log(hpc)` uncentered convention (point 1 above).
(c) The first-obs drop and the (2π) constant (points 2–3).
(d) The full-P decision rule in counterfactuals (point 9).

Anything else (Kalman initialization, optimizer choice, calibration of
δ/β/ψ/α) is second-order at the f-stat level.
