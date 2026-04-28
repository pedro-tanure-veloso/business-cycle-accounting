# BCKM Replication — Assessment & Path to a Verdict

*Last updated: 2026-04-28. Supersedes the earlier tiered divergence analysis,
which was written before the Step 9 calibration audit and contained an
incorrect "BCKM defaults" calibration table (α=0.35, ψ=2.24, δ=0.0464,
β=0.9722). Those values are **not** in the paper or in `datamine.m`; they
appear to have come from a prior session's misreading. The doc has been
rewritten end-to-end to reflect the current state of the project.*

---

## 1. Goal

Replicate BCKM (2016) US Business Cycle Accounting MLE estimates against
`BCA/BCA_info.md` Section 7 Tables 8–11 — in particular, the f-statistics:

| target            | fY[A] | fY[τ_l] | fY[τ_x] | fY[g] |
|-------------------|-------|---------|---------|-------|
| BCKM Table 11     | 0.16  | 0.46    | 0.32    | ~0    |

The **purpose of the replication is methodology validation**, not a
deliverable in itself. The user's actual research question is post-2014
US business-cycle decomposition. Replication is the gate we have to clear
before trusting forward analysis.

---

## 2. What is currently working

After Step 9's audit + calibration fix, the pipeline is BCKM-faithful in
every architectural respect we have checked:

| component                                | status                                         |
|------------------------------------------|------------------------------------------------|
| Calibration constants (α, ψ, δ, β)       | ✅ matches BCKM Table 1 + `datamine.m` exactly |
| Theta layout (Sbar(4) + P(16) + Q_tri(10))| ✅ 30 params, matches `mleqadj.m`              |
| Sbar parametrization, P₀ = (I-P)·Sbar    | ✅ matches `mleqadj.m:29` and `runmleadj.m`    |
| Per-iteration SS re-solve at Sbar (Path A)| ✅ matches `initmle.m` style                  |
| Steady-state Kalman, DARE per call       | ✅ matches `kfilter.m` / `mleqadj.m`           |
| Spectral-radius-only penalty (5e5, 0.995)| ✅ matches `mleqadj.m:134`                     |
| Analytical wedge extraction (`gwedges2`) | ✅ smoother and analytical agree to 4 decimals |
| BCKM `fstats3.m`-style f-stats           | ✅ Y₀-rebased levels, inverse-MSE             |
| Test suite                               | ✅ 65/65 pass                                 |

We have a complete BCA decomposition pipeline that matches the paper's
methodology. We can compute wedges, counterfactuals, and f-stats for any
window we feed it. The infrastructure is sound.

---

## 3. What we cannot do — the gap

We cannot reach BCKM's Table 8/9/10/11 numbers from any warm-start we have
tried. The two best basins we find are:

| run                          | LL   | fY[A] | fY[τ_l] | fY[τ_x] | fY[g] | story            |
|------------------------------|------|-------|---------|---------|-------|------------------|
| BCKM Table 11 target         | —    | 0.16  | 0.46    | 0.32    | ~0    | balanced         |
| pre-calibration-fix basin    | 1830 | 0.023 | 0.619   | 0.345   | 0.013 | labor-dominated  |
| post-calibration-fix basin   | 1819 | 0.030 | 0.035   | 0.917   | 0.018 | investment-dominated |

The two basins differ by **only 11 nats** in log-likelihood but tell
**economically opposite stories** about the Great Recession. That is the
core problem: the likelihood landscape is approximately flat between
plausible economic narratives.

Evaluating BCKM's *published* MLE parameters (from `worktemp.mat`) on our
data gives LL = −8110 — a 9930-nat gap from our optimum. So:

- We cannot reach BCKM's basin from our starting points.
- BCKM's basin is not even visible from where we land.
- Our basin is the global optimum *of our objective on our data*.

This rules out a pure optimizer-scope explanation. The remaining
candidates are data construction, sample window, or genuine
under-identification of the BCA likelihood.

---

## 4. What we ruled out in Step 9

A multi-phase numerical comparison against `worktemp.mat` (Phases A–F):

- **Phase B** (observables): RMSE 1–3% per series. Hours fixed by
  switching to PAYEMS × AWHNONAG. No remaining smoking gun in the
  observables themselves.
- **Phase C** (BCKM params → our smoother): smoothed wedges match BCKM's
  `worktemp.w` to within 1–16% RMSE. Smoother is **not structurally
  broken**.
- **Phase D** (wedge cross-correlations): xcorr(τ_l, τ_x) = 0.88 in our
  basin vs 0.45 in BCKM's. Indicates the basins differ in how they
  partition variance between labor and investment, not in the smoother
  mechanics.
- **Phase F** (Sbar coordinate translation): the LL surface is
  non-degenerate around BCKM's Sbar; the gap is not a parametrization
  artifact.

So what we have is two well-defined likelihood basins that the optimizer
chooses between based on warm-start, with a flat ridge in between.

---

## 5. The unresolved question

**Does running BCKM's MATLAB on their `data.mat` reproduce Tables 8–11
exactly, or does it land in one of several local optima?**

Three possibilities, each implying a different course of action:

- **(A) Their code → Tables 8–11 deterministically.** The published
  numbers are the unique global optimum on their dataset. Our gap is
  then purely data construction or optimizer scope.
- **(B) Their code → close-but-not-exact** (e.g., reproduces f-stats
  within ±0.05 but Table 8 P matrix differs slightly). The published
  numbers are *a* MLE optimum, the replication standard is "within ε".
- **(C) Their code → multiple basins like ours.** The published Tables
  are one draw of a methodologically under-identified problem. Wedge
  accounting is a flexible decomposition with weak likelihood
  identification (Cole–Ohanian critique).

We have not run BCKM's actual MATLAB code on their actual data. Until we
do, we cannot distinguish these three worlds.

---

## 6. Three-step plan to a verdict

| step | description                                                       | effort  |
|------|-------------------------------------------------------------------|---------|
| 1    | Run BCKM's MATLAB end-to-end in Octave on their `data.mat`.       | ~1 day  |
| 2    | Bootstrap our estimator (~100 random seeds, study f-stat distribution). | ~2 hr   |
| 3    | Sensitivity to calibration (vary α, ψ, δ, β by ±5%).              | ~1 hr   |

Step 1 resolves (A)/(B)/(C) directly. Step 2 tells us whether *our*
methodology is unimodal or bimodal independent of BCKM. Step 3 quantifies
how robust the conclusions are to normal calibration uncertainty — useful
even if Step 1 resolves (A).

---

## 7. Decision matrix for forward US analysis

| Step 1 outcome | Step 2 outcome | Verdict for forward analysis                        |
|----------------|----------------|-----------------------------------------------------|
| Reproduces (A) | unimodal       | Trust the framework. Chase residual data gap.       |
| Reproduces (A) | bimodal        | Optimizer issue specific to ours. Widen search.     |
| Approximate (B)| unimodal       | Trust with point-estimate caveat.                   |
| Approximate (B)| bimodal        | Trust with distributional caveat.                   |
| Multi-basin (C)| any            | Use Bayesian / report distribution, not point est.  |

The user's research goal — analyzing post-2014 US business cycles —
remains achievable in every cell of this matrix. What changes is what
we *say*: a single decomposition vs. a distribution over basins, and
how confidently we attribute slowdowns to specific wedges.

If outcome is (C), that itself is a publishable methodological finding:
that wedge accounting suffers from the same weak-identification problems
Cole and Ohanian (2007) flagged for the Great Depression analysis.

---

## 8. What changed in this iteration

What we modified in this round:

- **`bca_core/params.py`** — restored BCKM Table 1 calibration:
  α=1/3, ψ=2.5, δ_annual=0.05, ρ_annual=1/0.975−1 (so β_annual=0.975).
  Previously held the wrong values (α=0.35, ψ=2.24, δ=0.0464, β=0.9722)
  inherited from a prior session's misreading. Cross-checked against the
  paper (`BCA_info.md` Table 1) and against `matlab_reference/datamine.m`
  lines 11–15.
- **`tests/test_params.py`** — updated to assert the BCKM Table 1 values.
  All 65 tests pass.
- **`bca_core/var_estimation.py`** — extended the `eval_only`
  short-circuit to return smoothed states, enabling Phase C
  (smoothing-at-BCKM-params diagnostic).
- **`scripts/phase_c_bckm_smooth.py`** (new) — feeds BCKM's published
  MLE parameters through our Kalman smoother and compares wedges to
  `worktemp.w`. Confirmed our smoother is not structurally broken.
- **`scripts/eval_bckm_basin.py`** (new) — evaluates BCKM Table 8 P
  and Table 10 Q at our pipeline's data offsets and reports ΔLL plus
  Q_chol diagonal ratios. Used to bound how far our basin sits from
  BCKM's after the calibration fix.
- **`CLAUDE.md`** — added an explicit "do not change these calibration
  values" block and the rationale (so the next session does not regress
  to the old wrong values).
- **`Diary.md`** — Step 9 closure section with the full assessment,
  decision matrix, and three-step plan that this document mirrors.

---

## 9. Files in scope

| Matlab (BCKM ground truth)                                                  | Our counterpart                                |
|------------------------------------------------------------------------------|------------------------------------------------|
| [`mleqadj.m`](BCKM/Multicountry%20-%20End/USAN2/mleqadj.m)                  | [`bca_core/var_estimation.py`](bca_core/var_estimation.py) |
| [`initmle.m`](BCKM/Multicountry%20-%20End/USAN2/initmle.m)                  | (Sbar fsolve-init logic, in `var_estimation.py`) |
| [`runmleadj.m`](BCKM/Multicountry%20-%20End/USAN2/runmleadj.m)              | `estimate_var_mle` in `var_estimation.py`      |
| [`res_adjust.m`](BCKM/Multicountry%20-%20End/USAN2/res_adjust.m)            | `PrototypeModel` in [`bca_core/model.py`](bca_core/model.py)  |
| [`fixexpadj.m`](BCKM/Multicountry%20-%20End/USAN2/fixexpadj.m)              | `solve_counterfactual` in [`bca_core/counterfactuals.py`](bca_core/counterfactuals.py) |
| [`gwedges2.m`](BCKM/Multicountry%20-%20End/USAN2/gwedges2.m)                | analytical extraction in [`bca_core/wedges.py`](bca_core/wedges.py) |
| [`fstats3.m`](BCKM/Multicountry%20-%20End/USAN2/fstats3.m)                  | f-stats block in [`scripts/run_var_counterfactuals.py`](scripts/run_var_counterfactuals.py) |
| [`usdata.m`](BCKM/Multicountry%20-%20End/USAN2/usdata.m)                    | [`bca_core/data/pipeline.py`](bca_core/data/pipeline.py) + [`bca_core/data/adjustments.py`](bca_core/data/adjustments.py) |
| [`datamine.m`](BCKM/Multicountry%20-%20End/USAN2/datamine.m)                | calibration in [`bca_core/params.py`](bca_core/params.py) |
| [`runall.m`](BCKM/Multicountry%20-%20End/USAN2/runall.m)                    | [`scripts/run_var_counterfactuals.py`](scripts/run_var_counterfactuals.py) |
| [`kfilter.m`](BCKM/Multicountry%20-%20End/USAN2/kfilter.m)                  | `_steady_state_kalman` in `var_estimation.py`  |

---

## 10. What this doc replaces

The previous version of this file (commit `8e3f283` and earlier)
contained a tiered Tier-1/2/3/4 analysis written when the calibration
constants were still wrong. The largest single error in that doc was the
"Tier 1.5" calibration table, which listed α=0.35, ψ=2.24, δ=0.0464,
β=0.9722 as BCKM defaults — those values do **not** appear in the paper
or in `datamine.m`. They appear to have been a prior session's
misreading that propagated into `params.py` defaults. Step 9's audit
caught it, the paper and Matlab code were re-checked, and the values
have been restored to BCKM Table 1 / `datamine.m`.

Several Tier-1 items in the original doc were already addressed before
that audit (Sbar parametrization, BCKM-faithful penalty, asymptotic
Kalman, Y₀-rebased f-stats, hours via PAYEMS×AWHNONAG, calgz-style
detrending). Those are noted in section 2 of this rewrite as "currently
working." The Tier-1 items that remain genuinely open are folded into
the unresolved question in section 5 — pending Step 1 of the plan in
section 6, we cannot tell which of them is load-bearing.
