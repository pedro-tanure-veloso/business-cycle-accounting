# Step 9 — Two-Machine Handoff Plan

This document coordinates the two-machine experiment that closes Step 9.
Read it on whichever machine you are picking up the work, then follow
your column.

## Goal

Resolve the unresolved question in [DIVERGENCE_ANALYSIS.md](DIVERGENCE_ANALYSIS.md) §5:

> Does running BCKM's MATLAB on their `data.mat` reproduce Tables 8–11
> exactly, or does it land in one of several local optima?

Three possibilities — (A) deterministic reproduction, (B) approximate
reproduction, (C) multi-basin under-identification — and we walk the
decision matrix in §7 of that doc once we have evidence in hand.

The work is split across two machines so it can run in parallel:

- **Other machine** (origin of the project, where `BCKM/` and `data.mat` live): runs Step 1 (Octave on BCKM's MATLAB).
- **This machine** (where Step 9 closure was written): runs Steps 2–3 (bootstrap + calibration sensitivity in our Python pipeline).

We converge on this machine, synthesize results, write the verdict.

---

## Phase 1 — establish shared base (already done)

`pedro_macbookpro` has been fast-forward-merged to `main` and pushed,
so both machines branch off `main` going forward. Verify with
`git log --oneline main` — top commit should be the Step 9 closure.

---

## Phase 2 — each machine forks a task branch off `main`

### Other machine — Step 1 (Octave replication)

```bash
git fetch origin
git checkout main
git pull
git checkout -b step1-octave-replication
```

### This machine — Steps 2 & 3 (bootstrap + sensitivity)

```bash
git checkout main
git pull   # picks up Phase 1 merge
git checkout -b step23-bootstrap-sensitivity
```

---

## Phase 3 — work in parallel

**The single rule that prevents merge conflicts: don't both edit `Diary.md` or `DIVERGENCE_ANALYSIS.md` during this phase.** Each task writes to its own dedicated results file. Synthesis to Diary/Divergence happens at convergence.

### Other machine: Step 1

**Pre-flight check** before anything else. The repo already contains
the BCKM Matlab files we need — they live in **`matlab_reference/`**
(tracked in git; the larger gitignored `BCKM/` folder is the full
replication package and is not required for Step 1):

```bash
ls matlab_reference/{runall.m,mleqadj.m,runmleadj.m,kfilter.m,fstats3.m,initmle.m,data.mat,worktemp.mat}
```

All of those should resolve. If anything is missing, fall back to a
physical copy of `BCKM/Multicountry - End/USAN2/` from wherever the
project originated.

**What to do:**

1. Install Octave if not present (`brew install octave` on macOS).
2. `cd` into `matlab_reference/`.
3. Run BCKM's pipeline end-to-end:
   ```bash
   octave --no-gui runall.m   # or whichever script drives the full pipeline
   ```
   If `runall.m` doesn't take you all the way to f-stats, look for the
   script that calls `mleqadj` → estimates P/Sbar/Q → calls `fstats3`
   (likely `runmleadj.m`).
4. **Capture** (essential — this is the comparison data):
   - Final P matrix (4×4)
   - Final P₀ vector (4)
   - Final V (= Q·Q', 4×4 symmetric)
   - Final Sbar (4)
   - f-statistics (fY[A], fY[τ_l], fY[τ_x], fY[g])
   - Final log-likelihood
   - Save the raw outputs: `save -mat octave_output/run_default.mat`
5. **Probe multimodality** (this is what tells us A vs B vs C):
   - Re-run from at least 5 random warm-starts (perturb Sbar/P/Q in
     `runmleadj.m` by ~5%).
   - Capture each run's converged P, P₀, V, f-stats, LL.
   - If they all converge to the same Tables 8–11 values → outcome (A).
   - If most converge near Tables 8–11 with small spread → outcome (B).
   - If they converge to economically different basins → outcome (C).
6. Write findings to **`STEP1_OCTAVE_RESULTS.md`** at the repo root.
   Suggested structure:
   - § Run setup (Octave version, machine, what scripts were called).
   - § Default-start results: P, P₀, V, f-stats, LL — side-by-side
     with BCKM Tables 8/9/10/11.
   - § Multi-start results: table of N runs × {LL, fY[A], fY[τ_l],
     fY[τ_x], P diag} showing the spread.
   - § Verdict on A/B/C with a one-paragraph justification.
   - § Anomalies / things to flag for the synthesis on the other
     machine.
7. Save raw `.mat` outputs in `octave_output/` (gitignored OK if
   they're large; otherwise commit small summary `.mat`s).
8. Commit and push:
   ```bash
   git add STEP1_OCTAVE_RESULTS.md octave_output/  # adjust as needed
   git commit -m "Step 1: BCKM Octave replication — verdict (A|B|C)"
   git push -u origin step1-octave-replication
   ```

**Files this branch should NOT touch:** `Diary.md`,
`DIVERGENCE_ANALYSIS.md`, `bca_core/`, `tests/`, `scripts/`. If
something needs fixing in those, capture the observation in
`STEP1_OCTAVE_RESULTS.md` and surface it at convergence.

### This machine: Steps 2 & 3

**Step 2 — bootstrap our estimator.**

- New script `scripts/bootstrap_mle.py`:
  - For N ≥ 100 random seeds, perturb the Table-77 warm-start by
    `std=0.01` (the BCKM-faithful perturbation per CLAUDE.md), run
    full MLE, record the converged {P, P₀, Q_chol, LL, f-stats}.
  - Save to `data/bootstrap_results.npz`.
  - Plot histograms of fY[A], fY[τ_l], fY[τ_x], fY[g], LL.
- Question to answer: is *our* objective unimodal or multimodal? Are
  the labor-dominated (LL=1830) and investment-dominated (LL=1819)
  basins reachable from random restarts, or are there more?

**Step 3 — calibration sensitivity.**

- New script `scripts/sensitivity_calibration.py`:
  - Vary α, ψ, δ_annual, β_annual one at a time by ±5%.
  - For each setting, run a single MLE from the standard warm-start.
  - Record converged {P, P₀, f-stats, LL, ss_new}.
- Question to answer: how much do the f-statistic conclusions move
  under normal calibration uncertainty? If they swing wildly, the
  framework's economic conclusions are fragile even before we worry
  about the BCKM gap.

**Write findings to** `STEP23_BOOTSTRAP_SENSITIVITY.md` at the repo
root. Structure:
- § Bootstrap: histogram summary, basin count, distance between
  basins, how often each basin is hit from random starts.
- § Sensitivity: table of `(parameter, ±5%) → Δf-stats, ΔLL`.
- § Implications for forward US analysis (does the framework give
  stable conclusions or does it depend on warm-start luck?).
- § Pre-synthesis observations to combine with Step 1 findings.

**Files this branch should NOT touch:** `Diary.md`,
`DIVERGENCE_ANALYSIS.md`, `STEP1_OCTAVE_RESULTS.md`,
`bca_core/params.py` (don't change calibration defaults — sensitivity
runs override per-call). Capture observations and defer.

Commit and push:
```bash
git add STEP23_BOOTSTRAP_SENSITIVITY.md scripts/bootstrap_mle.py scripts/sensitivity_calibration.py
git commit -m "Steps 2-3: bootstrap + calibration sensitivity"
git push -u origin step23-bootstrap-sensitivity
```

---

## Phase 4 — converge on this machine

After both branches are pushed:

```bash
git checkout main
git pull
git merge step1-octave-replication
git merge step23-bootstrap-sensitivity
```

The two branches touched disjoint files (different result docs,
different scripts, no shared edits to Diary/Divergence/core code), so
no conflicts.

**Synthesis:**

1. Read `STEP1_OCTAVE_RESULTS.md` and `STEP23_BOOTSTRAP_SENSITIVITY.md`
   together.
2. Walk the decision matrix in [DIVERGENCE_ANALYSIS.md](DIVERGENCE_ANALYSIS.md) §7 with
   the actual outcomes: pick the row, read off the verdict.
3. Append a single **Step 9 verdict** entry to `Diary.md` (one entry,
   one machine, no contention).
4. Update `DIVERGENCE_ANALYSIS.md` with the resolved A/B/C and the
   implications for forward post-2014 US analysis.
5. Commit, push to `main`, delete the two task branches:
   ```bash
   git push origin :step1-octave-replication :step23-bootstrap-sensitivity
   git branch -D step1-octave-replication step23-bootstrap-sensitivity
   ```

---

## Boundaries summary

| concern                 | other machine (Step 1) | this machine (Steps 2–3)        |
|-------------------------|------------------------|---------------------------------|
| writes new file         | `STEP1_OCTAVE_RESULTS.md` + `octave_output/` | `STEP23_BOOTSTRAP_SENSITIVITY.md`, `scripts/bootstrap_mle.py`, `scripts/sensitivity_calibration.py` |
| edits `Diary.md`        | no                     | no (only at convergence here)   |
| edits `DIVERGENCE_ANALYSIS.md` | no              | no (only at convergence here)   |
| edits `bca_core/`       | no                     | no (sensitivity overrides per call) |
| edits `tests/`          | no                     | only if a new script needs a test |

If either side needs to break a boundary, stop and coordinate before
committing — the cost of a merge conflict on `Diary.md` is much
higher than the cost of a quick check.

---

## Quick reference for the other machine starting Claude Code fresh

When you launch Claude Code on the other machine after the pull:

1. It reads `CLAUDE.md` automatically — gets project conventions.
2. Show it this file: "read `STEP9_HANDOFF.md` and execute the
   'Other machine: Step 1' section."
3. It can consult `Diary.md` and `DIVERGENCE_ANALYSIS.md` for
   background but should not edit them.
4. The Octave run itself does not need Claude Code's help — Claude
   Code is most useful for capturing/structuring the results into
   `STEP1_OCTAVE_RESULTS.md` and walking the multi-start probe.
