# COVID Analysis — Session Diary

Append-only journal for Layer 2 work (COVID-era smoke test and beyond).
One entry per session. Format: date, what was done, what is in progress,
what is blocked, exact next step.

---

## 2026-05-01 — Session: COVID smoke test scaffolding

### What was done
- **Reframing confirmed** (user-validated): project goal is generalizability
  validation, not exact BCKM numerical match. The 1980-2014 replication is
  closed; all further work targets Layer 2 (arbitrary windows).
- **Methodology decisions** (user-decided):
  - Window: 2010Q1–2023Q4 (post-GR start; avoids LFWA extrapolation gap).
  - Bind: 2019Q4 (NBER cycle peak; natural COVID counterfactual anchor).
  - Trend: both full-window calgz and pre-COVID-fit (2010Q1–2019Q4 slope,
    applied forward) run side-by-side.
  - Data: FRED defaults (`y/x/g_source="fred"`). No BEA branches.
  - `labor_target_mean=None` (raw hours per capita, no BCKM 0.24279 anchor).
- **Infrastructure built**:
  - `covid_analysis/` folder created (mirrors `bckm_replication/` layout).
  - `covid_analysis/scripts/run_covid_analysis.py` — full driver: builds
    datasets, runs MLE (warm-started from BCKM θ), extracts wedges,
    computes counterfactuals + f-stats, prints wedge table, emits 5 figures.
  - `tests/test_covid_analysis.py` — three test classes: dataset shape/
    metadata, structural algebraic identities, smoke outputs.
  - `covid_analysis/REPORT.md` — placeholder (to be filled after first run).
  - `CLAUDE.md` updated: `covid_analysis/Diary.md` rule added under
    Progress Tracking.
- **Note on Step 1**: `mle_window` parameter already existed in
  `build_us_dataset()` — no pipeline changes needed.

### BLS TFP narrative priors (validation rubric)
- 2020Q2: τ_l should collapse (COVID lockdown labor shock).
- 2020Q2: A should be small/mechanical (output and hours collapse together).
- 2021Q4: A should spike above pre-COVID trend (BLS: +3.2% TFP in 2021,
  largest since 1983 — output recovered faster than capital/labor re-hiring).
- 2021Q4: g should be elevated (CARES Act + ARPA ≈ 10% of GDP cumulative).
- 2023Q4: τ_l ≈ recovered; A near 2019Q4 trend; g still elevated but cooling
  (BLS: TFP ~1.5% in 2024, back near the 1.0% cycle-average 2019–2025).

### What is in progress
- First actual run of `run_covid_analysis.py` (pending FRED fetch).
- Eyeball of the four figures against narrative priors.

### What is blocked
- Nothing structural. FRED_API_KEY is required for the first run.

### Exact next step
1. Run `FRED_API_KEY=... .venv/bin/python covid_analysis/scripts/run_covid_analysis.py`
2. Check the narrative-prior rubric printout (PASS/FAIL per row).
3. Eyeball `covid_analysis/figures/wedges_us_2010_2023.png` — τ_l should
   show the 2020Q2 collapse; A should show the 2021 spike.
4. Fill in `covid_analysis/REPORT.md` with actual results.
5. Run `pytest tests/ -v --tb=short` — all 79 existing + new COVID tests
   should pass.
