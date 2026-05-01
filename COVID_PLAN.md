# COVID Smoke Test: Validate BCA Pipeline Generalizability on US 2010Q1–2023Q4

> **Status (2026-05-01):** approved plan; work is queued, not started.
> All reorganization prerequisites (BCKM 1980-2014 closure, repo
> `bckm_replication/` move, BCA_info.md split) are done.
>
> **Setup on a fresh machine:**
> ```bash
> git pull
> pip install -e ".[dev]"
> # First run only — fetches FRED series into ~/.bca_cache/fred/
> export FRED_API_KEY=<your key>
> ```
> The 1980-2014 BCKM regression parquet
> (`bckm_replication/data/us_1980_2014_calgz.parquet`) is committed and
> needs no rebuilding. The two COVID parquets are built fresh in
> Step 2 below.

## Context

The previous session's work closed out the BCKM 1980Q1–2014Q4 replication
on the data layer: every BEA-NIPA toggle makes the LL gap monotonically
worse; FRED-default is the production setting; full forensics in
`bckm_replication/REPORT.md` and `bckm_replication/DATA_FORENSICS.md`. At
BCKM-θ on our df, f-stats already match Table 11 to ~0.01 in every
channel. The remaining 1980-2014 gap is a basin-escape problem (optimizer
drifts away from BCKM's Sbar attractor), intellectually interesting but
not on the critical path for the project's actual goal.

The reframed goal (user-validated 2026-05-01): **confidence that our code
produces credible BCA decompositions for arbitrary country/period
combinations**, not exact numerical match to BCKM 2014. The cheapest, most
informative validation is to apply the existing pipeline to a different US
window where we have strong narrative priors and see whether the wedge
decomposition agrees with them.

The COVID era is the obvious smoke test:
- **2020Q2 lockdowns**: should appear as a massive negative shock to
  **τ_l** (labor wedge). Output−hours collapse together by construction,
  so the efficiency channel is mechanical until the labor side resolves.
- **2020Q3–2021 recovery**: BLS reports record TFP growth — **+3.2% in
  2021, the largest since 1983** — driven by output recovering faster
  than capital re-investment / labor re-hiring. Our pipeline should show
  **A** spike up here, not down. This is the *strongest* prior on the
  whole exercise: if A shows persistently negative through 2021 in our
  decomposition, the pipeline has a labor/efficiency separation bug that
  the BCKM 1980-2014 window happened to mask.
- **2020–2021 fiscal**: CARES Act + ARPA expand transfers and federal
  spend by ~10% of GDP cumulatively. Should appear as a positive **g**
  wedge.
- **2023Q4 normalization**: BLS reports TFP back near trend (~1.5% in
  2024) by end of 2023. Wedges should be returning to pre-COVID
  baseline — τ_l ≈ 1, A near 2019Q4 trend, g elevated but cooling.

If the pipeline produces this qualitative pattern, we have *real* evidence
that the code generalizes. If it doesn't, we know exactly which channel
has a hidden bug that the BCKM window happened to mask.

## Methodology choices (user-decided 2026-05-01)

- **Window**: **2010Q1–2023Q4** (post-GR start; ends at 2023Q4 to avoid
  any data extrapolation — `LFWA64TTUSQ647N` working-age population
  publishes through 2023Q4 cleanly).
- **Bind**: **2019Q4** (NBER cycle peak; natural COVID counterfactual
  anchor).
- **Trend handling**: **both side-by-side** — full-window calgz as
  headline (BCKM-faithful methodology) + pre-COVID-fit (2010Q1–2019Q4
  then extrapolate forward) as robustness check.
- **Data sources**: FRED defaults
  (`y_source=x_source=g_source="fred"`). No BEA-branch parameterization
  (the 1980/2014 hardcodes in `bca_core/data/bea.py:195-273` and the BEA
  paths in `bca_core/data/pipeline.py:224,246-247,286-287` stay).

## Implementation steps

### 1. Add an `mle_window` parameter to the trend fit

`bca_core/data/adjustments.py:278-330` `remove_trend()` currently fits
the calgz slope on the full sample passed in. Add
`mle_window: tuple[str, str] | None = None` so a caller can specify "fit
slope on this sub-window, apply to full window." Default behavior
(`None` → full-window) preserves the BCKM 1980-2014 path bit-for-bit.

Plumb the same parameter through `bca_core/data/pipeline.py:26`
`build_us_dataset(mle_window=None, ...)` so callers can switch between
full-window and pre-COVID-fit detrending without copy-pasting.

### 2. Two cached datasets for the COVID window

Build and cache:
- `data/us_2010_2023_calgz.parquet` — full-window calgz (headline).
- `data/us_2010_2023_calgz_preCOVID.parquet` — pre-COVID-fit
  (robustness).

Reuse `build_us_dataset(start="2010Q1", end="2023Q4",
base_year_quarter="2019Q4", ...)` for both; the second pass adds
`mle_window=("2010Q1", "2019Q4")`. No data adapter changes — every FRED
ticker (including `LFWA64TTUSQ647N`) covers the full window without
extrapolation.

### 3. New driver script: `scripts/run_covid_smoke.py`

Single entry point that:
1. Loads both COVID parquets.
2. Runs `estimate_var_mle` on each (warm-start from BCKM Table 8 P/Q is
   fine — they're sensible priors even on a different window; if
   convergence fails, fall back to OLS-VAR seed).
3. Computes wedges via the existing path (no changes needed).
4. Calls `f_statistics_bckm` with explicit `anchor=idx_of_2019Q4` and a
   COVID counterfactual window (e.g. **2019Q4 → 2022Q4**, peak →
   recovery established). The hardcoded GR window indices in
   `bckm_replication/scripts/eval_bckm_fstats.py:72-74` should not be
   reused — compute window indices dynamically from the parquet's date
   index.
5. Emits four headline plots (`wedges_us_2010_2023.png`,
   `figure_2{B,C,D}_covid.png`) by adapting the patterns in
   `bckm_replication/scripts/plot_wedges.py:108` and
   `bckm_replication/scripts/plot_figure_2_at_bckm_theta.py:201-203`.
6. Prints a side-by-side table of wedge values at four reference
   quarters (2019Q4 = bind = 1.0, 2020Q2 = trough, 2021Q4 = recovery,
   2023Q4 = normalization) under both trend variants.

### 4. Validation rubric

Compare the printed table against the narrative priors:

| Reference | Prior on wedge direction | Tolerance |
|---|---|---|
| 2020Q2 | τ_l strongly negative; A small / mechanical | sign |
| 2021Q4 | A above pre-COVID trend (record TFP); g positive | sign |
| 2023Q4 | τ_l ≈ 1 (recovered); A near 2019Q4 trend; g still elevated | sign |

A "pass" is qualitative agreement on **all three** rows under **both**
trend variants. If the two trend variants disagree on direction at any
reference row, that itself is an interesting result worth documenting
(means the COVID anomaly is large enough that detrending methodology
matters — useful for the web app's vintage/methodology UX).

### 5. Diary entry (`Diary.md`)

Append a session entry capturing:
- The reframing from "match BCKM exactly" to "validate generalizability"
  (2026-05-01 user pivot).
- The decision rationale for COVID smoke test as the highest-leverage
  next step.
- The methodology choices (window 2010-2023, bind 2019Q4, dual trend
  variants).
- The TFP narrative priors (BLS: +3.2% in 2021, ~1.5% in 2024,
  cycle-average 1.0% 2019-2025 vs 0.6% 2007-2019).
- The window-end choice (2023Q4) to avoid `LFWA64TTUSQ647N`
  extrapolation.
- Reference forward to wherever results land (placeholder section to be
  filled in after the run).

### 6. Tests

Add `tests/test_covid_smoke.py` with:
- Shape / no-NaN checks on both parquets.
- A re-run of the algebraic-identity tests T1–T7 (already
  window-agnostic; see `tests/test_counterfactuals.py`) on the 2010-2023
  dataset to confirm structural correctness on the new window.
- A regression check that `wedges_us_2010_2023.png` and the f-stat
  printout exist after running `run_covid_smoke.py` (smoke-only; we
  don't pin numerical values yet).

Window-specific tests (`test_bckm_table12.py`, the 1980-2014 index
assertions in `test_bckm_reference.py`) are NOT touched — they remain
BCKM-replication regression tests and should keep passing on their
original window.

## Critical files

- **Edit**:
  - `bca_core/data/adjustments.py` — `mle_window` parameter on
    `remove_trend()`
  - `bca_core/data/pipeline.py` — plumb `mle_window` through
    `build_us_dataset`
  - `Diary.md` — session entry
- **New**:
  - `scripts/run_covid_smoke.py` — single driver
  - `tests/test_covid_smoke.py` — smoke + algebraic identities on new
    window
  - `data/us_2010_2023_calgz.parquet` (cached)
  - `data/us_2010_2023_calgz_preCOVID.parquet` (cached)
  - `wedges_us_2010_2023.png`, `figure_2{B,C,D}_covid.png`
- **Reuse unchanged**:
  - `bca_core/data/fred.py` — no changes needed (2023Q4 endpoint avoids
    the working-age extrapolation issue)
  - `bca_core/var_estimation.py` (window-agnostic)
  - `bca_core/counterfactuals.py` — `f_statistics_bckm` accepts `anchor`
    kwarg at line 300; compute `anchor=idx_of_2019Q4` from the date
    index
  - All Klein / state-space / wedge-extraction code

## Verification

```bash
# Build datasets and run the smoke test end-to-end
python scripts/run_covid_smoke.py 2>&1 | tee /tmp/covid_smoke.log

# Re-run all tests (window-agnostic + new COVID smoke + existing BCKM
# regression). Existing 79/79 fast tests must still pass; new COVID
# tests pass on first run.
pytest tests/ -v --tb=short

# Eyeball the four headline plots:
#   wedges_us_2010_2023.png   — 4-panel wedge time series
#   figure_2B_covid.png       — counterfactual y under each wedge
#   figure_2C_covid.png       — wedge contributions to peak-trough
#   figure_2D_covid.png       — f-stat decomposition over COVID window
```

The deliverable is the four plots + the wedge-table printout + a Diary
entry documenting whether the pipeline passed or failed the
narrative-prior rubric under both trend variants.

## Out of scope (explicitly)

- BEA-branch migration to non-1980-2014 windows. The 1980/2014 hardcodes
  in `bca_core/data/bea.py:195-273` and the BEA paths in
  `pipeline.py:224-287` stay. Re-evaluate after the FRED-default smoke
  test lands.
- Basin-escape / optimizer-stability fixes. Orthogonal to this work.
- Web-app productionization (FastAPI, React). Premature until
  generalizability is validated.
- Cross-country extension. Same rationale.
- Any change to BCKM-replication regression tests on the 1980-2014
  window.
- Any data extrapolation. The 2023Q4 endpoint is chosen specifically to
  avoid the `LFWA64TTUSQ647N` post-2023 gap; if the user later wants to
  push to 2024Q4+, that's a separate decision with its own methodology
  trade-offs.
