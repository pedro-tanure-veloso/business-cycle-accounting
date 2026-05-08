# BCA Web — Session Diary

Append-only journal for `bca_web` and the `bca_data_builder` pipeline that
feeds it (`scripts/build_quarterly_data.py`). Modeled on
`covid_analysis/Diary.md` and the Findings section in root `CLAUDE.md`:
durable memory across context resets. Newest entry at the top. Do not
delete or rewrite earlier entries unless a later finding directly
invalidates them — in that case mark the old one `SUPERSEDED`.

---

## 2026-05-08 — Eval window moved to 2023Q1; flat-data f-stat degeneracy surfaced; events-log generator scaffolded

Follow-on to the 2026-05-07 entry. Verified the previous session's fixes
end-to-end via `workflow_dispatch`, then iterated on what the new data
revealed.

### Workflow side: one more fix needed before the cron could run

**`scripts/build_quarterly_data.py` failed CI with
`ModuleNotFoundError: No module named 'scripts'`.** The script does
`from scripts.run_bca import build_us_dataset, run_pipeline`, which
requires the repo root on `sys.path`. When invoked as
`python scripts/build_quarterly_data.py` (the workflow's invocation),
Python puts `scripts/` on `sys.path` instead. Fixed by inserting
`sys.path.insert(0, str(Path(__file__).resolve().parent.parent))` at
the top of the script (commit `d39c856`). The Node.js 20 deprecation
warning in the workflow log is unrelated noise.

### What the post-fix data revealed: f-stat degeneracy in flat-data regimes

After the workflow re-ran cleanly, the dashboard reported
**phi=(eff:0.17, lab:0.02, inv:0.08, gov:0.73)** over 2024Q1–2025Q3
— Government dominant by a 4× margin. But the cf-time-series panel told
a very different story: data essentially flat (~98–100), labor cf
peaking at +50 (output would have boomed under labor wedge alone),
investment cf bottoming at −22 (output would have crashed under
investment alone), efficiency cf at −15. Government cf was the only
one that stayed near data, hugging 100 throughout.

**Interpretation finding.** When observed output is approximately flat
within the eval window, the inverse-SSR f-stat mechanically favors
whichever counterfactual stayed *closest to flat* — i.e., the wedge
that did the **least**, not the wedge that "explains" most. The
Government=0.73 headline is a degenerate-ranking artifact, not
evidence that fiscal/external channels drove output. The interesting
signal in this regime is the **offsetting-wedge story** visible in the
right panel: labor and investment moved aggressively in opposite
directions, leaving observed output flat.

This is methodological, not a bug. The f-stat formula is correct
(BCKM Table 11 idiom). The degeneracy is intrinsic to inverse-SSR when
the denominator (data variation within the window) is small. Two
implications:
1. The dashboard's *visual* hierarchy should privilege the cf-time-
   series panel over the f-stat bar chart in flat-data regimes.
2. The hypothesis-layer LLM needs this context or it will dutifully
   write "Government wedge dominated" narratives when the truth is
   "labor and investment offset."

### Decisions taken this session

**1. Eval window moved 2024Q1 → 2023Q1.** With ~7 quarters in the
2024Q1+ window and observed output nearly flat, both the f-stat and
the cf rebase were under-conditioned. 2023Q1 widens the window to ~11
quarters and includes the 2023 banking-stress / Fed-tightening period,
which has more output variation. User's prior on 2024Q1 (agentic-AI
inflection) is still defensible but with this window length the f-stat
fragility wins. The robustness sweep over alternative start dates
(2022Q4 / 2023Q1 / 2024Q1) remains on the to-do list.

**2. Unified the eval-window definition.** Previously `phi_start_idx`
was anchored at 2024Q1 for the f-stat, while `recent_idx = -7` slug-
sliced the cf time series independently — two different windows by
construction. Replaced both with a single `eval_window_start =
pd.Timestamp("2023-01-01")` and derived `win_start_idx`, `win_end_idx`
from it; both panels are now provably computed off the same window.
[`scripts/build_quarterly_data.py:184-225`].

**3. Gemini prompt enriched with cf peak swings + methodological
caveat.** The previous prompt fed only f-stats and SD-from-mean of
wedges, so Gemini had no way to see "labor cf to 150, data flat = the
offsetting story." Now the prompt also includes:
   - Signed peak deviation of each cf from baseline within the window
   - Observed data range across the window
   - Explicit "read carefully before interpreting" block: when data is
     near-flat, lead with offsetting-swings narrative; use f-stats to
     corroborate, not as the headline.
[`scripts/build_quarterly_data.py:52-110`].

**4. UI: methodology caveat surfaced in the Wedge Decomposition
section.** Above the plots, a one-line dynamic "Evaluation window:
2023Q1 → {data.quarter}" header (renders 2025Q3 today; updates
automatically). Below the plots, a "Caveat:" block explaining the
inverse-SSR / flat-data degeneracy and pointing readers to the right
panel for the offsetting-swings story. Chart title also updated to
"Output Components (2023Q1 = 100)". [`bca_web/src/App.tsx:247-253`,
`:268`, `:295-297`].

**5. UI: defensive `.slice(-8)` cap on the demand chart.**
Belt-and-suspenders so the right-side chart can never exceed 8
quarters regardless of JSON length. [`bca_web/src/App.tsx:199`].

**6. Standalone events-log generator: `scripts/generate_events.py`.**
New script that, given a quarter (e.g. `2023Q1`), calls Gemini 2.5
Flash with **Google Search grounding** (`tools=[GoogleSearch()]`) and
appends a structured event log to `data/events.md`. Idempotent —
scans existing `## YYYYQX` headers and skips quarters already present
unless `--force`. Each quarter is generated **once and frozen**, so
older entries don't drift between runs. The prompt enforces five
fixed buckets (Monetary / Fiscal / Regulatory / Tech / Geopolitical),
≤3 entries per bucket, dated to within a week, factual-only (no
commentary or forecasts). The intent is to feed this log into a
future revision of the hypothesis-generation prompt so candidate
mechanisms are anchored to dated events rather than the LLM's
narrative-of-the-moment training prior.

  *Status*: scaffolded but not yet wired. To be tested on a separate
  machine where `GEMINI_API_KEY` is exported, then integrated into the
  workflow as a second pre-build step, then ingested by
  `generate_hypotheses` in `build_quarterly_data.py`.

### Open / to-do

- **Test `generate_events.py` against 2023Q1** on a machine with
  `GEMINI_API_KEY` exported. Inspect output for hallucination,
  bucket discipline, date precision. If acceptable, backfill
  2022Q1–2025Q3.
- **Wire `data/events.md` into `generate_hypotheses`.** Pass the
  current and prior quarter's event entries into the prompt;
  instruct Gemini to tie each candidate mechanism to a specific
  dated event whose quarter aligns with a wedge move.
- **Add `generate_events.py` to the GH Actions workflow** as a
  pre-`build_quarterly_data.py` step. Idempotent semantics mean it's
  safe to run on every cron firing.
- **f-stat robustness sweep** over 2022Q4 / 2023Q1 / 2024Q1 start
  dates (carried over from 2026-05-07).
- **Imports KPI sign at App.tsx:185** still shows green ↑ on rising
  imports (carried over).
- **YAxis units / chart titles** for the Demand chart (carried over).
- **`pyproject.toml` `google-genai>=0.1` pin tightening**: verify after
  next CI run that the resolver picks a version exposing
  `client.models.generate_content` and `tools=[GoogleSearch()]`. If
  flaky, bump to `>=0.3`.

### Files touched

- `scripts/build_quarterly_data.py`
  - L9-19: `import sys`; insert repo root on `sys.path` before any
    `scripts.*` / `bca_core.*` import (workflow `ModuleNotFoundError` fix)
  - L52-110: prompt extended with cf peak swings, data range, and
    methodological caveat
  - L184-225: unified `eval_window_start = 2023-01-01` shared by the
    f-stat anchor and the cf rebase
- `bca_web/src/App.tsx`
  - L199: `.slice(-8)` defensive cap on demand chart
  - L247-253: dynamic Evaluation-window line
  - L268: chart title 2024Q1 → 2023Q1
  - L295-297: "Caveat:" block below the wedge plots
- `scripts/generate_events.py` (new): events-log generator with Google
  Search grounding, idempotent per quarter

---

## 2026-05-07 — Three live bugs in the dashboard, two fixed in this session

Triage of the live `latest_quarter.json` (auto-updated by GH Actions on
2026-05-07) surfaced three real problems and one cosmetic stale type.

### Fixed in this session

**1. Output-fluctuation decomposition was using `phi_statistics` on
absolute log-deviation values without anchoring.** The dashboard's
2024Q1+ window slices each cf path's accumulated drift since the MLE
sample start (2010Q1) into the SSR, so the inverse-SSR weights were
dominated by the level offset between data and each cf at the start of
the window — exactly the failure mode `f_statistics_bckm`'s docstring
warns about. Empirically: reported phi=(eff:0.00, lab:0.00, inv:0.13,
gov:0.87), recomputed-from-the-same-cf_time_series f=(eff:0.61,
lab:0.02, inv:0.07, gov:0.29). The "Government 87%" headline was an
artifact of the SSR being dominated by gov_hat's start-of-window level
matching data_hat's, not within-window shape.

  *Fix*: `scripts/build_quarterly_data.py` now imports and calls
  `f_statistics_bckm(r["data_hat"], r["cfs"], window=(i1, i2),
  anchor=i1)` with absolute window indices into the unsliced arrays.
  This rebases each path to the window start before computing SSR —
  matches what the cf_time_series chart visualizes (rebase-to-100 at
  2024Q1) and is the BCKM-faithful formula per vision-doc Panel 2B
  (which explicitly says "f-statistic," not phi).

  *Methodology note*: the 7-point eval window means the f-stat is
  fragile to start-date choice. User picked 2024Q1 on the prior that
  agentic-AI macro impact starts ~2024. Defensible, but a robustness
  sweep over 2022Q4 / 2023Q1 / 2024Q1 is on the to-do list (see
  Open below).

**2. Hypothesis Layer silently broken since 2026-05-07
(commit `a925593`).** That "Fix BCA pipeline data sources and UI
mapping" commit wholesale-rewrote `build_quarterly_data.py` and
clobbered the `generate_hypotheses` function: the previous working
version used `from google import genai` (the new `google-genai` SDK,
which is what `pyproject.toml` actually pins) with model
`gemini-2.5-flash` and structured-JSON output via
`response_mime_type="application/json"`. The post-rewrite version uses
the legacy `google.generativeai` SDK (not installed) with model
`gemini-2.0-flash`, which fails on `ImportError` in CI; the
`except Exception: return {}` swallows it silently and the dashboard
section is gated off. `pattern_identification` has been empty in the
production JSON for that whole window.

  *Fix*: restored the `f33aaa8`-era `generate_hypotheses` —
  `google-genai` SDK, `gemini-2.5-flash`, structured-JSON output.
  Prompt scoped to phi-stats + wedge SD levels (no full-payload JSON
  dump). `GEMINI_API_KEY` is already in GitHub Actions secrets per
  user; the next monthly cron (or manual `workflow_dispatch`) should
  populate the layer.

  *To verify after next CI run*: pin `google-genai` more tightly if
  the current `>=0.1` happens to resolve to a version older than
  whatever ships `client.models.generate_content`. Worth tightening
  to `>=0.3` if the next run still fails.

**3. Demand Contributions chart spans 2020Q2 → 2026Q1 (24 quarters)**,
which compresses everything post-COVID into a flat squiggle. User
wanted 8 quarters.

  *Fix*: `df_fred.tail(24)` → `df_fred.tail(8)` in
  `build_quarterly_data.py`.

### Cosmetic

**4. Stale TypeScript interface for `demand_contributions`.**
`bca_web/src/App.tsx:47` declared the row keys as
`Consumption | Investment | Government | "Net Exports"` while the JSON
keys are `A_Consumption | B_Investment | C_Government | D_Exports`
(A/B/C/D prefix forces Recharts stack order — commit `88033b1`). The
`<Bar>` elements in App.tsx already used the prefixed keys correctly, so
the chart was rendering fine; the type was just lying.

  *Fix*: interface updated to match JSON.

### Confirmed-not-a-bug

**5. Demand Contributions numbers are correct, just unlabeled.** User
initially read 2020Q3=+34.9 / 2020Q3-Consumption=+24.95 as implausible.
They're SAAR percent contributions to annualized GDP growth, fetched
directly from FRED's `*RY2Q224SBEA` and `*RL1Q225SBEA` tickers. Spot
check vs BEA: 2020Q3 SAAR GDP=+33.8, PCE contribution=+25.2pp — the
shown numbers are within revision noise. They look implausible because
the chart title is "Demand Contributions to Growth" with no SAAR
annotation while the KPIs above are labeled "QoQ annualized," forcing
the eye to compare two slightly different presentations of the same
unit. Recommend (not done): add `unit="%"` to the chart YAxis and clarify
the title to include "(SAAR, %)".

### Open / to-do

- **f-stat robustness sweep over start dates** — 2022Q4 (ChatGPT),
  2023Q1 (GPT-4), 2024Q1 (current default). User's prior is 2024Q1
  is right because agentic AI doing real economic work is a 2024+
  phenomenon, but with only 7 points the f-stat is fragile to the
  start-date choice. Worth checking whether the dominant wedge
  survives ±2 quarters of slack on the anchor. Possible
  implementation: a `--phi-start` CLI flag on
  `scripts/build_quarterly_data.py` (default 2024Q1) so the
  robustness sweep is a few `python scripts/build_quarterly_data.py
  --phi-start 2023Q1` runs from the command line — without touching
  the dashboard's fixed default. **NOT a UI control** — vision-doc
  Stage 3 explicitly forbids window pickers in `bca_web`.
- **Conceptual note on AI and the wedges.** AI-as-TFP-shock should
  land in the **efficiency** wedge, not investment. The investment
  wedge is a friction on the Euler equation — it lights up when
  investment is unusually low *given fundamentals*, not when capex
  booms. AI raising the marginal product of capital makes investment
  boom *without* moving τ_x; the f-stat decomposition should pick
  this up as efficiency, which the new f-stat formula (eff=0.61
  recomputed) does. Useful sanity check that the fix is producing
  economically sensible numbers.
- **Imports KPI at App.tsx:185** shows green ↑ on rising imports,
  but rising imports subtract from GDP. Cosmetic / semantic. Not
  fixed in this session.
- **YAxis units / chart titles** — see #5 above. Not done.
- **`pyproject.toml` `google-genai>=0.1` may need tightening** to
  whatever version exposes `client.models.generate_content` reliably.
  Verify after next CI run.

### Files touched

- `scripts/build_quarterly_data.py`
  - L19 import: `phi_statistics` → `f_statistics_bckm`
  - L29-78 `generate_hypotheses`: rewritten on `google-genai` SDK,
    `gemini-2.5-flash`, structured-JSON output
  - L148-161: f-stat with absolute window indices anchored at 2024Q1
  - L200: demand chart `tail(24)` → `tail(8)`
- `bca_web/src/App.tsx:47` — TS interface aligned with JSON keys
