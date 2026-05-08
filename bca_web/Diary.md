# BCA Web — Session Diary

Append-only journal for `bca_web` and the `bca_data_builder` pipeline that
feeds it (`scripts/build_quarterly_data.py`). Modeled on
`covid_analysis/Diary.md` and the Findings section in root `CLAUDE.md`:
durable memory across context resets. Newest entry at the top. Do not
delete or rewrite earlier entries unless a later finding directly
invalidates them — in that case mark the old one `SUPERSEDED`.

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
