"""Isolate counterfactual machinery from data construction.

Feed BCKM's own ``Y_raw`` through our ``extract_wedges_bckm_style`` +
``run_all_counterfactuals`` at BCKM-published θ, then compare the
resulting per-wedge level-ratio paths to BCKM's stored
``worktemp.components`` (``mzy, mly, mxy, mgy, mzh, mlh, ...``).

This is the missing test in the worktemp walkdown: Stage 5 (component
decomposition). ``diag_worktemp_compare.py`` only goes to Stage 3
(smoothed wedges) and ``diag_bckm_data_isolation.py`` stops at the LL
scalar — neither runs the counterfactual machinery on BCKM data.

Three possible outcomes:

  (1) MATCH — counterfactual machinery is clean. The only model-side
      residual at identical (data, θ) is the 16.6-nat L gap from
      ``diag_bckm_data_isolation.py``. Fix path is data construction.

  (2) DIFFER IN LEVEL — bug in ``run_counterfactual`` (the realized-wedge
      feed, or the ``(C_j − C_0)`` subtraction).

  (3) DIFFER IN SHAPE — bug in ``solve_counterfactual`` (P_y application,
      capital LOM in CF mode, or ``bckm_state_space_cf`` itself).

BCKM Y_raw column order is ``[y, x, h, g, c, c_implied]``; ours is
``[y, l, x, g]``. Permutation ``(0, 2, 1, 3)`` maps BCKM → ours.

BCKM gwedges2.m:170 builds the per-wedge level ratios as
``mzy(t) = 100·exp(YMz_y(t) − YMz_y(Y0))`` where Y0 = bind = 2008Q1.
We do the same to our hat-space CF outputs.

Read-only — writes nothing to disk. Run from repo root.
"""
from __future__ import annotations

import numpy as np

from bca_core.bckm_reference import load_bckm_reference
from bca_core.constants import (
    P_BCKM_TABLE8 as P_BCKM,
    QCHOL_BCKM_TABLE10 as QCHOL_BCKM,
    SBAR_BCKM_TABLE8 as SBAR_BCKM,
)
from bca_core.counterfactuals import run_all_counterfactuals
from bca_core.data.pipeline import build_us_dataset
from bca_core.model import PrototypeModel
from bca_core.params import CalibrationParams
from bca_core.var_estimation import estimate_var_mle
from bca_core.wedges import extract_wedges_bckm_style


# Mapping: (our wedge name, our variable) → BCKM components column name.
# BCKM convention: m{wedge_initial}{var_initial} where wedge ∈ {z,l,x,g}
# and var ∈ {y,h,x,c}. Our ``l`` (labor) is BCKM's ``h`` (hours).
_COMP_MAP = {
    ("efficiency", "y"): "mzy",
    ("efficiency", "l"): "mzh",
    ("efficiency", "x"): "mzx",
    ("labor",      "y"): "mly",
    ("labor",      "l"): "mlh",
    ("labor",      "x"): "mlx",
    ("investment", "y"): "mxy",
    ("investment", "l"): "mxh",
    ("investment", "x"): "mxx",
    ("government", "y"): "mgy",
    ("government", "l"): "mgh",
    ("government", "x"): "mgx",
}


def _hat_to_level_ratio(hat_series: np.ndarray, bind: int) -> np.ndarray:
    """BCKM gwedges2.m:170 — convert hat-deviations to bind-anchored level
    ratios, scaled to 100. ``hat_series`` is ``y_cf`` from
    ``run_counterfactual`` (log-deviation from SS, or for single-wedge CFs,
    the BCKM incremental log-deviation ``(H_j − H_0) @ state``).
    """
    return 100.0 * np.exp(hat_series - hat_series[bind])


def _compare(label: str, ours: np.ndarray, theirs: np.ndarray, bind: int,
             gr_window: tuple[int, int]) -> dict:
    """Print full-sample and GR-window diff statistics. Returns a dict of
    summary scalars for the verdict table.
    """
    diff = ours - theirs
    full_max = float(np.max(np.abs(diff)))
    full_mean = float(np.mean(np.abs(diff)))
    full_max_t = int(np.argmax(np.abs(diff)))

    i1, i2 = gr_window
    gr_diff = diff[i1:i2 + 1]
    gr_max = float(np.max(np.abs(gr_diff)))
    gr_mean = float(np.mean(np.abs(gr_diff)))
    gr_max_t = i1 + int(np.argmax(np.abs(gr_diff)))

    bound_full = "MATCH" if full_max < 0.5 else ("close" if full_max < 2.0 else "DIFFER")
    bound_gr = "MATCH" if gr_max < 0.5 else ("close" if gr_max < 2.0 else "DIFFER")

    print(f"  {label:8s}  full: max={full_max:6.3f}pp at t={full_max_t:3d} "
          f"mean={full_mean:5.3f}pp [{bound_full:6s}]   "
          f"GR: max={gr_max:6.3f}pp at t={gr_max_t:3d} "
          f"mean={gr_mean:5.3f}pp [{bound_gr:6s}]")
    return dict(full_max=full_max, full_mean=full_mean,
                gr_max=gr_max, gr_mean=gr_mean)


def main():
    print("=" * 84)
    print("BCKM components isolation: feed BCKM Y_raw → our CFs → compare to mzy/mly/...")
    print("=" * 84)

    bckm = load_bckm_reference()
    bind = bckm.bind
    print(f"\n  worktemp.mat: T={len(bckm.time)}  bind={bind} ({bckm.bdate})")
    print(f"  components shape = {bckm.components.shape}")

    # ── Build params + proto matching what diag_bckm_data_isolation.py uses
    # so the calibrated-α/δ/etc. constants are identical. g_share comes from
    # OUR data; this is fine because g_share only affects model.steady_state
    # which gets overwritten by ss_new at BCKM Sbar. ─────────────────────
    df, _ = build_us_dataset(
        start="1980Q1", end="2014Q4",
        data_path="bckm_replication/data/us_1980_2014_calgz.parquet",
        detrend_method="calgz", base_year_quarter="2008Q1",
    )
    g_share = float(df["g"].mean() / df["y"].mean())
    proto = PrototypeModel(
        CalibrationParams(gamma_annual=0.019, n_annual=0.0098, g_share=g_share)
    )

    # BCKM Y_raw cols [y, x, h, g, c, c_implied] → our (y, l, x, g)
    obs_hat_bckm = bckm.Y_raw[:, [0, 2, 1, 3]].copy()
    print(f"  obs_hat_bckm shape = {obs_hat_bckm.shape} (BCKM data, our column order)")

    # ── Run pipeline at BCKM-θ on BCKM Y_raw to get H, ss_new ────────────
    res = estimate_var_mle(
        obs_hat_bckm, proto, verbose=False,
        eval_only=(SBAR_BCKM, P_BCKM, QCHOL_BCKM),
    )
    print(f"  LL (our convention) = {res['log_likelihood']:+.4f} on BCKM Y_raw at BCKM-θ")
    print(f"  ss_new keys: {list(res['ss_new'].keys())}")

    # ── Extract analytical (gwedges2.m-faithful) wedges from BCKM Y_raw ──
    # This produces the realized hat-state path that gwedges2.m feeds into
    # the per-wedge CF construction.
    states_analytic = extract_wedges_bckm_style(
        obs_hat=obs_hat_bckm,
        obs_offset=res["obs_offset_wedge"],
        H=res["H"],
        ss=res["ss_new"],
        params=proto.p,
    )
    print(f"  states_analytic shape = {states_analytic.shape} "
          f"[log_k_hat, log_z_hat, taul_hat, taux_hat, log_g_hat]")

    # Sanity: print the bind-period state (should be small but non-zero —
    # the smoothed wedge at bind isn't exactly SS).
    print(f"  state at bind={bind}: "
          f"log_k_hat={states_analytic[bind, 0]:+.4f}, "
          f"log_z_hat={states_analytic[bind, 1]:+.4f}, "
          f"taul_hat={states_analytic[bind, 2]:+.4f}, "
          f"taux_hat={states_analytic[bind, 3]:+.4f}, "
          f"log_g_hat={states_analytic[bind, 4]:+.4f}")

    # ── Run all single-wedge counterfactuals at BCKM-θ ───────────────────
    cfs = run_all_counterfactuals(
        smoothed_states=states_analytic,
        model=proto,
        P_var=P_BCKM,
        P_0=None,  # P_0 isn't used in run_counterfactual any more
        ss=res["ss_new"],
        Sbar=SBAR_BCKM,
    )

    # ── Convert each CF to bind-anchored level ratio (×100) ──────────────
    # gwedges2.m:170: mz = 100 * exp(YMz(:, var) - YMz(Y0, var))
    print("\n" + "=" * 84)
    print(f"Per-wedge components: ours ({{wedge}}_{{var}}) vs BCKM (m{{w}}{{v}})")
    print("=" * 84)

    # GR window: 2008Q1 (=bind) through 2011Q4 = bind + 15 (BCKM fstats3.m)
    gr_window = (bind, bind + 15)

    # Pretty header
    print(f"\n  Comparison metric: |ours_level − bckm_level| in pp (level=100·ratio).")
    print(f"  GR window = quarters {gr_window[0]} ({bckm.time[gr_window[0]]}) "
          f"to {gr_window[1]} ({bckm.time[gr_window[1]]}) inclusive.\n")

    summary = {}
    for wedge in ["efficiency", "labor", "investment", "government"]:
        print(f"  ── {wedge} wedge ──")
        for var in ["y", "l", "x"]:
            comp_col = _COMP_MAP[(wedge, var)]
            ours_lvl = _hat_to_level_ratio(cfs[wedge][var], bind)
            bckm_lvl = bckm.components[comp_col].values
            label = f"{var} vs {comp_col}"
            summary[(wedge, var)] = _compare(label, ours_lvl, bckm_lvl, bind, gr_window)
        print()

    # ── GR-window side-by-side at key dates ──────────────────────────────
    print("=" * 84)
    print("GR-window side-by-side: efficiency-wedge contributions (ours / bckm)")
    print("=" * 84)
    print(f"  {'date':14s}  {'mzy_ours':>10s} {'mzy_bckm':>10s}    "
          f"{'mzh_ours':>10s} {'mzh_bckm':>10s}    "
          f"{'mzx_ours':>10s} {'mzx_bckm':>10s}")

    cf_eff_y = _hat_to_level_ratio(cfs["efficiency"]["y"], bind)
    cf_eff_l = _hat_to_level_ratio(cfs["efficiency"]["l"], bind)
    cf_eff_x = _hat_to_level_ratio(cfs["efficiency"]["x"], bind)
    mzy = bckm.components["mzy"].values
    mzh = bckm.components["mzh"].values
    mzx = bckm.components["mzx"].values

    sample_dates = [("2007Q4", bind - 1), ("2008Q1 bind", bind),
                    ("2008Q4", bind + 3), ("2009Q3 trough", bind + 6),
                    ("2010Q4", bind + 11), ("2011Q4 (GR end)", bind + 15)]
    for label, t in sample_dates:
        print(f"  {label:14s}  "
              f"{cf_eff_y[t]:>10.3f} {mzy[t]:>10.3f}    "
              f"{cf_eff_l[t]:>10.3f} {mzh[t]:>10.3f}    "
              f"{cf_eff_x[t]:>10.3f} {mzx[t]:>10.3f}")

    # ── Verdict ──────────────────────────────────────────────────────────
    print("\n" + "=" * 84)
    print("[Verdict]  Counterfactual machinery isolation")
    print("=" * 84)

    gr_maxes = [s["gr_max"] for s in summary.values()]
    full_maxes = [s["full_max"] for s in summary.values()]
    overall_gr = max(gr_maxes)
    overall_full = max(full_maxes)

    print(f"  Worst GR-window max|diff|   = {overall_gr:.3f}pp")
    print(f"  Worst full-sample max|diff| = {overall_full:.3f}pp")
    print()

    if overall_gr < 0.5:
        print("  → MATCH on GR window. Counterfactual machinery is clean at BCKM data + BCKM-θ.")
        print("    The model-side residual at identical (data, θ) is contained in the")
        print("    16.6-nat L gap from diag_bckm_data_isolation.py. Fix path is DATA.")
    elif overall_gr < 2.0:
        print("  → CLOSE on GR window. Sub-2pp residual — may be the same 16.6-nat residual")
        print("    propagating, or a small bug in the C0 subtraction. Worth a closer look,")
        print("    but the dominant signal is still data-side.")
    else:
        print("  → DIFFER on GR window. Counterfactual machinery has a bug INDEPENDENT of")
        print("    data construction. Investigate before refactoring data. Likely sites:")
        print("      - run_counterfactual: realized-wedge feed (zeroing inactive cols?)")
        print("      - solve_counterfactual: (C_j − C_0) subtraction algebra")
        print("      - bckm_state_space_cf: As-pinning or Gamma construction")


if __name__ == "__main__":
    main()
