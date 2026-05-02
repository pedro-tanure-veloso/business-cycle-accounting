"""COVID-era BCA smoke tests.

Three test classes:
  TestCovidDatasets          — shape / metadata checks on both parquets
  TestCovidStructuralIdentities — algebraic identities on 2010-2023 data
  TestCovidSmokeOutputs      — driver script exits cleanly + figures exist

Slow tests (those that load data or run the MLE) are marked @pytest.mark.slow.
All dataset tests are skipped if the parquets haven't been built yet.

To build the parquets and run all tests:
    FRED_API_KEY=... python covid_analysis/scripts/run_covid_analysis.py
    pytest tests/test_covid_analysis.py -v
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
FULL_PARQUET     = REPO_ROOT / "covid_analysis" / "data" / "us_2010_2023_calgz.parquet"
PRECOVID_PARQUET = REPO_ROOT / "covid_analysis" / "data" / "us_2010_2023_calgz_preCOVID.parquet"
DRIVER_SCRIPT    = REPO_ROOT / "covid_analysis" / "scripts" / "run_covid_analysis.py"


# ── fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def full_dataset():
    from bca_core.data.pipeline import build_us_dataset
    df, meta = build_us_dataset(
        start="2010Q1", end="2023Q4",
        detrend_method="calgz",
        base_year_quarter="2019Q4",
        labor_target_mean=None,
        data_path=str(FULL_PARQUET),
    )
    return df, meta


@pytest.fixture(scope="module")
def precovid_dataset():
    from bca_core.data.pipeline import build_us_dataset
    df, meta = build_us_dataset(
        start="2010Q1", end="2023Q4",
        detrend_method="calgz",
        base_year_quarter="2019Q4",
        labor_target_mean=None,
        mle_window=("2010Q1", "2019Q4"),
        data_path=str(PRECOVID_PARQUET),
    )
    return df, meta


@pytest.fixture(scope="module")
def full_pipeline(full_dataset):
    """Run the BCA pipeline on the full-window dataset."""
    import numpy as np
    from bca_core.params import CalibrationParams
    from bca_core.model import PrototypeModel
    from bca_core.var_estimation import estimate_var_mle, prepare_observables
    from bca_core.wedges import extract_wedges_bckm_style
    from bca_core.counterfactuals import run_all_counterfactuals
    from bca_core.constants import (
        P_BCKM_TABLE8 as P_BCKM,
        QCHOL_BCKM_TABLE10 as QCHOL_BCKM,
        SBAR_BCKM_TABLE8 as SBAR_BCKM,
    )

    df, meta = full_dataset
    g_share = float(df["g"].mean() / df["y"].mean())
    params = CalibrationParams(
        gamma_annual=meta.get("gamma_annual", 0.019),
        n_annual=meta.get("n_annual", 0.0098),
        g_share=g_share,
    )
    proto = PrototypeModel(params)
    ss = proto.steady_state()
    obs_hat, _ = prepare_observables(df, ss, center=False)
    data_means = np.array([
        df["y"].mean(),
        (df["x"] / df["y"]).mean(),
        df["l"].mean(),
        (df["g"] / df["y"]).mean(),
    ])
    res = estimate_var_mle(
        obs_hat, proto, verbose=False, data_means=data_means,
        warm_start=(SBAR_BCKM, P_BCKM, QCHOL_BCKM),
    )
    states = extract_wedges_bckm_style(
        obs_hat=obs_hat, obs_offset=res["obs_offset"],
        H=res["H"], ss=res["ss_new"], params=params,
    )
    obs_dev = obs_hat - res["obs_offset"]
    data_hat = {"y": obs_dev[:, 0], "l": obs_dev[:, 1], "x": obs_dev[:, 2]}
    P0_implied = (np.eye(4) - res["P"]) @ res["Sbar"]
    cfs = run_all_counterfactuals(
        states, proto, res["P"], P_0=P0_implied,
        ss=res["ss_new"], Sbar=res["Sbar"],
    )
    return df, meta, params, proto, ss, res, states, obs_hat, data_hat, cfs


# ── TestCovidDatasets ─────────────────────────────────────────────────────────

@pytest.mark.slow
@pytest.mark.skipif(not FULL_PARQUET.exists(), reason=f"{FULL_PARQUET} not built yet")
class TestCovidDatasets:

    def test_shape_full(self, full_dataset):
        df, _ = full_dataset
        # 2010Q1 to 2023Q4 = 56 quarters
        assert len(df) == 56, f"Expected 56 rows, got {len(df)}"
        assert df.shape[1] >= 5
        assert not df[["y", "l", "x", "g"]].isnull().any().any()

    @pytest.mark.skipif(not PRECOVID_PARQUET.exists(),
                        reason=f"{PRECOVID_PARQUET} not built yet")
    def test_shape_precovid(self, precovid_dataset):
        df, _ = precovid_dataset
        assert len(df) == 56
        assert not df[["y", "l", "x", "g"]].isnull().any().any()

    def test_bind_is_2019Q4_full(self, full_dataset):
        _, meta = full_dataset
        assert meta.get("base_year_quarter") == "2019Q4", \
            f"Expected base_year_quarter='2019Q4', got {meta.get('base_year_quarter')}"

    @pytest.mark.skipif(not PRECOVID_PARQUET.exists(),
                        reason=f"{PRECOVID_PARQUET} not built yet")
    def test_bind_is_2019Q4_precovid(self, precovid_dataset):
        _, meta = precovid_dataset
        assert meta.get("base_year_quarter") == "2019Q4"

    def test_y_at_bind_is_one_full(self, full_dataset):
        import pandas as pd
        df, meta = full_dataset
        bind_dt = pd.Period("2019Q4", freq="Q").start_time
        try:
            bind_idx = df.index.get_loc(bind_dt)
        except KeyError:
            bind_idx = meta.get("base_idx")
        assert bind_idx is not None
        assert abs(float(df["y"].iloc[bind_idx]) - 1.0) < 0.02, \
            f"y at bind = {float(df['y'].iloc[bind_idx]):.4f}, expected ≈1.0"


# ── TestCovidStructuralIdentities ─────────────────────────────────────────────

@pytest.mark.slow
@pytest.mark.skipif(not FULL_PARQUET.exists(), reason=f"{FULL_PARQUET} not built yet")
class TestCovidStructuralIdentities:
    """Window-agnostic algebraic identities on 2010-2023 data."""

    def test_all_wedges_reproduce_data(self, full_pipeline):
        """All-active CF ≈ data path (incremental identity)."""
        from bca_core.counterfactuals import run_all_counterfactuals
        df, meta, params, proto, ss, res, states, obs_hat, data_hat, cfs = full_pipeline

        # Sum of four incremental CFs should ≈ all-active CF ≈ data
        for var in ["y", "l", "x"]:
            incremental_sum = np.zeros(len(df))
            for wname in ["efficiency", "labor", "investment", "government"]:
                incremental_sum += cfs[wname][var]
            # All-active CF from run_all_counterfactuals is already incremental
            # so sum ≈ data_hat[var]; check correlation > 0.9
            corr = float(np.corrcoef(incremental_sum, data_hat[var])[0, 1])
            assert corr > 0.8, \
                f"CF sum vs data correlation for {var} = {corr:.3f} (expected > 0.8)"

    def test_cf_additivity(self, full_pipeline):
        """Four incremental CFs sum to ≈ all-active CF (additive decomposition)."""
        df, meta, params, proto, ss, res, states, obs_hat, data_hat, cfs = full_pipeline
        from bca_core.counterfactuals import solve_counterfactual, run_counterfactual

        P0 = (np.eye(4) - res["P"]) @ res["Sbar"]
        cf_all = run_counterfactual(
            states,
            solve_counterfactual(
                proto, res["P"], active_wedges=[0, 1, 2, 3],
                ss=res["ss_new"], Sbar=res["Sbar"],
            ),
        )

        for var in ["y", "l", "x"]:
            total = sum(cfs[wname][var] for wname in
                        ["efficiency", "labor", "investment", "government"])
            max_diff = float(np.abs(total - cf_all[var]).max())
            assert max_diff < 0.05, \
                f"CF additivity for {var}: max|sum - all-active| = {max_diff:.4f}"

    def test_fstat_sums_to_one(self, full_pipeline):
        """f-statistics rows sum to 1.0 for each variable."""
        from bca_core.counterfactuals import f_statistics_bckm
        df, meta, params, proto, ss, res, states, obs_hat, data_hat, cfs = full_pipeline

        def find_idx(dates, year, quarter):
            qmap = {1: ("01", "Q1"), 2: ("04", "Q2"),
                    3: ("07", "Q3"), 4: ("10", "Q4")}
            month, qstr = qmap[quarter]
            for i, d in enumerate(dates):
                s = str(d)
                if str(year) in s and (month in s or qstr in s):
                    return i
            return None

        anchor = find_idx(df.index, 2019, 4)
        end    = find_idx(df.index, 2022, 4)
        assert anchor is not None and end is not None

        f_df = f_statistics_bckm(data_hat, cfs, window=(anchor, end), anchor=anchor)

        for var in ["y", "l", "x"]:
            try:
                row_sum = float(f_df.loc[:, var].sum()
                                if var in f_df.columns else f_df.loc[var].sum())
            except Exception:
                row_sum = float(f_df[var].sum())
            assert abs(row_sum - 1.0) < 1e-6, \
                f"f-stat sum for {var} = {row_sum:.8f} (expected 1.0)"


# ── TestCovidSmokeOutputs ─────────────────────────────────────────────────────

class TestCovidSmokeOutputs:
    """Smoke-only: driver exits cleanly and figure files are produced."""

    @pytest.mark.slow
    def test_driver_exits_cleanly(self):
        """run_covid_analysis.py exits with code 0."""
        result = subprocess.run(
            [sys.executable, str(DRIVER_SCRIPT)],
            capture_output=True,
            text=True,
            cwd=str(REPO_ROOT),
        )
        if result.returncode != 0:
            pytest.fail(
                f"Driver script exited {result.returncode}.\n"
                f"STDOUT:\n{result.stdout[-2000:]}\n"
                f"STDERR:\n{result.stderr[-2000:]}"
            )

    @pytest.mark.slow
    @pytest.mark.skipif(not FULL_PARQUET.exists(), reason="run driver first")
    def test_figures_exist(self):
        fig_dir = REPO_ROOT / "covid_analysis" / "figures"
        expected = [
            "wedges_us_2010_2023.png",
            "figure_2C_covid.png",
            "figure_2D_covid.png",
            "figure_2E_covid.png",
        ]
        missing = [f for f in expected if not (fig_dir / f).exists()]
        assert not missing, f"Missing figures: {missing}"
