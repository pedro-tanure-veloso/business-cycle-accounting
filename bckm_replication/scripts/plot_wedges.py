"""
Plot extracted wedges for the US — compare against BCKM (2016) Figure 2.

Usage:
    FRED_API_KEY=... python scripts/plot_wedges.py
"""

import matplotlib.pyplot as plt
import numpy as np

from bca_core.data.pipeline import build_us_dataset
from bca_core.params import CalibrationParams
from bca_core.wedges import extract_static_wedges


def main():
    # Fetch data and extract wedges
    df, meta = build_us_dataset(start="1980Q1", end="2014Q4")

    params = CalibrationParams(
        gamma_annual=meta["gamma_annual"],
        n_annual=meta["n_annual"],
    )

    wedges = extract_static_wedges(df, params)

    # Normalize: index to 100 at first observation (BCKM convention)
    dates = df.index

    A_norm = 100 * wedges["A"] / wedges["A"].iloc[0]
    taul_norm = 100 * wedges["one_minus_tau_l"] / wedges["one_minus_tau_l"].iloc[0]
    g_norm = 100 * wedges["g"] / wedges["g"].iloc[0]

    # Also plot data series for reference
    y_norm = 100 * df["y"] / df["y"].iloc[0]
    l_norm = 100 * df["l"] / df["l"].iloc[0]
    x_norm = 100 * df["x"] / df["x"].iloc[0]

    # NBER recession dates (approximate quarters)
    recessions = [
        ("1980-01-01", "1980-07-01"),
        ("1981-07-01", "1982-10-01"),
        ("1990-07-01", "1991-01-01"),
        ("2001-01-01", "2001-10-01"),
        ("2007-10-01", "2009-04-01"),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    fig.suptitle(
        "US Business Cycle Accounting — Wedges & Data (1980Q1–2014Q4)",
        fontsize=14,
        fontweight="bold",
    )

    def add_recessions(ax):
        for start, end in recessions:
            ax.axvspan(start, end, alpha=0.15, color="gray")

    # Panel 1: Output
    ax = axes[0, 0]
    ax.plot(dates, y_norm, "b-", linewidth=1.5)
    ax.set_title("Output (detrended)")
    ax.set_ylabel("Index (1980Q1 = 100)")
    add_recessions(ax)
    ax.grid(True, alpha=0.3)

    # Panel 2: Labor
    ax = axes[0, 1]
    ax.plot(dates, l_norm, "b-", linewidth=1.5)
    ax.set_title("Hours worked")
    ax.set_ylabel("Index (1980Q1 = 100)")
    add_recessions(ax)
    ax.grid(True, alpha=0.3)

    # Panel 3: Investment
    ax = axes[0, 2]
    ax.plot(dates, x_norm, "b-", linewidth=1.5)
    ax.set_title("Investment (detrended)")
    ax.set_ylabel("Index (1980Q1 = 100)")
    add_recessions(ax)
    ax.grid(True, alpha=0.3)

    # Panel 4: Efficiency wedge
    ax = axes[1, 0]
    ax.plot(dates, A_norm, "r-", linewidth=1.5)
    ax.set_title("Efficiency wedge (A)")
    ax.set_ylabel("Index (1980Q1 = 100)")
    add_recessions(ax)
    ax.grid(True, alpha=0.3)

    # Panel 5: Labor wedge
    ax = axes[1, 1]
    ax.plot(dates, taul_norm, "r-", linewidth=1.5)
    ax.set_title("Labor wedge (1 − τₗ)")
    ax.set_ylabel("Index (1980Q1 = 100)")
    add_recessions(ax)
    ax.grid(True, alpha=0.3)

    # Panel 6: Government wedge
    ax = axes[1, 2]
    ax.plot(dates, g_norm, "r-", linewidth=1.5)
    ax.set_title("Government wedge (g)")
    ax.set_ylabel("Index (1980Q1 = 100)")
    add_recessions(ax)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("wedges_us_1980_2014.png", dpi=150, bbox_inches="tight")
    print("Saved: wedges_us_1980_2014.png")
    plt.show()


if __name__ == "__main__":
    main()
