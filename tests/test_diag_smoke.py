"""Smoke test for ``bckm_replication/scripts/diag_tx_counterfactual.py``.

Skipped when ``bckm_replication/data/us_1980_2014_calgz.parquet`` is
missing — the script needs the full US dataset to evaluate at BCKM-θ.
When the data is present, running ``main()`` is the high-level
regression check that future refactors of ``solve_counterfactual`` /
``run_counterfactual`` do not reintroduce the 2026-04-29 cf-fix bugs
(the script asserts ``max diff`` between all-active CF policies and H
is below 1e-8 in its own output).
"""
from __future__ import annotations

import importlib.util
import io
import os
import sys
from contextlib import redirect_stdout
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_PATH = REPO_ROOT / "bckm_replication" / "data" / "us_1980_2014_calgz.parquet"
SCRIPT_PATH = (
    REPO_ROOT / "bckm_replication" / "scripts" / "diag_tx_counterfactual.py"
)


@pytest.mark.bckm
@pytest.mark.skipif(
    not DATA_PATH.exists(),
    reason="US dataset not present (run scripts to build it first)",
)
def test_diag_tx_counterfactual_runs_without_error():
    spec = importlib.util.spec_from_file_location(
        "diag_tx_counterfactual", SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(REPO_ROOT))
    try:
        spec.loader.exec_module(module)
        main = module.main
    finally:
        sys.path.pop(0)

    cwd = os.getcwd()
    os.chdir(REPO_ROOT)
    try:
        buf = io.StringIO()
        with redirect_stdout(buf):
            main()
        out = buf.getvalue()
    finally:
        os.chdir(cwd)

    # The script self-reports the all-active CF–H invariant. If the cf-fix
    # is still in place, every reported diff must be tiny.
    for line in out.splitlines():
        if line.lstrip().startswith("diff P_") and "vs H[" in line:
            # Format: "  diff P_y vs H[0]: 1.23e-15"
            val_str = line.rsplit(":", 1)[-1].strip()
            val = float(val_str)
            assert val < 1e-8, f"cf-H invariant broken on line: {line.strip()}"
