"""Canonical BCKM constants in CODE convention.

Convention note (read this before changing anything in this file):
================================================================

BCKM publishes the wedge VAR transition matrix in `BCA_info.md` Section 7
Table 8 in a "rows = drivers, columns = receivers" convention — i.e.,
``Table8[i, j]`` is the coefficient of state ``i_{t-1}`` in the equation
for state ``j_t``. Reading row 0 of the table tells you "what z does":
its self-persistence in column 0, its outgoing spillover to taul in
column 1, etc. This is convenient for narrative interpretation but is
the **transpose** of the standard textbook VAR.

BCKM's matlab code (``mleqadj.m:222``) uses the textbook convention
``state_{t+1} = P · state_t``, so ``P[i, j]`` is the coefficient of
state ``j_t`` in the equation for state ``i_{t+1}``. The matrix stored
in ``worktemp.params`` is in this convention. Reading row 0: "what
determines z's update".

The two conventions are transposes of each other. From 2026-04 through
2026-04-30 the codebase had Table 8 hardcoded in nine independent places
in the **paper convention** but used everywhere as if it were in the
**code convention**, producing a transposed P at every BCKM-θ
evaluation. This corrupted the warm-start, every probe, and the
counterfactual decomposition. Quantified impact at the published θ:
LL = +1195 (wrong) vs +1697 (correct), a 501-nat error.

**This module exports P in the CODE convention** — the same matrix
``mleqadj.m`` stores in ``worktemp.params``. If you need the table as
printed in the paper, take the transpose. Do **not** re-transcribe
Table 8 anywhere else; import from here.

Verified element-wise against ``octave_output/P_bckm.csv`` (dumped from
``mleqadj.m`` on 2026-04-30) to within 4.3e-5 (Table 8's 4-decimal
publication precision).
"""

from __future__ import annotations

import numpy as np

# Sbar in BCKM's CODE convention: [log(z), tauls, tauxs, log(g)] at the
# US converged MLE. Source: ``BCA_info.md`` §7 Table 8 (4-decimal
# precision); octave dump in ``octave_output/Sbar_bckm.csv`` carries the
# 16-decimal version (differences are <1e-4 and add ~0.01 nats to LL).
SBAR_BCKM_TABLE8 = np.array([0.1336, 0.3691, -0.0460, -1.9355])

# P in CODE convention (row i = equation for state i, col j = lag of
# state j). State order ``[z, taul, taux, g]``. This is the TRANSPOSE
# of the matrix as it appears in ``BCA_info.md`` §7 Table 8 — see the
# module docstring above.
P_BCKM_TABLE8 = np.array(
    [
        [0.9887, -0.0012, -0.0045,  0.0063],
        [0.0307,  1.0011,  0.0449,  0.0017],
        [-0.0089, -0.0275, 0.9675,  0.0016],
        [-0.0407, 0.0175, -0.0426,  0.9945],
    ]
)

# Q (lower-triangular Cholesky factor of V = QQ′) in BCKM Table 10. This
# matrix has no row/col convention switch — the off-diagonals are zero
# above the diagonal in both conventions. Source: ``BCA_info.md`` §7
# Table 10. Verified element-wise against ``octave_output/Qchol_bckm.csv``.
QCHOL_BCKM_TABLE10 = np.array(
    [
        [0.0077, 0.0000, 0.0000, 0.0000],
        [0.0024, 0.0043, 0.0000, 0.0000],
        [-0.0041, 0.0023, 0.0088, 0.0000],
        [0.0003, 0.0153, 0.0121, 0.0139],
    ]
)
