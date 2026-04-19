# Solve LP min c'x s.t. A_eq x = b_eq, x >= 0 via scipy HiGHS.
# Inputs  (Stata matrices):  __tm_c, __tm_Aeq, __tm_beq
# Outputs (Stata scalars):   __tm_lp_fun, __tm_lp_ok

from sfi import Matrix, Scalar
import numpy as np
from scipy.optimize import linprog

_c   = np.asarray(Matrix.get("__tm_c"),   dtype=float).flatten()
_Aeq = np.asarray(Matrix.get("__tm_Aeq"), dtype=float)
_beq = np.asarray(Matrix.get("__tm_beq"), dtype=float).flatten()
if _Aeq.ndim == 1:
    _Aeq = _Aeq.reshape(1, -1)

_bounds = [(0, None)] * _c.size
_res = linprog(c=_c, A_eq=_Aeq, b_eq=_beq, bounds=_bounds, method="highs")

Scalar.setValue("__tm_lp_fun", float(_res.fun) if _res.fun is not None else float("nan"))
Scalar.setValue("__tm_lp_ok",  1.0 if _res.success else 0.0)
