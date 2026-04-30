# ============================================================================
# _testmechs_lp_call.py
#
# Production LP solver for testmechs. Uses GLPK via swiglpk.
#
# Why GLPK and not HiGHS:
#   R's TestMechs package uses GLPK (via Rglpk_solve_LP). Using the same
#   solver in Stata makes the Stata implementation faithful to the R one and
#   eliminates LP-solver-induced numerical differences.
#
# Inputs  (Stata matrices):  __tm_c, __tm_Aeq, __tm_beq
# Outputs (Stata scalars):   __tm_lp_fun, __tm_lp_ok
#
# Requires the Python package swiglpk (pip install swiglpk) and the GLPK
# system library (e.g. brew install glpk on macOS, apt install libglpk-dev on
# Debian/Ubuntu).
# ============================================================================

from sfi import Matrix, Scalar
import numpy as np
import swiglpk as glpk

# ---------- Read LP inputs from Stata matrices -----------------------------

_c   = np.asarray(Matrix.get("__tm_c"),   dtype=float).flatten()
_Aeq = np.asarray(Matrix.get("__tm_Aeq"), dtype=float)
_beq = np.asarray(Matrix.get("__tm_beq"), dtype=float).flatten()
if _Aeq.ndim == 1:
    _Aeq = _Aeq.reshape(1, -1)

n_var  = int(_c.size)
n_cons = int(_beq.size)

# ---------- Build and solve the LP via GLPK --------------------------------

lp = glpk.glp_create_prob()
glpk.glp_set_obj_dir(lp, glpk.GLP_MIN)

# Variables: x_j >= 0, no upper bound, objective coefficient _c[j-1]
glpk.glp_add_cols(lp, n_var)
for j in range(1, n_var + 1):
    glpk.glp_set_col_bnds(lp, j, glpk.GLP_LO, 0.0, 0.0)
    glpk.glp_set_obj_coef(lp, j, float(_c[j - 1]))

# Equality constraints: A_eq x = b_eq
glpk.glp_add_rows(lp, n_cons)
for i in range(1, n_cons + 1):
    glpk.glp_set_row_bnds(lp, i, glpk.GLP_FX, float(_beq[i - 1]), float(_beq[i - 1]))

# Load the constraint matrix (1-indexed triplet form)
nnz = n_cons * n_var
ia = glpk.intArray(nnz + 1)
ja = glpk.intArray(nnz + 1)
ar = glpk.doubleArray(nnz + 1)
k = 1
for i in range(1, n_cons + 1):
    for j in range(1, n_var + 1):
        ia[k] = i
        ja[k] = j
        ar[k] = float(_Aeq[i - 1, j - 1])
        k += 1
glpk.glp_load_matrix(lp, nnz, ia, ja, ar)

# Solve using simplex (matches Rglpk default)
parm = glpk.glp_smcp()
glpk.glp_init_smcp(parm)
parm.msg_lev = glpk.GLP_MSG_OFF
ret = glpk.glp_simplex(lp, parm)

status = glpk.glp_get_status(lp)
ok  = 1.0 if status == glpk.GLP_OPT else 0.0
fun = float(glpk.glp_get_obj_val(lp)) if ok == 1.0 else float("nan")

glpk.glp_delete_prob(lp)

# ---------- Write results back to Stata ------------------------------------

Scalar.setValue("__tm_lp_fun", fun)
Scalar.setValue("__tm_lp_ok",  ok)
