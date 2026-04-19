from sfi import Matrix, Scalar
import numpy as np
from scipy.sparse import csc_matrix
import osqp

# Read QP problem:
#   min 0.5 x' P x + q' x
#   s.t. l <= A x <= u
# Inputs are Stata matrices __tm_qp_P, __tm_qp_q, __tm_qp_A, __tm_qp_l, __tm_qp_u.
# In our code P = Dmat and we pass Dmat directly (no factor-of-2 flip).

P = np.asarray(Matrix.get("__tm_qp_P"), dtype=float)
q = np.asarray(Matrix.get("__tm_qp_q"), dtype=float).flatten()
A = np.asarray(Matrix.get("__tm_qp_A"), dtype=float)
l = np.asarray(Matrix.get("__tm_qp_l"), dtype=float).flatten()
u = np.asarray(Matrix.get("__tm_qp_u"), dtype=float).flatten()

# OSQP requires infinity replaced with large finite values (or OSQP_INFTY);
# OSQP 0.6 accepts np.inf directly as well.
# Replace our Mata-side "-1e20" sentinel with -np.inf for clarity.
l[l <= -1e19] = -np.inf
u[u >=  1e19] =  np.inf

# OSQP requires P and A in sparse (CSC) format.
P_sp = csc_matrix(P)
A_sp = csc_matrix(A)

prob = osqp.OSQP()
# Match R TestMechs settings: eps_abs = eps_rel = 1e-8, verbose off.
prob.setup(P=P_sp, q=q, A=A_sp, l=l, u=u,
           eps_abs=1e-8, eps_rel=1e-8, verbose=False, polish=True,
           max_iter=20000)

res = prob.solve()

# status_val: OSQP uses string statuses; treat "solved" and "solved inaccurate" as ok.
ok = 1.0 if res.info.status in ("solved", "solved inaccurate") else 0.0

# Write solution x back. Pad with NaN if unsuccessful.
nvar = P.shape[0]
if ok == 1.0 and res.x is not None:
    x_out = np.asarray(res.x, dtype=float).reshape(-1, 1).tolist()
else:
    x_out = np.full((nvar, 1), np.nan).tolist()

# Create (and overwrite) the matrix with the correct dimensions
Matrix.create("__tm_qp_x", nvar, 1, 0)
Matrix.store("__tm_qp_x", x_out)

Scalar.setValue("__tm_qp_ok", ok)
