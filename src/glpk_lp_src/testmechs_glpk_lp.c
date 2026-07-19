// ============================================================================
// testmechs_glpk_lp.c
//
// Stata plugin that solves a linear program using GLPK (statically linked).
// Replaces the Python + swiglpk path in testmechs_test_sharpnull.
//
// Stata invocation:
//   plugin call _testmechs_glpk_lp, "__tm_c" "__tm_Aeq" "__tm_beq" "__tm_lp_fun" "__tm_lp_ok"
//
// Arguments:
//   argv[0] = name of Stata matrix with objective coefficients (n_var x 1)
//   argv[1] = name of Stata matrix with equality constraints (n_cons x n_var)
//   argv[2] = name of Stata matrix with RHS (n_cons x 1)
//   argv[3] = name of Stata scalar to receive the optimal objective value
//   argv[4] = name of Stata scalar to receive the success flag (1.0 or 0.0)
//
// The LP solved is:
//   min c' x
//   s.t. A_eq x = b_eq
//        x >= 0    (no upper bounds)
//
// On success: writes optimal objective to __tm_lp_fun, and 1.0 to __tm_lp_ok.
// On failure: writes NaN to __tm_lp_fun, and 0.0 to __tm_lp_ok.
// ============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "stplugin.h"
#include "glpk.h"

STDLL stata_call(int argc, char *argv[])
{
    ST_retcode rc;
    char msg[256];

    // ---- 1. Validate arguments ----
    if (argc != 5) {
        SF_error("testmechs_glpk_lp: expected 5 args (c Aeq beq fun_scalar ok_scalar)\n");
        return 198;
    }

    const char *c_name    = argv[0];
    const char *Aeq_name  = argv[1];
    const char *beq_name  = argv[2];
    const char *fun_name  = argv[3];
    const char *ok_name   = argv[4];

    // ---- 2. Get matrix shapes ----
    int n_var  = (int)SF_row((char*)c_name);   // c is n_var x 1 (or 1 x n_var)
    int c_cols = (int)SF_col((char*)c_name);

    // c should be n_var x 1 OR 1 x n_var — accept either
    if (c_cols != 1 && n_var == 1) {
        int tmp = n_var;
        n_var = c_cols;
        c_cols = tmp;
    }
    if (c_cols != 1) {
        snprintf(msg, sizeof(msg),
                 "testmechs_glpk_lp: c matrix must be n_var x 1 or 1 x n_var (got %d x %d)\n",
                 (int)SF_row((char*)c_name), (int)SF_col((char*)c_name));
        SF_error(msg);
        return 198;
    }

    int n_cons     = (int)SF_row((char*)Aeq_name);
    int Aeq_ncols  = (int)SF_col((char*)Aeq_name);
    if (Aeq_ncols != n_var) {
        snprintf(msg, sizeof(msg),
                 "testmechs_glpk_lp: Aeq has %d cols but c has %d entries\n",
                 Aeq_ncols, n_var);
        SF_error(msg);
        return 198;
    }

    int beq_nrows = (int)SF_row((char*)beq_name);
    int beq_ncols = (int)SF_col((char*)beq_name);
    if (beq_ncols != 1 && beq_nrows == 1) {
        int tmp = beq_nrows;
        beq_nrows = beq_ncols;
        beq_ncols = tmp;
    }
    if (beq_ncols != 1 || beq_nrows != n_cons) {
        snprintf(msg, sizeof(msg),
                 "testmechs_glpk_lp: beq must be n_cons x 1 (got %d x %d, expected %d x 1)\n",
                 (int)SF_row((char*)beq_name), (int)SF_col((char*)beq_name), n_cons);
        SF_error(msg);
        return 198;
    }

    if (n_var <= 0 || n_cons <= 0) {
        snprintf(msg, sizeof(msg),
                 "testmechs_glpk_lp: invalid LP shape (n_var=%d, n_cons=%d)\n",
                 n_var, n_cons);
        SF_error(msg);
        return 198;
    }

    // ---- 3. Read c into a C array ----
    // SF_mat_el uses 1-indexed row and column
    double *c = (double*)malloc(n_var * sizeof(double));
    if (!c) {
        SF_error("testmechs_glpk_lp: out of memory (c)\n");
        return 909;
    }
    for (int j = 0; j < n_var; j++) {
        double val;
        // Handle both (n_var x 1) and (1 x n_var) — we normalized n_var above
        // but SF_mat_el needs the actual original shape.
        int orig_rows = (int)SF_row((char*)c_name);
        int orig_cols = (int)SF_col((char*)c_name);
        if (orig_cols == 1) {
            rc = SF_mat_el((char*)c_name, j + 1, 1, &val);
        } else {
            rc = SF_mat_el((char*)c_name, 1, j + 1, &val);
        }
        if (rc) {
            free(c);
            snprintf(msg, sizeof(msg),
                     "testmechs_glpk_lp: SF_mat_el failed on c at index %d\n", j + 1);
            SF_error(msg);
            return rc;
        }
        c[j] = val;
    }

    // ---- 4. Read Aeq into a C array (row-major) ----
    double *Aeq = (double*)malloc((size_t)n_cons * (size_t)n_var * sizeof(double));
    if (!Aeq) {
        free(c);
        SF_error("testmechs_glpk_lp: out of memory (Aeq)\n");
        return 909;
    }
    for (int i = 0; i < n_cons; i++) {
        for (int j = 0; j < n_var; j++) {
            double val;
            rc = SF_mat_el((char*)Aeq_name, i + 1, j + 1, &val);
            if (rc) {
                free(c); free(Aeq);
                snprintf(msg, sizeof(msg),
                         "testmechs_glpk_lp: SF_mat_el failed on Aeq at (%d,%d)\n",
                         i + 1, j + 1);
                SF_error(msg);
                return rc;
            }
            Aeq[i * n_var + j] = val;
        }
    }

    // ---- 5. Read beq ----
    double *beq = (double*)malloc(n_cons * sizeof(double));
    if (!beq) {
        free(c); free(Aeq);
        SF_error("testmechs_glpk_lp: out of memory (beq)\n");
        return 909;
    }
    for (int i = 0; i < n_cons; i++) {
        double val;
        int orig_rows = (int)SF_row((char*)beq_name);
        int orig_cols = (int)SF_col((char*)beq_name);
        if (orig_cols == 1) {
            rc = SF_mat_el((char*)beq_name, i + 1, 1, &val);
        } else {
            rc = SF_mat_el((char*)beq_name, 1, i + 1, &val);
        }
        if (rc) {
            free(c); free(Aeq); free(beq);
            snprintf(msg, sizeof(msg),
                     "testmechs_glpk_lp: SF_mat_el failed on beq at index %d\n", i + 1);
            SF_error(msg);
            return rc;
        }
        beq[i] = val;
    }

    // ---- 6. Build GLPK problem ----
    glp_prob* lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MIN);

    // Variables: x_j >= 0, no upper bound, obj coef c[j]
    glp_add_cols(lp, n_var);
    for (int j = 1; j <= n_var; j++) {
        glp_set_col_bnds(lp, j, GLP_LO, 0.0, 0.0);
        glp_set_obj_coef(lp, j, c[j - 1]);
    }

    // Equality constraints: A_eq x = b_eq
    glp_add_rows(lp, n_cons);
    for (int i = 1; i <= n_cons; i++) {
        glp_set_row_bnds(lp, i, GLP_FX, beq[i - 1], beq[i - 1]);
    }

    // Load constraint matrix (1-indexed triplet)
    int nnz = n_cons * n_var;
    int *ia = (int*)malloc((nnz + 1) * sizeof(int));
    int *ja = (int*)malloc((nnz + 1) * sizeof(int));
    double *ar = (double*)malloc((nnz + 1) * sizeof(double));
    if (!ia || !ja || !ar) {
        free(c); free(Aeq); free(beq);
        free(ia); free(ja); free(ar);
        glp_delete_prob(lp);
        SF_error("testmechs_glpk_lp: out of memory (GLPK triplet)\n");
        return 909;
    }
    int k = 1;
    for (int i = 1; i <= n_cons; i++) {
        for (int j = 1; j <= n_var; j++) {
            ia[k] = i;
            ja[k] = j;
            ar[k] = Aeq[(i - 1) * n_var + (j - 1)];
            k++;
        }
    }
    glp_load_matrix(lp, nnz, ia, ja, ar);

    // ---- 7. Solve with simplex (matches R Rglpk and swiglpk default) ----
    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_OFF;
    int ret = glp_simplex(lp, &parm);
    int status = glp_get_status(lp);

    double ok  = (status == GLP_OPT) ? 1.0 : 0.0;
    double fun = (ok == 1.0) ? glp_get_obj_val(lp) : nan("");

    // ---- 8. Write results back to Stata ----
    rc = SF_scal_save((char*)fun_name, fun);
    if (rc == 0) {
        rc = SF_scal_save((char*)ok_name, ok);
    }

    // ---- 9. Free everything ----
    free(c); free(Aeq); free(beq);
    free(ia); free(ja); free(ar);
    glp_delete_prob(lp);

    if (rc) {
        snprintf(msg, sizeof(msg),
                 "testmechs_glpk_lp: failed to write result scalar\n");
        SF_error(msg);
        return rc;
    }

    return 0;
}
