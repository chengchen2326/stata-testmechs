// ============================================================================
// testmechs_dqrdc2_rank.c
//
// Stata plugin that computes the rank of a Stata matrix using R's dqrdc2
// routine (LINPACK QR with limited pivoting strategy).
//
// Stata invocation:
//   plugin call testmechs_dqrdc2_rank, "MatName" "tol_str" "result_scalar"
//
// Arguments:
//   argv[0] = name of input Stata matrix (must be square or rectangular)
//   argv[1] = tolerance as a string (e.g. "1e-5" or "0.00001")
//   argv[2] = name of Stata scalar to receive the rank
//
// On success: writes rank to the scalar named in argv[2]; returns 0.
// ============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stplugin.h"

// Fortran subroutine signature (compiled from R's dqrdc2.f).
// gfortran adds a trailing underscore to the symbol name.
extern void dqrdc2_(double* x, int* ldx, int* n, int* p, double* tol,
                    int* k, double* qraux, int* jpvt, double* work);

STDLL stata_call(int argc, char *argv[])
{
    ST_retcode rc;
    char msg[256];

    // ---- 1. Validate arguments ----
    if (argc != 3) {
        SF_error("testmechs_dqrdc2_rank: expected 3 args (matrix tol scalar)\n");
        return 198;
    }

    const char *mat_name    = argv[0];
    const char *tol_str     = argv[1];
    const char *result_name = argv[2];

    // ---- 2. Parse tolerance ----
    double tol = atof(tol_str);
    if (tol <= 0.0) {
        snprintf(msg, sizeof(msg),
                 "testmechs_dqrdc2_rank: invalid tol '%s'\n", tol_str);
        SF_error(msg);
        return 198;
    }

    // ---- 3. Get matrix shape ----
    // SF_row / SF_col return the number of rows/cols of the named matrix.
    int n = (int)SF_row((char*)mat_name);
    int p = (int)SF_col((char*)mat_name);

    if (n <= 0 || p <= 0) {
        snprintf(msg, sizeof(msg),
                 "testmechs_dqrdc2_rank: matrix '%s' has invalid shape %d x %d\n",
                 mat_name, n, p);
        SF_error(msg);
        return 198;
    }

    // ---- 4. Read matrix into column-major C array ----
    // dqrdc2 expects x stored column-major with leading dimension ldx >= n.
    // x[i + j*ldx] is the (i, j) entry, both 0-indexed in C.
    int ldx = n;
    double *x = (double*)malloc((size_t)ldx * (size_t)p * sizeof(double));
    if (!x) {
        SF_error("testmechs_dqrdc2_rank: out of memory (x)\n");
        return 909;
    }

    for (int j = 0; j < p; j++) {
        for (int i = 0; i < n; i++) {
            double val;
            // SF_mat_el is 1-indexed
            rc = SF_mat_el((char*)mat_name, i + 1, j + 1, &val);
            if (rc) {
                free(x);
                snprintf(msg, sizeof(msg),
                         "testmechs_dqrdc2_rank: SF_mat_el failed at (%d,%d)\n",
                         i + 1, j + 1);
                SF_error(msg);
                return rc;
            }
            x[i + (size_t)j * ldx] = val;
        }
    }

    // ---- 5. Allocate work arrays for dqrdc2 ----
    int    k = 0;
    double *qraux = (double*)calloc((size_t)p,     sizeof(double));
    int    *jpvt  = (int*)   calloc((size_t)p,     sizeof(int));
    double *work  = (double*)calloc((size_t)p * 2, sizeof(double));

    if (!qraux || !jpvt || !work) {
        free(x); free(qraux); free(jpvt); free(work);
        SF_error("testmechs_dqrdc2_rank: out of memory (work arrays)\n");
        return 909;
    }

    // jpvt(j) on entry: column index (1-based for Fortran)
    for (int j = 0; j < p; j++) jpvt[j] = j + 1;

    // ---- 6. Call dqrdc2 ----
    dqrdc2_(x, &ldx, &n, &p, &tol, &k, qraux, jpvt, work);

    // ---- 7. Write rank back to Stata scalar ----
    rc = SF_scal_save((char*)result_name, (double)k);

    // ---- 8. Free everything ----
    free(x);
    free(qraux);
    free(jpvt);
    free(work);

    if (rc) {
        snprintf(msg, sizeof(msg),
                 "testmechs_dqrdc2_rank: SF_scal_save('%s') failed\n",
                 result_name);
        SF_error(msg);
        return rc;
    }

    return 0;
}
