# Known Limitations

This document describes known limitations and numerical differences between the Stata `testmechs` package and the R `TestMechs` package it is based on. It is intended for researchers who care about exact numerical reproducibility across implementations.

## Summary

The Stata port is a faithful translation of the R reference implementation. All algorithmic logic — partial-order construction for multivariate mediators, the LP-based degree-of-freedom computation, the quadratic programming step for the test statistic, and the rank-based computations on the constraint matrices — mirrors the R code.

However, because the Stata and R ecosystems rely on different underlying numerical libraries for linear programming, quadratic programming, and matrix rank computations, the Stata outputs will not always be bit-exact identical to the R outputs. The discrepancies we have observed are small in all single-mediator cases we have tested, and grow in multi-mediator cases when the underlying optimization problem is near-degenerate.

All tests that we have run yield the same qualitative statistical conclusion (reject or fail to reject the sharp null at standard significance levels) in both the Stata and R implementations.

## Numerical agreement on our validation cases

Using the `mother_data.dta` file (the Baranov et al. experimental sample filtered to `THP_sample == 1`), we compared the `testmechs_test_sharpnull` command to its R equivalent under the default settings (`method = "CS"`, `num_Ybins = 5`, clustered on `uc`).

| Mediator(s)                              | Stata pval   | R pval       | Agreement |
|------------------------------------------|--------------|--------------|-----------|
| `grandmother` only                       | 0.02283916   | 0.02283916   | Identical to 7 decimal places |
| `relationship_husb` only                 | 0.02838328   | 0.02838332   | Agrees to 7 decimal places (~4e-8 apart) |
| `grandmother` and `relationship_husb`    | 0.17107925   | 0.65408650   | Same conclusion at α=0.05 (fail to reject); numerical values differ |

In all three cases, the test statistic itself is identical between Stata and R to at least 6 decimal places. The differences shown above arise only in the degree-of-freedom computation for the default `new_dof_CS = FALSE` algorithm.

## Source of the multi-mediator discrepancy

The default degree-of-freedom algorithm for the Cox-Shi (CS) test requires solving a sequence of linear programs to identify which inequality constraints are binding at the optimum ("implicit equalities"). The number of such binding constraints directly determines the degrees of freedom of the chi-squared reference distribution and hence the p-value.

The R package uses GLPK (via the `Rglpk` R package) to solve these LPs. The Stata package uses HiGHS (via `scipy.optimize.linprog`). When the LP is near-degenerate (multiple vertices achieve the same optimal value), different solvers select different vertices, producing slightly different values for the variables we use to identify binding constraints. In single-mediator problems this difference is typically below the threshold used to declare a constraint binding and has no effect. In multi-mediator problems with moderately large support, it can change the count of binding constraints by several, changing the reported p-value meaningfully even though the underlying test statistic is unchanged.

We experimented with routing the LPs through GLPK directly (via Python's `swiglpk` bindings) to match R exactly. This improved the multi-mediator agreement but degraded the single-mediator cases, because the QR-based rank computations downstream of the LP are themselves sensitive to the LP solver's choices. We therefore retain HiGHS as the default, which matches R exactly in the simpler cases and disagrees numerically but not qualitatively in the more complex cases.

## The `newdofcs` option

The R package exposes an alternative, simpler degree-of-freedom algorithm via the argument `new_dof_CS = TRUE`. Our Stata package exposes the same alternative via the `newdofcs` option. This algorithm does not require LP solves and is entirely numerical linear algebra; it agrees with R in all cases we have tested. If exact reproducibility across implementations matters for your application, using `newdofcs` will give you that at the cost of using a slightly different dof formula than the R default.

## Dependencies

The default code path requires Python (via Stata's Python integration) with the following packages:
- `numpy`
- `scipy` (for LP via `scipy.optimize.linprog`)
- `osqp` (version 0.6.x; the quadratic programming solver that R's `TestMechs` also uses under the hood)

If any of these are unavailable, users can add the `newdofcs` option, which uses only Mata-native linear algebra and has no Python dependency.

## What we have NOT tested

- The `ARP` and `FSST` methods are not yet supported in the Stata port. Only `method(CS)` is implemented.
- The `reg_formula` argument for non-experimental designs is not yet supported.
- The `frac_ATs_affected` argument is not yet supported.
- The `use_binary` optimization branch is not yet ported.

Users requiring any of these features should use the R package until the Stata port extends to cover them.

## Reporting a discrepancy

If you run into a case where the Stata `testmechs` and the R `TestMechs` produce materially different statistical conclusions (not just numerical differences at the 3rd or 4th decimal place), please open an issue with:
- The data (or a reproducible subset).
- The exact command line used in both packages.
- Both output p-values and test statistics.
- Your Python, scipy, osqp versions (`python: import scipy, osqp; print(scipy.__version__, osqp.__version__)`).
