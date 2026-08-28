# External validation scripts

These scripts verify Stata port output against R's TestMechs numerically.
They require R with `TestMechs`, `fixest`, `sandwich`, and `osqp` installed,
and the private data file `data/mother_data_extended.dta`
(mother_data + noise_rand + treat2/3/4, generated locally with
R `rnorm(seed=123)` — not committed to the repo).

Unlike `tests/*.do`, these are **not** unit tests — they compare against
an external tool (R) and are meant for hand-run verification during
package development, not for CI.


## IV Influence Function parity (2026-08-28)

Verifies that our Stata IV IF implementation (Phase 2b, commit 65bb166)
produces the same per-observation influence function as R's
`sandwich::estfun(feols) %*% t(sandwich::bread(feols))`.

### How to re-run

```bash
# R side: dumps IF for 3 instruments (treat2/3/4) to /tmp/r_iv_if_*.csv
Rscript tests/verify/dump_iv_if_R_all.R

# Stata side: dumps same to /tmp/stata_iv_if_*.csv
# (must run after R because it borrows R's y_bin for exact alignment)
stata -e "do tests/verify/dump_iv_if_stata_all.do"

# Compare
python3 tests/verify/diff_iv_if.py
```

### Expected result (as of 2026-08-28)

All 3 instruments match to double-precision (max|diff| < 3e-14):

```
treat2: max|diff|=2.7978e-14  mean|diff|=2.7837e-15  [PASS]
treat3: max|diff|=1.5987e-14  mean|diff|=1.0300e-15  [PASS]
treat4: max|diff|=2.3981e-14  mean|diff|=1.2516e-15  [PASS]
```


## IV end-to-end pval parity (2026-08-28)

Verifies `test_sharp_null` with IV `reg_formula` produces the same pval
as R for K=5 and K=10 cases.

### How to re-run

```bash
Rscript tests/verify/check_iv_e2e_R.R
# then run the corresponding Stata commands manually and compare
```

### Expected result

| Case | R T_CC | Stata T | R pval | Stata pval | Match |
|---|---|---|---|---|---|
| grandmother K=2 | (via binary_m path) | 6.719854 | 0.0239976 | 0.0347378 | **NO — architectural difference** |
| relationship_husb K=5 | 10.5492844 | 10.549284 | 0.1595230 | 0.15952302 | YES |
| combined K=10 | 9.2446004 | 9.244610 | 0.2355686 | 0.23556796 | YES |

**K=2 discrepancy is expected:** R auto-dispatches binary M to a
specialised `test_sharp_null_binary_m` path (with different Cox-Shi
variant). The Stata port only implements the general CS code path,
so K=2 uses the same math as K=5+ instead of the R binary-M
optimisation. Users get valid CS test output for K=2, just not
identical to R's binary-M path.

To port the R binary-M path to Stata would require translating
`r_reference/TestMechs-master/R/test_sharp_null_binary_m.R` and
`test_sharp_null_coxandshi_binary_m.R` — deferred pending Jon's
input on priority.
