clear all
set more off

* ============================================================
* Test translation for R testthat: test-regression-equivalence.R
*
* Translates the subset that maps to Stata port functionality.
* Skipped (not implemented in Stata port):
*   - augmented mother_data_plus_fake with `| fake` FE syntax
*   - fixest `| interviewer` FE syntax
*   - fixest `| ... | treat ~ treativ` IV-with-FE syntax
*   - bounds_ade_ats tests
*   - partial_density_plot tests
*
* Core invariant tested here: baseline (no reg_formula) must equal
* reg_formula("~ treat") to numerical precision. `~ treat` is
* mathematically the same as sample-mean estimation per (d, m, y)
* cell, so any regression-branch bug shows up here immediately.
* ============================================================

adopath ++ .
adopath ++ src

local tol 1e-6

* ============================================================
* Test 1: test_sharp_null baseline == "~ treat"  (binary M)
* R block: "test_sharp_null regression adjustments match baseline (binary M)"
* ============================================================
use "data/mother_data.dta", clear

testmechs_test_sharpnull treat grandmother motherfinancial, ///
    method(CS) numybins(5) cluster(uc)
scalar base_pval_gm    = r(pval)
scalar base_tstat_gm   = r(test_stat)

testmechs_test_sharpnull treat grandmother motherfinancial, ///
    method(CS) numybins(5) cluster(uc) reg_formula("~ treat")
scalar trivial_pval_gm  = r(pval)
scalar trivial_tstat_gm = r(test_stat)

assert abs(base_pval_gm  - trivial_pval_gm)  < `tol'
assert abs(base_tstat_gm - trivial_tstat_gm) < `tol'


* ============================================================
* Test 2: test_sharp_null baseline == "~ treat"  (non-binary M, K=5)
* R block: "test_sharp_null regression adjustments match baseline (nonbinary M)"
*
* THIS EXERCISES YESTERDAY'S FIX: reg_prob was corrupting IF
* vectors for rare-cell m-values with K>=5. If that fix regresses,
* this test fails.
* ============================================================
use "data/mother_data.dta", clear

testmechs_test_sharpnull treat relationship_husb motherfinancial, ///
    method(CS) numybins(5) cluster(uc)
scalar base_pval_rel    = r(pval)
scalar base_tstat_rel   = r(test_stat)

testmechs_test_sharpnull treat relationship_husb motherfinancial, ///
    method(CS) numybins(5) cluster(uc) reg_formula("~ treat")
scalar trivial_pval_rel  = r(pval)
scalar trivial_tstat_rel = r(test_stat)

assert abs(base_pval_rel  - trivial_pval_rel)  < `tol'
assert abs(base_tstat_rel - trivial_tstat_rel) < `tol'


* ============================================================
* Test 3: test_sharp_null baseline == "~ treat"  (combined K=10)
* R block: "lb_frac_affected matches across mediator vector regressions"
* (adapted here for test_sharp_null to cover the K=10 case)
* ============================================================
use "data/mother_data.dta", clear

testmechs_test_sharpnull treat relationship_husb grandmother motherfinancial, ///
    method(CS) numybins(5) cluster(uc)
scalar base_pval_comb    = r(pval)
scalar base_tstat_comb   = r(test_stat)

testmechs_test_sharpnull treat relationship_husb grandmother motherfinancial, ///
    method(CS) numybins(5) cluster(uc) reg_formula("~ treat")
scalar trivial_pval_comb  = r(pval)
scalar trivial_tstat_comb = r(test_stat)

assert abs(base_pval_comb  - trivial_pval_comb)  < `tol'
assert abs(base_tstat_comb - trivial_tstat_comb) < `tol'


* ============================================================
* Test 4: lb_frac_affected baseline == "~ treat"  (binary M)
* R block: "lb_frac_affected respects regression choices" (partial)
* ============================================================
use "data/mother_data.dta", clear

testmechs_lb_fracaffected treat grandmother motherfinancial, ///
    numybins(5) atgroup(0) allowmindefiers
scalar base_lb_gm = r(lb)

testmechs_lb_fracaffected treat grandmother motherfinancial, ///
    numybins(5) atgroup(0) reg_formula("~ treat") allowmindefiers
scalar trivial_lb_gm = r(lb)

assert abs(base_lb_gm - trivial_lb_gm) < `tol'


* ============================================================
* Test 5: lb_frac_affected baseline == "~ treat"  (non-binary M)
* R block: "lb_frac_affected respects regression choices" (partial)
* ============================================================
use "data/mother_data.dta", clear

testmechs_lb_fracaffected treat relationship_husb motherfinancial, ///
    numybins(5) allowmindefiers
scalar base_lb_rel = r(lb)

testmechs_lb_fracaffected treat relationship_husb motherfinancial, ///
    numybins(5) reg_formula("~ treat") allowmindefiers
scalar trivial_lb_rel = r(lb)

assert abs(base_lb_rel - trivial_lb_rel) < `tol'


noi di as result "test_regression_equivalence.do passed"
