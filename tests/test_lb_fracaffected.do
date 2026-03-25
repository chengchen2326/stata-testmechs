clear all
set more off

* Test translation for R testthat: test-lb-frac-affected.R
* Runs only features currently supported by Stata MVP.

adopath ++ .
adopath ++ src

local tol 1e-6

* ============================================================
* Test 1: lb is ~0 when (Y,M) are independent of D
* ============================================================
use "data/mother_data.dta", clear
preserve
    keep if treat == 1
    tempfile t1
    save `t1', replace
restore

preserve
    keep if treat == 1
    replace treat = 0
    append using `t1'

    testmechs_lb_fracaffected treat grandmother motherfinancial, numybins(5)
    scalar lb_zero_gm = r(lb)
    assert abs(lb_zero_gm - 0) < `tol'

    testmechs_lb_fracaffected treat relationship_husb motherfinancial, numybins(5)
    scalar lb_zero_rel = r(lb)
    assert abs(lb_zero_rel - 0) < `tol'

    testmechs_lb_fracaffected treat grandmother relationship_husb motherfinancial, numybins(5)
    scalar lb_zero_pool2m = r(lb)
    assert abs(lb_zero_pool2m - 0) < `tol'
restore


* ============================================================
* Test 2: lb is ~1 when Y is colinear with D
* ============================================================
use "data/mother_data.dta", clear
preserve
    keep if treat == 1
    replace motherfinancial = 1
    tempfile t2
    save `t2', replace
restore

preserve
    keep if treat == 1
    replace treat = 0
    replace motherfinancial = 2
    append using `t2'

    testmechs_lb_fracaffected treat grandmother motherfinancial, numybins(5)
    scalar lb_one_gm = r(lb)
    assert abs(lb_one_gm - 1) < `tol'

    testmechs_lb_fracaffected treat relationship_husb motherfinancial, numybins(5)
    scalar lb_one_rel = r(lb)
    assert abs(lb_one_rel - 1) < `tol'
restore


* ============================================================
* Test 3: subgroup bounds (atgroup) with no effect for M=0 and full effect for M=1
* ============================================================
use "data/mother_data.dta", clear
preserve
    keep if treat == 1
    replace motherfinancial = 1
    tempfile t3
    save `t3', replace
restore

preserve
    keep if treat == 1
    replace treat = 0
    replace motherfinancial = 1 if grandmother == 0
    replace motherfinancial = 2 if grandmother == 1
    append using `t3'

    testmechs_lb_fracaffected treat grandmother motherfinancial, numybins(5) atgroup(0)
    scalar lb_group0 = r(lb)
    assert abs(lb_group0 - 0) < `tol'

    testmechs_lb_fracaffected treat grandmother motherfinancial, numybins(5) atgroup(1)
    scalar lb_group1 = r(lb)
    assert abs(lb_group1 - 1) < `tol'
restore


* ============================================================
* Test 4: binary analytic formula match (burstzyn_data)
* ============================================================
use "data/burstzyn_data.dta", clear
keep if !missing(condition2, signed_up_number, applied_out_fl)

quietly count if condition2 == 1
scalar N1 = r(N)
quietly count if condition2 == 0
scalar N0 = r(N)

quietly summarize signed_up_number if condition2 == 0, meanonly
scalar frac_ats = r(mean)
quietly summarize signed_up_number if condition2 == 1, meanonly
scalar frac_nts = 1 - r(mean)
scalar frac_cs = 1 - frac_ats - frac_nts

quietly count if condition2 == 1 & applied_out_fl == 1 & signed_up_number == 1
scalar p111 = r(N) / N1
quietly count if condition2 == 0 & applied_out_fl == 1 & signed_up_number == 1
scalar p011 = r(N) / N0
quietly count if condition2 == 1 & applied_out_fl == 0 & signed_up_number == 1
scalar p101 = r(N) / N1
quietly count if condition2 == 0 & applied_out_fl == 0 & signed_up_number == 1
scalar p001 = r(N) / N0
quietly count if condition2 == 1 & signed_up_number == 1
scalar pm11 = r(N) / N1
quietly count if condition2 == 0 & signed_up_number == 1
scalar pm01 = r(N) / N0

scalar diff_pd_1 = max(0, p111 - p011, p101 - p001, pm11 - pm01)
scalar lb_ats_formula = max(0, (diff_pd_1 - frac_cs) / frac_ats)

testmechs_lb_fracaffected condition2 signed_up_number applied_out_fl, numybins(2) atgroup(1)
scalar lb_ats_pkg = r(lb)
assert abs(lb_ats_pkg - lb_ats_formula) < 1e-3

quietly count if condition2 == 1 & applied_out_fl == 1 & signed_up_number == 0
scalar p110 = r(N) / N1
quietly count if condition2 == 0 & applied_out_fl == 1 & signed_up_number == 0
scalar p010 = r(N) / N0
quietly count if condition2 == 1 & applied_out_fl == 0 & signed_up_number == 0
scalar p100 = r(N) / N1
quietly count if condition2 == 0 & applied_out_fl == 0 & signed_up_number == 0
scalar p000 = r(N) / N0
quietly count if condition2 == 1 & signed_up_number == 0
scalar pm10 = r(N) / N1
quietly count if condition2 == 0 & signed_up_number == 0
scalar pm00 = r(N) / N0

scalar diff_pd_0 = max(0, p110 - p010, p100 - p000, pm10 - pm00)
scalar lb_nts_formula = max(0, diff_pd_0 / frac_nts)

testmechs_lb_fracaffected condition2 signed_up_number applied_out_fl, numybins(2) atgroup(0)
scalar lb_nts_pkg = r(lb)
assert abs(lb_nts_pkg - lb_nts_formula) < 1e-3

scalar lb_pooled_formula = (frac_ats / (frac_ats + frac_nts)) * lb_ats_formula + ///
                           (frac_nts / (frac_ats + frac_nts)) * lb_nts_formula

testmechs_lb_fracaffected condition2 signed_up_number applied_out_fl, numybins(2)
scalar lb_pooled_pkg = r(lb)
assert abs(lb_pooled_pkg - lb_pooled_formula) < 1e-3


* ============================================================
* Test 5: multi-valued M extremes agree with binarized versions
* ============================================================
use "data/mother_data.dta", clear
keep if !missing(treat, relationship_husb, motherfinancial)

generate byte relationship5 = (relationship_husb == 5)
generate byte relationshipnonzero = (relationship_husb > 1)

testmechs_lb_fracaffected treat relationship5 motherfinancial, ///
    numybins(5) atgroup(1) maxdefiersshare(0.02)
scalar lb_rel5_bin = r(lb)

testmechs_lb_fracaffected treat relationship_husb motherfinancial, ///
    numybins(5) atgroup(5) maxdefiersshare(0.02)
scalar lb_rel5_multi = r(lb)

assert abs(lb_rel5_bin - lb_rel5_multi) < 1e-3

testmechs_lb_fracaffected treat relationshipnonzero motherfinancial, ///
    numybins(5) atgroup(0) maxdefiersshare(0.02)
scalar lb_relnz_bin = r(lb)

testmechs_lb_fracaffected treat relationship_husb motherfinancial, ///
    numybins(5) atgroup(0) maxdefiersshare(0.02)
scalar lb_relnz_multi = r(lb)

assert abs(lb_relnz_bin - lb_relnz_multi) < 1e-3


* ============================================================
* Test 6: defier-share sensitivity slope (binary M)
* ============================================================
use "data/mother_data.dta", clear
keep if !missing(treat, grandmother, motherfinancial)

quietly summarize grandmother if treat == 1, meanonly
scalar frac_nts_md = 1 - r(mean)
scalar eps = 0.0001

testmechs_lb_fracaffected treat grandmother motherfinancial, numybins(5) atgroup(0)
scalar lb_base = r(lb)

testmechs_lb_fracaffected treat grandmother motherfinancial, numybins(5) atgroup(0) maxdefiersshare(`=eps')
scalar lb_eps = r(lb)

scalar slope_num = (lb_base - lb_eps) / eps
scalar slope_theory = (1 / frac_nts_md) * (1 - lb_base)
assert abs(slope_num - slope_theory) < 1e-3


* ============================================================
* Unsupported in current Stata MVP
* ============================================================
* TODO: R test "matches Lee bounds lower bound w binary M,Y" compares against
* bounds_ade_ats(). That comparison command is not part of the current Stata MVP,
* so the direct cross-command equivalence test is intentionally skipped here.

noi di as result "test_lb_fracaffected.do passed"
