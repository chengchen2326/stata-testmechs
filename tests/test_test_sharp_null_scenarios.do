clear all
set more off

* ============================================================
* Test translation for R testthat: test-test_sharp_null.R
*
* Translates the CS-method subset (Stata port doesn't support
* ARP or FSST). Also skips the use_binary check because Stata
* port has no use_binary option.
*
* All tests check correctness (reject / not reject), not
* numerical parity with R.
*
* Uses fixed RNG seed for reproducibility.
* ============================================================

adopath ++ .
adopath ++ src

capture program drop testmechs_test_sharpnull
capture program drop testmechs__reg_prob
discard

local alpha 0.05


* ============================================================
* Scenario A: burstzyn data, treated==control -> do NOT reject
* R block: "Do not reject when treated and control outcomes are the same"
* ============================================================
use "data/burstzyn_data.dta", clear
keep if !missing(condition2, signed_up_number, index)

* Build df_test: bind rows where condition2==1 with a copy where
* condition2 has been set to 0. Result: two identical outcome
* distributions labeled as different treatments.
preserve
    keep if condition2 == 1
    tempfile treated_only
    save `treated_only', replace
restore

preserve
    keep if condition2 == 1
    replace condition2 = 0
    append using `treated_only'

    testmechs_test_sharpnull condition2 signed_up_number index, method(CS)
    scalar pval_A = r(pval)
restore

assert pval_A >= `alpha'
noi di "Scenario A: pval = " %6.4f pval_A " >= " %4.2f `alpha' " (not rejected as expected)"


* ============================================================
* Scenario B: Y independent of D, binary M -> do NOT reject
* R block: "Do not reject when Y does not depends on D, binary M"
*
* N=1000, d~Bernoulli(0.5), m~Bernoulli(0.5), m[d==1]=1,
* y = 2*m + N(0, 0.5)   (Y depends on M but not D)
* ============================================================
clear
set seed 20260828
set obs 1000
gen byte treated  = (runiform() < 0.5)
gen byte mediator = (runiform() < 0.5)
replace mediator = 1 if treated == 1
gen double outcome = 2*mediator + rnormal(0, 0.5)

testmechs_test_sharpnull treated mediator outcome, method(CS) numybins(5)
scalar pval_B = r(pval)

assert pval_B >= `alpha'
noi di "Scenario B: pval = " %6.4f pval_B " >= " %4.2f `alpha' " (not rejected as expected)"


* ============================================================
* Scenario C: Y depends on both D and M, binary M -> REJECT
* R block: "Should reject when Y depends on both M and D, binary M"
*
* Same setup as B but y = 2*m + 50*d + N(0, 0.5)
* ============================================================
clear
set seed 20260828
set obs 1000
gen byte treated  = (runiform() < 0.5)
gen byte mediator = (runiform() < 0.5)
replace mediator = 1 if treated == 1
gen double outcome = 2*mediator + 50*treated + rnormal(0, 0.5)

testmechs_test_sharpnull treated mediator outcome, method(CS) numybins(5)
scalar pval_C = r(pval)

assert pval_C < `alpha'
noi di "Scenario C: pval = " %6.4f pval_C " < " %4.2f `alpha' " (rejected as expected)"


* ============================================================
* Scenario D: Y independent of D, non-binary M -> do NOT reject
* R block: "Do not reject when Y does not depends on D, non-binary M"
*
* N=1000, m ~ {0,1,2}, m[d==1] ~ {1,2}, y = 2*m + N(0, 0.5)
* ============================================================
clear
set seed 20260828
set obs 1000
gen byte treated  = (runiform() < 0.5)
* m ~ uniform on {0,1,2}
gen double _u = runiform()
gen byte mediator = 0
replace mediator = 1 if _u >= 1/3 & _u < 2/3
replace mediator = 2 if _u >= 2/3
drop _u
* For treated: m ~ uniform on {1,2}
replace mediator = 1 + (runiform() < 0.5) if treated == 1
gen double outcome = 2*mediator + rnormal(0, 0.5)

testmechs_test_sharpnull treated mediator outcome, method(CS) numybins(5)
scalar pval_D = r(pval)

assert pval_D >= `alpha'
noi di "Scenario D: pval = " %6.4f pval_D " >= " %4.2f `alpha' " (not rejected as expected)"


* ============================================================
* Scenario E: Y depends on both D and M, non-binary M -> REJECT
* R block: "Do not reject when Y does not depends on D, non-binary M"
*   (R comment mislabels this - the code is a reject-case)
*
* Same setup as D but y = 2*m + 50*d + N(0, 0.5)
* ============================================================
clear
set seed 20260828
set obs 1000
gen byte treated  = (runiform() < 0.5)
gen double _u = runiform()
gen byte mediator = 0
replace mediator = 1 if _u >= 1/3 & _u < 2/3
replace mediator = 2 if _u >= 2/3
drop _u
replace mediator = 1 + (runiform() < 0.5) if treated == 1
gen double outcome = 2*mediator + 50*treated + rnormal(0, 0.5)

testmechs_test_sharpnull treated mediator outcome, method(CS) numybins(5)
scalar pval_E = r(pval)

assert pval_E < `alpha'
noi di "Scenario E: pval = " %6.4f pval_E " < " %4.2f `alpha' " (rejected as expected)"


noi di as result "test_test_sharp_null_scenarios.do passed"
