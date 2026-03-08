clear all
set more off

adopath ++ .
adopath ++ src

use "data/mother_data.dta", clear

testmechs_test_sharpnull treat grandmother motherfinancial, method(CS) numybins(5) cluster(uc)
return list

assert abs(r(pval) - 0.02283916) < 1e-4

noi di as result "test_test_sharpnull_cs.do passed"
