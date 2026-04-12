clear all
set more off

adopath ++ .
adopath ++ src

use "data/mother_data.dta", clear

testmechs_test_sharpnull treat grandmother motherfinancial, method(CS) numybins(5) cluster(uc)
return list

assert abs(r(pval) - 0.02283916) < 1e-4

testmechs_test_sharpnull treat relationship_husb grandmother motherfinancial, method(CS) numybins(5) cluster(uc)
scalar p_multi_12 = r(pval)

testmechs_test_sharpnull treat grandmother relationship_husb motherfinancial, method(CS) numybins(5) cluster(uc)
scalar p_multi_21 = r(pval)

assert abs(p_multi_12 - p_multi_21) < 1e-10
assert abs(p_multi_12 - 0.6540863) < 1e-4

noi di as result "test_test_sharpnull_cs.do passed"
