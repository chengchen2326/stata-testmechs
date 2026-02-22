clear all
set more off

* Minimal validation dataset (binary D, M, Y).
input byte(d m y)
1 1 1
1 1 1
1 0 1
1 0 1
0 1 0
0 1 0
0 0 1
0 0 1
end

* Ensure local ado is on adopath when running from repo root.
adopath ++ .

noi di "Running pooled bound"
testmechs_lb_fracaffected d m y
scalar pooled = r(lb)
return list

noi di "Running AT-group bound (m=1)"
testmechs_lb_fracaffected d m y, atgroup(1)
scalar at_lb = r(lb)
return list

noi di "Running NT-group bound (m=0)"
testmechs_lb_fracaffected d m y, atgroup(0)
scalar nt_lb = r(lb)
return list

* Expected in this toy setup:
* - AT-group bound = 1
* - NT-group bound = 0
* - pooled bound   = 0.5
assert abs(at_lb - 1) < 1e-8
assert abs(nt_lb - 0) < 1e-8
assert abs(pooled - 0.5) < 1e-8

noi di as result "test_lb_fracaffected.do passed"
