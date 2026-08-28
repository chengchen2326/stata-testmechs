clear all
set more off

* ============================================================
* Test translation for R testthat: test-discretize_y.R
*
* Tests inline copy of testmechs__discretize_y (quantile binning).
* KEEP IN SYNC with src/testmechs_test_sharpnull.ado line ~300.
* ============================================================

adopath ++ .
adopath ++ src

mata:
real colvector dy_test(real colvector y, real scalar numBins)
{
    real colvector ys, probs, cuts, yout
    real scalar n, i, p, h, h0, g

    y = vec(y)
    ys = uniqrows(sort(y,1))
    if (rows(ys) <= numBins) return(y)

    n = rows(y)
    ys = sort(y,1)
    probs = (1::(numBins-1)):/numBins
    cuts = J(rows(probs),1,.)
    for (i=1; i<=rows(probs); i++) {
        p = probs[i]
        h = (n-1)*p + 1
        h0 = floor(h)
        g = h - h0
        if (h0 >= n) cuts[i] = ys[n]
        else cuts[i] = (1-g)*ys[h0] + g*ys[h0+1]
    }
    cuts = uniqrows(sort(cuts,1))

    yout = J(n,1,1)
    for (i=1; i<=rows(cuts); i++) {
        yout = yout :+ (y :> cuts[i])
    }
    return(yout)
}

// Test 1: 5 values into 2 bins -> c(1,1,1,2,2)
r1 = dy_test((1\2\3\4\5), 2)
e1 = (1\1\1\2\2)
if (max(abs(r1-e1)) > 1e-12) _error(9)

// Test 2: 6 values into 3 bins -> c(1,1,2,2,3,3)
r2 = dy_test((1\2\3\4\5\6), 3)
e2 = (1\1\2\2\3\3)
if (max(abs(r2-e2)) > 1e-12) _error(9)

// Test 3: ties at low end -> c(1,1,1,1,1,2,2,3,3,3)
r3 = dy_test((1\1\1\1\1\6\7\8\9\10), 3)
e3 = (1\1\1\1\1\2\2\3\3\3)
if (max(abs(r3-e3)) > 1e-12) _error(9)

printf("All 3 discretize_y tests passed\n")
end

noi di as result "test_discretize_y.do passed"
