program define testmechs_test_sharpnull, rclass
    version 16.0
    syntax varlist(min=3 max=3 numeric) [if] [in] , METHOD(string) [NUMYBINS(integer 5) CLUSTER(varname)]

    if ("`method'" != "CS") {
        di as err "Only method(CS) is supported in this MVP implementation"
        exit 198
    }

    marksample touse

    local d : word 1 of `varlist'
    local m : word 2 of `varlist'
    local y : word 3 of `varlist'

    tempvar touse2 ywork
    quietly gen byte `touse2' = `touse'
    quietly replace `touse2' = 0 if missing(`d', `m', `y')
    if ("`cluster'" != "") quietly replace `touse2' = 0 if missing(`cluster')

    quietly count if `touse2'
    if (r(N) == 0) {
        di as err "No non-missing observations in estimation sample"
        exit 2000
    }

    quietly levelsof `m' if `touse2', local(mlevels)
    local nm : word count `mlevels'
    if (`nm' != 2) {
        di as err "Only binary mediator is supported; m must have exactly two observed levels"
        exit 198
    }

    quietly count if `touse2' & !inlist(`d',0,1)
    if (r(N) > 0) {
        di as err "Treatment variable d must be coded 0/1"
        exit 198
    }

    quietly count if `touse2' & `d' == 0
    local n0 = r(N)
    quietly count if `touse2' & `d' == 1
    local n1 = r(N)
    if (`n0' == 0 | `n1' == 0) {
        di as err "Both treatment groups (d==0 and d==1) must be present"
        exit 2000
    }

    quietly levelsof `y' if `touse2', local(yuniq)
    local nyuniq : word count `yuniq'

    if (`nyuniq' <= `numybins') {
        quietly gen double `ywork' = `y' if `touse2'
    }
    else {
        local plist
        forvalues i = 1/`=`numybins'-1' {
            local p = 100*`i'/`numybins'
            local plist `plist' `p'
        }
        tempname qmat
        quietly _pctile `y' if `touse2', p(`plist')
        matrix `qmat' = r(r1)'

        local cuts
        local prev
        forvalues i = 1/`=rowsof(`qmat')' {
            scalar __q = `qmat'[`i',1]
            local q : display %24.16g scalar(__q)
            if ("`prev'" == "" | `q' != `prev') {
                local cuts `cuts' `q'
                local prev `q'
            }
        }

        quietly gen int `ywork' = . if `touse2'
        local j = 0
        foreach c of local cuts {
            local ++j
            quietly replace `ywork' = `j' if `touse2' & missing(`ywork') & `y' <= `c'
        }
        quietly replace `ywork' = `j' + 1 if `touse2' & missing(`ywork')
    }

    if ("`cluster'" == "") {
        tempvar clusterid
        quietly gen long `clusterid' = _n if `touse2'
        local cluster `clusterid'
    }

    mata: testmechs__sharpnull_cs("`d'","`m'","`ywork'","`cluster'","`touse2'",0.05)

    return scalar pval = scalar(__tm_pval)
    di as result "p-value = " %10.8f r(pval)
end

mata:
real scalar testmechs__subset_min_qp(real colvector Y, real matrix Sigma, real matrix A, real colvector b, real scalar tol, real colvector best_mu)
{
    real scalar q, nmask, mask, j, stat, best
    real matrix Q, As, G
    real colvector bs, lambda, mu, feas, idx
    real rowvector pick

    q = rows(A)
    Q = invsym(Sigma)
    nmask = 2^q
    best = .
    best_mu = J(rows(Y),1,.)

    for (mask = 0; mask < nmask; mask++) {
        pick = J(1, q, 0)
        for (j = 1; j <= q; j++) {
            if (bitand(mask, 2^(j-1)) != 0) pick[1,j] = 1
        }
        idx = selectindex(pick')

        if (rows(idx) == 0) {
            mu = Y
        }
        else {
            As = A[idx,.]
            bs = b[idx,.]
            G = As * Sigma * As'
            if (min(symeigenvalues(G)) <= tol) continue
            lambda = invsym(G) * (As*Y - bs)
            mu = Y - Sigma * As' * lambda
        }

        feas = A*mu - b
        if (max(feas) > 1e-7) continue

        stat = quadcross(Y-mu, Q, Y-mu)
        if (missing(best) | stat < best) {
            best = stat
            best_mu = mu
        }
    }

    return(best)
}

void testmechs__sharpnull_cs(string scalar dvar, string scalar mvar, string scalar yvar,
    string scalar clvar, string scalar tousevar, real scalar alpha)
{
    real colvector d, m, y, cl, yvals, mvals
    real scalar n, j, k, g, pd0, pd1, tol, test_stat, df, cv, pval, tauhat, betahat, bindnorm
    real colvector beta, Yorig, Y, muhat, eval, idx, p00, p01, p10, p11, b, resid, isbind, numerator, denominator, valid, clu
    real matrix ifs, ifs_cl, Sigma, S, V, M, A, AsigA

    d = st_data(., dvar, tousevar)
    m = st_data(., mvar, tousevar)
    y = st_data(., yvar, tousevar)
    cl = st_data(., clvar, tousevar)

    yvals = uniqrows(sort(y,1))
    mvals = uniqrows(sort(m,1))

    n = rows(y)
    pd0 = mean(d :== 0)
    pd1 = mean(d :== 1)

    j = rows(yvals)
    p00 = J(j,1,.)
    p01 = J(j,1,.)
    p10 = J(j,1,.)
    p11 = J(j,1,.)

    for (k = 1; k <= j; k++) {
        p00[k] = mean((y :== yvals[k]) :& (m :== mvals[1]) :& (d :== 0)) / pd0
        p01[k] = mean((y :== yvals[k]) :& (m :== mvals[1]) :& (d :== 1)) / pd1
        p10[k] = mean((y :== yvals[k]) :& (m :== mvals[2]) :& (d :== 0)) / pd0
        p11[k] = mean((y :== yvals[k]) :& (m :== mvals[2]) :& (d :== 1)) / pd1
    }

    beta = (p00 - p01) \ (p11 - p10)
    Yorig = -beta

    ifs = J(n, 2*j, .)
    for (k = 1; k <= j; k++) {
        ifs[,k] = (d :== 0) :* (((y :== yvals[k]) :& (m :== mvals[1])) :- p00[k]) / pd0 -
                   (d :== 1) :* (((y :== yvals[k]) :& (m :== mvals[1])) :- p01[k]) / pd1

        ifs[,j+k] = (d :== 1) :* (((y :== yvals[k]) :& (m :== mvals[2])) :- p11[k]) / pd1 -
                     (d :== 0) :* (((y :== yvals[k]) :& (m :== mvals[2])) :- p10[k]) / pd0
    }

    clu = uniqrows(sort(cl,1))
    g = rows(clu)
    ifs_cl = J(g, cols(ifs), 0)
    for (k = 1; k <= g; k++) {
        ifs_cl[k,.] = colsum(select(ifs, cl :== clu[k]))
    }

    Sigma = variance(ifs_cl :* (g/n)) / g

    tol = 1e-8
    eval = symeigenvalues(Sigma)
    if (min(eval) < tol) {
        symeigensystem(Sigma, V, eval)
        idx = selectindex(eval :> tol)
        M = I(rows(eval))[idx,.]
        Y = M * V' * Yorig
        S = diag(eval[idx])
        A = V * M'
        b = A * Y - Yorig
    }
    else {
        Y = Yorig
        S = Sigma
        A = I(rows(Sigma))
        b = J(rows(A),1,0)
    }

    test_stat = testmechs__subset_min_qp(Y, S, A, b, tol, muhat)

    resid = A*muhat - b
    isbind = abs(resid) :< 1e-5
    df = sum(isbind)

    if (df == 0) {
        cv = 0
        pval = 1
    }
    else if (df != 1) {
        cv = invchi2(df, 1-alpha)
        pval = chi2tail(df, test_stat)
    }
    else {
        AsigA = A * S * A'
        bindnorm = sqrt(quadcross(isbind, AsigA, isbind))
        numerator = -bindnorm :* (A*muhat - b)
        denominator = bindnorm :* sqrt(diagonal(AsigA)) - AsigA * isbind
        valid = selectindex((denominator :> 0) :& (isbind :== 0))

        if (rows(valid) > 0) {
            tauhat = min(numerator[valid] :/ denominator[valid])
            betahat = 2*alpha*normal(tauhat)
        }
        else {
            tauhat = 0
            betahat = alpha
        }

        cv = invchi2(1, 1-betahat)
        pval = chi2tail(1, test_stat) / (2*normal(tauhat))
    }

    st_numscalar("__tm_pval", pval)
    st_numscalar("__tm_test_stat", test_stat)
    st_numscalar("__tm_cv", cv)
}
end
