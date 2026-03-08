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
	local m0 : word 1 of `mlevels'
    local m1 : word 2 of `mlevels'

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

* --- R-style cutpoints: compute on ALL non-missing y (not just touse2) ---
tempname q20 q40 q60 q80
quietly _pctile `y' if !missing(`y'), p(20 40 60 80)

scalar `q20' = r(r1)
scalar `q40' = r(r2)
scalar `q60' = r(r3)
scalar `q80' = r(r4)

scalar __tm_q20 = scalar(`q20')
scalar __tm_q40 = scalar(`q40')
scalar __tm_q60 = scalar(`q60')
scalar __tm_q80 = scalar(`q80')

* --- Apply fixed cutpoints to estimation sample (touse2) ---
quietly gen byte `ywork' = . if `touse2'
quietly replace `ywork' = 1 if `touse2' & `y' <= scalar(`q20')
quietly replace `ywork' = 2 if `touse2' & `y' >  scalar(`q20') & `y' <= scalar(`q40')
quietly replace `ywork' = 3 if `touse2' & `y' >  scalar(`q40') & `y' <= scalar(`q60')
quietly replace `ywork' = 4 if `touse2' & `y' >  scalar(`q60') & `y' <= scalar(`q80')
quietly replace `ywork' = 5 if `touse2' & `y' >  scalar(`q80')

* Debug: check how many bins in estimation sample
quietly levelsof `ywork' if `touse2', local(ybins)
scalar __tm_j_stata = wordcount("`ybins'")

	
   tempvar clusterid

if ("`cluster'" == "") {
    quietly gen long `clusterid' = _n if `touse2'
}
else {
    capture confirm variable `cluster'
    if _rc {
        di as err "cluster() variable not found: `cluster'"
        exit 111
    }

    * Robust: if cluster is string, use encode (keeps distinct strings distinct).
    * If cluster is numeric, use egen group() on the estimation sample.
    capture confirm numeric variable `cluster'
    if _rc {
        quietly encode `cluster' if `touse2', gen(`clusterid')
    }
    else {
        quietly egen long `clusterid' = group(`cluster') if `touse2'
    }

    * sanity check: must have at least 2 clusters
    quietly levelsof `clusterid' if `touse2', local(clu)
    local G : word count `clu'
	scalar __tm_G_stata = `G'
    if (`G' < 2) {
        di as err "Need at least 2 clusters in estimation sample (cluster(`cluster') has only `G')"
        exit 198
    }
}

mata: testmechs__sharpnull_cs("`d'","`m'","`ywork'","`clusterid'","`touse2'",0.05, `m0', `m1')

    return scalar pval      = scalar(__tm_pval)
return scalar test_stat = scalar(__tm_test_stat)
return scalar cv        = scalar(__tm_cv)

di as result "test_stat = " %10.6f scalar(__tm_test_stat)
di as result "cv        = " %10.6f scalar(__tm_cv)
di as result "p-value   = " %10.8f scalar(__tm_pval)
end



mata:

real scalar __bitand_scalar(real scalar a, real scalar b)
{
    real scalar res, bit, aa, bb

    aa  = floor(a)
    bb  = floor(b)
    res = 0
    bit = 1

    // works for nonnegative integer masks
    while (aa > 0 | bb > 0) {
        if (mod(aa,2)==1 & mod(bb,2)==1) res = res + bit
        aa  = floor(aa/2)
        bb  = floor(bb/2)
        bit = bit*2
    }
    return(res)
}

real matrix bitand(real matrix A, real matrix B)
{
    real matrix AA, BB, R
    real scalar i, j

    // expand scalars to match the other argument
    if (rows(A)==1 & cols(A)==1) AA = J(rows(B), cols(B), A[1,1])
    else AA = A
    if (rows(B)==1 & cols(B)==1) BB = J(rows(AA), cols(AA), B[1,1])
    else BB = B

    if (rows(AA)!=rows(BB) | cols(AA)!=cols(BB)) _error(3200, "bitand(): nonconformable arguments")

    R = J(rows(AA), cols(AA), 0)
    for (i=1; i<=rows(AA); i++) {
        for (j=1; j<=cols(AA); j++) {
            R[i,j] = __bitand_scalar(AA[i,j], BB[i,j])
        }
    }
    return(R)
}

real scalar testmechs__subset_min_qp(real colvector Y, real matrix Sigma, real matrix A, real colvector b, real scalar tol, real colvector best_mu)
{
    real scalar q, nmask, mask, j, stat, best
    real matrix Q, As, G
    real colvector bs, lambda, mu, feas, idx, v, pick

    Y = vec(Y)
    b = vec(b)

    q = rows(A)
    Q = pinv(Sigma)

    nmask = 2^q
    best  = .
	real scalar bestmask, bestminfeas
    bestmask = .
    bestminfeas = .
	real scalar n_total, n_feas, n_best
    n_total = 0
    n_feas  = 0
    n_best  = 0
    best_mu = J(rows(Y),1,.)

    for (mask=0; mask<nmask; mask++) {

        // build subset indicator without bitand()
        pick = J(q,1,0)
		n_total = n_total + 1
        for (j=1; j<=q; j++) {
            if (mod(floor(mask/2^(j-1)),2)==1) pick[j] = 1
        }
        idx = selectindex(pick)

        if (rows(idx)==0) {
            mu = Y
        }
        else {
            As = A[idx,.]
            bs = b[idx]

            G = As*Sigma*As'
            // small ridge to stabilize pseudo-inverse
            G = G + (1e-12)*I(rows(G))

            lambda = pinv(G) * (As*Y - bs)
            mu = Y - Sigma*As'*lambda
        }

        feas = A*mu - b
        real scalar feas_tol
        feas_tol = 1e-8
        if (min(feas) < -feas_tol) continue
        n_feas = n_feas + 1
        v = vec(Y - mu)
        stat = (v' * Q * v)[1,1]

        if (missing(best) | stat < best - 1e-12) {
            best = stat
            best_mu = mu
			bestmask = mask
            bestminfeas = min(feas)
        }
    }
    st_numscalar("__tm_qp_total", n_total)
    st_numscalar("__tm_qp_feas",  n_feas)
    st_numscalar("__tm_qp_best",  n_best)
	st_numscalar("__tm_bestmask", bestmask)
    st_numscalar("__tm_bestminfeas", bestminfeas)
    st_numscalar("__tm_beststat", best)
    return(best)
}

void testmechs__sharpnull_cs(string scalar dvar, string scalar mvar, string scalar yvar,
    string scalar clvar, string scalar tousevar, real scalar alpha, real scalar m0, real scalar m1)
{
    real colvector d, m, y, cl, yvals, mvals
    real scalar n, j, k, g, pd0, pd1, tol, test_stat, df, cv, pval, tauhat, betahat, bindnorm
    real colvector beta, Yorig, Y, muhat, eval, idx, p00, p01, p10, p11, b, resid, isbind, numerator, denominator, valid, clu
    real matrix ifs, ifs_cl, Sigma, S, V, M, A, AsigA

       real colvector keep

    // Hard filter: build an explicit index from the 0/1 touse variable
    keep = selectindex(st_data(., tousevar) :!= 0)
    if (rows(keep) == 0) _error(2000, "No observations in estimation sample")

    // Read only kept observations
    d  = st_data(keep, dvar)
    m  = st_data(keep, mvar)
    y  = st_data(keep, yvar)
    cl = st_data(keep, clvar)
	
	// ---- DEBUG A: right after reading cl (before missing filter) ----
d  = vec(d)
m  = vec(m)
y  = vec(y)
cl = vec(cl)

real colvector cluA
real scalar GA
cluA = uniqrows(sort(cl,1))
GA   = rows(cluA)

st_numscalar("__tm_G_A", GA)
st_numscalar("__tm_clA_min", min(cl))
st_numscalar("__tm_clA_max", max(cl))
st_numscalar("__tm_n_A", rows(cl))

    // Drop any remaining missings defensively
    keep = selectindex((d :< .) :& (m :< .) :& (y :< .) :& (cl :< .))
    d  = d[keep]
    m  = m[keep]
    y  = y[keep]
    cl = cl[keep]
	
	// ---- DEBUG B: after missing filter / subsetting ----
d  = vec(d)
m  = vec(m)
y  = vec(y)
cl = vec(cl)

real colvector cluB
real scalar GB
cluB = uniqrows(sort(cl,1))
GB   = rows(cluB)

st_numscalar("__tm_G_B", GB)
st_numscalar("__tm_clB_min", min(cl))
st_numscalar("__tm_clB_max", max(cl))
st_numscalar("__tm_n_B", rows(cl))
// debug C 
st_numscalar("__tm_step", 10)

    if (rows(d) == 0) _error(2000, "No complete-case observations after filtering")
	
yvals = uniqrows(sort(y,1))
mvals = (m0 \ m1)
if (m0 == m1) _error(198, "Binary mediator required (m levels are not distinct)")

// ---- TEST2 DEBUG: what does Mata see for clusters? ----
cl = vec(cl)

real colvector clu_tmp
real scalar G_mata

clu_tmp = uniqrows(sort(cl,1))
G_mata  = rows(clu_tmp)

st_numscalar("__tm_G_mata_fn", G_mata)
st_numscalar("__tm_cl_min_fn", min(cl))
st_numscalar("__tm_cl_max_fn", max(cl))

if (G_mata < 2) _error(198, "Need at least 2 clusters for cluster-robust inference")

    n = rows(y)
    pd0 = mean(d :== 0)
    pd1 = mean(d :== 1)

    j = rows(yvals)
	st_numscalar("__tm_j", j)
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
	st_matrix("__tm_beta", beta)
    st_matrix("__tm_Yorig", -beta)
    Yorig = -beta
	
	//debugD
	st_numscalar("__tm_step", 20)

    ifs = J(n, 2*j, .)
    for (k = 1; k <= j; k++) {
        ifs[,k] = (d :== 0) :* (((y :== yvals[k]) :& (m :== mvals[1])) :- p00[k]) / pd0 -
                   (d :== 1) :* (((y :== yvals[k]) :& (m :== mvals[1])) :- p01[k]) / pd1

        ifs[,j+k] = (d :== 1) :* (((y :== yvals[k]) :& (m :== mvals[2])) :- p11[k]) / pd1 -
                     (d :== 0) :* (((y :== yvals[k]) :& (m :== mvals[2])) :- p10[k]) / pd0
    }
	
	//debug E 
	st_numscalar("__tm_step", 30)

    clu = uniqrows(sort(cl,1))
    g = rows(clu)
    ifs_cl = J(g, cols(ifs), 0)
    for (k = 1; k <= g; k++) {
        ifs_cl[k,.] = colsum(select(ifs, cl :== clu[k]))
    }

    Sigma = (ifs_cl' * ifs_cl) / (n^2)  // since ifs_cl are cluster sums
	
	real scalar g_corr
    g_corr = (g > 1 ? g/(g-1) : 1)
    Sigma = g_corr * Sigma
	
	real matrix Sigma_var, Sigma_meat

Sigma_var  = variance(ifs_cl :* (g/n)) / g
Sigma_meat = (ifs_cl' * ifs_cl) / g   // 经典 1/g * sum S_g S_g'

st_numscalar("__tm_sig_trace_var",  trace(Sigma_var))
st_numscalar("__tm_sig_trace_meat", trace(Sigma_meat))
	
	//debug F
	st_numscalar("__tm_step", 40)

// DIAGNOSTIC: disable eigen-reduction, stay in full space
Y = Yorig
S = Sigma
A = -I(rows(Sigma))
b = J(rows(A),1,0)

// rank debug
st_numscalar("__tm_rank", rows(Sigma))
	
	//debug G
	st_numscalar("__tm_step", 50)

    test_stat = testmechs__subset_min_qp(Y, S, A, b, tol, muhat)
    //debug H
	st_numscalar("__tm_step", 60)
	
   resid = A*muhat - b
bind_idx = selectindex(abs(resid) :< 1e-8)

if (rows(bind_idx)==0) {
    df = 0
    cv = 0
    pval = 1
}
else {
    real matrix Sb
    Sb = S[bind_idx, bind_idx]
    df = rank(Sb)
    cv = invchi2(df, 1-alpha)
    pval = chi2tail(df, test_stat)
} 

    st_numscalar("__tm_pval", pval)
    st_numscalar("__tm_test_stat", test_stat)
    st_numscalar("__tm_cv", cv)
}
end
