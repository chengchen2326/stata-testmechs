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
    if (`nm' < 2) {
        di as err "Mediator variable m must have at least two observed levels"
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

mata: testmechs__sharpnull_cs("`d'","`m'","`ywork'","`clusterid'","`touse2'",0.05)

    return scalar pval      = scalar(__tm_pval)
return scalar test_stat = scalar(__tm_test_stat)
return scalar cv        = scalar(__tm_cv)
    return scalar df        = scalar(__tm_df)
    return scalar rank_sigma = scalar(__tm_rank_sigma)
    return scalar rank_pos_sigma = scalar(__tm_rank_pos_sigma)
    return scalar Aobs_rows = scalar(__tm_Aobs_rows)
    return scalar Aobs_cols = scalar(__tm_Aobs_cols)
    return scalar Ashp_rows = scalar(__tm_Ashp_rows)
    return scalar Ashp_cols = scalar(__tm_Ashp_cols)
    return matrix beta_obs = __tm_beta_obs
    return matrix sigma_eigenvalues = __tm_sigma_evals

di as result "test_stat = " %10.6f scalar(__tm_test_stat)
di as result "df        = " %10.0f scalar(__tm_df)
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
    string scalar clvar, string scalar tousevar, real scalar alpha)
{
    real colvector d, m, y, cl, keep, yvals, mvals, beta_obs, beta_shp, Y
    real colvector clu, binding_index, resid, evals
    real matrix A_obs, A_shp, A_full, ifs, ifs_cl, sigma_obs
    real matrix eigvecs, M, Atilde, sigma_work, A_work, sigma_inv
    real colvector btilde, b_work, Xstar, muhat
    real scalar n, j, K, g, pd0, pd1, tol, test_stat, cv, pval, df
    real scalar i, k, min_eval, rank_sigma, rank_pos
    real colvector p_m_d0, p_m_d1, p_ym_0_vec, p_ym_1_vec

    tol = 1e-8

    keep = selectindex(st_data(., tousevar) :!= 0)
    if (rows(keep) == 0) _error(2000, "No observations in estimation sample")

    d  = vec(st_data(keep, dvar))
    m  = vec(st_data(keep, mvar))
    y  = vec(st_data(keep, yvar))
    cl = vec(st_data(keep, clvar))

    keep = selectindex((d:<.) :& (m:<.) :& (y:<.) :& (cl:<.))
    d  = d[keep]
    m  = m[keep]
    y  = y[keep]
    cl = cl[keep]

    n = rows(y)
    if (n==0) _error(2000, "No complete observations after filtering")

    yvals = uniqrows(sort(y,1))
    mvals = uniqrows(sort(m,1))
    j = rows(yvals)
    K = rows(mvals)
    if (K < 2) _error(198, "Mediator must have at least 2 values")

    clu = uniqrows(sort(cl,1))
    g = rows(clu)
    if (g < 2) _error(198, "Need at least 2 clusters for cluster-robust inference")

    pd0 = mean(d:==0)
    pd1 = mean(d:==1)

    // ---- construct A.obs, A.shp, beta.shp for inequalities_only=TRUE, default ordering ----
    A_shp = J(K + (K^2 + j*K), K^2 + j*K, 0)
    // First K rows from shape equations after dropping nuisance and extraneous theta(l>k)
    for (k=1; k<=K; k++) {
        for (i=1; i<=K; i++) {
            if (i <= k & i != k) A_shp[k, (k-1)*K + i] = 1
        }
        A_shp[k, (K^2 + (k-1)*j + 1)..(K^2 + k*j)] = -J(1,j,1)
    }
    // nonnegativity constraints
    A_shp[(K+1)..rows(A_shp), .] = I(cols(A_shp))
    beta_shp = J(rows(A_shp),1,0)

    // A.obs for inequalities_only=TRUE
    A_obs = J(2*K + 2*K + j*K, K^2 + j*K, 0)
    for (k=1; k<=K; k++) {
        // P(M=k|D=0)
        for (i=1; i<=K; i++) if (i<=k) A_obs[k, (i-1)*K + k] = 1
        // P(M=k|D=1)
        for (i=1; i<=K; i++) if (k<=i) A_obs[K+k, (k-1)*K + i] = 1
        // opposite-sign duplicated equalities
        A_obs[2*K + k, .] = -A_obs[k,.]
        A_obs[3*K + k, .] = -A_obs[K+k,.]
        // sup Delta(A) constraints
        for (i=1; i<=j; i++) A_obs[4*K + (k-1)*j + i, K^2 + (k-1)*j + i] = 1
    }

    // ---- construct beta.obs ----
    p_ym_0_vec = J(j*K,1,0)
    p_ym_1_vec = J(j*K,1,0)
    p_m_d0 = J(K,1,0)
    p_m_d1 = J(K,1,0)
    for (k=1; k<=K; k++) {
        for (i=1; i<=j; i++) {
            p_ym_0_vec[(k-1)*j + i] = mean((y:==yvals[i]) :& (m:==mvals[k]) :& (d:==0)) / pd0
            p_ym_1_vec[(k-1)*j + i] = mean((y:==yvals[i]) :& (m:==mvals[k]) :& (d:==1)) / pd1
        }
        p_m_d0[k] = sum(p_ym_0_vec[((k-1)*j+1)..(k*j)])
        p_m_d1[k] = sum(p_ym_1_vec[((k-1)*j+1)..(k*j)])
    }
    beta_obs = p_m_d0 \ p_m_d1 \ (-p_m_d0) \ (-p_m_d1) \ (p_ym_1_vec - p_ym_0_vec)

    // ---- clustered analytic variance for beta.obs ----
    ifs = J(n, rows(beta_obs), 0)
    for (k=1; k<=K; k++) {
        // p_m_d0 and p_m_d1 IFs
        ifs[,k] = (d:==0) :* (((m:==mvals[k])) :- p_m_d0[k]) / pd0
        ifs[,K+k] = (d:==1) :* (((m:==mvals[k])) :- p_m_d1[k]) / pd1
        ifs[,2*K+k] = -ifs[,k]
        ifs[,3*K+k] = -ifs[,K+k]
        for (i=1; i<=j; i++) {
            real scalar col
            col = 4*K + (k-1)*j + i
            ifs[,col] = (d:==1) :* (((y:==yvals[i]) :& (m:==mvals[k])) :- p_ym_1_vec[(k-1)*j+i]) / pd1 -
                        (d:==0) :* (((y:==yvals[i]) :& (m:==mvals[k])) :- p_ym_0_vec[(k-1)*j+i]) / pd0
        }
    }

    ifs_cl = J(g, cols(ifs), 0)
    for (k=1; k<=g; k++) ifs_cl[k,.] = colsum(select(ifs, cl:==clu[k]))
    sigma_obs = variance(ifs_cl :* (g/n)) / g

    // ---- cox_shi_nonuisance input objects ----
    Y = A_obs * (beta_shp :- J(rows(beta_shp),1,0)) :- beta_obs
    A_full = A_obs * A_shp

    // rank-deficiency handling via eigen decomposition
    symeigensystem(sigma_obs, eigvecs, evals)
    min_eval = min(evals)
    rank_sigma = rank(sigma_obs)
    if (min_eval < tol) {
        real colvector pos
        pos = selectindex(evals :> tol)
        rank_pos = rows(pos)
        M = I(rows(evals))[pos,.]
        if (rank_pos==1) M = M'
        Xstar = M * eigvecs' * Y
        sigma_work = diag(evals[pos])
        Atilde = eigvecs * M'
        btilde = Atilde * Xstar - Y
        Y = Xstar
        A_work = A_full * Atilde
        b_work = A_full * btilde
    }
    else {
        rank_pos = rows(sigma_obs)
        sigma_work = sigma_obs
        A_work = A_full
        b_work = J(rows(A_work),1,0)
    }

    // our QP solver expects A*mu >= b. convert from A*mu <= b.
    test_stat = testmechs__subset_min_qp(Y, sigma_work, -A_work, -b_work, tol, muhat)
    resid = A_work*muhat - b_work
    binding_index = selectindex(abs(resid) :< 1e-5)
    df = rows(binding_index)

    if (df==0) {
        cv = 0
        pval = 1
    }
    else {
        cv = invchi2(df, 1-alpha)
        pval = chi2tail(df, test_stat)
    }

    st_matrix("__tm_beta_obs", beta_obs)
    st_matrix("__tm_sigma_evals", sort(evals,-1)')
    st_numscalar("__tm_Aobs_rows", rows(A_obs))
    st_numscalar("__tm_Aobs_cols", cols(A_obs))
    st_numscalar("__tm_Ashp_rows", rows(A_shp))
    st_numscalar("__tm_Ashp_cols", cols(A_shp))
    st_numscalar("__tm_rank_sigma", rank_sigma)
    st_numscalar("__tm_rank_pos_sigma", rank_pos)
    st_numscalar("__tm_test_stat", test_stat)
    st_numscalar("__tm_df", df)
    st_numscalar("__tm_cv", cv)
    st_numscalar("__tm_pval", pval)
}
end
