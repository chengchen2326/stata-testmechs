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

    tempvar touse2 clusterid
    quietly gen byte `touse2' = `touse'
    quietly replace `touse2' = 0 if missing(`d', `m', `y')
    if ("`cluster'" != "") quietly replace `touse2' = 0 if missing(`cluster')

    quietly count if `touse2'
    if (r(N) == 0) {
        di as err "No non-missing observations in estimation sample"
        exit 2000
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

    if ("`cluster'" == "") {
        quietly gen long `clusterid' = _n if `touse2'
    }
    else {
        capture confirm variable `cluster'
        if _rc {
            di as err "cluster() variable not found: `cluster'"
            exit 111
        }

        capture confirm numeric variable `cluster'
        if _rc {
            quietly encode `cluster' if `touse2', gen(`clusterid')
        }
        else {
            quietly egen long `clusterid' = group(`cluster') if `touse2'
        }

        quietly levelsof `clusterid' if `touse2', local(clu)
        local G : word count `clu'
        if (`G' < 2) {
            di as err "Need at least 2 clusters in estimation sample (cluster(`cluster') has only `G')"
            exit 198
        }
    }

    mata: testmechs__sharpnull_cs("`d'","`m'","`y'","`clusterid'","`touse2'",0.05,`numybins',1)

    return scalar pval      = scalar(__tm_pval)
    return scalar test_stat = scalar(__tm_test_stat)
    return scalar cv        = scalar(__tm_cv)

    di as result "test_stat = " %10.6f scalar(__tm_test_stat)
    di as result "cv        = " %10.6f scalar(__tm_cv)
    di as result "p-value   = " %10.8f scalar(__tm_pval)
end

mata:

real colvector testmechs__discretize_y(real colvector y, real scalar numBins)
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

real matrix testmechs__build_Aobs(real scalar K, real scalar Jy)
{
    real scalar ntheta, ndelta, p, r, k, l, yidx, oldidx
    real colvector keep_theta, map_theta
    real matrix A

    ntheta = K*(K+1)/2
    ndelta = K*Jy
    p = ntheta + ndelta

    keep_theta = J(K*K,1,0)
    for (k=1; k<=K; k++) {
        for (l=1; l<=K; l++) {
            oldidx = (k-1)*K + l
            if (l<=k) keep_theta[oldidx] = 1
        }
    }
    map_theta = J(K*K,1,0)
    r = 0
    for (oldidx=1; oldidx<=K*K; oldidx++) {
        if (keep_theta[oldidx]) {
            r = r + 1
            map_theta[oldidx] = r
        }
    }

    A = J(4*K + K*Jy, p, 0)

    // +P(M=k|D=0): indices k + K*(0:(K-1))
    for (k=1; k<=K; k++) {
        for (l=1; l<=K; l++) {
            oldidx = k + K*(l-1)
            if (map_theta[oldidx] > 0) A[k, map_theta[oldidx]] = 1
        }
    }

    // +P(M=k|D=1): indices ((k-1)*K+1):(k*K)
    for (k=1; k<=K; k++) {
        for (l=1; l<=K; l++) {
            oldidx = (k-1)*K + l
            if (map_theta[oldidx] > 0) A[K+k, map_theta[oldidx]] = 1
        }
    }

    A[(2*K+1)..(3*K),.] = -A[1..K,.]
    A[(3*K+1)..(4*K),.] = -A[(K+1)..(2*K),.]

    // p_ym(1)-p_ym(0) <= delta
    r = 4*K
    for (k=1; k<=K; k++) {
        for (yidx=1; yidx<=Jy; yidx++) {
            r = r + 1
            A[r, ntheta + (k-1)*Jy + yidx] = 1
        }
    }

    return(A)
}

real matrix testmechs__build_Ashp(real scalar K, real scalar Jy)
{
    real scalar ntheta, ndelta, p, r, k, l, yidx, oldidx
    real colvector keep_theta, map_theta
    real matrix A

    ntheta = K*(K+1)/2
    ndelta = K*Jy
    p = ntheta + ndelta

    keep_theta = J(K*K,1,0)
    for (k=1; k<=K; k++) {
        for (l=1; l<=K; l++) {
            oldidx = (k-1)*K + l
            if (l<=k) keep_theta[oldidx] = 1
        }
    }
    map_theta = J(K*K,1,0)
    r = 0
    for (oldidx=1; oldidx<=K*K; oldidx++) {
        if (keep_theta[oldidx]) {
            r = r + 1
            map_theta[oldidx] = r
        }
    }

    A = J(K + p, p, 0)

    // For each k: sum_{l != k} theta_lk - sum_y delta_yk <= 0
    for (k=1; k<=K; k++) {
        for (l=1; l<=K; l++) {
            if (l==k) continue
            oldidx = (k-1)*K + l
            if (map_theta[oldidx] > 0) A[k, map_theta[oldidx]] = 1
        }
        for (yidx=1; yidx<=Jy; yidx++) {
            A[k, ntheta + (k-1)*Jy + yidx] = -1
        }
    }

    // nonnegativity
    A[(K+1)..(K+p),.] = I(p)

    return(A)
}

real colvector testmechs__qp_active_set(real matrix H, real colvector f, real matrix G, real colvector h, real scalar tol)
{
    real scalar p, m, iter, maxiter, dropi, worst, j, found
    real colvector x, v, lambda, active, ia, rhs, sol, cand
    real matrix KKT, GW

    p = cols(H)
    m = rows(G)
    maxiter = 500

    H = H + 1e-8*I(p)
    x = cholsolve(cholesky(H), f)
    active = J(m,1,0)

    for (iter=1; iter<=maxiter; iter++) {
        v = G*x - h
        worst = max(v)

        if (worst <= tol) {
            ia = selectindex(active)
            if (rows(ia)==0) return(x)

            GW = G[ia,.]
            if (rank(GW) < rows(GW)) {
                _error(498, "Active set became linearly dependent before multiplier step")
            }

            KKT = (H, GW' \ GW, J(rows(ia), rows(ia), 0))
            rhs = (f \ h[ia])
            sol = lusolve(KKT, rhs)
            if (any(sol :>= .)) {
                _error(498, "KKT solve failed in multiplier step")
            }

            lambda = sol[(p+1)..rows(rhs)]

            if (min(lambda) >= -tol) return(x)

            dropi = ia[selectindex(lambda :== min(lambda))[1]]
            active[dropi] = 0
            continue
        }

        cand = selectindex((v :== worst) :& (active :== 0))
        if (rows(cand)==0) cand = selectindex((v :> tol) :& (active :== 0))
        if (rows(cand)==0) return(x)

        found = 0
     for (j=1; j<=rows(cand); j++) {
         active[cand[j]] = 1
         ia = selectindex(active)
         GW = G[ia,.]
 
         if (rank(GW) == rows(GW)) {
             KKT = (H, GW' \ GW, J(rows(ia), rows(ia), 0))
             rhs = (f \ h[ia])
             sol = lusolve(KKT, rhs)
 
             if (all(sol :< .)) {
                 found = 1
                 break
             }
         }
 
         active[cand[j]] = 0
     }
 
     if (!found) {
         _error(498, "No violated constraint can be added without making KKT numerically singular")
     }
 
     x = sol[1..p]
    }

    _error(498, "Active-set QP did not converge within maxiter")
    return(x)
}

void testmechs__sharpnull_cs(string scalar dvar, string scalar mvar, string scalar yvar,
    string scalar clvar, string scalar tousevar, real scalar alpha, real scalar numybins, real scalar new_dof_CS)
{
    real colvector d, m, y, cl, yvals, mvals, beta_obs, beta_shp, beta, d_Z
    real colvector p_m0, p_m1, p_ym0, p_ym1, beta_red, xhat, lvec, uvec, dvec, khat, clu, hqp
    real colvector if0, if1, ifm0, ifm1
    real scalar n, k, jy, i, ii, g, pd0, pd1, test_stat, dof_n, cv, pval, r, q, d_nuis, tol
    real matrix ifs, ifs_cl, sigma_obs, A_obs, A_shp, A, sigma, eval, evec, B_Z, C_Z
    real matrix sigmaInv, Amat, Amat_aug, Dmat, Gqp

    real colvector keep
    keep = selectindex(st_data(., tousevar) :!= 0)
    if (rows(keep) == 0) _error(2000, "No observations in estimation sample")

    d  = vec(st_data(keep, dvar))
    m  = vec(st_data(keep, mvar))
    y  = vec(st_data(keep, yvar))
    cl = vec(st_data(keep, clvar))

    keep = selectindex((d :< .) :& (m :< .) :& (y :< .) :& (cl :< .))
    d = d[keep]; m = m[keep]; y = y[keep]; cl = cl[keep]
    if (rows(d) == 0) _error(2000, "No complete-case observations after filtering")

    y = testmechs__discretize_y(y, numybins)

    yvals = uniqrows(sort(y,1))
    mvals = uniqrows(sort(m,1))
    jy = rows(yvals)
    k = rows(mvals)

    n = rows(y)
    pd0 = mean(d :== 0)
    pd1 = mean(d :== 1)

    ifs = J(n, 4*k + k*jy, 0)
    p_ym0 = J(k*jy,1,.)
    p_ym1 = J(k*jy,1,.)
    p_m0  = J(k,1,.)
    p_m1  = J(k,1,.)

    // my_values ordering: arrange(m, y)
    ii = 0
    for (i=1; i<=k; i++) {
        ifm0 = J(n,1,0)
        ifm1 = J(n,1,0)
        for (r=1; r<=jy; r++) {
            ii = ii + 1
            p_ym0[ii] = mean((y :== yvals[r]) :& (m :== mvals[i]) :& (d :== 0)) / pd0
            p_ym1[ii] = mean((y :== yvals[r]) :& (m :== mvals[i]) :& (d :== 1)) / pd1

            if0 = (d :== 0) :* (((y :== yvals[r]) :& (m :== mvals[i]) :& (d :== 0)) :- p_ym0[ii]) / pd0
            if1 = (d :== 1) :* (((y :== yvals[r]) :& (m :== mvals[i]) :& (d :== 1)) :- p_ym1[ii]) / pd1

            ifs[,4*k + ii] = if1 - if0
            ifm0 = ifm0 + if0
            ifm1 = ifm1 + if1
        }
        p_m0[i] = sum(p_ym0[((i-1)*jy+1)..(i*jy)])
        p_m1[i] = sum(p_ym1[((i-1)*jy+1)..(i*jy)])

        ifs[, i] = ifm0
        ifs[, k+i] = ifm1
        ifs[, 2*k+i] = -ifm0
        ifs[, 3*k+i] = -ifm1
    }

    beta_obs = p_m0 \ p_m1 \ (-p_m0) \ (-p_m1) \ (p_ym1 - p_ym0)

    clu = uniqrows(sort(cl,1))
    g = rows(clu)
    ifs_cl = J(g, cols(ifs), 0)
    for (i=1; i<=g; i++) {
        ifs_cl[i,.] = colsum(ifs[selectindex(cl :== clu[i]), .])
    }
    sigma_obs = (g > 1 ? g/(g-1) : 1) * (ifs_cl' * ifs_cl) / (n^2)

    A_obs = testmechs__build_Aobs(k, jy)
    A_shp = testmechs__build_Ashp(k, jy)
    beta_shp = J(rows(A_shp),1,0)

    beta = beta_obs \ beta_shp
    A = A_obs \ A_shp
    sigma = (sigma_obs, J(rows(sigma_obs), rows(A_shp), 0) \
             J(rows(A_shp), cols(sigma_obs), 0), J(rows(A_shp), rows(A_shp), 0))

    d_Z = J(rows(A),1,0)
    C_Z = A

        symeigensystem((sigma + sigma')/2, evec, eval)
    eval = Re(eval)

    if (min(eval) < 1e-8) {
        keep = selectindex(eval :> 1e-8)
        if (rows(keep) == 0) _error(498, "No positive eigenvalues in sigma")

        B_Z = evec[, keep]
        beta_red = B_Z' * beta
        d_Z = d_Z + B_Z * beta_red - beta
        sigma = diag(eval[keep])
    }
    else {
        B_Z = I(rows(beta))
        beta_red = beta
    }

    sigmaInv = invsym(sigma)
    r = rows(sigma)
    q = rows(C_Z)
    d_nuis = cols(C_Z)

    Amat = (B_Z, -C_Z \ J(d_nuis, cols(B_Z), 0), I(d_nuis))
    uvec = d_Z \ J(d_nuis,1,1)
    lvec = J(rows(d_Z),1,-1e20) \ J(d_nuis,1,0)

    Dmat = J(r + d_nuis, r + d_nuis, 0)
    Dmat[1..r,1..r] = 2*sigmaInv
    dvec = (2*sigmaInv*beta_red) \ J(d_nuis,1,0)

    // convert box constraints l <= Amat x <= u into Gx <= h
    Gqp = Amat
    hqp = uvec
    for (i=1; i<=rows(lvec); i++) {
        if (lvec[i] > -1e19) {
            Gqp = Gqp \ (-Amat[i,.])
            hqp = hqp \ (-lvec[i])
        }
    }

    tol = 1e-8
    xhat = testmechs__qp_active_set(Dmat, dvec, Gqp, hqp, tol)
	if (any(xhat :>= .)) _error(498, "xhat contains missing values")

    test_stat = ((beta_red - xhat[1..r])' * sigmaInv * (beta_red - xhat[1..r]))[1,1]

    if (new_dof_CS) {
        Amat_aug = Amat \ -(J(d_nuis, cols(B_Z), 0), I(d_nuis))
        uvec = uvec \ J(d_nuis,1,0)
        khat = selectindex(abs(uvec - Amat_aug*xhat) :< sqrt(tol))
        if (rows(khat)==0) {
            dof_n = 0
        }
        else {
            dof_n = rank(Amat_aug[khat,.]) - rank(Amat_aug[khat,(cols(B_Z)+1)..cols(Amat_aug)])
            if (dof_n < 0) dof_n = 0
        }
    }
    else {
        dof_n = 0
    }
    
	if (dof_n <= 0) {
    cv = 0
    pval = (test_stat <= 0 ? 1 : 0)
}
else {
    cv = invchi2(dof_n, 1-alpha)
    pval = chi2tail(dof_n, test_stat)
}

    st_numscalar("__tm_pval", pval)
    st_numscalar("__tm_test_stat", test_stat)
    st_numscalar("__tm_cv", cv)
}
end
