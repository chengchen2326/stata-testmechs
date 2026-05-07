program _testmechs_dqrdc2_rank, plugin
program define testmechs_test_sharpnull, rclass
    // Syntax: d m1 [m2 m3 ...] y , options
    // First word = treatment (d), last word = outcome (y), middle words = mediator(s)
    version 16.0
    syntax varlist(min=3 max=100 numeric) [if] [in] , METHOD(string) ///
    [NUMYBINS(integer 5) CLUSTER(varname) MAXDEFIERSSHARE(real 0) NEWDOFCS]

    if ("`method'" != "CS") {
        di as err "Only method(CS) is supported in this MVP implementation"
        exit 198
    }

    if (`maxdefiersshare' < 0 | `maxdefiersshare' > 1) {
        di as err "maxdefiersshare() must be between 0 and 1"
        exit 198
    }
	
	// Python availability check (only needed when newdofcs is NOT specified)
    // Checks each dependency independently and gives a specific install hint.
    if ("`newdofcs'" == "") {
        // Step 1: Stata's Python integration must be configured
        capture python query
        if (_rc) {
            di as err "testmechs_test_sharpnull requires Stata 16+ with Python integration."
            di as err "Configure Python in Stata with:  python set exec /path/to/python3"
            di as err "See {help python} for details."
            exit 199
        }

        // Step 2: numpy (used by both the LP and QP wrappers)
        capture python: import numpy
        if (_rc) {
            di as err "testmechs_test_sharpnull requires the Python package 'numpy'."
            di as err "Install it with:  pip install numpy"
            exit 199
        }

        // Step 3: swiglpk (GLPK bindings; LP step uses the same solver as R)
        capture python: import swiglpk
        if (_rc) {
            di as err "testmechs_test_sharpnull requires the Python package 'swiglpk'."
            di as err "Install it with:  pip install swiglpk"
            di as err "On macOS you may also need:  brew install glpk"
            di as err "On Debian/Ubuntu:           sudo apt install libglpk-dev"
            exit 199
        }

        // Step 4: osqp 0.6.x (used by the QP step)
        capture python: import osqp
        if (_rc) {
            di as err "testmechs_test_sharpnull requires the Python package 'osqp' (version 0.6.x)."
            di as err "Install it with:  pip install 'osqp<1'"
            di as err "(osqp version 1.x is NOT compatible)"
            exit 199
        }
    }

    marksample touse

    // Parse varlist: first word = d, last word = y, middle word(s) = mediators
    local nwords : word count `varlist'
    local d : word 1 of `varlist'
    local y : word `nwords' of `varlist'
    local mvars ""
    forvalues wi = 2/`= `nwords' - 1' {
        local mvars `mvars' `: word `wi' of `varlist''
    }
    local mvars = strtrim("`mvars'")

    tempvar touse2 clusterid
    quietly gen byte `touse2' = `touse'
    quietly replace `touse2' = 0 if missing(`d', `y')
    foreach mv of local mvars {
        quietly replace `touse2' = 0 if missing(`mv')
    }
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
    
    // Pass space-separated mediator names to Mata; Mata builds uni_m and PO internally
	local ndofcs = cond("`newdofcs'" == "", 0, 1)
    mata: testmechs__sharpnull_cs("`d'","`mvars'","`y'","`clusterid'","`touse2'",0.05,`numybins',`ndofcs',`maxdefiersshare')

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

// Build A.obs matrix.
//
// PO is a K x K indicator matrix where PO[l,k] = 1 iff mediator-tuple label l
// is componentwise <= label k.  For a univariate mediator with labels sorted
// ascending, PO[l,k] = (l <= k), recovering the original total-order behavior.
//
// theta[l,k] (prob that M(0)=l given M(1)=k) is stored at linear index
// (k-1)*K + l (column-major).  We retain theta[l,k] iff PO[l,k] == 1.
real matrix testmechs__build_Aobs(real scalar K, real scalar Jy,
                                   real scalar max_defiers_share, real matrix PO)
{
    real scalar ntheta, ndelta, p, r, k, l, yidx, oldidx
    real colvector keep_theta, map_theta
    real matrix A

    // keep_theta[(k-1)*K+l] = 1 iff theta[l,k] is retained under partial order
    keep_theta = J(K*K, 1, (max_defiers_share==0 ? 0 : 1))
    if (max_defiers_share==0) {
        for (k=1; k<=K; k++) {
            for (l=1; l<=K; l++) {
                oldidx = (k-1)*K + l
                // Replace hard-coded l<=k with the partial-order condition PO[l,k]
                if (PO[l,k]) keep_theta[oldidx] = 1
            }
        }
    }

    ntheta = sum(keep_theta)   // = sum(PO) when max_defiers_share==0, else K^2
    ndelta = K*Jy
    p = ntheta + ndelta

    map_theta = J(K*K, 1, 0)
    r = 0
    for (oldidx=1; oldidx<=K*K; oldidx++) {
        if (keep_theta[oldidx]) {
            r = r + 1
            map_theta[oldidx] = r
        }
    }

    A = J(4*K + K*Jy, p, 0)

    // P(M=k|D=0) = sum_l theta[k,l]; index of theta[k,l] = (l-1)*K + k
    for (k=1; k<=K; k++) {
        for (l=1; l<=K; l++) {
            oldidx = k + K*(l-1)
            if (map_theta[oldidx] > 0) A[k, map_theta[oldidx]] = 1
        }
    }

    // P(M=k|D=1) = sum_l theta[l,k]; index of theta[l,k] = (k-1)*K + l
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

// Build A.shp matrix.
//
// Uses the same PO partial-order matrix as testmechs__build_Aobs.
real matrix testmechs__build_Ashp(real scalar K, real scalar Jy,
                                   real scalar max_defiers_share, real matrix PO)
{
    real scalar ntheta, ndelta, p, r, k, l, yidx, oldidx, def_row
    real colvector keep_theta, map_theta
    real matrix A

    keep_theta = J(K*K, 1, (max_defiers_share==0 ? 0 : 1))
    if (max_defiers_share==0) {
        for (k=1; k<=K; k++) {
            for (l=1; l<=K; l++) {
                oldidx = (k-1)*K + l
                if (PO[l,k]) keep_theta[oldidx] = 1
            }
        }
    }

    ntheta = sum(keep_theta)
    ndelta = K*Jy
    p = ntheta + ndelta

    map_theta = J(K*K, 1, 0)
    r = 0
    for (oldidx=1; oldidx<=K*K; oldidx++) {
        if (keep_theta[oldidx]) {
            r = r + 1
            map_theta[oldidx] = r
        }
    }

    A = J(K + (max_defiers_share!=0) + p, p, 0)

    // Row k: sum_{l != k, PO[l,k]==1} theta[l,k] - sum_y delta[y,k] <= 0
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

    // Defiers constraint: sum of theta[l,k] where l is NOT <= k in partial order
    if (max_defiers_share!=0) {
        def_row = K + 1
        for (k=1; k<=K; k++) {
            for (l=1; l<=K; l++) {
                oldidx = (k-1)*K + l
                // Replace l>k (total order) with !PO[l,k] (partial order)
                if (!PO[l,k] & map_theta[oldidx] > 0) A[def_row, map_theta[oldidx]] = -1
            }
        }
    }

    // nonnegativity
    A[(K + (max_defiers_share!=0) + 1)..(K + (max_defiers_share!=0) + p),.] = I(p)

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

// -----------------------------------------------------------------------------
// Helpers for R's default (new_dof_CS = FALSE) algorithm.
// -----------------------------------------------------------------------------

// Rank of A with R-style relative tolerance (replaces R's qr(M, tol=qr_tol)$rank).
// Uses SVD: count singular values exceeding tol * sigma_max.
real scalar testmechs__rank_tol(real matrix A, real scalar tol)
{
    real scalar k
    
    if (rows(A) == 0 | cols(A) == 0) return(0)
    
    // Push matrix to Stata, call plugin, read result back
    st_matrix("__tm_rank_M", A)
    stata("scalar __tm_rank_result = .")
    stata("plugin call _testmechs_dqrdc2_rank, __tm_rank_M " + char(34) + strofreal(tol, "%21.15e") + char(34) + " __tm_rank_result")
    k = st_numscalar("__tm_rank_result")
    
    return(k)
}

// Rank-revealing column pivot via pivoted modified Gram-Schmidt.
// Mirrors R's qr(M, tol=qr_tol)$pivot[1:target_rank]: picks the first
// `target_rank` linearly-independent columns of M (relative tol = tol).
// NOTE: R uses LAPACK QR-with-column-pivoting; this is a simple greedy
// equivalent, sufficient for our use (only for producing G1_ind).
// QR column pivoting: mimics R's qr(M, tol=qr_tol)$pivot[1:target_rank].
// At each step, selects the remaining column with the largest residual norm
// (orthogonal to already-chosen columns).  This matches LAPACK DGEQP3.
real colvector testmechs__rank_revealing_pivot(real matrix M, real scalar tol,
                                                real scalar target_rank)
{
    real colvector pivot, colnorms, sel_idx, r_vec
    real matrix Q, R_mat, candidates
    real scalar j, k, n, m, ref_norm, best_j, best_norm, cand_norm, nv, cur_idx
    real colvector v, chosen_mask
    
    n = rows(M)
    m = cols(M)
    if (target_rank <= 0 | n == 0 | m == 0) return(J(0, 1, 0))
    
    // Initial column norms of M
    colnorms = J(m, 1, 0)
    for (j = 1; j <= m; j++) colnorms[j] = sqrt(M[., j]' * M[., j])
    ref_norm = max(colnorms)
    if (ref_norm <= 0) return(J(0, 1, 0))
    
    pivot = J(target_rank, 1, 0)
    Q = J(n, 0, 0)
    chosen_mask = J(m, 1, 0)
    k = 0
    
    while (k < target_rank) {
        // Find the remaining column with largest residual norm
        best_j = 0
        best_norm = -1
        for (j = 1; j <= m; j++) {
            if (chosen_mask[j]) continue
            v = M[., j]
            if (cols(Q) > 0) {
                v = v - Q * (Q' * v)
                // re-orthogonalize for stability
                v = v - Q * (Q' * v)
            }
            cand_norm = sqrt(v' * v)
            if (cand_norm > best_norm) {
                best_norm = cand_norm
                best_j = j
            }
        }
        
        // Check if best column is still numerically independent
        if (best_norm <= tol * ref_norm | best_j == 0) break
        
        // Accept this column
        k = k + 1
        pivot[k] = best_j
        chosen_mask[best_j] = 1
        
        // Orthogonalize and add to Q
        v = M[., best_j]
        if (cols(Q) > 0) {
            v = v - Q * (Q' * v)
            v = v - Q * (Q' * v)
        }
        nv = sqrt(v' * v)
        Q = Q, (v / nv)
    }
    
    if (k < target_rank) pivot = pivot[1..k]
    return(pivot)
}

// Solve min c'x s.t. A_eq x = b_eq, x >= 0  via scipy.optimize.linprog (HiGHS).
// Pushes data through Stata matrices, calls the `_testmechs_lp_call` ado
// subroutine (which runs a python: block), then reads scalars back.
// Result objective -> Stata scalar __tm_lp_fun; success flag -> __tm_lp_ok.
void testmechs__lp_solve(real colvector c, real matrix A_eq, real colvector b_eq)
{
    st_matrix("__tm_c",   c)
    st_matrix("__tm_Aeq", A_eq)
    st_matrix("__tm_beq", b_eq)
    stata("python script " + char(34) + "`c(sysdir_plus)'py/_testmechs_lp_call.py" + char(34))
}

// Solve QP min 0.5 x' P x + q' x s.t. l <= A x <= u via OSQP 0.6.x.
// Pushes data through Stata matrices, calls _testmechs_qp_call.py.
// Returns the solution vector (p x 1).  Aborts with _error on failure.
real colvector testmechs__qp_solve(real matrix P, real colvector q, real matrix A,
                                    real colvector l, real colvector u)
{
    real colvector xhat
    real scalar ok
    
    st_matrix("__tm_qp_P", P)
    st_matrix("__tm_qp_q", q)
    st_matrix("__tm_qp_A", A)
    st_matrix("__tm_qp_l", l)
    st_matrix("__tm_qp_u", u)
    
    stata("python script " + char(34) + "`c(sysdir_plus)'py/_testmechs_qp_call.py" + char(34))
    
    ok = st_numscalar("__tm_qp_ok")
    if (ok != 1) _error(498, "OSQP failed to solve the QP")
    
    xhat = st_matrix("__tm_qp_x")
    return(xhat)
}



// Main CS test function.
//
// mvars_str is a space-separated string of mediator variable names (one or more).
// When J==1 this reduces exactly to the previous single-mediator behavior because
// uniqrows() sorts the single column ascending, giving PO[l,k] = (l<=k).
void testmechs__sharpnull_cs(string scalar dvar, string scalar mvars_str, string scalar yvar,
    string scalar clvar, string scalar tousevar, real scalar alpha, real scalar numybins,
    real scalar newdofcs, real scalar max_defiers_share)
{
    real colvector d, m, y, cl, yvals, mvals, beta_obs, beta_shp, beta, d_Z
    real colvector p_m0, p_m1, p_ym0, p_ym1, beta_red, xhat, lvec, uvec, dvec, khat, clu, hqp
    real colvector if0, if1, ifm0, ifm1
    real scalar n, k, jy, i, ii, g, pd0, pd1, test_stat, dof_n, cv, pval, r, q, d_nuis, tol
    real matrix ifs, ifs_cl, sigma_obs, A_obs, A_shp, A, sigma, eval, evec, B_Z, C_Z
    real matrix sigmaInv, Amat, Amat_aug, Dmat, Gqp

    // For multivariate mediator support
    string rowvector mvars_tok
    real matrix M_mat, m_supp, PO
    real scalar J_med, ki, li, ji
    real colvector uni_m, match_vec, complete_mask, idxs

    real colvector keep
    keep = selectindex(st_data(., tousevar) :!= 0)
    if (rows(keep) == 0) _error(2000, "No observations in estimation sample")

    // Parse mediator variable names
    mvars_tok = tokens(mvars_str)
    J_med = cols(mvars_tok)

    d     = vec(st_data(keep, dvar))
    y     = vec(st_data(keep, yvar))
    cl    = vec(st_data(keep, clvar))
    M_mat = st_data(keep, mvars_tok)    // n x J_med matrix of raw mediator values

    // Filter to complete cases across all variables
    complete_mask = (d :< .) :& (y :< .) :& (cl :< .)
    for (ji=1; ji<=J_med; ji++) {
        complete_mask = complete_mask :& (M_mat[,ji] :< .)
    }
    keep = selectindex(complete_mask)
    d = d[keep]; y = y[keep]; cl = cl[keep]; M_mat = M_mat[keep,]
    if (rows(d) == 0) _error(2000, "No complete-case observations after filtering")

    y = testmechs__discretize_y(y, numybins)

    // --- Build univariate integer label uni_m and partial-order matrix PO ---
    //
    // Unique mediator tuples are sorted lexicographically by uniqrows() (ascending
    // by column 1, then column 2, ...).  This gives a deterministic labeling that
    // is invariant to the order of observations in the dataset.
    //
    // For J==1, m_supp is just the sorted unique values, so PO[l,k] = (l<=k),
    // which is the standard total order -- identical to the previous behavior.
    //
    // TODO: if the numerical result for the multivariate case is off by more than
    // 1e-3 from the R benchmark, try first-occurrence labeling (R's unique() order)
    // instead of lexicographic.  The pval should be invariant to relabeling in
    // theory, but numerical coincidences may matter.
    m_supp = uniqrows(M_mat)    // K x J_med, lexicographically sorted
    k = rows(m_supp)

    // Label each observation with its tuple's integer index (1..K)
    uni_m = J(rows(d), 1, 0)
    for (ki=1; ki<=k; ki++) {
        match_vec = J(rows(d), 1, 1)
        for (ji=1; ji<=J_med; ji++) {
            match_vec = match_vec :& (M_mat[,ji] :== m_supp[ki,ji])
        }
        idxs = selectindex(match_vec)
        if (rows(idxs) > 0) uni_m[idxs] = J(rows(idxs), 1, ki)
    }

    // PO[l,k] = 1 iff every component of tuple l is <= the corresponding component of tuple k
    PO = J(k, k, 0)
    for (li=1; li<=k; li++) {
        for (ki=1; ki<=k; ki++) {
            // m_supp[li,]' and m_supp[ki,]' are J_med x 1 column vectors
            if (all(m_supp[li,]' :<= m_supp[ki,]')) PO[li,ki] = 1
        }
    }

    // From here on, use integer-labeled univariate m = uni_m
    m     = uni_m
    mvals = (1::k)

    yvals = uniqrows(sort(y,1))
    jy    = rows(yvals)

    n   = rows(y)
    pd0 = mean(d :== 0)
    pd1 = mean(d :== 1)

    ifs   = J(n, 4*k + k*jy, 0)
    p_ym0 = J(k*jy, 1, .)
    p_ym1 = J(k*jy, 1, .)
    p_m0  = J(k, 1, .)
    p_m1  = J(k, 1, .)

    // Build IFs; my_values ordering: arrange(m, y) i.e. outer loop over m, inner over y
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

        ifs[, i]       = ifm0
        ifs[, k+i]     = ifm1
        ifs[, 2*k+i]   = -ifm0
        ifs[, 3*k+i]   = -ifm1
    }

    beta_obs = p_m0 \ p_m1 \ (-p_m0) \ (-p_m1) \ (p_ym1 - p_ym0)

    clu    = uniqrows(sort(cl,1))
    g      = rows(clu)
    ifs_cl = J(g, cols(ifs), 0)
    for (i=1; i<=g; i++) {
        ifs_cl[i,.] = colsum(ifs[selectindex(cl :== clu[i]), .])
    }
    sigma_obs = (g > 1 ? g/(g-1) : 1) * (ifs_cl' * ifs_cl) / (n^2)

	
	


    // Pass PO to both builders so they can apply the correct partial order
    A_obs    = testmechs__build_Aobs(k, jy, max_defiers_share, PO)
    A_shp    = testmechs__build_Ashp(k, jy, max_defiers_share, PO)
    beta_shp = J(rows(A_shp), 1, 0)
    if (max_defiers_share != 0) beta_shp[k+1] = -max_defiers_share

    beta  = beta_obs \ beta_shp
    A     = A_obs \ A_shp
    sigma = (sigma_obs, J(rows(sigma_obs), rows(A_shp), 0) \
             J(rows(A_shp), cols(sigma_obs), 0), J(rows(A_shp), rows(A_shp), 0))

    d_Z = J(rows(A), 1, 0)
    C_Z = A

    symeigensystem((sigma + sigma')/2, evec, eval)
    eval = Re(eval)

    if (min(eval) < 1e-8) {
        keep = selectindex(eval :> 1e-8)
        if (rows(keep) == 0) _error(498, "No positive eigenvalues in sigma")

        B_Z      = evec[, keep]
        beta_red = B_Z' * beta
        d_Z      = d_Z + B_Z * beta_red - beta
        sigma    = diag(eval[keep])
    }
    else {
        B_Z      = I(rows(beta))
        beta_red = beta
    }

    sigmaInv = invsym(sigma)
    r        = rows(sigma)
    q        = rows(C_Z)
    d_nuis   = cols(C_Z)

    Amat = (B_Z, -C_Z \ J(d_nuis, cols(B_Z), 0), I(d_nuis))
    uvec = d_Z \ J(d_nuis, 1, 1)
    lvec = J(rows(d_Z), 1, -1e20) \ J(d_nuis, 1, 0)

    Dmat = J(r + d_nuis, r + d_nuis, 0)
    Dmat[1..r,1..r] = 2*sigmaInv
    dvec = (2*sigmaInv*beta_red) \ J(d_nuis, 1, 0)

    // Solve the QP via OSQP (match R TestMechs which uses osqp 0.6.x):
    //   min  0.5 x' Dmat x + dvec' x
    //   s.t. lvec <= Amat x <= uvec
    // Note: Dmat here includes the factor of 2 (Dmat = 2*sigmaInv block); dvec is
    // the same sign as in R (not negated).  R's OSQP call uses P=Dmat, q=-dvec.
    xhat = testmechs__qp_solve(Dmat, -dvec, Amat, lvec, uvec)
    if (any(xhat :>= .)) _error(498, "xhat contains missing values")

    test_stat = ((beta_red - xhat[1..r])' * sigmaInv * (beta_red - xhat[1..r]))[1,1]

    if (newdofcs) {
        Amat_aug = Amat \ -(J(d_nuis, cols(B_Z), 0), I(d_nuis))
        uvec     = uvec \ J(d_nuis, 1, 0)
        khat     = selectindex(abs(uvec - Amat_aug*xhat) :< sqrt(tol))
        if (rows(khat)==0) {
            dof_n = 0
        }
        else {
            dof_n = rank(Amat_aug[khat,.]) - rank(Amat_aug[khat,(cols(B_Z)+1)..cols(Amat_aug)])
            if (dof_n < 0) dof_n = 0
        }
    }
    else {
        // R's default new_dof_CS = FALSE branch (test_sharp_null.R L468-696):
        // Solve LPs via scipy.linprog to detect implicit equality constraints,
        // then compute dof_n from a QR rank analysis.
        real scalar j_lp, V_min, flag_ok, qr_tol, psi_tol, rkA, rank_BZ
        real colvector c_lp, beq_lp, psis, beta_red_star, bz_bstar, diff_vec
        real colvector impeq_mask, G1_ind, G2_ind, is_in_G1
        real matrix I_dineq, Aeq_lp, A_impeq, A_full_eq, G, G1, G2, B_Z1, B_Z2, Gamma

        I_dineq       = I(q)
        beta_red_star = xhat[1..r]
        bz_bstar      = B_Z * beta_red_star
        diff_vec      = d_Z - bz_bstar

        // --- 3.1 First LP: feasibility check (R L493-520) -----------------
        //   min (d_Z - B_Z beta_red_star)' x
        //   s.t. C_Z' x = 0, sum(x) = 1, x >= 0
        c_lp   = diff_vec
        Aeq_lp = C_Z' \ J(1, q, 1)
        beq_lp = J(d_nuis, 1, 0) \ 1
        testmechs__lp_solve(c_lp, Aeq_lp, beq_lp)
        V_min   = st_numscalar("__tm_lp_fun")
        flag_ok = st_numscalar("__tm_lp_ok")

        if (flag_ok != 1 | V_min >= 0.00005) {
            // Feasibility LP infeasible/non-optimal or V_min not ~ 0: dof_n = 0
            dof_n = 0
        }
        else {
            // --- 3.2 K+1 LPs: find implicit equalities (R L522-549) ------
            //   For each j = 1..q:
            //     min -I_q[j,:] x
            //     s.t. C_Z' x = 0, diff_vec' x = V_min, sum(x) = 1, x >= 0
            psis   = J(q, 1, .)
            Aeq_lp = C_Z' \ diff_vec' \ J(1, q, 1)
            beq_lp = J(d_nuis, 1, 0) \ V_min \ 1
            for (j_lp = 1; j_lp <= q; j_lp++) {
                c_lp = -I_dineq[j_lp, .]'
                testmechs__lp_solve(c_lp, Aeq_lp, beq_lp)
                if (st_numscalar("__tm_lp_ok") == 1) {
                    psis[j_lp] = st_numscalar("__tm_lp_fun")
                }
                else {
                    psis[j_lp] = -1
                }
            }

            // --- 3.3 Compute dof_n from QR rank analysis (R L623-685) ----
            psi_tol = 1e-6
            qr_tol  = 1e-5


            // A_impeq = rows of I_q where psis[j] > -psi_tol
            impeq_mask = (psis :> -psi_tol)
            if (sum(impeq_mask) > 0) {
                A_impeq   = I_dineq[selectindex(impeq_mask), .]
                A_full_eq = A_impeq \ -(C_Z') \ (bz_bstar - d_Z)'
            }
            else {
                A_full_eq = -(C_Z') \ (bz_bstar - d_Z)'
            }

            // rkA ~ rank(A_full_eq' A_full_eq, tol = qr_tol)
            rkA     = testmechs__rank_tol(A_full_eq' * A_full_eq, qr_tol)
            rank_BZ = testmechs__rank_tol(B_Z, qr_tol)

            if (rank_BZ == q) {
                // B_Z full-rank branch (R L633-634)
                dof_n = q - rkA
            }
            else {
                // Rank-deficient B_Z branch (R L636-685)
                G = A_full_eq'
                if (rkA >= rows(G)) {
                    dof_n = 0
                }
                else {
                    // QR pivot on t(G) == A_full_eq: pick rkA linearly independent columns
                    G1_ind = testmechs__rank_revealing_pivot(A_full_eq, qr_tol, rkA)
                    is_in_G1 = J(rows(G), 1, 0)
                    for (j_lp = 1; j_lp <= rows(G1_ind); j_lp++) {
                        is_in_G1[G1_ind[j_lp]] = 1
                    }
                    G2_ind = selectindex(is_in_G1 :== 0)

                    G1   = G[G1_ind, .]
                    G2   = G[G2_ind, .]
                    B_Z1 = B_Z[G1_ind, .]
                    B_Z2 = B_Z[G2_ind, .]

				
                    // Gamma solves G1' Gamma = -G2' (R: qr.coef(qr(t(G1)), t(G2)))
                    Gamma = -qrsolve(G1', G2')
                    dof_n = testmechs__rank_tol(Gamma' * B_Z1 + B_Z2, qr_tol)
					
                }
            }
            if (dof_n < 0) dof_n = 0
        }
    }

    if (dof_n <= 0) {
        cv   = 0
        pval = (test_stat <= 0 ? 1 : 0)
    }
    else {
        cv   = invchi2(dof_n, 1-alpha)
        pval = chi2tail(dof_n, test_stat)
    }

    st_numscalar("__tm_pval",      pval)
    st_numscalar("__tm_test_stat", test_stat)
    st_numscalar("__tm_cv",        cv)
}
end

// -----------------------------------------------------------------------------
// Ado subroutine called by testmechs__lp_solve() via stata("_testmechs_lp_call").
// Reads Stata matrices __tm_c, __tm_Aeq, __tm_beq, runs scipy.optimize.linprog
// (HiGHS), and writes Stata scalars __tm_lp_fun and __tm_lp_ok.
// -----------------------------------------------------------------------------
program define _testmechs_lp_call
    version 16.0
    python:
from sfi import Matrix, Scalar
import numpy as np
from scipy.optimize import linprog

_c   = np.asarray(Matrix.get("__tm_c"),   dtype=float).flatten()
_Aeq = np.asarray(Matrix.get("__tm_Aeq"), dtype=float)
_beq = np.asarray(Matrix.get("__tm_beq"), dtype=float).flatten()
if _Aeq.ndim == 1:
    _Aeq = _Aeq.reshape(1, -1)

_bounds = [(0, None)] * _c.size
_res = linprog(c=_c, A_eq=_Aeq, b_eq=_beq, bounds=_bounds, method="highs")

Scalar.setValue("__tm_lp_fun", float(_res.fun) if _res.fun is not None else float("nan"))
Scalar.setValue("__tm_lp_ok",  1.0 if _res.success else 0.0)
end
