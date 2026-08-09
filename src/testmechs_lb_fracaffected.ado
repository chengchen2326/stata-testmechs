capture program drop _testmechs_glpk_lp
capture program _testmechs_glpk_lp, plugin

program define testmechs_lb_fracaffected, rclass
    version 16.0

    syntax varlist(min=3 numeric) [if] [in] [, atgroup(string) numybins(string) maxdefiersshare(string) allowmindefiers reg_formula(string)]

    if ("`numybins'" == "") local numybins 5
    if ("`maxdefiersshare'" == "") local maxdefiersshare 0
    local numybins = real("`numybins'")
    local maxdefiersshare = real("`maxdefiersshare'")

    marksample touse
    local nvars : word count `varlist'
    local d : word 1 of `varlist'
    local y : word `nvars' of `varlist'
    local mvars
forvalues i = 2/`=`nvars'-1' {
    local mi : word `i' of `varlist'
    local mvars `mvars' `mi'
}

local nmvars : word count `mvars'

* Canonicalize mediator order when two mediators are supplied,
* so results do not depend on the order in which users type m1 m2.
if `nmvars' == 2 {
    local m1 : word 1 of `mvars'
    local m2 : word 2 of `mvars'

    * Canonicalize to descending lexical order.
    * This makes "relationship_husb grandmother" the internal order.
    if "`m1'" < "`m2'" {
        local mvars `m2' `m1'
    }
}

	
    if `maxdefiersshare' < 0 {
        di as err "maxdefiersshare() must be nonnegative"
        exit 198
    }

    tempvar wt touse2 ywork missm mgroup
quietly gen double `wt' = 1 if `touse'
quietly egen byte `missm' = rowmiss(`mvars') if `touse'
quietly gen byte `touse2' = `touse' & !missing(`d', `y', `wt') & `missm' == 0

* Preserve original mediator values if there is only one mediator.
* Only create grouped joint mediator values when there are two mediators.
if `nmvars' == 1 {
    local m : word 1 of `mvars'
    quietly gen double `mgroup' = `m' if `touse2'
}
else {
    quietly egen long `mgroup' = group(`mvars') if `touse2'
}

    quietly count if `touse2' & !inlist(`d', 0, 1)
    if r(N) > 0 {
        di as err "treatment variable must be coded 0/1"
        exit 198
    }

    quietly count if `touse2' & `d' == 1
    local n1 = r(N)
    quietly count if `touse2' & `d' == 0
    local n0 = r(N)
    if (`n1' == 0 | `n0' == 0) {
        di as err "both treatment groups (d==0 and d==1) must be present"
        exit 2000
    }

    if `numybins' > 0 quietly xtile `ywork' = `y' if `touse2', nq(`numybins')
    else quietly gen double `ywork' = `y' if `touse2'

    quietly levelsof `mgroup' if `touse2', local(mlevels)
    local K : word count `mlevels'
    if `K' < 2 {
        di as err "mediator must have at least 2 observed levels"
        exit 198
    }

    quietly summarize `wt' if `touse2' & `d' == 1, meanonly
    scalar W1 = r(sum)
    quietly summarize `wt' if `touse2' & `d' == 0, meanonly
    scalar W0 = r(sum)

    * Full set of y bin values (needed when reg_formula iterates over all (y,m))
    quietly levelsof `ywork' if `touse2', local(yfulllevels)

    * Guard: reg_formula only sensible with discrete y
    if "`reg_formula'" != "" {
        local Ky : word count `yfulllevels'
        if `Ky' > `numybins' + 5 {
            di as err "reg_formula requires y to be discrete (currently `Ky' distinct values)."
            di as err "Set numybins() to discretize y first."
            exit 198
        }
    }

    tempname p1 p0 maxdiff
    matrix `p1' = J(`K',1,0)
    matrix `p0' = J(`K',1,0)
    matrix `maxdiff' = J(`K',1,0)

    local k = 0
    foreach mv of local mlevels {
        local ++k
        quietly summarize `wt' if `touse2' & `d' == 1 & `mgroup' == `mv', meanonly
        matrix `p1'[`k',1] = r(sum) / W1
        quietly summarize `wt' if `touse2' & `d' == 0 & `mgroup' == `mv', meanonly
        matrix `p0'[`k',1] = r(sum) / W0

        if "`reg_formula'" == "" {
            quietly levelsof `ywork' if `touse2' & `mgroup' == `mv', local(yvals)
            scalar thisdiff = 0
            foreach yv of local yvals {
                quietly summarize `wt' if `touse2' & `d' == 1 & `mgroup' == `mv' & `ywork' == `yv', meanonly
                scalar p1y = r(sum) / W1
                quietly summarize `wt' if `touse2' & `d' == 0 & `mgroup' == `mv' & `ywork' == `yv', meanonly
                scalar p0y = r(sum) / W0
                scalar thisdiff = thisdiff + max(p1y - p0y, 0)
            }
            matrix `maxdiff'[`k',1] = thisdiff
        }
        else {
            * Regression branch: iterate over the FULL set of y values
            * (not just those observed at this m), because a regression
            * can produce a nonzero coefficient for an unobserved (y,m).
            * Strip the leading "~" if user included it.
            local reg_rhs = strtrim("`reg_formula'")
            if substr("`reg_rhs'", 1, 1) == "~" {
                local reg_rhs = strtrim(substr("`reg_rhs'", 2, .))
            }
            scalar thisdiff = 0
            foreach yv of local yfulllevels {
                testmechs__reg_prob, y_val(`yv') m_val(`mv') ywork(`ywork') mgroup(`mgroup') dvar(`d') touse2(`touse2') trans(treated) rhs("`reg_rhs'")
                scalar p1y = r(p)
                testmechs__reg_prob, y_val(`yv') m_val(`mv') ywork(`ywork') mgroup(`mgroup') dvar(`d') touse2(`touse2') trans(control) rhs("`reg_rhs'")
                scalar p0y = r(p)
                scalar thisdiff = thisdiff + max(p1y - p0y, 0)
            }
            matrix `maxdiff'[`k',1] = thisdiff
        }
    }

    local atindex 0
    if "`atgroup'" != "" {
        local atgroup_num = real("`atgroup'")
        local k = 0
        foreach mv of local mlevels {
            local ++k
            if (`mv' == `atgroup_num') local atindex `k'
        }
        if `atindex' == 0 {
            di as err "atgroup() must equal one observed mediator level"
            exit 198
        }
    }

    * ==================== GLPK-based LP solver ====================
    * Replaces the previous Python simplex with two GLPK plugin calls.
    * Same math: LP1 finds the min defiers share; LP2 gives the fractional
    * lower bound. GLPK is called via testmechs_lbfrac__glpk() Mata helper.

    * Pack marginals and maxdiffs into row vectors readable by Mata
    tempname pp1 pp0 mdd
    matrix `pp1' = `p1'
    matrix `pp0' = `p0'
    matrix `mdd' = `maxdiff'

    local allowmin_lp = cond("`allowmindefiers'" != "", 1, 0)

    scalar __tm_lp_K        = `K'
    scalar __tm_lp_atindex  = `atindex'
    scalar __tm_lp_maxdef   = `maxdefiersshare'
    scalar __tm_lp_allowmin = `allowmin_lp'

    mata: testmechs_lbfrac__solve_lps("`pp1'", "`pp0'", "`mdd'")

    tempname lb min_def maxdef_used
    scalar `lb'          = scalar(__tm_lp_lb)
    scalar `min_def'     = scalar(__tm_lp_min_def)
    scalar `maxdef_used' = scalar(__tm_lp_maxdef_used)

    return scalar lb = `lb'
    return scalar min_defier_share = `min_def'
    return scalar maxdefiersshare_used = `maxdef_used'

    di as txt "testmechs_lb_fracaffected"
    di as res "  lower bound = " %9.6f `lb'
end




*==============================================================
* Mata helpers for the GLPK-based LP solver.
* Implements LP1 (feasibility, min defiers share) and LP2 (fractional
* lower bound). Called from the ado via testmechs_lbfrac__solve_lps().
*==============================================================
mata:

// Solve GLPK LP: min c'x s.t. Aeq x = beq, x >= 0.
// Returns the objective on success, missing on failure.
real scalar testmechs_lbfrac__glpk(real colvector c, real matrix Aeq, real colvector beq)
{
    st_matrix("__tm_c",   c)
    st_matrix("__tm_Aeq", Aeq)
    st_matrix("__tm_beq", beq)
    stata(`"plugin call _testmechs_glpk_lp, "__tm_c" "__tm_Aeq" "__tm_beq" "__tm_lp_fun" "__tm_lp_ok""')
    if (st_numscalar("__tm_lp_ok") != 1) return(.)
    return(st_numscalar("__tm_lp_fun"))
}

// LP1: feasibility. Find the minimum sum_{i>j} X[i,j] subject to the
// marginal constraints. Returns the minimum defiers share.
real scalar testmechs_lbfrac__lp1(real colvector p1, real colvector p0)
{
    real scalar K, nvar, i, j
    real colvector c, beq
    real matrix Aeq

    K    = length(p1)
    nvar = K * K

    // Objective: 1 on X[i,j] for i>j, else 0. Flatten with idx(i,j)=(i-1)*K+j
    c = J(nvar, 1, 0)
    for (i = 1; i <= K; i++) {
        for (j = 1; j <= K; j++) {
            if (i > j) c[(i-1) * K + j] = 1
        }
    }

    // 2K equality constraints: K for j-marginals, K for i-marginals
    Aeq = J(2 * K, nvar, 0)
    beq = J(2 * K, 1, 0)

    // Row j (1..K): sum_i X[i,j] = p1[j]
    for (j = 1; j <= K; j++) {
        for (i = 1; i <= K; i++) {
            Aeq[j, (i-1) * K + j] = 1
        }
        beq[j] = p1[j]
    }
    // Row K+i (i=1..K): sum_j X[i,j] = p0[i]
    for (i = 1; i <= K; i++) {
        for (j = 1; j <= K; j++) {
            Aeq[K + i, (i-1) * K + j] = 1
        }
        beq[K + i] = p0[i]
    }

    return(testmechs_lbfrac__glpk(c, Aeq, beq))
}

// LP2: fractional lower bound.
//
// Variables (in this order):
//   X[i,j] for i,j in 1..K              (K^2 vars, idx(i,j) = (i-1)*K + j)
//   s[k] for k in 1..K                  (K vars,   at positions K^2+1..K^2+K)
//   t                                   (1 var,    at position K^2+K+1)
//   slack_k for k in 1..K               (K slacks for the "defier cap" rows)
//   slack_maxdef                        (1 slack for the maxdef row)
// Total variables: K^2 + K + 1 + K + 1 = K^2 + 2K + 2
//
// Constraints (all as equalities after adding slacks):
//   K rows:  sum_i X[i,j] - p1[j]*t = 0                     (j-marginal, equality)
//   K rows:  sum_j X[i,j] - p0[i]*t = 0                     (i-marginal, equality)
//   1 row:   sum_{g in groups} X[g,g] = 1                   (equality)
//   K rows:  -sum_{i!=k} X[i,k] - s[k] + md[k]*t + slack_k = 0   (defier cap)
//   1 row:   sum_{i>j} X[i,j] - maxdef*t + slack_maxdef = 0      (maxdef bound)
// Total rows: 3K + 2
real scalar testmechs_lbfrac__lp2(real colvector p1, real colvector p0,
                                   real colvector md, real scalar atindex,
                                   real scalar maxdef)
{
    real scalar K, nX, nS, nT, nSlack, nvar, nrows, i, j, k, g, rowc, tix, slack_start, gi
    real colvector c, beq, groups
    real matrix Aeq

    K = length(p1)
    nX = K * K
    nS = K
    nT = 1
    nSlack = K + 1
    nvar = nX + nS + nT + nSlack

    tix = nX + nS + nT                    // 1-based index of t
    slack_start = nX + nS + nT + 1        // 1-based index of first slack

    // Which g contribute to the "sum X[g,g] = 1" row and to the objective
    if (atindex == 0) groups = (1::K)
    else              groups = J(1, 1, atindex)

    // Objective: 1 on s[g] for g in groups, 0 otherwise. s[k] is at nX + k.
    c = J(nvar, 1, 0)
    for (gi = 1; gi <= length(groups); gi++) {
        g = groups[gi]
        c[nX + g] = 1
    }

    nrows = 3 * K + 2
    Aeq = J(nrows, nvar, 0)
    beq = J(nrows, 1, 0)

    rowc = 0

    // Rows 1..K: sum_i X[i,j] - p1[j]*t = 0
    for (j = 1; j <= K; j++) {
        rowc = rowc + 1
        for (i = 1; i <= K; i++) {
            Aeq[rowc, (i-1)*K + j] = 1
        }
        Aeq[rowc, tix] = -p1[j]
    }

    // Rows K+1..2K: sum_j X[i,j] - p0[i]*t = 0
    for (i = 1; i <= K; i++) {
        rowc = rowc + 1
        for (j = 1; j <= K; j++) {
            Aeq[rowc, (i-1)*K + j] = 1
        }
        Aeq[rowc, tix] = -p0[i]
    }

    // Row 2K+1: sum_{g in groups} X[g,g] = 1
    rowc = rowc + 1
    for (gi = 1; gi <= length(groups); gi++) {
        g = groups[gi]
        Aeq[rowc, (g-1)*K + g] = 1
    }
    beq[rowc] = 1

    // Rows 2K+2..3K+1: for each k, -sum_{i!=k} X[i,k] - s[k] + md[k]*t + slack_k = 0
    for (k = 1; k <= K; k++) {
        rowc = rowc + 1
        for (i = 1; i <= K; i++) {
            if (i != k) Aeq[rowc, (i-1)*K + k] = -1
        }
        Aeq[rowc, nX + k] = -1
        Aeq[rowc, tix] = md[k]
        Aeq[rowc, slack_start + (k - 1)] = 1
    }

    // Row 3K+2: sum_{i>j} X[i,j] - maxdef*t + slack_maxdef = 0
    rowc = rowc + 1
    for (i = 1; i <= K; i++) {
        for (j = 1; j <= K; j++) {
            if (i > j) Aeq[rowc, (i-1)*K + j] = 1
        }
    }
    Aeq[rowc, tix] = -maxdef
    Aeq[rowc, slack_start + K] = 1

    return(testmechs_lbfrac__glpk(c, Aeq, beq))
}

// Driver called from the ado. Reads Stata state, computes both LPs,
// writes results back to Stata scalars.
void testmechs_lbfrac__solve_lps(string scalar p1name, string scalar p0name, string scalar mdname)
{
    real colvector p1, p0, md
    real scalar K, atindex, maxdef, allowmin
    real scalar min_def, maxdef_used, lb

    K        = st_numscalar("__tm_lp_K")
    atindex  = st_numscalar("__tm_lp_atindex")
    maxdef   = st_numscalar("__tm_lp_maxdef")
    allowmin = st_numscalar("__tm_lp_allowmin")

    p1 = st_matrix(p1name)
    p0 = st_matrix(p0name)
    md = st_matrix(mdname)

    // Ensure column vectors
    if (cols(p1) > 1) p1 = p1'
    if (cols(p0) > 1) p0 = p0'
    if (cols(md) > 1) md = md'

    min_def = testmechs_lbfrac__lp1(p1, p0)
    if (min_def == .) _error(498, "GLPK feasibility LP failed")

    if (min_def > maxdef) {
        if (allowmin == 1) maxdef_used = min_def + 1e-6
        else _error(2000, "data incompatible with maxdefiersshare when allowmindefiers is off")
    }
    else {
        maxdef_used = maxdef
    }

    lb = testmechs_lbfrac__lp2(p1, p0, md, atindex, maxdef_used)
    if (lb == .) _error(498, "GLPK fractional LP failed")

    st_numscalar("__tm_lp_lb",          lb)
    st_numscalar("__tm_lp_min_def",     min_def)
    st_numscalar("__tm_lp_maxdef_used", maxdef_used)
}

end

*==============================================================
* testmechs__reg_prob: helper used by the reg_formula branch.
* For a single (y_val, m_val) cell, run a regression under one of
* R's control_transform / treated_transform and return the
* coefficient on the treatment variable in r(p).
*
* Mirrors R's compute_regression_probs() (test_sharp_null.R).
*==============================================================
program define testmechs__reg_prob, rclass
    version 16.0
    syntax , y_val(real) m_val(real) ywork(varname) mgroup(varname) dvar(varname) touse2(varname) trans(string) rhs(string)

    tempvar lhs
    quietly gen double `lhs' = (`ywork' == `y_val' & `mgroup' == `m_val') if `touse2'

    if "`trans'" == "treated" {
        quietly replace `lhs' = `dvar' * `lhs' if `touse2'
    }
    else if "`trans'" == "control" {
        quietly replace `lhs' = (`dvar' - 1) * `lhs' if `touse2'
    }
    else {
        di as err "trans must be 'treated' or 'control'"
        exit 198
    }

    * Guard: if lhs is constant, regression fails; return 0
    quietly summarize `lhs' if `touse2'
    if r(sd) == 0 {
        return scalar p = 0
        exit
    }

    * Detect IV syntax: does rhs contain "(<var> = <var>)"?
    local iv_pattern 0
    if regexm("`rhs'", "\((.+)=(.+)\)") local iv_pattern 1

    if `iv_pattern' {
        local endog = strtrim(regexs(1))
        local instr = strtrim(regexs(2))

        * Remove the IV portion from rhs to get controls
        local controls : subinstr local rhs "(`endog' = `instr')" "", all
        local controls : subinstr local controls "(`endog'= `instr')" "", all
        local controls : subinstr local controls "(`endog' =`instr')" "", all
        local controls : subinstr local controls "(`endog'=`instr')" "", all

        * Strip "+" and collapse whitespace
        local controls : subinstr local controls "+" " ", all
        local controls = strtrim(stritrim("`controls'"))

        * Try IV first; on r(481) collinearity (e.g. perfect instrument),
        * fall back to OLS — R's fixest silently produces this result.
        capture quietly ivregress 2sls `lhs' `controls' (`endog' = `instr') if `touse2'
        if _rc == 481 {
            quietly regress `lhs' `endog' `controls' if `touse2'
            return scalar p = _b[`endog']
        }
        else if _rc == 0 {
            return scalar p = _b[`endog']
        }
        else {
            ivregress 2sls `lhs' `controls' (`endog' = `instr') if `touse2'
        }
    }
    else {
        * Strip "+" and collapse whitespace
        local clean_rhs : subinstr local rhs "+" " ", all
        local clean_rhs = strtrim(stritrim("`clean_rhs'"))

        quietly regress `lhs' `clean_rhs' if `touse2'
        return scalar p = _b[`dvar']
    }
end
