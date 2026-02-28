program define testmechs_lb_fracaffected, rclass
    version 16.0

    // MVP translation of TestMechs::lb_frac_affected default path.
    // R args mapped to Stata positional varlist: d m y.
    syntax varlist(min=3 max=3 numeric) [if] [in] ///
	[, ATGroup(string) NUMYBins(string) MAXDefiersShare(string)]
	
	if ("`atgroup'" == "") local atgroup 0
	if ("`numybins'" == "") local numybins 5
    if ("`maxdefiersshare'" == "") local maxdefiersshare 0
    * convert to numeric (still stored as locals, but validated)
	local atgroup = real("`atgroup'")
	local numybins = real("`numybins'")
    local maxdefiersshare = real("`maxdefiersshare'")

    marksample touse
    gettoken d rest : varlist
    gettoken m y : rest

    // MVP scope guardrails: implement only default/recommended R path.
    if "`continuousy'" != "" {
        di as err "continuous_Y=TRUE path is not implemented in MVP"
        exit 198
    }
    if "`regformula'" != "" {
        di as err "reg_formula path is not implemented in MVP"
        exit 198
    }
    if `maxdefiersshare' < 0 {
        di as err "maxdefiersshare() must be nonnegative"
        exit 198
    }
    if "`allowmindefiers'" != "" {
        di as err "allow_min_defiers option is not implemented in MVP"
        exit 198
    }
    if "`returnmindefiers'" != "" {
        di as err "return_min_defiers option is not implemented in MVP"
        exit 198
    }

    tempvar wt touse2 ywork
    quietly gen double `wt' = 1 if `touse'
    if "`weight'" != "" {
        quietly replace `wt' = `exp' if `touse'
    }

    // Match remove_missing_from_df behavior for core variables.
    quietly gen byte `touse2' = `touse' & !missing(`d', `m', `y', `wt')

    // Enforce binary treatment as in the main use-case/tests.
    quietly count if `touse2' & !inlist(`d', 0, 1)
    if r(N) > 0 {
        di as err "treatment variable must be coded 0/1 in MVP"
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

    // Optional discretization of Y (R: num_Ybins).
    if `numybins' > 0 {
        quietly xtile `ywork' = `y' if `touse2', nq(`numybins')
    }
    else {
        quietly gen double `ywork' = `y' if `touse2'
    }

    // MVP currently supports a binary mediator only.
    quietly levelsof `m' if `touse2', local(mlevels)
    local K : word count `mlevels'
    if `K' != 2 {
        di as err "MVP currently supports binary mediator only; found `K' levels"
        exit 198
    }

    local m_lo : word 1 of `mlevels'
    local m_hi : word 2 of `mlevels'

    // Weighted totals within treatment arms.
    quietly summarize `wt' if `touse2' & `d' == 1, meanonly
    scalar W1 = r(sum)
    quietly summarize `wt' if `touse2' & `d' == 0, meanonly
    scalar W0 = r(sum)

    // R: p_m_1 and p_m_0.
    quietly summarize `wt' if `touse2' & `d' == 1 & `m' == `m_lo', meanonly
    scalar p_m1_lo = r(sum) / W1
    quietly summarize `wt' if `touse2' & `d' == 1 & `m' == `m_hi', meanonly
    scalar p_m1_hi = r(sum) / W1

    quietly summarize `wt' if `touse2' & `d' == 0 & `m' == `m_lo', meanonly
    scalar p_m0_lo = r(sum) / W0
    quietly summarize `wt' if `touse2' & `d' == 0 & `m' == `m_hi', meanonly
    scalar p_m0_hi = r(sum) / W0

    // R: max_p_diffs[m] = sum_y max( p(y,m|d=1)-p(y,m|d=0), 0 ).
    scalar maxdiff_lo = 0
    scalar maxdiff_hi = 0

    foreach mv in `m_lo' `m_hi' {
        quietly levelsof `ywork' if `touse2' & `m' == `mv', local(yvals)
        scalar thisdiff = 0
        foreach yv of local yvals {
            quietly summarize `wt' if `touse2' & `d' == 1 & `m' == `mv' & `ywork' == `yv', meanonly
            scalar p1 = r(sum) / W1
            quietly summarize `wt' if `touse2' & `d' == 0 & `m' == `mv' & `ywork' == `yv', meanonly
            scalar p0 = r(sum) / W0
            scalar thisdiff = thisdiff + max(p1 - p0, 0)
        }
        if `mv' == `m_lo' scalar maxdiff_lo = thisdiff
        if `mv' == `m_hi' scalar maxdiff_hi = thisdiff
    }

    // Binary-M LP with bounded defier share (matches R default path constraints).
    // Type shares as function of defier share d:
    //   theta_nt(d) = p_m1_lo - d
    //   theta_at(d) = p_m0_hi - d
    //   theta_c(d)  = p_m1_hi - p_m0_hi + d
    scalar theta_c0 = p_m1_hi - p_m0_hi
    scalar d_min = max(0, p_m0_hi - p_m1_hi)
    scalar d_data_max = min(p_m1_lo, p_m0_hi)
    if d_min > d_data_max + 1e-10 {
        di as err "binary mediator marginals are internally inconsistent"
        exit 498
    }

    scalar min_defier_share = d_min
    scalar maxdef_eff = `maxdefiersshare'
    if maxdef_eff < d_min {
        scalar maxdef_eff = d_min + 1e-6
        di as txt "note: maxdefiersshare() below data-feasible minimum; using " %9.6f maxdef_eff
    }
    scalar d_max = min(d_data_max, maxdef_eff)

    if (`atgroup' != `m_lo') & (`atgroup' != `m_hi') {
        di as err "at_group must equal one of mediator levels: `m_lo' or `m_hi'"
        exit 198
    }

    scalar best_lb = .
    scalar best_d = .

    // Candidate d values: endpoints + kinks from positive-part terms.
    local d_candidates `=d_min' `=d_max' `=maxdiff_lo' `=maxdiff_hi-theta_c0'
    foreach dc of local d_candidates {
        scalar d_try = `dc'
        if d_try < d_min - 1e-10 continue
        if d_try > d_max + 1e-10 continue

        scalar theta_nt_try = p_m1_lo - d_try
        scalar theta_at_try = p_m0_hi - d_try
        scalar theta_c_try = theta_c0 + d_try

        scalar v_lo_try = max(maxdiff_lo - d_try, 0)
        scalar v_hi_try = max(maxdiff_hi - theta_c_try, 0)

        scalar lb_try = .
        if `atgroup' == `m_lo' {
            if theta_nt_try > 0 scalar lb_try = v_lo_try / theta_nt_try
        }
        else if `atgroup' == `m_hi' {
            if theta_at_try > 0 scalar lb_try = v_hi_try / theta_at_try
        }
        else {
            scalar denom_try = theta_nt_try + theta_at_try
            if denom_try > 0 scalar lb_try = (v_lo_try + v_hi_try) / denom_try
        }

        if missing(best_lb) | (!missing(lb_try) & lb_try < best_lb) {
            scalar best_lb = lb_try
            scalar best_d = d_try
        }
    }

    // Evaluate all reported quantities at the optimizing defier share.
    scalar theta_nt = p_m1_lo - best_d
    scalar theta_at = p_m0_hi - best_d
    scalar theta_c  = theta_c0 + best_d

    scalar v_lo = max(maxdiff_lo - best_d, 0)
    scalar v_hi = max(maxdiff_hi - theta_c, 0)

    scalar lb_lo = .
    scalar lb_hi = .
    if theta_nt > 0 scalar lb_lo = v_lo / theta_nt
    if theta_at > 0 scalar lb_hi = v_hi / theta_at

    scalar lb = best_lb

    return scalar lb = lb
    return scalar lb_group_lo = lb_lo
    return scalar lb_group_hi = lb_hi
    return scalar theta_nt = theta_nt
    return scalar theta_at = theta_at
    return scalar theta_c = theta_c
    return scalar defier_share = best_d
    return scalar min_defier_share = min_defier_share
    return scalar maxdefiersshare_used = maxdef_eff
    return scalar max_p_diff_lo = maxdiff_lo
    return scalar max_p_diff_hi = maxdiff_hi
    return scalar m_lo = `m_lo'
    return scalar m_hi = `m_hi'

    di as txt "testmechs_lb_fracaffected (MVP)"
    di as res "  lower bound = " %9.6f lb
end
