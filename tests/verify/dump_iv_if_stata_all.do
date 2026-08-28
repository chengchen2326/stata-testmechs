clear all
set more off

* Dump IV IF for all 3 instruments (treat2, treat3, treat4).
* Cell: y_bin=1, m_val=0, treated transform.

cd "~/Library/CloudStorage/OneDrive-Personal/Brown/26 RA with Jon/stata-testmechs"

* ---- Define a mata helper once (batch mode, not interactive) ----
mata
void build_iv_if(string scalar Mrow_name, string scalar resid_name,
                 string scalar xhat_endog_name, string rowvector control_names,
                 real scalar nreg, string scalar out_matname)
{
    real rowvector Mrow
    real colvector resid, IF
    real matrix xhat
    real scalar i, n

    Mrow  = st_matrix(Mrow_name)
    resid = st_data(., resid_name)
    xhat  = st_data(., xhat_endog_name), st_data(., control_names)
    xhat  = xhat, J(rows(xhat), 1, 1)

    n  = rows(xhat)
    IF = J(n, 1, 0)
    for (i=1; i<=n; i++) {
        if (resid[i] < .) IF[i] = (Mrow * xhat[i,]') * resid[i] * nreg
    }
    st_matrix(out_matname, IF)
}
end

use "data/mother_data_extended.dta", clear

gen byte touse2 = !missing(treat, grandmother, motherfinancial, uc, ///
                            age_baseline, edu_mo_baseline, wealth_baseline, ///
                            treat2, treat3, treat4)
keep if touse2

gen row = _n
preserve
    import delimited using "/tmp/r_iv_if_treat3.csv", clear
    keep row y_bin
    tempfile r_ybin
    save `r_ybin', replace
restore
merge 1:1 row using `r_ybin', keep(master matched) nogen

gen double lhs = treat * (y_bin == 1 & grandmother == 0)

* Loop over 3 instruments
foreach instr in treat2 treat3 treat4 {
    di ""
    di "========== IV with `instr' =========="

    capture ivregress 2sls lhs age_baseline edu_mo_baseline wealth_baseline (treat = `instr')
    local rc_iv = _rc

    if `rc_iv' == 481 {
        di "  -> r(481): fall back to OLS"
        regress lhs treat age_baseline edu_mo_baseline wealth_baseline
    }
    else if `rc_iv' != 0 {
        di as err "  Unexpected rc: `rc_iv'"
        continue
    }

    di "Coef treat = " %10.8f _b[treat]

    tempname V_iv M_iv Mrow_iv
    matrix V_iv = e(V)
    scalar sigma2_iv = e(rmse)^2
    scalar nreg_iv   = e(N)
    matrix M_iv = V_iv / sigma2_iv
    local rn : rownames M_iv
    local j_iv = 0
    local ii = 0
    foreach nm of local rn {
        local ++ii
        if "`nm'" == "treat" local j_iv = `ii'
    }
    matrix Mrow_iv = M_iv[`j_iv', 1..colsof(M_iv)]

    capture drop resid_iv
    predict double resid_iv, residuals

    capture drop xhat_treat
    if `rc_iv' == 481 {
        gen double xhat_treat = treat
    }
    else {
        regress treat age_baseline edu_mo_baseline wealth_baseline `instr'
        predict double xhat_treat, xb
    }

    mata: build_iv_if("Mrow_iv", "resid_iv", "xhat_treat", ///
                      ("age_baseline", "edu_mo_baseline", "wealth_baseline"), ///
                      st_numscalar("nreg_iv"), "if_stata_mat")

    preserve
        clear
        svmat double if_stata_mat, names(if_treat)
        gen row = _n
        keep row if_treat
        order row if_treat
        export delimited using "/tmp/stata_iv_if_`instr'.csv", replace
    restore

    di "IF first 5:"
    mata: st_matrix("if_stata_mat")[1..5, 1]
}

di ""
di "All 3 Stata IF files dumped to /tmp/stata_iv_if_*.csv"
