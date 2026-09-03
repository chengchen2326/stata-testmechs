clear all
set more off

* ============================================================
* Generate a partial density plot to include in the README as a
* static example figure, illustrating what test_sharp_null tests.
*
* Reproduces R README's nt_plot for the Baranov et al. (2020) data:
* P(Y=y, M=0 | D=1) vs P(Y=y, M=0 | D=0) for the "never-takers"
* (grandmother NOT present, M=0).
*
* Under monotonicity + full mediation, we should have
*   P(Y=y, M=0 | D=1) <= P(Y=y, M=0 | D=0)  for all y.
* A treated bar taller than the control bar violates the sharp null.
* ============================================================

cd "~/Library/CloudStorage/OneDrive-Personal/Brown/26 RA with Jon/stata-testmechs"
use "data/mother_data.dta", clear

* Discretize Y into 5 quantile bins (matching numybins(5) in test_sharpnull)
xtile y_bin = motherfinancial, nq(5)

* Compute P(Y=y, M=0 | D=1) and P(Y=y, M=0 | D=0) for y = 1..5
* i.e. jointly conditional on treatment arm and grandmother==0.
tempname tab_treated tab_control
preserve
    keep if treat == 1
    * denom = N in treated arm
    scalar n_treated = _N
    * for each y bin, count where M=0 AND y=y_bin, divide by n_treated
    forval y = 1/5 {
        count if grandmother == 0 & y_bin == `y'
        scalar p_treated_`y' = r(N) / n_treated
    }
restore

preserve
    keep if treat == 0
    scalar n_control = _N
    forval y = 1/5 {
        count if grandmother == 0 & y_bin == `y'
        scalar p_control_`y' = r(N) / n_control
    }
restore

* Build a small dataset for plotting
clear
set obs 10
gen y_bin  = mod(_n-1, 5) + 1
gen treat_str = cond(_n <= 5, "Treated (D=1)", "Control (D=0)")
gen double prob = .

forval y = 1/5 {
    replace prob = scalar(p_treated_`y') if y_bin == `y' & treat_str == "Treated (D=1)"
    replace prob = scalar(p_control_`y') if y_bin == `y' & treat_str == "Control (D=0)"
}

* Save data for reproducibility
list y_bin treat_str prob, sepby(y_bin)

* Bar chart: y_bin on x-axis, side-by-side bars for treat vs control
graph bar (asis) prob, ///
    over(treat_str, gap(0)) ///
    over(y_bin, gap(50)) ///
    asyvars ///
    bar(1, color(navy%80))    ///
    bar(2, color(cranberry%80)) ///
    ytitle("P(Y=y, M=0 | D)", size(medsmall)) ///
    b1title("Y bin (motherfinancial quintiles)", size(medsmall)) ///
    title("Partial density: no grandmother present (M=0)", size(medium)) ///
    subtitle("Sharp null + monotonicity => treated bars <= control bars for all Y", size(small)) ///
    legend(order(1 "P(Y=y, M=0 | D=1) - Treated" 2 "P(Y=y, M=0 | D=0) - Control") ///
           rows(1) size(small)) ///
    graphregion(color(white))

* Export as PNG for README embedding
graph export "figures/partial_density_grandmother_nt.png", ///
    width(1200) replace

di ""
di "Figure saved to figures/partial_density_grandmother_nt.png"
