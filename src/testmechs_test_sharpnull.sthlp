{smcl}
{* *! version 0.2.0}{...}
{title:Title}

{p 4 8}
{cmd:testmechs_test_sharpnull} {hline 2} Cox and Shi test of the sharp null of full mediation


{title:Syntax}

{p 8 12}
{cmd:testmechs_test_sharpnull} {it:d} {it:m1} [{it:m2} ...] {it:y}{cmd:,} {opt method(CS)}
   [{opt numybins(#)} {opt cluster(varname)} {opt maxdefiersshare(#)} {opt allowmindefiers} {opt reg_formula(string)}]


{title:Description}

{pstd}
Stata translation of R function {cmd:test_sharp_null()} for {cmd:method = "CS"}
(Cox and Shi, 2023). Tests the sharp null of full mediation: the hypothesis
that the treatment {it:d} affects the outcome {it:y} only through the mediator
{it:m1} (or the combination {it:m1 m2}). Rejection means there is a direct
effect of the treatment that is not routed through the specified mediators.

{pstd}
Supports single or two mediators (multi-valued discrete), continuous outcomes
discretized via {opt numybins()}, cluster-robust inference, relaxed monotonicity
via {opt maxdefiersshare()}, and regression adjustment (OLS with controls or
2SLS IV) via {opt reg_formula()}.


{title:Options}

{phang}{opt method(CS)} must be {cmd:CS}. Other R-package methods (ARP, FSST,
toru) are not implemented.

{phang}{opt numybins(#)} number of bins for discretizing {it:y}. Default is 5.
The test relies on a central limit theorem approximation within cells defined
by {it:(y, m, d)}, so the number of bins should be small enough that the CLT
is reasonable within each cell.

{phang}{opt cluster(varname)} cluster variable for cluster-robust variance. If
omitted, standard (i.i.d.) variance is used.

{phang}{opt maxdefiersshare(#)} upper bound on the population share of defiers,
relaxing the default monotonicity assumption that treatment can only increase
{it:m}. Default is 0 (strict monotonicity).

{phang}{opt allowmindefiers} if the requested {opt maxdefiersshare()} is
infeasible against the empirical distribution, raise it to the minimum feasible
value plus a tiny tolerance. Matches R's {cmd:allow_min_defiers = TRUE}.

{phang}{opt reg_formula(string)} regression formula (with a leading tilde) for
computing regression-adjusted partial probabilities under conditional
random assignment. Two forms are accepted:

{p 12 12 2}
- OLS with controls: {cmd:"~ treat + age + educ + wealth"}. Runs
{cmd:regress lhs treat age educ wealth} in each {it:(y, m)} cell and takes the
coefficient on {it:d}.

{p 12 12 2}
- 2SLS IV: {cmd:"~ age + educ + wealth + (treat = instrument)"}. Runs
{cmd:ivregress 2sls lhs age educ wealth (treat = instrument)}. If the instrument
is perfectly collinear with the treatment (Stata error {cmd:r(481)}), the helper
silently falls back to OLS with the endogenous variable as a regular regressor,
matching R's {cmd:fixest::feols} behavior in the same degenerate case.

{p 8 12 2}
Cluster-robust inference under {opt reg_formula()} uses analytic influence
functions of the sandwich form: for OLS, {cmd:IF_i = M * x_i * e_i * n_reg} with
{cmd:M = (X'X)^-1}; for 2SLS, {cmd:IF_i = M * xhat_i * e_iv,i * n_reg} where
{it:xhat} substitutes first-stage fitted values for the endogenous variable
and {cmd:M = e(V)/e(rmse)^2}. This reproduces R's
{cmd:sandwich::estfun(feols) %*% t(sandwich::bread(feols))} to machine
precision on the same data (see {cmd:tests/verify/} in the source repo).


{title:Returned results}

{pstd}
Scalars in {cmd:r()}:

{synoptset 24 tabbed}{...}
{synopt:{cmd:r(pval)}}p-value of the sharp-null test{p_end}
{synopt:{cmd:r(test_stat)}}Cox-Shi test statistic{p_end}
{synopt:{cmd:r(cv)}}chi-squared critical value at 5%{p_end}


{title:Examples}

{pstd}Baseline test (randomly-assigned setting):{p_end}
{phang2}{cmd:. use "data/mother_data.dta", clear}{p_end}
{phang2}{cmd:. testmechs_test_sharpnull treat grandmother motherfinancial, method(CS) numybins(5) cluster(uc)}{p_end}
{p 8 8 2}Expected: {cmd:p-value = 0.02283916} (matches R).

{pstd}Multi-valued mediator (relationship quality, 1-5 scale):{p_end}
{phang2}{cmd:. testmechs_test_sharpnull treat relationship_husb motherfinancial, method(CS) numybins(5) cluster(uc)}{p_end}
{p 8 8 2}Expected: {cmd:p-value = 0.02838332} (matches R).

{pstd}Combination of two mediators:{p_end}
{phang2}{cmd:. testmechs_test_sharpnull treat relationship_husb grandmother motherfinancial, method(CS) numybins(5) cluster(uc)}{p_end}
{p 8 8 2}Expected: {cmd:p-value = 0.65408650} (matches R).

{pstd}With OLS baseline controls (non-experimental setting):{p_end}
{phang2}{cmd:. testmechs_test_sharpnull treat grandmother motherfinancial, ///}{p_end}
{phang2}{cmd:      method(CS) numybins(5) cluster(uc) ///}{p_end}
{phang2}{cmd:      reg_formula("~ treat + age_baseline + edu_mo_baseline + wealth_baseline")}{p_end}
{p 8 8 2}Expected: {cmd:p-value = 0.03105106} (matches R exactly).

{pstd}With 2SLS IV:{p_end}
{phang2}{cmd:. set seed 0}{p_end}
{phang2}{cmd:. gen double iv = treat + rnormal(0, 0.1)}{p_end}
{phang2}{cmd:. testmechs_test_sharpnull treat grandmother motherfinancial, ///}{p_end}
{phang2}{cmd:      method(CS) numybins(5) cluster(uc) ///}{p_end}
{phang2}{cmd:      reg_formula("~ age_baseline + edu_mo_baseline + wealth_baseline + (treat = iv)")}{p_end}
{p 8 8 2}Note: Stata's {cmd:rnormal()} and R's {cmd:rnorm()} use different random
streams, so the constructed {it:iv} will differ from R's example even with the same
seed. Only the analytical procedure is shared.


{title:Notes on numerical agreement with R}

{pstd}
For the default (no {opt reg_formula()}) and OLS-with-controls paths, p-values
match R to at least 7 decimal places. For 2SLS IV, the per-observation influence
function matches R's sandwich IF to machine precision (~1e-14) and end-to-end
p-values match R exactly for K >= 3 (multi-valued or combined mediators).

{pstd}
For K = 2 (binary M) IV cases, R auto-dispatches to a specialised
{cmd:test_sharp_null_binary_m} code path with a different Cox-Shi variant
tuned for the binary-M no-nuisance case. The Stata port only implements the
general CS code path, so binary-M IV p-values may differ from R's binary-M
p-values. The general CS output remains a valid CS test.


{title:Not implemented}

{pstd}
The following R package features are not implemented in the Stata port and
will return an error:

{phang}{cmd:continuousy} option{p_end}
{phang}{cmd:method(ARP)}, {cmd:method(FSST)}, {cmd:method(toru)}{p_end}
{phang}Specialised binary-M code path (used by R for K=2){p_end}
{phang}{cmd:fixest}-style fixed effect syntax ({cmd:| fe_var}){p_end}


{title:Reference}

{pstd}
Kwon, S. and J. Roth (2024). Testing Mechanisms.
{browse "https://www.jonathandroth.com/assets/files/TestingMechanisms_Draft.pdf"}
{p_end}

{pstd}
R package source and this Stata port on GitHub:
{browse "https://github.com/chengchen2326/stata-testmechs"}
{p_end}
