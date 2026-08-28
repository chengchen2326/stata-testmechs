{smcl}
{* *! version 0.3.0}{...}
{title:Title}

{p 4 8}
{cmd:testmechs_lb_fracaffected} {hline 2} lower bound on the fraction of always-takers affected by treatment


{title:Syntax}

{p 8 12}
{cmd:testmechs_lb_fracaffected} {it:d} {it:m1} [{it:m2} ...] {it:y} {ifin}
   [{cmd:,} {opt atgroup(#)} {opt numybins(#)} {opt maxdefiersshare(#)}
    {opt allowmindefiers} {opt reg_formula(string)}]


{title:Description}

{pstd}
Stata translation of R function {cmd:lb_frac_affected()}. Computes a lower
bound on the fraction of "always-takers" (or "never-takers", or the pooled
group) whose outcome is affected by the treatment despite having the same
value of the mediator {it:m} under both treatments. This gives a sense of the
strength of mechanisms other than {it:m}: it tells us what fraction of the
always-takers (or never-takers) have a direct effect of the treatment not
routed through {it:m}.

{pstd}
Supports single or two mediators (multi-valued discrete), continuous outcomes
discretized via {opt numybins()}, relaxed monotonicity via
{opt maxdefiersshare()}, and regression adjustment (OLS with controls or
2SLS IV) via {opt reg_formula()}.


{title:Options}

{phang}{opt atgroup(#)} requests the lower bound for one specific mediator
level {it:#} (e.g., {opt atgroup(0)} for never-takers when {it:m} is 0/1).
If omitted, the command reports the pooled always-taker weighted-average
bound (matching R's {cmd:at_group = NULL}).

{phang}{opt numybins(#)} discretizes {it:y} into # quantile bins before
computation. Default is 5.

{phang}{opt maxdefiersshare(#)} upper bound on the population share of defiers,
relaxing the default monotonicity assumption. Default is 0.

{phang}{opt allowmindefiers} if the requested {opt maxdefiersshare()} is
infeasible against the empirical distribution, raise it to the minimum feasible
value plus a tiny tolerance. Matches R's {cmd:allow_min_defiers = TRUE}. If
omitted, infeasible {opt maxdefiersshare()} causes an error.

{phang}{opt reg_formula(string)} regression formula (with a leading tilde) for
computing regression-adjusted partial probabilities. Two forms:

{p 12 12 2}
- OLS with controls: {cmd:"~ treat + age + educ + wealth"}. Runs
{cmd:regress lhs treat age educ wealth} in each {it:(y, m)} cell and takes the
coefficient on {it:d}.

{p 12 12 2}
- 2SLS IV: {cmd:"~ age + educ + wealth + (treat = instrument)"}. Runs
{cmd:ivregress 2sls}. Perfect-instrument cases (Stata error {cmd:r(481)})
fall back to OLS, matching R's {cmd:fixest} behavior.


{title:Returned results}

{pstd}
Scalars in {cmd:r()}:

{synoptset 28 tabbed}{...}
{synopt:{cmd:r(lb)}}the requested lower bound (group-specific or pooled){p_end}
{synopt:{cmd:r(min_defier_share)}}minimum data-feasible defier share{p_end}
{synopt:{cmd:r(maxdefiersshare_used)}}effective defier cap used by the solver{p_end}


{title:Examples}

{pstd}Lower bound for never-takers ({it:m} = 0), binary mediator:{p_end}
{phang2}{cmd:. use "data/mother_data.dta", clear}{p_end}
{phang2}{cmd:. testmechs_lb_fracaffected treat grandmother motherfinancial, numybins(5) atgroup(0)}{p_end}
{p 8 8 2}Expected: {cmd:lower bound = 0.185891} (matches R).

{pstd}Same, allowing 1 percent defiers:{p_end}
{phang2}{cmd:. testmechs_lb_fracaffected treat grandmother motherfinancial, ///}{p_end}
{phang2}{cmd:      numybins(5) atgroup(0) maxdefiersshare(0.01)}{p_end}
{p 8 8 2}Expected: {cmd:lower bound = 0.171642} (matches R).

{pstd}Multi-valued mediator, pooled across always-taker groups:{p_end}
{phang2}{cmd:. testmechs_lb_fracaffected treat relationship_husb motherfinancial, ///}{p_end}
{phang2}{cmd:      numybins(5) allowmindefiers}{p_end}
{p 8 8 2}Expected: {cmd:lower bound = 0.100221} (matches R).

{pstd}Combination of two mediators:{p_end}
{phang2}{cmd:. testmechs_lb_fracaffected treat relationship_husb grandmother motherfinancial, ///}{p_end}
{phang2}{cmd:      numybins(5) allowmindefiers}{p_end}
{p 8 8 2}Expected: {cmd:lower bound = 0.072513} (matches R).

{pstd}With OLS baseline controls (non-experimental setting):{p_end}
{phang2}{cmd:. testmechs_lb_fracaffected treat grandmother motherfinancial, ///}{p_end}
{phang2}{cmd:      numybins(5) atgroup(0) ///}{p_end}
{phang2}{cmd:      reg_formula("~ treat + age_baseline + edu_mo_baseline + wealth_baseline")}{p_end}


{title:Not implemented}

{pstd}
The following R package features are not implemented in the Stata port and
will return an error:

{phang}{cmd:continuousy} option{p_end}
{phang}{cmd:returnmindefiers} option{p_end}
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
