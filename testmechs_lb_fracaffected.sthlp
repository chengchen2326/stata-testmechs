{smcl}
{* *! version 0.1.0}
{title:Title}

{p 4 8}
{cmd:testmechs_lb_fracaffected} {hline 2} MVP lower bound on fraction of always-takers affected

{title:Syntax}

{p 8 12}
{cmd:testmechs_lb_fracaffected} {it:d m y} {ifin} [{it:weight}] [{cmd:,} {opt atgroup(#)} {opt numybins(#)}]

{title:Description}

{pstd}
MVP translation of R function {cmd:lb_frac_affected} (default path only).
This command currently supports:

{pmore}
- binary treatment {it:d} coded 0/1
{break}
- binary mediator {it:m}
{break}
- discrete outcome {it:y} (or discretized via {opt numybins()})
{break}
- no defiers ({cmd:max_defiers_share = 0} only)

{pstd}
Unsupported R options/paths are intentionally not implemented in MVP and return a clear error.

{title:Options}

{phang}
{opt atgroup(#)} requests the subgroup bound for one mediator level. Must equal one of the two observed mediator values.

{phang}
{opt numybins(#)} discretizes {it:y} into # quantile bins before computation.

{title:Returned results}

{pstd}
Scalars in {cmd:r()}:

{synoptset 24 tabbed}
{synopt:{cmd:r(lb)}}requested lower bound (group-specific or pooled){p_end}
{synopt:{cmd:r(lb_group_lo)}}bound for low mediator level group{p_end}
{synopt:{cmd:r(lb_group_hi)}}bound for high mediator level group{p_end}
{synopt:{cmd:r(theta_nt)}}share of never-takers (binary-M mapping){p_end}
{synopt:{cmd:r(theta_at)}}share of always-takers (binary-M mapping){p_end}
{synopt:{cmd:r(theta_c)}}share of compliers (binary-M mapping){p_end}
{synopt:{cmd:r(max_p_diff_lo)}}max partial-density difference for low mediator level{p_end}
{synopt:{cmd:r(max_p_diff_hi)}}max partial-density difference for high mediator level{p_end}
{synopt:{cmd:r(m_lo)}}numeric value of low mediator level{p_end}
{synopt:{cmd:r(m_hi)}}numeric value of high mediator level{p_end}

{title:Example}

{phang2}{cmd:. do tests/test_lb_fracaffected.do}
