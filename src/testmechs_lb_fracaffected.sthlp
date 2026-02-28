{smcl}
{* *! version 0.2.0}
{title:Title}

{p 4 8}
{cmd:testmechs_lb_fracaffected} {hline 2} lower bound on fraction of always-takers affected

{title:Syntax}

{p 8 12}
{cmd:testmechs_lb_fracaffected} {it:d m y} {ifin} [{it:weight}] [{cmd:,} {opt atgroup(#)} {opt numybins(#)} {opt maxdefiersshare(#)} {opt allowmindefiers}]

{title:Description}

{pstd}
Translation of R function {cmd:lb_frac_affected} for the default discrete/binned-Y path.
This command supports:

{pmore}
- binary treatment {it:d} coded 0/1
{break}
- multi-valued discrete mediator {it:m}
{break}
- discrete outcome {it:y} (or discretized via {opt numybins()})
{break}
- bounded defiers via {opt maxdefiersshare(#)} with allow-min-defiers behavior

{pstd}
Not implemented: {cmd:regformula()} and continuous-Y path.

{title:Options}

{phang}
{opt atgroup(#)} requests the subgroup bound for one mediator level. If omitted,
the command reports the pooled always-taker weighted-average bound.

{phang}
{opt numybins(#)} discretizes {it:y} into # quantile bins before computation.
Default is 5.

{phang}
{opt maxdefiersshare(#)} sets an upper bound on the share of defiers. Default is 0.

{phang}
{opt allowmindefiers} accepted for compatibility with the R interface. Current
implementation follows the R default behavior ({cmd:allow_min_defiers = TRUE}) and
will automatically increase {cmd:maxdefiersshare()} to the data-feasible minimum
(plus a tiny tolerance) when needed.

{title:Returned results}

{pstd}
Scalars in {cmd:r()}:

{synoptset 24 tabbed}
{synopt:{cmd:r(lb)}}requested lower bound (group-specific or pooled){p_end}
{synopt:{cmd:r(min_defier_share)}}minimum data-feasible defier share{p_end}
{synopt:{cmd:r(maxdefiersshare_used)}}effective defier cap used by the solver{p_end}

{title:Example}

{phang2}{cmd:. do tests/test_lb_fracaffected.do}
