{smcl}
{* *! version 0.1}
{title:testmechs_test_sharpnull}

{p 4 4 2}
{cmd:testmechs_test_sharpnull} {it:d m1 [m2 ...] y}, {cmd:method(CS)} [{cmd:numybins(}{it:#}{cmd:)} {cmd:cluster(}{it:varname}{cmd:)}]

{title:Description}
{pstd}
Minimal MVP implementation of TestMechs sharp-null test for the Cox & Shi method with one or more numeric mediators.

{title:Options}
{phang}{cmd:method(string)} must be {cmd:CS}.
{phang}{cmd:numybins(integer)} number of bins for discretizing y (default 5).
{phang}{cmd:cluster(varname)} cluster variable for cluster-robust variance.

{title:Returns}
{phang}{cmd:r(pval)} p-value.
