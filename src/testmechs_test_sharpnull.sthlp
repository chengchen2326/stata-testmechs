{smcl}
{* *! version 0.1}
{title:testmechs_test_sharpnull}

{p 4 4 2}
{cmd:testmechs_test_sharpnull} {it:d m y}, {cmd:method(CS)} [{cmd:numybins(}{it:#}{cmd:)} {cmd:cluster(}{it:varname}{cmd:)}]

{title:Description}
{pstd}
Minimal MVP implementation of TestMechs sharp-null test for the Cox & Shi method with a binary mediator.

{title:Options}
{phang}{cmd:method(string)} must be {cmd:CS}.
{phang}{cmd:numybins(integer)} number of bins for discretizing y (default 5).
{phang}{cmd:cluster(varname)} cluster variable for cluster-robust variance.

{title:Returns}
{phang}{cmd:r(pval)} p-value.

{title:Example}

{pstd}
The example dataset {cmd:mother_data.dta} is not shipped with the package
(following Stata convention that packages provide code, not data). Download
it from the repository into the current directory:

{phang2}{cmd:. copy "https://raw.githubusercontent.com/chengchen2326/stata-testmechs/Stata%2BGLPK%2Bdqrdc2/data/mother_data.dta" mother_data.dta}{p_end}
{phang2}{cmd:. use "mother_data.dta", clear}{p_end}
{phang2}{cmd:. testmechs_test_sharpnull treat relationship_husb motherfinancial, method(CS) numybins(5) cluster(uc)}{p_end}

{pstd}
Expected output: {cmd:p-value = 0.02838332} (matches the R TestMechs baseline).
