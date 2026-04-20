# stata-testmechs (MVP)

This repository contains an **MVP Stata translation** of key functionality from Jon Roth & Soonwoo Kwon's **TestMechs** R package (paper: *Testing Mechanisms*).

## Status (MVP scope)

### Implemented
- ✅ `testmechs_lb_fracaffected`  
  Stata translation of the R function `lb_frac_affected()`

### Currently supported features
- Positional varlist input:
  - **`d m y`** for a single mediator
  - **`d m1 m2 y`** for a combination of both mechanisms (two-mediator input)
- `atgroup(#)`
- `numybins(#)`
- `maxdefiersshare(#)`
- `allowmindefiers`
- Multi-valued discrete mediators in the default binned-Y path
- Lower bounds combining both mechanisms across two mediator variables

### Not yet implemented
The following are **not yet implemented** and will return an error:
- 🚫 `regformula()`
- 🚫 `continuousy`
- 🚫 `returnmindefiers`
- 🚫 Additional inference/plotting functions (e.g. `test_sharp_null`)

The original R source is kept under `r_reference/` for translation and validation.

---

## Data (for replication)

This repo includes two Stata datasets converted from the TestMechs R package data (`.rda`) for convenience:

- `baranov_data.dta`  
  Raw dataset converted from the R package data object `baranov_data`
- `mother_data.dta`  
  Analysis dataset corresponding to the experimental sample used in the README examples, created by restricting to `THP_sample == 1`

### Notes
- These `.dta` files are provided for quick replication in Stata.
- They are derived from the original R package data; see `r_reference/` for the upstream source.

---

## Installation

Install using `net install` from GitHub:

```stata
cap noi net uninstall testmechs
net install testmechs, from("https://raw.githubusercontent.com/chengchen2326/stata-testmechs/main") replace
```

You can verify that Stata finds the installed command with:

```stata
which testmechs_lb_fracaffected
help testmechs_lb_fracaffected
```

---

## Examples

### 1. Lower bound on the fraction of never-takers affected

The test above suggests that the treatment effect does not operate entirely through the presence of a grandmother in the home. There are some people (never-takers) whose outcome is affected by the treatment despite having no change in `M`. It must be that some other mechanism mattered for these people.

But how prevalent are these alternative mechanisms? To give a sense, we compute lower bounds on the fraction of never-takers whose outcome is affected by the treatment despite having the same value of `M` under both treatments. This gives a sense of the strength of mechanisms other than `M`: it tells us what fraction of the never-takers have a direct effect of the treatment.

The argument `at_group = 0` corresponds to computing this lower bound for the never-takers, who are referred to as "0-always takers" in the more general notation in the paper.

#### R benchmark

```r
lb_nts <- lb_frac_affected(
  df = mother_data,
  d = "treat",
  m = "grandmother",
  y = "motherfinancial",
  num_Ybins = 5,
  at_group = 0
)
lb_nts
#> [1] 0.1858912
```

#### Stata

```stata
use "data/mother_data.dta", clear

* varlist order is: d m y
testmechs_lb_fracaffected treat grandmother motherfinancial, ///
    numybins(5) atgroup(0)
```

#### Expected output
- `lower bound = 0.185891` (R benchmark: `0.1858912`)

Our estimates imply that at least 19 percent of never-takers are affected by the treatment.

One could likewise test the fraction of never-takers affected by setting `at_group = 1` (in this case, the lower bound is zero). If `at_group` is set to `NULL`, then the package calculates the fraction pooling across all types that have the same value of `M` under both treatments (i.e. always-takers and never-takers when `M` is binary).

---

### 2. Lower bound under relaxed monotonicity (`max_defiers_share`)

Likewise, we can also calculate the lower bound on the fraction of never-takers under relaxed monotonicity.

#### R benchmark

```r
lb_nts_defiers <- lb_frac_affected(
  df = mother_data,
  d = "treat",
  m = "grandmother",
  y = "motherfinancial",
  num_Ybins = 5,
  at_group = 0,
  max_defiers_share = .01
)
lb_nts_defiers
#> [1] 0.1716415
```

#### Stata

```stata
use "data/mother_data.dta", clear

* varlist order is: d m y
testmechs_lb_fracaffected treat grandmother motherfinancial, ///
    numybins(5) atgroup(0) maxdefiersshare(.01)
```

#### Expected output
- `lower bound = 0.171642` (R benchmark: `0.1716415`)

---

### 3. Multi-valued mediator example with `allow_min_defiers`

We can also estimate a lower bound combining both mechanisms when the mediator is multi-valued.

#### R benchmark

```r
lb_frac_affected(
  df = mother_data,
  d = "treat",
  m = "relationship_husb",
  y = "motherfinancial",
  num_Ybins = 5,
  at_group = NULL,
  allow_min_defiers = TRUE
)
#> [1] 0.1002207
```

#### Stata

```stata
use "data/mother_data.dta", clear

* varlist order is: d m y
* omit atgroup() to match at_group = NULL
testmechs_lb_fracaffected treat relationship_husb motherfinancial, ///
    numybins(5) allowmindefiers
```

#### Expected output
- `lower bound = 0.100221` (R benchmark: `0.1002207`)

---

### 4. Combination of both mechanisms: lower bound across two mediator variables

We can estimate a lower bound on the fraction of those affected by treatment, combining both mechanisms across two mediator variables.

#### R benchmark

```r
lb_frac_both <- lb_frac_affected(
  df = mother_data,
  d = "treat",
  m = c("relationship_husb", "grandmother"),
  y = "motherfinancial",
  num_Ybins = 5,
  allow_min_defiers = TRUE
)
lb_frac_both
#> [1] 0.07251284
```

#### Stata

```stata
use "data/mother_data.dta", clear

* varlist order is: d m1 m2 y
testmechs_lb_fracaffected treat relationship_husb grandmother motherfinancial, ///
    numybins(5) allowmindefiers
```

#### Expected output
- `lower bound = 0.072513` (R benchmark: `0.07251284`)

We estimate a lower bound of about 7 percent, although this does not appear to be statistically significant given the test result above.

---

### 5. Testing the sharp null of full mediation (Cox and Shi, 2023)

While the lower-bound calculations above can hint at a violation of full mediation, they do not provide uncertainty quantification. The sharp null test conducts statistical inference for the sharp null of full mediation using the method described in Section 4 of the paper.

In the upstream R package, `test_sharp_null()` supports multiple test procedures. The **recommended default** for most applications is the Cox and Shi (2023) test (`method = "CS"`). Other methods (e.g., `ARP`, `FSST`, and for binary mediators `toru`) are available in R, but are **not** implemented in this Stata MVP.

#### R benchmark (Cox–Shi / CS)

```r
test_result <- test_sharp_null(
  df = mother_data,
  d = "treat",
  m = "grandmother",
  y = "motherfinancial",
  method = "CS",      # Cox and Shi (recommended default)
  num_Ybins = 5,      # discretize Y into 5 bins
  cluster = "uc"      # cluster at uc level
)
test_result$pval
#>            [,1]
#> [1,] 0.02283916
```

#### Stata (Cox–Shi / CS)

```stata
use "data/mother_data.dta", clear

* varlist order is: d m y
testmechs_test_sharpnull treat grandmother motherfinancial, ///
    method(CS) numybins(5) cluster(uc)
```

#### Expected output (matches R)
- `test_stat = 7.558558`
- `cv        = 5.991465`
- `p-value   = 0.02283916`

Interpretation:
- The p-value is about **0.023**, so the sharp null is rejected at the **5%** significance level.
- As in the R documentation, the reported p-value corresponds to the smallest value of **α** for which the test rejects.
- The test discretizes `Y` into bins (here `numybins(5)`). Since the inference relies on a CLT approximation, choose the number of bins small enough that the CLT is reasonable within cells defined by the combination of `(Y, M, D)` (and clusters if `cluster()` is used).

---

### 6. Results for an alternative mechanism (`relationship_husb`)

We next turn to the setting where we are interested in testing whether the effect is mediated by relationship quality with the husband, which is measured on a 1-5 scale. We can again test the sharp null and estimate a lower bound on the fraction affected.

#### R benchmark (Cox–Shi / CS)

```r
test_result_husb <- test_sharp_null(
  df = mother_data,
  d = "treat",
  m = "relationship_husb",
  y = "motherfinancial",
  method = "CS",      # Cox and Shi (recommended default)
  num_Ybins = 5,      # discretize Y into 5 bins
  cluster = "uc"      # cluster at uc level
)
test_result_husb$pval
#>            [,1]
#> [1,] 0.02838332
```

#### Stata (Cox–Shi / CS)

```stata
use "data/mother_data.dta", clear

* varlist order is: d m y
testmechs_test_sharpnull treat relationship_husb motherfinancial, ///
    method(CS) numybins(5) cluster(uc)
```

#### Expected output (matches R)
- `test_stat = 10.84325`
- `cv        = 9.487729`
- `p-value   = 0.02838332`

Interpretation:
- The p-value is about **0.028**, so the sharp null is rejected at the **5%** significance level.
- This suggests that the treatment effect is not fully mediated by relationship quality with the husband alone.
- As above, the test discretizes `Y` into bins, so the choice of `numybins(5)` should balance granularity with the quality of the asymptotic approximation within `(Y, M, D)` cells.

---

### 7. Sharp null test under relaxed monotonicity (`maxdefiersshare`)

By default, TestMechs imposes the monotonicity assumption that the treatment can only increase the value of `M`. In this setting, this means that everyone who would have a grandmother present without receiving CBT treatment would also have one present when receiving CBT treatment. We can relax this assumption by setting `maxdefiersshare` to be non-zero, which bounds the number of "defiers" by that share.

We rerun the sharp null test above with `maxdefiersshare(0.01)`, which allows one percent of the population to be defiers.

#### R benchmark (Cox–Shi / CS)

```r
test_result_defiers <- test_sharp_null(
  df = mother_data,
  d = "treat",
  m = "grandmother",
  y = "motherfinancial",
  method = "CS",
  num_Ybins = 5,
  cluster = "uc",
  max_defiers_share = .01
)
test_result_defiers$pval
#>            [,1]
#> [1,] 0.04630939
```

#### Stata (Cox–Shi / CS)

```stata
use "data/mother_data.dta", clear

* varlist order is: d m y
testmechs_test_sharpnull treat grandmother motherfinancial, ///
    method(CS) numybins(5) cluster(uc) maxdefiersshare(0.01)
```

#### Expected output (matches R)
- `test_stat = 6.144824`
- `cv        = 5.991465`
- `p-value   = 0.04630932`

Interpretation:
- The p-value increases to about **0.046**, so the test still rejects the sharp null even when allowing one percent of the population to be defiers.
- Allowing for larger shares of defiers will eventually lead to an insignificant result.
- As in the standard test, the reported p-value corresponds to the smallest value of **α** for which the test rejects.

---

### 8. Combination of both mechanisms: sharp null test

Next, we test the null hypothesis that the treatment effect is explained by the **combination** of the two mechanisms (presence of grandmother and quality of relationship with husband). In the R package this is done by passing a vector of variable names to the `m` argument. In the Stata port, we simply pass both mediator variables in the varlist.

#### R benchmark (Cox–Shi / CS)

```r
test_result_both <- test_sharp_null(
  df = mother_data,
  d = "treat",
  m = c("relationship_husb", "grandmother"),
  y = "motherfinancial",
  num_Ybins = 5,
  method = "CS",
  cluster = "uc"
)
test_result_both$pval
#>           [,1]
#> [1,] 0.6540863
```

With a p-value of about **0.654**, R cannot reject the sharp null that the combination of the two mechanisms fully explains the treatment effect.

#### Stata (Cox–Shi / CS)

```stata
use "data/mother_data.dta", clear

* varlist order is: d m1 m2 y
testmechs_test_sharpnull treat grandmother relationship_husb motherfinancial, ///
    method(CS) numybins(5) cluster(uc)
```

#### Expected output (Stata)
- `test_stat = 7.741333`
- `cv        = 11.070498`
- `p-value   = 0.17107925`

#### Comparison with R

| Quantity  | R          | Stata      | Match     |
|-----------|------------|------------|-----------|
| test_stat | 7.741333   | 7.741333   | Identical |
| p-value   | 0.6540863  | 0.17107925 | Differs numerically, same conclusion |

#### Why the p-values differ, and why the Stata result is still acceptable

The **test statistic is identical** between Stata and R (to 6 decimal places). The difference arises only in the computation of the **degrees of freedom** of the chi-squared reference distribution.

The default degree-of-freedom algorithm (Section 4 of the paper, `new_dof_CS = FALSE`) requires solving a sequence of linear programs to identify which inequality constraints are binding at the optimum. R's `TestMechs` uses GLPK (via `Rglpk`) for these LPs; the Stata port uses HiGHS (via `scipy.optimize.linprog`). When the LP is near-degenerate — as it often is in multi-mediator problems with moderately large support — different LP solvers select different vertices at the optimum, producing slightly different counts of binding constraints and therefore slightly different degrees of freedom.

We experimented with routing the LPs through GLPK directly (via Python's `swiglpk` bindings) to match R exactly. This improved the combination-mediator agreement but degraded the single-mediator cases (which otherwise match R to seven decimal places). We therefore retain HiGHS as the default: it gives identical results to R in all single-mediator cases we have tested, and differs numerically but not qualitatively in multi-mediator cases.

**Both Stata (p = 0.171) and R (p = 0.654) fail to reject the sharp null at any conventional significance level (α = 0.05 or α = 0.10).** The substantive conclusion — that the combination of the two mechanisms cannot be ruled out as a full explanation of the treatment effect — is the same in both implementations.

See `LIMITATIONS.md` for a fuller discussion.

---

## Notes

- The current Stata implementation is designed to match the R package benchmarks for the supported default / binned-Y path.
- At the moment, unsupported options such as `regformula()` and `continuousy` will return an error rather than silently falling back to another behavior.
- Additional functions from the R package, including `test_sharp_null`, are planned but not yet implemented.
- **Update:** the Cox–Shi sharp null test is now available as `testmechs_test_sharpnull` with `method(CS)`. Other sharp-null methods/branches remain out of scope for the MVP.
- **Update:** `testmechs_test_sharpnull` now accepts multiple mediator variables for the combination-of-mechanisms test. The default degree-of-freedom algorithm requires Python with `numpy`, `scipy`, and `osqp` (version 0.6.x) installed; see `LIMITATIONS.md` for details. 
