# stata-testmechs (MVP)

This repository contains an **MVP Stata translation** of key functionality from Jon Roth & Soonwoo Kwon's **TestMechs** R package (paper: *Testing Mechanisms*).

## Status (MVP scope)

### Implemented
- ✅ `testmechs_lb_fracaffected`  
  Stata translation of the R function `lb_frac_affected()`
- ✅ `testmechs_test_sharpnull`  
  Stata translation of the R function `test_sharp_null()` for `method = "CS"` (Cox and Shi, 2023). Supports single mediators and combinations of mediators.

### Currently supported features
- Positional varlist input:
  - **`d m y`** for a single mediator
  - **`d m1 m2 y`** for a combination of both mechanisms (two-mediator input)
- `atgroup(#)`
- `numybins(#)`
- `maxdefiersshare(#)`
- `allowmindefiers`
- `reg_formula("...")` — OLS with controls and 2SLS IV syntax `(endog = instr)`. Perfect-instrument cases (instrument collinear with treatment) fall back to OLS to match R's `fixest` behavior. Cluster-robust inference via analytic influence functions (sandwich-form `M · x · e · n_reg` for OLS, `M · xhat · e_iv · n_reg` for 2SLS).
- Multi-valued discrete mediators in the default binned-Y path
- Lower bounds combining both mechanisms across two mediator variables
- Cox–Shi sharp null test with cluster-robust standard errors

### Not yet implemented
The following are **not yet implemented** and will return an error:
- 🚫 `continuousy`
- 🚫 `returnmindefiers`
- 🚫 Sharp-null methods other than `CS` (e.g., `ARP`, `FSST`, `toru`)
- 🚫 Specialised binary-M code path (see "Notes" below)
- 🚫 `fixest`-style fixed-effect syntax (`| interviewer`) and IV-with-FE combined syntax
- 🚫 `bounds_ade_ats` and `partial_density_plot`

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

### Dependencies

`testmechs_test_sharpnull` uses three bundled Stata plugins for all its numerical work. It requires:

- **Stata 16+**

That's it — no Python installation, no external solver libraries. Rank computation (via `dqrdc2`), the linear programs (via GLPK 5.0), and the quadratic program (via OSQP 0.6.3) all run inside statically-linked Stata plugins.

### Bundled plugins

`testmechs_test_sharpnull` ships with three precompiled Stata plugins.

**`_testmechs_dqrdc2_rank.plugin`** calls R's `dqrdc2` Fortran routine for computing matrix rank. This is necessary because R's `qr(M, tol)$rank` uses a 1995 modification of LINPACK's `dqrdc` written by Ross Ihaka specifically for R, with a custom column-pivoting strategy that gives different results from generic SVD-based rank on near-rank-deficient matrices. Without this plugin, the Stata p-values would differ from R for several test cases. See `src/dqrdc2_src/` for the Fortran/C sources.

**`_testmechs_glpk_lp.plugin`** solves the linear programs used by the Cox–Shi test. It statically links GLPK 5.0 (the same LP solver used by R's `Rglpk_solve_LP`), so users do not need to install GLPK or `swiglpk`. See `src/glpk_lp_src/` for the C source and `src/glpk_src/` for the bundled GLPK tarball.

**`honestosqp_plugin.plugin`** solves the quadratic program used by the Cox–Shi test. It statically links OSQP 0.6.3. The plugin is a light modification of Mauricio Cáceres Bravo's OSQP plugin from HonestDiD, extended to expose the `eps_abs` and `eps_rel` convergence tolerances so we can set them to 1e-8 (to match R's TestMechs) instead of HonestDiD's default 1e-5. See `src/osqp_qp_src/` for the C source and modification details.

All three plugins are pre-shipped for **macOS Apple Silicon**, **macOS Intel (x86_64)**, **Linux x86_64**, and **Windows x86_64**. Linux and Windows binaries are built automatically on GitHub Actions from `.github/workflows/build-plugins-linux.yml` and `.github/workflows/build-plugins-windows.yml`; macOS binaries are built locally per the build scripts under `src/dqrdc2_src/`, `src/glpk_lp_src/`, and `src/osqp_qp_src/`. Only the macOS Apple Silicon binaries have been run end-to-end against the R baseline; the other three platforms are architecture-verified but not yet Stata-verified.

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
- `test_stat = 10.843249`
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

#### Expected output (matches R)
- `test_stat = 7.741333`
- `cv        = 18.307038`
- `p-value   = 0.65408650`

Interpretation:
- With a p-value of about **0.654**, the test cannot reject the sharp null that the combination of the two mechanisms fully explains the treatment effect.
- This matches the R benchmark exactly (R: 0.6540863, Stata: 0.65408650).

---

## Non-experimental setting

The examples above focus on data from a randomized controlled trial, where treatment `D` is randomly assigned. TestMechs assumes by default that we have an RCT, and treatment effects are estimated by comparing means for the treated and control group. However, TestMechs can also be applied in settings where we have conditional randomization given covariates or an instrumental variable for the treatment. Both `testmechs_test_sharpnull` and `testmechs_lb_fracaffected` accept the `reg_formula("...")` option, which allows the researcher to provide a regression formula (OLS or 2SLS IV) to estimate treatment effects after adjusting linearly for observable characteristics.

The formula string uses R-style syntax (with a leading tilde), mirroring the R package's convention:
- **OLS with controls:** `reg_formula("~ treat + age_baseline + edu_mo_baseline + wealth_baseline")`
- **2SLS IV:** `reg_formula("~ age_baseline + edu_mo_baseline + wealth_baseline + (treat = iv)")` — the parenthesised `(endog = instr)` clause specifies the endogenous variable and its instrument.

While adjusting for covariates is not necessary in our running example (which is an RCT), we can still adjust for covariates to increase precision. Below we illustrate using the baseline covariates `age_baseline`, `edu_mo_baseline`, and `wealth_baseline`. For brevity we focus on `testmechs_test_sharpnull`; the same `reg_formula` argument works with `testmechs_lb_fracaffected`.

### 9. Sharp null with regression-adjusted probabilities (OLS controls)

The key new argument is `reg_formula("~ treat + age_baseline + edu_mo_baseline + wealth_baseline")`. This says that we should estimate the effect of the treatment (`treat`) on `Y` and `M` (or functions thereof) by running a regression with the treatment variable and baseline controls on the right-hand side.

#### R benchmark (Cox–Shi / CS with controls)

```r
test_result_gm_ols <- test_sharp_null(
  df = mother_data,
  d = "treat",
  m = "grandmother",
  y = "motherfinancial",
  reg_formula = "~ treat + age_baseline + edu_mo_baseline + wealth_baseline",
  method = "CS",
  num_Ybins = 5,
  cluster = "uc"
)
test_result_gm_ols$pval
#>            [,1]
#> [1,] 0.03105106
```

#### Stata (Cox–Shi / CS with controls)

```stata
use "data/mother_data.dta", clear

* varlist order is: d m y
testmechs_test_sharpnull treat grandmother motherfinancial, ///
    method(CS) numybins(5) cluster(uc) ///
    reg_formula("~ treat + age_baseline + edu_mo_baseline + wealth_baseline")
```

#### Expected output (matches R)
- `test_stat = 6.944244`
- `cv        = 5.991465`
- `p-value   = 0.03105106`

Interpretation:
- The p-value is about **0.031**, so the sharp null is rejected at the **5%** significance level after adjusting for baseline covariates.
- Compared to the unadjusted result in example 5 (p = 0.023), the adjusted p-value is somewhat larger — the covariates absorb some of the variation that made the unadjusted test more powerful here.
- Inference uses cluster-robust analytic influence functions: `IF_i = M[j,:] · x_i · e_i · n_reg` where `M = (X'X)^{-1}`, mirroring R's `sandwich::estfun · t(bread)` computation.

---

### 10. Sharp null with an instrumental variable (2SLS IV)

We can also use `reg_formula` to estimate treatment effects using instrumental variables. To illustrate, we construct a noisy instrument for `treat` (analogous to the R README example) and then pass IV syntax to `reg_formula`.

The IV portion is written as `(endog = instr)` inside the formula, matching Stata's `ivregress` convention: everything outside the parentheses is a control, and the parenthesised clause specifies the endogenous variable and its instrument.

#### R benchmark (Cox–Shi / CS with IV)

```r
set.seed(0)
mother_data$iv <- mother_data$treat + rnorm(n = length(mother_data$treat), sd = 0.1)  # iv = treat + noise

test_result_gm_iv <- test_sharp_null(
  df = mother_data,
  d = "treat",
  m = "grandmother",
  y = "motherfinancial",
  reg_formula = "~ age_baseline + edu_mo_baseline + wealth_baseline | treat ~ iv",
  method = "CS",
  num_Ybins = 5,
  cluster = "uc"
)
test_result_gm_iv$pval
#>           [,1]
#> [1,] 0.0311524
```

#### Stata (Cox–Shi / CS with IV)

```stata
use "data/mother_data.dta", clear

* Construct a noisy instrument (analogous to R's set.seed(0) + rnorm(sd=0.1))
set seed 0
gen double iv = treat + rnormal(0, 0.1)

* varlist order is: d m y
* IV syntax inside reg_formula uses (endog = instr)
testmechs_test_sharpnull treat grandmother motherfinancial, ///
    method(CS) numybins(5) cluster(uc) ///
    reg_formula("~ age_baseline + edu_mo_baseline + wealth_baseline + (treat = iv)")
```

#### Notes on expected output

- **Numerical result depends on random-number generator.** Stata's `rnormal(0, 0.1)` and R's `rnorm(sd = 0.1)` use different random-number streams, so the constructed `iv` variable will differ between the two runs even with the same seed. The Stata p-value will therefore not equal the R benchmark exactly; only the analytical procedure is the same. To get an exact numerical match against R, users would need to construct `iv` in R and export it (e.g. via `write_dta()`) or use a common noise source.
- The 2SLS IV point estimate is computed via `ivregress 2sls`, and the per-observation influence function is built via the sandwich form `IF_i = M[j,:] · x̂_i · e_iv,i · n_reg` where `M = e(V)/rmse²` is the IV bread, `x̂_i` is the first-stage fitted regressor row, and `e_iv,i` is the IV residual. This exactly reproduces R's `sandwich::estfun(feols) · t(sandwich::bread(feols))` to machine precision on the same data — see `tests/verify/` for the element-wise verification against R.
- **Perfect-instrument case:** if the instrument is collinear with the endogenous variable (Stata returns `r(481)`), the helper falls back to OLS treating the endog as a regular regressor. This mirrors R's `fixest` behavior in the same degenerate case, so the two produce identical point estimates and IFs when the "instrument" is a copy of the treatment.

---

## Notes on numerical agreement with R

All Cox–Shi p-values reported above match the R benchmark to at least 7 decimal places for the non-IV examples. Achieving this agreement required two design choices:

1. **GLPK as the LP solver.** Earlier versions of this package used HiGHS (via `scipy.optimize.linprog`). HiGHS gave correct test statistics but selected slightly different LP vertices in near-degenerate cases, which propagated to small differences in the constraint-binding count and therefore the chi-squared degrees of freedom. Switching to GLPK (the same solver R uses) eliminates this source of disagreement. GLPK is now embedded directly in the bundled plugin `_testmechs_glpk_lp.plugin` — no `swiglpk` or system GLPK required.

2. **R's `dqrdc2` for rank computation.** R's `qr(M, tol)$rank` uses `dqrdc2.f` — a 1995 modification by Ross Ihaka of LINPACK's `dqrdc`, available only in R's source tree. It applies a custom column-pivoting strategy that diverges from generic SVD-based or LAPACK QR rank computations on near-rank-deficient matrices. The difference is rare but consequential: for the `relationship_husb` test, SVD-based rank gave dof = 5 where `dqrdc2` gives dof = 4, leading to p-values of 0.0546 vs 0.0284. We therefore ship a precompiled Stata plugin (`src/_testmechs_dqrdc2_rank.plugin`) that calls `dqrdc2` directly from Mata, matching R's rank exactly.

For **`reg_formula` cases**, the OLS-with-controls path also matches R exactly (see example 9). The IV path reproduces R's sandwich-form influence function per observation to machine precision (verified in `tests/verify/`), so end-to-end p-values match R exactly for `K ≥ 3` (multi-valued or combined mediators). For `K = 2` (binary M) IV cases, however, R auto-dispatches to a specialised `test_sharp_null_binary_m` code path (with a different Cox–Shi variant tuned for the binary-M no-nuisance case). The Stata port only implements the general CS code path, so **binary-M IV p-values may differ from R's binary-M-path p-values**. The general CS output remains a valid CS test; porting the binary-M specialisation is on the roadmap pending priority.

---

## Notes

- The current Stata implementation is designed to match the R package benchmarks for the supported default / binned-Y path and the `reg_formula` OLS / 2SLS IV paths.
- Unsupported options such as `continuousy` will return an error rather than silently falling back to another behavior.
- Additional sharp-null methods from the R package (`ARP`, `FSST`, `toru`) are not yet implemented and will return an error.
- The specialised binary-M code path (used by R for K = 2 IV cases) is not ported; users get valid general-CS output on those cases but the exact p-value will differ from R.
- The default degree-of-freedom algorithm (`new_dof_CS = FALSE`) is implemented; the alternative (`new_dof_CS = TRUE`) Cox–Shi formulas remain on the roadmap.
- `fixest`-style fixed-effect syntax (`| interviewer`) and IV-with-FE combined syntax are not supported; use `i.interviewer` explicitly on the RHS as a workaround for fixed effects.
