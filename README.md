# stata-testmechs (MVP)

This repository contains an **MVP Stata translation** of key functionality from Jon Roth & Soonwoo Kwon’s **TestMechs** R package (paper: *Testing Mechanisms*).

## Status (MVP scope)
- ✅ Implemented: `testmechs_lb_fracaffected` (MVP version of R `lb_frac_affected()`), **default path only**
- ✅ Supported inputs:
  - Positional varlist: **`d m y`** (all numeric)
  - Options: `atgroup(#)`, `numybins(#)`, `maxdefiersshare(#)`, `allowmindefiers`
  - Multi-valued discrete mediators (default binned-Y path)
- 🚫 Not yet implemented (will error):
  - `regformula()`, `continuousy`, `returnmindefiers`
  - Additional inference/plotting functions (e.g., `test_sharp_null`)

The original R source is kept under `r_reference/` for translation and validation.

## Data (for replication)

This repo includes two Stata datasets converted from the TestMechs R package data (`.rda`) for convenience:

- `baranov_data.dta`: raw dataset converted from the R package data object `baranov_data`.
- `mother_data.dta`: analysis dataset corresponding to the experimental sample used in the README example, created by restricting to `THP_sample == 1`.

Notes:
- These `.dta` files are provided for quick replication in Stata.
- They are derived from the original R package data; see `r_reference/` for the upstream source.

## Installation

Install using `net install` from GitHub:

```stata
cap noi net uninstall testmechs
net install testmechs, from("https://raw.githubusercontent.com/chengchen2326/stata-testmechs/main") replace
```
## Examples

### Lower bound on fraction affected (default/recommended path)

This replicates the README example from the TestMechs R package for the experimental sample.

```stata
use "mother_data.dta", clear

* varlist order is: d m y
testmechs_lb_fracaffected treat relationship_husb motherfinancial, numybins(5)
```
Expected output (MVP):

lower bound = 0.100221 (R benchmark: 0.1002207)