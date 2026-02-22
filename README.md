# stata-testmechs (MVP)

This repository contains an **MVP Stata translation** of key functionality from Jon Roth & Soonwoo Kwon’s **TestMechs** R package (paper: *Testing Mechanisms*).

## Status (MVP scope)
- ✅ Implemented: `testmechs_lb_fracaffected` (MVP version of R `lb_frac_affected()`), **default path only**
- ✅ Supported inputs:
  - Positional varlist: **`d m y`** (all numeric)
  - Options: `atgroup(#)` and `numybins(#)`
- 🚫 Not yet implemented (will error in MVP):
  - `regformula()`, `continuousy`
  - `maxdefiersshare()>0`, `allowmindefiers`, `returnmindefiers`
  - Additional inference/plotting functions (e.g., `test_sharp_null`)

The original R source is kept under `r_reference/` for translation and validation.
