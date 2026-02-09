# stata-testmechs

Stata package (MVP) translating key functionality from Jon Roth's TestMechs R package.

## Status
- MVP in progress: default-path implementation only (reduced feature set vs R).
- Added `testmechs_lb_fracaffected` as an MVP translation of `lb_frac_affected.R` for the binary-mediator default path.

## Reference
The original R source is stored in `r_reference/` for translation and validation.

## Usage
From Stata in repo root:

```stata
adopath ++ .
* var order is: treatment mediator outcome

testmechs_lb_fracaffected d m y
testmechs_lb_fracaffected d m y, atgroup(1)
testmechs_lb_fracaffected d m y, numybins(5)
```

## Testing
Run the dedicated MVP test do-file from the repo root:

```stata
do tests/test_lb_fracaffected.do
```
