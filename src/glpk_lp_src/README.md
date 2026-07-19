# GLPK LP plugin

Source code for `_testmechs_glpk_lp.plugin`, the Stata plugin that solves the linear programs used by `testmechs_test_sharpnull` (Cox–Shi test, default `new_dof_CS=FALSE` branch). It replaces the previous Python + swiglpk path.

## Files

- `testmechs_glpk_lp.c` — plugin entry point; reads three Stata matrices, calls GLPK's C API, writes two Stata scalars
- `stplugin.c`, `stplugin.h` — StataCorp's plugin glue (identical to the copy in `../dqrdc2_src/`)
- `build_macarm64.sh` — build for macOS Apple Silicon (tested)

## Interface

The plugin expects the following Stata objects to exist and be readable:

- Input matrices: `__tm_c` (n×1 or 1×n), `__tm_Aeq` (m×n), `__tm_beq` (m×1 or 1×m)
- Output scalars: `__tm_lp_fun`, `__tm_lp_ok`

The Mata caller (in `../testmechs_test_sharpnull.ado`) does:

```stata
plugin call _testmechs_glpk_lp, "__tm_c" "__tm_Aeq" "__tm_beq" "__tm_lp_fun" "__tm_lp_ok"
```

The plugin solves

    min c' x   s.t.  A_eq x = b_eq,  x >= 0

using GLPK's simplex method (`glp_simplex`), which matches the solver R's `Rglpk_solve_LP` uses by default. On success `__tm_lp_ok` is set to 1 and `__tm_lp_fun` to the optimal objective value; on failure `__tm_lp_ok` is 0 and `__tm_lp_fun` is missing.

## Building

The plugin statically links GLPK 5.0, so you need to build the GLPK static library first. See `../glpk_src/README.md`.

Once `../glpk_src/glpk-5.0/src/.libs/libglpk.a` exists, run the platform-specific build script from this directory:

```bash
./build_macarm64.sh
```

The resulting plugin is written to `../_testmechs_glpk_lp_<platform>.plugin`.
