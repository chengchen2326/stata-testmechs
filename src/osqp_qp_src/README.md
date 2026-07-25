# OSQP quadratic programming plugin

The QP solver used by `testmechs_test_sharpnull` is [OSQP](https://osqp.org/) — the Operator Splitting Quadratic Program solver. It is bundled here as a Stata plugin so the user does not need Python or the OSQP Python package.

## Files in this directory

Source (adapted from [HonestDiD](https://github.com/mcaceresb/stata-honestdid)):

- `honestosqp.c`, `honestosqp.h` — plugin entry point (adapted from HonestDiD; see modifications below)
- `sf_printf.c`, `sf_printf.h` — Stata-friendly printf helpers (from HonestDiD, unchanged)
- `stplugin.c`, `stplugin.h` — StataCorp's plugin glue (from HonestDiD, unchanged)
- `osqp.mata` — Mata interface to the plugin (adapted from HonestDiD; see modifications below)

## Files elsewhere in the repo

- `../_honestosqp_macarm64.plugin` — precompiled plugin for macOS Apple Silicon (verified)
- `../_honestosqp_macintel64.plugin` — precompiled for macOS Intel (untested; upstream binary)
- `../_honestosqp_linux64.plugin` — precompiled for Linux x86_64 (untested; upstream binary)
- `../_honestosqp_win64.plugin` — precompiled for Windows x86_64 (untested; upstream binary)

## Provenance and modifications

The starting point for these files is the OSQP-plugin machinery in Mauricio Cáceres Bravo, Ashesh Rambachan, and Jonathan Roth's [HonestDiD Stata package](https://github.com/mcaceresb/stata-honestdid) (2025). Mauricio pointed us to it and authorised the reuse; Jon Roth is a co-author of both packages.

We modified two files to expose OSQP's convergence tolerances (`eps_abs`, `eps_rel`) to the Mata layer, because our port of TestMechs requires tighter tolerances (1e-8) than HonestDiD's default (1e-5) to match the results of R's `TestMechs` package.

### `honestosqp.c` (three modifications)

1. Added `eps_abs` and `eps_rel` to the local variable declarations.
2. Added two `SF_scal_use` calls next to the existing `__honestosqp_max_iter` read, populating `__honestosqp_eps_abs` and `__honestosqp_eps_rel` from Stata scalars.
3. Changed the hard-coded `settings->eps_abs = 1e-5` and `settings->eps_rel = 1e-5` to fall back to `1e-5` only if the corresponding Stata scalar is missing.

### `osqp.mata` (two modifications)

1. In `OSQP_setup`, initialise `__honestosqp_eps_abs` and `__honestosqp_eps_rel` the same way `__honestosqp_max_iter` is initialised (keep the user's value if set, otherwise leave missing).
2. In `OSQP_cleanup`, clear those two scalars alongside `__honestosqp_max_iter`.

Callers can now do

```mata
st_numscalar("__honestosqp_eps_abs", 1e-8)
st_numscalar("__honestosqp_eps_rel", 1e-8)
result = OSQP(P, q, A, u, l, 1)
```

to obtain 1e-8 tolerance. `testmechs__lp_solve` (in `../testmechs_test_sharpnull.ado`) does this so that the QP objective values match R's `TestMechs` to eight decimal places.

Upstream users of HonestDiD are unaffected: the modifications are strictly additive (new scalars, missing → old default).

## Public API (used by the ado)

We solve

    min_x   0.5 x' P x + q' x
    s.t.    l <= A x <= u

by calling

```mata
result = OSQP(P, q, A, u, l, 1)
```

`result` is a `struct OSQP_workspace_abridged` with fields `rc`, `info_status`, `info_obj_val`, `solution_x`.

The plugin binary is called `honestosqp_plugin`. Our ado registers it under that exact name (a `capture program drop honestosqp_plugin` guards against conflicts if the user has both `honestdid` and `testmechs` installed).

## Build

The plugin statically links OSQP 0.6.3. Follow Mauricio's `stata-honestdid/src/compile.sh` to build the OSQP static library — the short version is:

```bash
git clone --recursive https://github.com/osqp/osqp
cd osqp
git reset --hard v0.6.3
git submodule update --init --recursive
mkdir -p build && cd build
cmake -DCMAKE_OSX_ARCHITECTURES=arm64 -DBUILD_SHARED_LIBS=OFF ..
make -j$(sysctl -n hw.ncpu)
```

Then compile the plugin:

```bash
clang -Wall -O3 -bundle -DSYSTEM=APPLEMAC -target arm64-apple-macos11 \
    -I../osqp_src/osqp/include \
    -o ../_honestosqp_macarm64.plugin \
    honestosqp.c stplugin.c \
    ../osqp_src/osqp/build/out/libosqp.a
```

## Acknowledgement

We are grateful to Mauricio Cáceres Bravo for writing the OSQP plugin and Mata wrapper, and for suggesting we reuse it here.
