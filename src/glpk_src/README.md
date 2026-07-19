# GLPK 5.0 (bundled for the LP plugin)

`_testmechs_glpk_lp.plugin` is a Stata plugin that solves linear programs by calling GLPK's C API directly. GLPK is statically linked into the plugin so end users do not need to install GLPK on their machines.

This directory holds the GLPK source and (once you build it) the static library that the plugin links against.

## Files in this directory (in the repo)

- `glpk-5.0.tar.gz` — GLPK 5.0 source distribution (from https://ftp.gnu.org/gnu/glpk/), bundled for reproducibility

## Files created by the build (not committed)

- `glpk-5.0/` — extracted source tree
- `glpk-5.0/src/.libs/libglpk.a` — the static library the plugin links against

## Build instructions (macOS or Linux)

From this directory:

```bash
# Extract the source
tar -xzf glpk-5.0.tar.gz
cd glpk-5.0

# Configure: static library only, minimal dependencies
./configure --enable-static --disable-shared --without-gmp --disable-dl

# Build
make -j$(sysctl -n hw.ncpu 2>/dev/null || nproc)

# Verify
ls -la src/.libs/libglpk.a
```

Then go into `../glpk_lp_src/` and run the appropriate `build_*.sh` for your platform.

## Why GLPK 5.0

R's `Rglpk` package is the reference implementation for TestMechs' LP solver. `Rglpk` links dynamically against whichever GLPK the user has installed, but GLPK 5.0 (released 2020) is the current version and by far the most common one distributed by package managers (Homebrew, Debian, etc.). Bundling GLPK 5.0 gives us the closest possible match to what R users experience.

## License

GLPK is licensed under GPL-3. Because the plugin links GLPK statically, the resulting `.plugin` binary is a combined work distributed under GPL-3. The rest of the `testmechs` Stata package remains under its own license (see repository root).
