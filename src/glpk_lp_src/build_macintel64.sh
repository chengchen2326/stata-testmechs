#!/bin/bash
# Build _testmechs_glpk_lp_macintel64.plugin for macOS Intel (x86_64).
# Output: ../_testmechs_glpk_lp_macintel64.plugin
#
# Prerequisites: the GLPK static library must have been built first for
# x86_64. On an Intel Mac, running the standard configure + make from
# ../glpk_src/glpk-5.0/ should produce the right library automatically.
# On Apple Silicon you would need to cross-compile — pass
# CFLAGS="-arch x86_64" LDFLAGS="-arch x86_64" to configure, or run
# under Rosetta.
#
# Requires: clang (Xcode command-line tools).
#
# UNTESTED: this build script has not been verified end-to-end. Please
# test the resulting plugin gives the same p-values as the macarm64
# version on real data before shipping.

set -e
cd "$(dirname "$0")"

OUTPUT_NAME="_testmechs_glpk_lp_macintel64.plugin"
GLPK_ROOT="../glpk_src/glpk-5.0"

if [ ! -f "${GLPK_ROOT}/src/.libs/libglpk.a" ]; then
    echo "ERROR: GLPK static library not found at ${GLPK_ROOT}/src/.libs/libglpk.a"
    echo "Build it first — see ../glpk_src/README.md."
    exit 1
fi

echo "[1/1] Compiling and linking plugin..."
clang -bundle -DSYSTEM=APPLEMAC -target x86_64-apple-macos10.13 \
    -I"${GLPK_ROOT}/src" \
    testmechs_glpk_lp.c stplugin.c \
    "${GLPK_ROOT}/src/.libs/libglpk.a" \
    -lm \
    -o "../${OUTPUT_NAME}"

echo ""
echo "Build complete: ../${OUTPUT_NAME}"
ls -la "../${OUTPUT_NAME}"
