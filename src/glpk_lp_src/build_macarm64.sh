#!/bin/bash
# Build _testmechs_glpk_lp_macarm64.plugin for macOS Apple Silicon.
# Output: ../_testmechs_glpk_lp_macarm64.plugin
#
# Prerequisites: the GLPK static library must have been built first.
# See ../glpk_src/README.md for instructions.
#
# Requires: clang (Xcode command-line tools).

set -e
cd "$(dirname "$0")"

OUTPUT_NAME="_testmechs_glpk_lp_macarm64.plugin"
GLPK_ROOT="../glpk_src/glpk-5.0"

# Sanity check: the GLPK static library must exist.
if [ ! -f "${GLPK_ROOT}/src/.libs/libglpk.a" ]; then
    echo "ERROR: GLPK static library not found at ${GLPK_ROOT}/src/.libs/libglpk.a"
    echo "Build it first — see ../glpk_src/README.md."
    exit 1
fi

echo "[1/1] Compiling and linking plugin..."
clang -bundle -DSYSTEM=APPLEMAC -target arm64-apple-macos11 \
    -I"${GLPK_ROOT}/src" \
    testmechs_glpk_lp.c stplugin.c \
    "${GLPK_ROOT}/src/.libs/libglpk.a" \
    -lm \
    -o "../${OUTPUT_NAME}"

echo ""
echo "Build complete: ../${OUTPUT_NAME}"
ls -la "../${OUTPUT_NAME}"
