#!/bin/bash
# Build _testmechs_glpk_lp_linux64.plugin for 64-bit Linux on x86_64.
# Output: ../_testmechs_glpk_lp_linux64.plugin
#
# Prerequisites: the GLPK static library must have been built first.
# On Debian/Ubuntu:
#   sudo apt install gcc make
# Then from ../glpk_src/glpk-5.0/, run:
#   ./configure --enable-static --disable-shared --without-gmp --disable-dl
#   make -j$(nproc)
# See ../glpk_src/README.md.
#
# UNTESTED: this build script has not been verified end-to-end. Please
# test the resulting plugin gives the same p-values as the macarm64
# version on real data before shipping.

set -e
cd "$(dirname "$0")"

OUTPUT_NAME="_testmechs_glpk_lp_linux64.plugin"
GLPK_ROOT="../glpk_src/glpk-5.0"

if [ ! -f "${GLPK_ROOT}/src/.libs/libglpk.a" ]; then
    echo "ERROR: GLPK static library not found at ${GLPK_ROOT}/src/.libs/libglpk.a"
    echo "Build it first — see ../glpk_src/README.md."
    exit 1
fi

echo "[1/1] Compiling and linking plugin..."
gcc -shared -DSYSTEM=OPUNIX -fPIC \
    -I"${GLPK_ROOT}/src" \
    testmechs_glpk_lp.c stplugin.c \
    "${GLPK_ROOT}/src/.libs/libglpk.a" \
    -lm \
    -o "../${OUTPUT_NAME}"

echo ""
echo "Build complete: ../${OUTPUT_NAME}"
ls -la "../${OUTPUT_NAME}"
