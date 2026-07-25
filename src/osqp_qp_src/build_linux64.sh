#!/bin/bash
# Build _honestosqp_linux64.plugin for 64-bit Linux on x86_64.
# Output: ../_honestosqp_linux64.plugin
#
# Prerequisites: OSQP static library built for Linux x86_64.
# On Debian/Ubuntu:
#   sudo apt install gcc make cmake
# Then follow ../osqp_qp_src/README.md to build libosqp.a.
#
# UNTESTED: this build script has not been verified end-to-end. Please
# test the resulting plugin gives the same p-values as the macarm64
# version on real data before shipping.

set -e
cd "$(dirname "$0")"

OUTPUT_NAME="_honestosqp_linux64.plugin"
OSQP_ROOT="../osqp_src/osqp"

if [ ! -f "${OSQP_ROOT}/build/out/libosqp.a" ]; then
    echo "ERROR: OSQP static library not found at ${OSQP_ROOT}/build/out/libosqp.a"
    echo "Build it first — see ../osqp_qp_src/README.md."
    exit 1
fi

echo "[1/1] Compiling and linking plugin..."
gcc -Wall -O3 -shared -fPIC -DSYSTEM=OPUNIX \
    -I"${OSQP_ROOT}/include" \
    -o "../${OUTPUT_NAME}" \
    honestosqp.c stplugin.c \
    "${OSQP_ROOT}/build/out/libosqp.a"

echo ""
echo "Build complete: ../${OUTPUT_NAME}"
ls -la "../${OUTPUT_NAME}"
