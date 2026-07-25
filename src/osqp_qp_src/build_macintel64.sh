#!/bin/bash
# Build _honestosqp_macintel64.plugin for macOS Intel (x86_64).
# Output: ../_honestosqp_macintel64.plugin
#
# Prerequisites: OSQP static library built for x86_64. On an Intel Mac,
# CMake will pick x86_64 automatically. On Apple Silicon you must cross-compile
# by adding -DCMAKE_OSX_ARCHITECTURES=x86_64 to the cmake command.
# See ../osqp_qp_src/README.md for full build instructions.
#
# Requires: clang (Xcode command-line tools).
#
# UNTESTED: this build script has not been verified end-to-end. Please
# test the resulting plugin gives the same p-values as the macarm64
# version on real data before shipping.

set -e
cd "$(dirname "$0")"

OUTPUT_NAME="_honestosqp_macintel64.plugin"
OSQP_ROOT="../osqp_src/osqp"

if [ ! -f "${OSQP_ROOT}/build/out/libosqp.a" ]; then
    echo "ERROR: OSQP static library not found at ${OSQP_ROOT}/build/out/libosqp.a"
    echo "Build it first — see ../osqp_qp_src/README.md."
    exit 1
fi

echo "[1/1] Compiling and linking plugin..."
clang -Wall -O3 -bundle -DSYSTEM=APPLEMAC -target x86_64-apple-macos10.13 \
    -I"${OSQP_ROOT}/include" \
    -o "../${OUTPUT_NAME}" \
    honestosqp.c stplugin.c \
    "${OSQP_ROOT}/build/out/libosqp.a"

echo ""
echo "Build complete: ../${OUTPUT_NAME}"
ls -la "../${OUTPUT_NAME}"
