#!/bin/bash
# Build _honestosqp_macarm64.plugin for macOS Apple Silicon.
# Output: ../_honestosqp_macarm64.plugin
#
# Prerequisites: the OSQP static library must have been built first.
# See ../osqp_qp_src/README.md for full instructions. Short version:
#   git clone --recursive https://github.com/osqp/osqp ../osqp_src/osqp
#   cd ../osqp_src/osqp
#   git reset --hard v0.6.3
#   git submodule update --init --recursive
#   mkdir build && cd build
#   cmake -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
#         -DCMAKE_OSX_ARCHITECTURES=arm64 \
#         -DBUILD_SHARED_LIBS=OFF ..
#   make -j$(sysctl -n hw.ncpu)
#
# Requires: clang (Xcode command-line tools).

set -e
cd "$(dirname "$0")"

OUTPUT_NAME="_honestosqp_macarm64.plugin"
OSQP_ROOT="../osqp_src/osqp"

if [ ! -f "${OSQP_ROOT}/build/out/libosqp.a" ]; then
    echo "ERROR: OSQP static library not found at ${OSQP_ROOT}/build/out/libosqp.a"
    echo "Build it first — see ../osqp_qp_src/README.md."
    exit 1
fi

echo "[1/1] Compiling and linking plugin..."
clang -Wall -O3 -bundle -DSYSTEM=APPLEMAC -target arm64-apple-macos11 \
    -I"${OSQP_ROOT}/include" \
    -o "../${OUTPUT_NAME}" \
    honestosqp.c stplugin.c \
    "${OSQP_ROOT}/build/out/libosqp.a"

echo ""
echo "Build complete: ../${OUTPUT_NAME}"
ls -la "../${OUTPUT_NAME}"
