#!/bin/bash
# Build _testmechs_dqrdc2_rank_macarm64.plugin for macOS Apple Silicon.
# Output: ../_testmechs_dqrdc2_rank_macarm64.plugin
#
# Requires: clang (Xcode command-line tools), gfortran (e.g. brew install gcc).
# Works on macOS 11+ on Apple Silicon (arm64).

set -e
cd "$(dirname "$0")"

OUTPUT_NAME="_testmechs_dqrdc2_rank_macarm64.plugin"

echo "[1/3] Compiling Fortran sources..."
gfortran -c -O2 dqrdc2.f       -o dqrdc2.o
gfortran -c -O2 dnrm2.f        -o dnrm2.o
gfortran -c -O2 blas.f         -o blas.o
gfortran -c -O2 xerbla_stub.f  -o xerbla_stub.o

echo "[2/3] Compiling and linking plugin..."
clang -bundle -DSYSTEM=APPLEMAC -target arm64-apple-macos11 \
    testmechs_dqrdc2_rank.c stplugin.c \
    dqrdc2.o dnrm2.o blas.o xerbla_stub.o \
    -lgfortran -L/opt/homebrew/opt/gcc/lib/gcc/current \
    -o "../${OUTPUT_NAME}"

echo "[3/3] Cleaning up..."
rm -f dqrdc2.o dnrm2.o blas.o xerbla_stub.o

echo ""
echo "Build complete: ../${OUTPUT_NAME}"
ls -la "../${OUTPUT_NAME}"
