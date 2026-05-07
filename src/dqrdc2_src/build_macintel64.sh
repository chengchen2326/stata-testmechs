#!/bin/bash
# Build _testmechs_dqrdc2_rank_macintel64.plugin for macOS Intel.
# Output: ../_testmechs_dqrdc2_rank_macintel64.plugin
#
# Requires: clang (Xcode command-line tools), gfortran (e.g. brew install gcc).
# Works on macOS 10.13+ on Intel (x86_64).
#
# UNTESTED: this build script has not been verified end-to-end. Once you
# have access to a macOS Intel machine, please test that the resulting
# plugin produces the same p-values as the Apple Silicon version (see
# README and section "Testing the sharp null...").

set -e
cd "$(dirname "$0")"

OUTPUT_NAME="_testmechs_dqrdc2_rank_macintel64.plugin"

echo "[1/3] Compiling Fortran sources..."
gfortran -c -O2 -arch x86_64 dqrdc2.f       -o dqrdc2.o
gfortran -c -O2 -arch x86_64 dnrm2.f        -o dnrm2.o
gfortran -c -O2 -arch x86_64 blas.f         -o blas.o
gfortran -c -O2 -arch x86_64 xerbla_stub.f  -o xerbla_stub.o

echo "[2/3] Compiling and linking plugin..."
clang -bundle -DSYSTEM=APPLEMAC -target x86_64-apple-macos10.13 \
    testmechs_dqrdc2_rank.c stplugin.c \
    dqrdc2.o dnrm2.o blas.o xerbla_stub.o \
    -lgfortran \
    -o "../${OUTPUT_NAME}"

echo "[3/3] Cleaning up..."
rm -f dqrdc2.o dnrm2.o blas.o xerbla_stub.o

echo ""
echo "Build complete: ../${OUTPUT_NAME}"
ls -la "../${OUTPUT_NAME}"
