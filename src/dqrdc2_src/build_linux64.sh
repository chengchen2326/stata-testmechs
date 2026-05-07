#!/bin/bash
# Build _testmechs_dqrdc2_rank_linux64.plugin for 64-bit Linux on x86_64.
# Output: ../_testmechs_dqrdc2_rank_linux64.plugin
#
# Requires: gcc and gfortran. On Debian/Ubuntu:
#   sudo apt install gcc gfortran
#
# UNTESTED: this build script has not been verified end-to-end. Once you
# have access to a Linux x86_64 machine, please test that the resulting
# plugin produces the same p-values as the Apple Silicon version (see
# README and section "Testing the sharp null...").

set -e
cd "$(dirname "$0")"

OUTPUT_NAME="_testmechs_dqrdc2_rank_linux64.plugin"

echo "[1/3] Compiling Fortran sources..."
gfortran -c -O2 -fPIC dqrdc2.f       -o dqrdc2.o
gfortran -c -O2 -fPIC dnrm2.f        -o dnrm2.o
gfortran -c -O2 -fPIC blas.f         -o blas.o
gfortran -c -O2 -fPIC xerbla_stub.f  -o xerbla_stub.o

echo "[2/3] Compiling and linking plugin..."
gcc -shared -DSYSTEM=OPUNIX -fPIC \
    testmechs_dqrdc2_rank.c stplugin.c \
    dqrdc2.o dnrm2.o blas.o xerbla_stub.o \
    -lgfortran \
    -o "../${OUTPUT_NAME}"

echo "[3/3] Cleaning up..."
rm -f dqrdc2.o dnrm2.o blas.o xerbla_stub.o

echo ""
echo "Build complete: ../${OUTPUT_NAME}"
ls -la "../${OUTPUT_NAME}"
