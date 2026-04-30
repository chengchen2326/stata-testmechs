#!/bin/bash
# Build _testmechs_dqrdc2_rank.plugin for macOS Apple Silicon.
# Output: ../_testmechs_dqrdc2_rank.plugin

set -e
cd "$(dirname "$0")"

echo "Compiling Fortran sources..."
gfortran -c -O2 dqrdc2.f       -o dqrdc2.o
gfortran -c -O2 dnrm2.f        -o dnrm2.o
gfortran -c -O2 blas.f         -o blas.o
gfortran -c -O2 xerbla_stub.f  -o xerbla_stub.o

echo "Compiling and linking plugin..."
clang -bundle -DSYSTEM=APPLEMAC -target arm64-apple-macos11 \
    testmechs_dqrdc2_rank.c stplugin.c \
    dqrdc2.o dnrm2.o blas.o xerbla_stub.o \
    -lgfortran -L/opt/homebrew/opt/gcc/lib/gcc/current \
    -o ../_testmechs_dqrdc2_rank.plugin

echo "Cleaning up..."
rm -f dqrdc2.o dnrm2.o blas.o xerbla_stub.o

echo "Build complete: ../_testmechs_dqrdc2_rank.plugin"
ls -la ../_testmechs_dqrdc2_rank.plugin
