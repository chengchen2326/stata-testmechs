@echo off
REM Build _testmechs_dqrdc2_rank_win64.plugin for 64-bit Windows on x86_64.
REM Output: ..\_testmechs_dqrdc2_rank_win64.plugin
REM
REM Requires MSYS2 with MinGW-w64 toolchain (https://www.msys2.org/).
REM Run this script from an MSYS2 MinGW64 shell or a Windows command prompt
REM with MinGW64 binaries on PATH. Need: gcc, gfortran (from mingw-w64-x86_64-gcc-fortran package).
REM
REM UNTESTED: this build script has not been verified end-to-end. Once you
REM have access to a Windows x86_64 machine, please test that the resulting
REM plugin produces the same p-values as the Apple Silicon version (see
REM README and section "Testing the sharp null...").

setlocal enabledelayedexpansion
cd /d "%~dp0"

set OUTPUT_NAME=_testmechs_dqrdc2_rank_win64.plugin

echo [1/3] Compiling Fortran sources...
gfortran -c -O2 dqrdc2.f       -o dqrdc2.o      || goto :error
gfortran -c -O2 dnrm2.f        -o dnrm2.o       || goto :error
gfortran -c -O2 blas.f         -o blas.o        || goto :error
gfortran -c -O2 xerbla_stub.f  -o xerbla_stub.o || goto :error

echo [2/3] Compiling and linking plugin...
gcc -shared -DSYSTEM=STWIN32 ^
    testmechs_dqrdc2_rank.c stplugin.c ^
    dqrdc2.o dnrm2.o blas.o xerbla_stub.o ^
    -lgfortran -static-libgcc -static-libgfortran ^
    -o "..\%OUTPUT_NAME%" || goto :error

echo [3/3] Cleaning up...
del /q dqrdc2.o dnrm2.o blas.o xerbla_stub.o

echo.
echo Build complete: ..\%OUTPUT_NAME%
dir "..\%OUTPUT_NAME%"
goto :end

:error
echo BUILD FAILED at step.
exit /b 1

:end
endlocal
