@echo off
REM Build _testmechs_glpk_lp_win64.plugin for 64-bit Windows on x86_64.
REM Output: ..\_testmechs_glpk_lp_win64.plugin
REM
REM Prerequisites:
REM   - MSYS2 with the MinGW-w64 toolchain (https://www.msys2.org/)
REM   - GLPK static library built for MinGW-w64. Extract
REM     ..\glpk_src\glpk-5.0.tar.gz inside an MSYS2 MinGW64 shell and run:
REM       ./configure --enable-static --disable-shared --without-gmp --disable-dl
REM       make -j$(nproc)
REM     This produces ..\glpk_src\glpk-5.0\src\.libs\libglpk.a linkable by MinGW gcc.
REM
REM Run this script from a Windows command prompt with MinGW64 binaries on PATH.
REM
REM UNTESTED: this build script has not been verified end-to-end. Please
REM test the resulting plugin gives the same p-values as the macarm64
REM version on real data before shipping.

setlocal enabledelayedexpansion
cd /d "%~dp0"

set OUTPUT_NAME=_testmechs_glpk_lp_win64.plugin
set GLPK_ROOT=..\glpk_src\glpk-5.0

if not exist "%GLPK_ROOT%\src\.libs\libglpk.a" (
    echo ERROR: GLPK static library not found at %GLPK_ROOT%\src\.libs\libglpk.a
    echo Build it first ^-- see ..\glpk_src\README.md.
    exit /b 1
)

echo [1/1] Compiling and linking plugin...
gcc -shared -DSYSTEM=STWIN32 ^
    -I"%GLPK_ROOT%\src" ^
    testmechs_glpk_lp.c stplugin.c ^
    "%GLPK_ROOT%\src\.libs\libglpk.a" ^
    -lm -static-libgcc ^
    -o "..\%OUTPUT_NAME%" || goto :error

echo.
echo Build complete: ..\%OUTPUT_NAME%
dir "..\%OUTPUT_NAME%"
goto :end

:error
echo BUILD FAILED.
exit /b 1

:end
endlocal
