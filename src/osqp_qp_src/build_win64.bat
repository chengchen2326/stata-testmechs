@echo off
REM Build _honestosqp_win64.plugin for 64-bit Windows on x86_64.
REM Output: ..\_honestosqp_win64.plugin
REM
REM Prerequisites:
REM   - MSYS2 with the MinGW-w64 toolchain (https://www.msys2.org/)
REM   - OSQP static library built for MinGW-w64. Extract and build inside
REM     an MSYS2 MinGW64 shell following the instructions in
REM     ..\osqp_qp_src\README.md. The resulting libosqp.a should be at
REM     ..\osqp_src\osqp\build\out\libosqp.a and linkable by MinGW gcc.
REM
REM Run this script from a Windows command prompt with MinGW64 binaries on PATH.
REM
REM UNTESTED: this build script has not been verified end-to-end. Please
REM test the resulting plugin gives the same p-values as the macarm64
REM version on real data before shipping.

setlocal enabledelayedexpansion
cd /d "%~dp0"

set OUTPUT_NAME=_honestosqp_win64.plugin
set OSQP_ROOT=..\osqp_src\osqp

if not exist "%OSQP_ROOT%\build\out\libosqp.a" (
    echo ERROR: OSQP static library not found at %OSQP_ROOT%\build\out\libosqp.a
    echo Build it first ^-- see ..\osqp_qp_src\README.md.
    exit /b 1
)

echo [1/1] Compiling and linking plugin...
gcc -Wall -O3 -shared -DSYSTEM=STWIN32 ^
    -I"%OSQP_ROOT%\include" ^
    -o "..\%OUTPUT_NAME%" ^
    honestosqp.c stplugin.c ^
    "%OSQP_ROOT%\build\out\libosqp.a" ^
    -static-libgcc || goto :error

echo.
echo Build complete: ..\%OUTPUT_NAME%
dir "..\%OUTPUT_NAME%"
goto :end

:error
echo BUILD FAILED.
exit /b 1

:end
endlocal
