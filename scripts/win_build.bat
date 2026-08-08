@echo off
REM Developer helper: configure + build chemcad with MSVC + Ninja on Windows.
REM Usage: scripts\win_build.bat [prefix]   (prefix = conda-style RDKit root)
setlocal
call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat" >nul || exit /b 1
set "PATH=C:\cmake\cmake-3.30.5-windows-x86_64\bin;%PATH%"
set "PREFIX=%~1"
if "%PREFIX%"=="" set "PREFIX=%LOCALAPPDATA%\chemcad-deps\rdkit\Library"
set "NINJA=%VSINSTALLDIR%Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja\ninja.exe"
if not exist "%NINJA%" for /f "delims=" %%i in ('where ninja') do set "NINJA=%%i"
cd /d "%~dp0.."
cmake -S . -B build-win -G Ninja -DCMAKE_BUILD_TYPE=Release ^
      -DCMAKE_MAKE_PROGRAM="%NINJA%" -DCMAKE_PREFIX_PATH="%PREFIX%" || exit /b 1
cmake --build build-win || exit /b 1
