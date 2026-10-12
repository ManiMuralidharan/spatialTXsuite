@echo off
setlocal
rem Run from the folder this file is in (the project root)
cd /d "%~dp0"

rem Find Rscript: on PATH, else the newest R under Program Files
set "RSCRIPT="
where Rscript >nul 2>nul && set "RSCRIPT=Rscript"
if not defined RSCRIPT for /d %%D in ("%ProgramFiles%\R\R-*") do if exist "%%D\bin\Rscript.exe" set "RSCRIPT=%%D\bin\Rscript.exe"
if not defined RSCRIPT for /d %%D in ("%LOCALAPPDATA%\Programs\R\R-*") do if exist "%%D\bin\Rscript.exe" set "RSCRIPT=%%D\bin\Rscript.exe"
if not defined RSCRIPT (
  echo Rscript was not found. Install R from https://cran.r-project.org/bin/windows/base/
  echo or add the folder containing Rscript.exe to PATH.
  pause
  exit /b 1
)
echo Using: %RSCRIPT%
"%RSCRIPT%" tests\smoke_test.R
pause
