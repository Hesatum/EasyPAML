@echo off
REM EasyPAML launcher for Windows. Uses the .venv created by install.bat;
REM without it, uses the system Python.
cd /d "%~dp0"
if exist ".venv\Scripts\python.exe" (
    ".venv\Scripts\python.exe" EasyPAML.py
    goto :check
)
where python >nul 2>&1
if not errorlevel 1 (
    python EasyPAML.py
    goto :check
)
py EasyPAML.py

:check
if errorlevel 1 (
    echo.
    echo  EasyPAML could not start.
    echo  Run install.bat if you have not installed it yet.
    pause
)
