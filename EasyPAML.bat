@echo off
REM EasyPAML - lancador Windows. Usa o .venv criado pelo install.bat;
REM sem ele, usa o Python do sistema (instalacao antiga com pip --user).
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
    echo  Erro ao iniciar EasyPAML.
    echo  Execute install.bat se ainda nao instalou.
    pause
)
