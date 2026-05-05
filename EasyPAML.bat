@echo off
cd /d "%~dp0"
python EasyPAML.py
if errorlevel 1 (
    echo.
    echo  Erro ao iniciar EasyPAML.
    echo  Execute install.bat se ainda nao instalou.
    pause
)
