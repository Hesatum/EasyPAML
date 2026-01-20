@echo off
chcp 65001 >nul
setlocal enabledelayedexpansion

REM ============================================================================
REM               EasyPAML - LAUNCHER
REM ============================================================================

cd /d "%~dp0"

REM Verificar se está instalado
if not exist "requirements.txt" (
    echo.
    echo ✗ ERRO: requirements.txt não encontrado!
    echo Execute SETUP.bat primeiro para instalar.
    echo.
    pause
    exit /b 1
)

REM Verificar Python
python --version >nul 2>&1
if errorlevel 1 (
    echo.
    echo ✗ ERRO: Python não encontrado!
    echo Execute SETUP.bat para configurar tudo.
    echo.
    pause
    exit /b 1
)

REM Verificar CODEML
codeml.exe --help >nul 2>&1
if errorlevel 1 (
    echo.
    echo ⚠ AVISO: CODEML não encontrado no PATH
    echo O programa pode não funcionar corretamente.
    echo.
    echo Solução: Execute SETUP.bat novamente como Administrador
    echo.
)

REM Iniciar EasyPML
echo.
echo Iniciando EasyPAML...
echo.

python EasyPAML.py

if errorlevel 1 (
    echo.
    echo ✗ Erro ao executar EasyPAML
    echo Verifique os logs acima
    echo.
    pause
)

endlocal
