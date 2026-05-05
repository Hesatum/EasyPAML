@echo off
setlocal enabledelayedexpansion
chcp 65001 >nul 2>&1
cd /d "%~dp0"

echo.
echo  ============================================================
echo   EasyPAML - Instalador Automatico para Windows
echo   Analise de Selecao Positiva com CODEML/PAML
echo  ============================================================
echo.

REM ── 1. Encontrar Python ─────────────────────────────────────────────────────
set PYTHON=
echo [1/4] Verificando Python...

python --version >nul 2>&1
if not errorlevel 1 (
    set PYTHON=python
    goto :python_found
)

py --version >nul 2>&1
if not errorlevel 1 (
    set PYTHON=py
    goto :python_found
)

python3 --version >nul 2>&1
if not errorlevel 1 (
    set PYTHON=python3
    goto :python_found
)

echo.
echo  ERRO: Python nao encontrado!
echo.
echo  Solucao passo a passo:
echo   1. Abra o navegador e acesse: https://www.python.org/downloads/
echo   2. Clique no botao amarelo "Download Python 3.x.x"
echo   3. Abra o instalador baixado
echo   4. IMPORTANTE: marque a caixa "Add Python to PATH"
echo   5. Clique em "Install Now"
echo   6. Feche este janela e abra o install.bat novamente
echo.
pause
exit /b 1

:python_found
for /f "tokens=*" %%V in ('!PYTHON! --version 2^>^&1') do set PY_VER=%%V
echo  OK: !PY_VER! encontrado

REM ── 2. Verificar versao minima (3.8) ───────────────────────────────────────
for /f "tokens=2 delims= " %%v in ('!PYTHON! --version 2^>^&1') do set PY_FULL=%%v
for /f "tokens=1,2 delims=." %%a in ("!PY_FULL!") do (
    set PY_MAJ=%%a
    set PY_MIN=%%b
)
if !PY_MAJ! LSS 3 (
    echo  ERRO: Python !PY_FULL! e muito antigo. Instale Python 3.8 ou superior.
    pause
    exit /b 1
)
if !PY_MAJ! EQU 3 if !PY_MIN! LSS 8 (
    echo  ERRO: Python !PY_FULL! e muito antigo. Instale Python 3.8 ou superior.
    pause
    exit /b 1
)

REM ── 3. Instalar dependencias ────────────────────────────────────────────────
echo.
echo [2/4] Instalando dependencias Python...
echo  (Isso pode levar alguns minutos na primeira vez)
echo.

!PYTHON! -m pip install --upgrade pip --quiet --user
if errorlevel 1 (
    echo  Aviso: nao foi possivel atualizar pip, tentando continuar...
)

!PYTHON! -m pip install -r requirements.txt --user
if errorlevel 1 (
    echo.
    echo  ERRO ao instalar dependencias!
    echo.
    echo  Possiveis causas:
    echo   - Sem conexao com internet
    echo   - Antivirus bloqueando pip
    echo   - requirements.txt corrompido
    echo.
    echo  Tente executar manualmente:
    echo   !PYTHON! -m pip install -r requirements.txt --user
    echo.
    pause
    exit /b 1
)
echo  OK: Todas as dependencias instaladas

REM ── 4. Criar launcher EasyPAML.bat ─────────────────────────────────────────
echo.
echo [3/4] Criando launcher...

set APP_DIR=%~dp0
set APP_DIR=!APP_DIR:~0,-1!

(
    echo @echo off
    echo cd /d "!APP_DIR!"
    echo !PYTHON! EasyPAML.py
    echo if errorlevel 1 pause
) > "!APP_DIR!\EasyPAML.bat"

echo  OK: Launcher criado (EasyPAML.bat)

REM ── 5. Atalho na Area de Trabalho ───────────────────────────────────────────
echo.
echo [4/4] Criando atalho na area de trabalho...

powershell -NoProfile -ExecutionPolicy Bypass -Command ^
    "$s = (New-Object -ComObject WScript.Shell).CreateShortcut([Environment]::GetFolderPath('Desktop') + '\EasyPAML.lnk');" ^
    "$s.TargetPath = '!APP_DIR!\EasyPAML.bat';" ^
    "$s.WorkingDirectory = '!APP_DIR!';" ^
    "$s.Description = 'EasyPAML - Analise de Selecao Positiva';" ^
    "$s.Save()" >nul 2>&1

if errorlevel 1 (
    echo  Aviso: atalho nao criado (permissao negada). Use EasyPAML.bat diretamente.
) else (
    echo  OK: Atalho "EasyPAML" criado na area de trabalho
)

echo.
echo  ============================================================
echo   INSTALACAO CONCLUIDA!
echo  ============================================================
echo.
echo  Para usar o EasyPAML:
echo    - Duplo-clique em "EasyPAML" na area de trabalho
echo    - OU duplo-clique em EasyPAML.bat nesta pasta
echo    - OU execute: !PYTHON! EasyPAML.py
echo.
echo  Dados de exemplo em: exemplos_teste\
echo.
pause
endlocal
