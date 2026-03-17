@echo off
REM ============================================================================
REM         TESTE DE INSTALAÇÃO - Simula novo usuário
REM ============================================================================
REM Este arquivo verifica se tudo está pronto para funcionar

chcp 65001 >nul
cls

echo.
echo ╔════════════════════════════════════════════════════════════════╗
echo ║         TESTE DE INSTALAÇÃO - EASYPML                         ║
echo ║         (Simula novo usuário verificando ambiente)            ║
echo ╚════════════════════════════════════════════════════════════════╝
echo.

setlocal enabledelayedexpansion
set PASS=0
set FAIL=0

REM Teste 1: SETUP.bat existe
echo [1/7] Verificando SETUP.bat...
if exist "SETUP.bat" (
    echo ✓ SETUP.bat encontrado
    set /a PASS+=1
) else (
    echo ✗ SETUP.bat NÃO encontrado!
    set /a FAIL+=1
)

REM Teste 2: requirements.txt existe
echo [2/7] Verificando requirements.txt...
if exist "requirements.txt" (
    echo ✓ requirements.txt encontrado
    set /a PASS+=1
) else (
    echo ✗ requirements.txt NÃO encontrado!
    set /a FAIL+=1
)

REM Teste 3: EasyPAML.py existe
echo [3/7] Verificando EasyPAML.py...
if exist "EasyPAML.py" (
    echo ✓ EasyPAML.py encontrado
    set /a PASS+=1
) else (
    echo ✗ EasyPAML.py NÃO encontrado!
    set /a FAIL+=1
)

REM Teste 4: src/ existe
echo [4/7] Verificando pasta src/...
if exist "src\" (
    echo ✓ Pasta src/ encontrada
    set /a PASS+=1
) else (
    echo ✗ Pasta src/ NÃO encontrada!
    set /a FAIL+=1
)

REM Teste 5: codeml.exe existe
echo [5/7] Verificando bin/codeml.exe...
if exist "bin\codeml.exe" (
    echo ✓ CODEML encontrado
    set /a PASS+=1
) else (
    echo ✗ CODEML NÃO encontrado em bin/!
    set /a FAIL+=1
)

REM Teste 6: exemplos_teste/ existe
echo [6/7] Verificando exemplos_teste/...
if exist "exemplos_teste\" (
    echo ✓ Pasta exemplos_teste/ encontrada
    set /a PASS+=1
) else (
    echo ✗ Pasta exemplos_teste/ NÃO encontrada!
    set /a FAIL+=1
)

REM Teste 7: Documentação existe
echo [7/7] Verificando documentação...
if exist "COMECE_AQUI.txt" (
    echo ✓ Documentação encontrada (COMECE_AQUI.txt)
    set /a PASS+=1
) else (
    echo ✗ Documentação NÃO encontrada!
    set /a FAIL+=1
)

echo.
echo ╔════════════════════════════════════════════════════════════════╗
if !FAIL! equ 0 (
    echo ║            ✓ TODOS OS TESTES PASSARAM!                      ║
    echo ║         Projeto pronto para distribuição!                   ║
) else (
    echo ║           ✗ %FAIL% teste(s) falharam!                       ║
    echo ║         Verifique os arquivos acima                         ║
)
echo ╚════════════════════════════════════════════════════════════════╝
echo.

echo Resumo:
echo   ✓ Testes passados: %PASS%/7
echo   ✗ Testes falhados: %FAIL%/7
echo.

if !FAIL! equ 0 (
    echo Próximos passos:
    echo   1. Comprima a pasta EasyPML para distribuição
    echo   2. Novo usuário executa SETUP.bat
    echo   3. Tudo funciona automaticamente!
    echo.
) else (
    echo ATENÇÃO: Corrija os problemas antes de distribuir!
    echo.
)

pause
endlocal
