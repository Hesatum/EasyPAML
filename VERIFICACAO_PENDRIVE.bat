@echo off
setlocal enabledelayedexpansion
chcp 65001 >nul
cls

REM ============================================================================
REM         VERIFICAÇÃO FINAL - EasyPML antes de copiar para pendrive
REM ============================================================================

echo.
echo ╔════════════════════════════════════════════════════════════════╗
echo ║       VERIFICAÇÃO FINAL - EASYPML                             ║
echo ║  (Execute ANTES de copiar para pendrive)                      ║
echo ╚════════════════════════════════════════════════════════════════╝
echo.

setlocal enabledelayedexpansion
set TOTAL=0
set PASS=0

REM ============= TESTE 1: Arquivos principais =============
echo [TESTE 1] Verificando arquivos principais...
set /a TOTAL+=1

if exist "SETUP.bat" (
    set /a PASS+=1
    echo   ✓ SETUP.bat encontrado
) else (
    echo   ✗ ERRO: SETUP.bat NÃO encontrado!
)

if exist "RUN_EASYPML.bat" (
    set /a PASS+=1
    echo   ✓ RUN_EASYPML.bat encontrado
) else (
    echo   ✗ ERRO: RUN_EASYPML.bat NÃO encontrado!
)

if exist "EasyPAML.py" (
    set /a PASS+=1
    echo   ✓ EasyPAML.py encontrado
) else (
    echo   ✗ ERRO: EasyPAML.py NÃO encontrado!
)

if exist "requirements.txt" (
    set /a PASS+=1
    echo   ✓ requirements.txt encontrado
) else (
    echo   ✗ ERRO: requirements.txt NÃO encontrado!
)

REM ============= TESTE 2: Pastas essenciais =============
echo.
echo [TESTE 2] Verificando pastas essenciais...
set /a TOTAL+=1

if exist "src\" (
    set /a PASS+=1
    echo   ✓ Pasta src/ encontrada
) else (
    echo   ✗ ERRO: Pasta src/ NÃO encontrada!
)

if exist "bin\" (
    set /a PASS+=1
    echo   ✓ Pasta bin/ encontrada
) else (
    echo   ✗ ERRO: Pasta bin/ NÃO encontrada!
)

if exist "exemplos_teste\" (
    set /a PASS+=1
    echo   ✓ Pasta exemplos_teste/ encontrada
) else (
    echo   ✗ ERRO: Pasta exemplos_teste/ NÃO encontrada!
)

REM ============= TESTE 3: CODEML =============
echo.
echo [TESTE 3] Verificando CODEML...
set /a TOTAL+=1

if exist "bin\codeml.exe" (
    set /a PASS+=1
    echo   ✓ CODEML (codeml.exe) encontrado
    for /f "tokens=*" %%i in ('powershell -Command ""{0:N0}"" -f (Get-Item bin\codeml.exe).Length') do (
        echo   └─ Tamanho: %%i bytes
    )
) else (
    echo   ✗ ERRO: CODEML NÃO encontrado em bin/codeml.exe!
)

REM ============= TESTE 4: Documentação =============
echo.
echo [TESTE 4] Verificando documentação...
set /a TOTAL+=1

if exist "COMECE_AQUI.txt" (
    set /a PASS+=1
    echo   ✓ COMECE_AQUI.txt encontrado
) else (
    echo   ✗ AVISO: COMECE_AQUI.txt NÃO encontrado
)

if exist "INSTALACAO.txt" (
    set /a PASS+=1
    echo   ✓ INSTALACAO.txt encontrado
) else (
    echo   ✗ AVISO: INSTALACAO.txt NÃO encontrado
)

if exist "GUIA_USUARIO.txt" (
    set /a PASS+=1
    echo   ✓ GUIA_USUARIO.txt encontrado
) else (
    echo   ✗ AVISO: GUIA_USUARIO.txt NÃO encontrado
)

if exist "README.md" (
    set /a PASS+=1
    echo   ✓ README.md encontrado
) else (
    echo   ✗ AVISO: README.md NÃO encontrado
)

REM ============= TESTE 5: Python =============
echo.
echo [TESTE 5] Verificando Python...
set /a TOTAL+=1

python --version >nul 2>&1
if !errorlevel! equ 0 (
    set /a PASS+=1
    for /f "tokens=*" %%i in ('python --version') do (
        echo   ✓ %%i
    )
) else (
    echo   ✗ AVISO: Python não encontrado no PATH
    echo   └─ Novo usuário precisará instalar Python
)

REM ============= TESTE 6: Conteúdo de arquivos críticos =============
echo.
echo [TESTE 6] Verificando conteúdo dos arquivos...
set /a TOTAL+=1

findstr /c:"SETUP" SETUP.bat >nul 2>&1
if !errorlevel! equ 0 (
    set /a PASS+=1
    echo   ✓ SETUP.bat tem conteúdo válido
) else (
    echo   ✗ ERRO: SETUP.bat vazio ou corrompido!
)

REM ============= TESTE 7: Estrutura src/ =============
echo.
echo [TESTE 7] Verificando estrutura src/...
set /a TOTAL+=1

if exist "src\backend\" (
    set /a PASS+=1
    echo   ✓ src/backend/ encontrado
) else (
    echo   ✗ AVISO: src/backend/ NÃO encontrado
)

if exist "src\gui\" (
    set /a PASS+=1
    echo   ✓ src/gui/ encontrado
) else (
    echo   ✗ AVISO: src/gui/ NÃO encontrado
)

REM ============= TESTE 8: Arquivos de dados =============
echo.
echo [TESTE 8] Verificando exemplos_teste/...
set /a TOTAL+=1

setlocal disabledelayedexpansion
for /d %%i in (exemplos_teste\*) do (
    set /a PASS+=1
    echo   ✓ Subpasta encontrada: %%~nxi
)
endlocal enabledelayedexpansion

REM ============= RESULTADO FINAL =============
echo.
echo.
echo ╔════════════════════════════════════════════════════════════════╗

if !PASS! geq 16 (
    echo ║        ✓ VERIFICAÇÃO COMPLETA - PRONTO PARA PENDRIVE!       ║
    echo ║                   Todos os testes passaram                  ║
) else if !PASS! geq 12 (
    echo ║        ⚠ VERIFICAÇÃO PARCIAL - AVISOS ENCONTRADOS           ║
    echo ║         Revise os itens acima antes de copiar               ║
) else (
    echo ║         ✗ VERIFICAÇÃO FALHOU - NÃO COPIE AINDA!            ║
    echo ║         Corrija os erros acima antes de prosseguir          ║
)

echo ╚════════════════════════════════════════════════════════════════╝
echo.

echo Resultado: %PASS% de 19 testes passaram
echo.

if !PASS! geq 16 (
    echo ✅ PRÓXIMAS ETAPAS:
    echo    1. Copie a pasta EasyPML para o pendrive
    echo    2. Novo usuário com Python instalado executa SETUP.bat
    echo    3. Tudo funciona automaticamente!
    echo.
    echo ✅ IMPORTANTE no novo computador:
    echo    • Python 3.8+ DEVE estar instalado
    echo    • Python DEVE ter "Add Python to PATH" marcado
    echo    • Após SETUP.bat, FECHAR e ABRIR novo terminal
    echo.
) else (
    echo ⚠ ATENÇÃO: Corrija os problemas acima antes de distribuir!
    echo.
)

pause
endlocal
