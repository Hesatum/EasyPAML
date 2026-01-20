@echo off
setlocal enabledelayedexpansion
chcp 65001 >nul
cls

REM ============================================================================
REM               EasyPAML - INSTALADOR AUTOMÁTICO
REM ============================================================================

echo.
echo ╔════════════════════════════════════════════════════════════════╗
echo ║           EasyPAML - INSTALADOR AUTOMÁTICO                     ║
echo ║        Seleção Positiva em Sequências Genômicas               ║
echo ╚════════════════════════════════════════════════════════════════╝
echo.

REM Verificar se Python está instalado
echo [1/6] Verificando Python...
python --version >nul 2>&1
if errorlevel 1 (
    echo.
    echo ✗ ERRO: Python não encontrado!
    echo.
    echo Solução: Acesse https://www.python.org/downloads/
    echo   1. Baixe Python 3.8 ou superior
    echo   2. Durante instalação, MARQUE "Add Python to PATH"
    echo   3. Execute este arquivo novamente
    echo.
    pause
    exit /b 1
) else (
    for /f "tokens=*" %%i in ('python --version') do set PYTHON_VERSION=%%i
    echo ✓ !PYTHON_VERSION! encontrado
)

echo.
echo [2/6] Atualizando pip...
python -m pip install --upgrade pip
if errorlevel 1 (
    echo ✗ Erro ao atualizar pip
    pause
    exit /b 1
) else (
    echo ✓ Pip atualizado com sucesso
)

echo.
echo [3/6] Instalando dependências...
python -m pip install -r requirements.txt
if errorlevel 1 (
    echo.
    echo ✗ Erro ao instalar dependências
    echo Tente executar novamente ou verifique sua conexão de internet
    pause
    exit /b 1
) else (
    echo ✓ Dependências instaladas com sucesso
)

echo.
echo [4/6] Adicionando CODEML ao PATH do Windows...
REM Adicionar bin/ ao PATH permanentemente
for /f "tokens=2*" %%A in ('reg query "HKLM\SYSTEM\CurrentControlSet\Control\Session Manager\Environment" /v PATH ^| findstr /i path') do (
    set CURRENT_PATH=%%B
)

setx PATH "%CD%\bin;!CURRENT_PATH!"
echo ✓ CODEML adicionado ao PATH

echo.
echo [5/6] Configurando variáveis de ambiente...
REM Criar variáveis de ambiente
setx EASYPML_HOME "%CD%"
setx EASYPML_LANG PT
setx PYTHONIOENCODING utf-8
echo ✓ Variáveis de ambiente configuradas

echo.
echo [6/6] Criando atalho de inicialização...
REM Criar atalho no Desktop
powershell -Command ^
    "$WshShell = New-Object -ComObject WScript.Shell; " ^
    "$Shortcut = $WshShell.CreateShortCut('%USERPROFILE%\Desktop\EasyPAML.lnk'); " ^
    "$Shortcut.TargetPath = '%CD%\RUN_EASYPML.bat'; " ^
    "$Shortcut.WorkingDirectory = '%CD%'; " ^
    "$Shortcut.Description = 'EasyPML - Análise de Seleção Positiva'; " ^
    "$Shortcut.Save()"

echo ✓ Atalho criado no Desktop

echo.
echo ╔════════════════════════════════════════════════════════════════╗
echo ║                  INSTALAÇÃO COMPLETA!                         ║
echo ╚════════════════════════════════════════════════════════════════╝
echo.
echo ✓ Python configurado
echo ✓ CODEML adicionado ao PATH
echo ✓ Dependências instaladas
echo ✓ Variáveis de ambiente definidas
echo ✓ Atalho criado no Desktop
echo.
echo Para iniciar EasyPAML:
echo   • Duplo-clique em "EasyPAML" no Desktop
echo   OU
echo   • Execute: python EasyPAML.py
echo.
echo Próximos passos:
echo   1. Feche todos os terminais e abra novamente
echo   2. Prepare seus dados (arquivos FASTA)
echo   3. Consulte o arquivo GUIA_USUARIO.txt para instruções
echo   4. Teste com exemplos em pasta "exemplos_teste"
echo.
pause
endlocal
