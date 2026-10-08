@echo off
setlocal enabledelayedexpansion
chcp 65001 >nul 2>&1
cd /d "%~dp0"

echo.
echo  ============================================================
echo   EasyPAML installer for Windows
echo   Positive selection analysis with PAML/codeml
echo  ============================================================
echo.

REM ── 1. Find Python ──────────────────────────────────────────────────────────
set PYTHON=
echo [1/4] Checking Python...

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
echo  ERROR: Python not found.
echo.
echo  To fix it:
echo   1. Open https://www.python.org/downloads/
echo   2. Click "Download Python 3.x.x"
echo   3. Run the installer
echo   4. Tick "Add Python to PATH"
echo   5. Click "Install Now"
echo   6. Close this window and run INSTALL_WINDOWS.bat again
echo.
pause
exit /b 1

:python_found
for /f "tokens=*" %%V in ('!PYTHON! --version 2^>^&1') do set PY_VER=%%V
echo  OK: !PY_VER! found

REM ── 2. Minimum version (3.8) ───────────────────────────────────────────────
for /f "tokens=2 delims= " %%v in ('!PYTHON! --version 2^>^&1') do set PY_FULL=%%v
for /f "tokens=1,2 delims=." %%a in ("!PY_FULL!") do (
    set PY_MAJ=%%a
    set PY_MIN=%%b
)
if !PY_MAJ! LSS 3 (
    echo  ERROR: Python !PY_FULL! is too old. Install Python 3.8 or newer.
    pause
    exit /b 1
)
if !PY_MAJ! EQU 3 if !PY_MIN! LSS 8 (
    echo  ERROR: Python !PY_FULL! is too old. Install Python 3.8 or newer.
    pause
    exit /b 1
)

REM ── 3. .venv and dependencies ──────────────────────────────────────────────
REM Creates an isolated Python environment in .venv inside this folder. If venv
REM is unavailable, falls back to pip --user.
echo.
echo [2/4] Installing Python dependencies in .venv ...
echo  (a few minutes the first time)
echo.

if not exist ".venv\Scripts\python.exe" (
    !PYTHON! -m venv .venv
)
if not exist ".venv\Scripts\python.exe" goto :install_user

".venv\Scripts\python.exe" -m pip install --upgrade pip --quiet --disable-pip-version-check
".venv\Scripts\python.exe" -m pip install -r tools\requirements-lock.txt --disable-pip-version-check --quiet
if not errorlevel 1 (
    echo  OK: dependencies installed in .venv ^(tested versions, requirements-lock.txt^)
    goto :deps_done
)
echo  Warning: the tested versions are not available for this Python; installing newer ones
".venv\Scripts\python.exe" -m pip install -r tools\requirements.txt --disable-pip-version-check
if errorlevel 1 goto :pip_error
echo  OK: dependencies installed in .venv
goto :deps_done

:install_user
echo  Warning: could not create .venv; installing for the user (--user)
!PYTHON! -m pip install --upgrade pip --quiet --user
!PYTHON! -m pip install -r tools\requirements.txt --user
if errorlevel 1 goto :pip_error
echo  OK: dependencies installed (--user)
goto :deps_done

:pip_error
echo.
echo  ERROR installing the dependencies.
echo.
echo  Possible causes:
echo   - no internet connection
echo   - an antivirus blocking pip
echo.
echo  Try running these in this folder:
echo   !PYTHON! -m venv .venv
echo   .venv\Scripts\python.exe -m pip install -r tools\requirements.txt
echo.
pause
exit /b 1

:deps_done

REM ── 4. CODEML ───────────────────────────────────────────────────────────────
echo.
echo [3/4] Checking codeml...
if exist "bin\codeml.exe" (
    echo  OK: bin\codeml.exe found
) else (
    echo  Warning: bin\codeml.exe not found. Download PAML from
    echo  https://github.com/abacus-gene/paml/releases and copy codeml.exe to the bin folder
)

set APP_DIR=%~dp0
set APP_DIR=!APP_DIR:~0,-1!

REM ── 5. Launchers and desktop shortcut ──────────────────────────────────────
REM Written here, so the folder shows only the installer until EasyPAML is installed.
echo.
echo [4/4] Creating EasyPAML.bat, easypaml-cli.bat and a desktop shortcut...

setlocal disabledelayedexpansion
> "EasyPAML.bat" (
    echo @echo off
    echo setlocal
    echo REM Opens the EasyPAML window. Created by INSTALL_WINDOWS.bat.
    echo cd /d "%%~dp0"
    echo if exist ".venv\Scripts\python.exe" ^(
    echo     ".venv\Scripts\python.exe" src\easypaml_window.py
    echo ^) else ^(
    echo     %PYTHON% src\easypaml_window.py
    echo ^)
    echo set "EASYPAML_RC=%%errorlevel%%"
    echo REM "if errorlevel 1" misses negative codes ^(native crashes, e.g. -1073741819^),
    echo REM so compare with 0 and keep the console open to show the error.
    echo if not "%%EASYPAML_RC%%"=="0" ^(
    echo     echo.
    echo     echo  EasyPAML closed with an error ^^^(exit code %%EASYPAML_RC%%^^^).
    echo     echo  Read the messages above. If it does not start at all, run INSTALL_WINDOWS.bat again.
    echo     pause
    echo ^)
)
> "easypaml-cli.bat" (
    echo @echo off
    echo REM EasyPAML command line: easypaml-cli.bat --help. Created by INSTALL_WINDOWS.bat.
    echo if exist "%%~dp0.venv\Scripts\python.exe" ^(
    echo     "%%~dp0.venv\Scripts\python.exe" "%%~dp0src\easypaml_cli.py" %%*
    echo ^) else ^(
    echo     %PYTHON% "%%~dp0src\easypaml_cli.py" %%*
    echo ^)
)
endlocal
echo  OK: EasyPAML.bat ^(window^) and easypaml-cli.bat ^(command line^) created

REM a plain command, no script file, so no execution policy needs changing
powershell -NoProfile -Command ^
    "$s = (New-Object -ComObject WScript.Shell).CreateShortcut([Environment]::GetFolderPath('Desktop') + '\EasyPAML.lnk');" ^
    "$s.TargetPath = '!APP_DIR!\EasyPAML.bat';" ^
    "$s.WorkingDirectory = '!APP_DIR!';" ^
    "$s.Description = 'EasyPAML - positive selection analysis';" ^
    "$s.Save()" >nul 2>&1

if errorlevel 1 (
    echo  Warning: shortcut not created ^(permission denied^). Use EasyPAML.bat.
) else (
    echo  OK: "EasyPAML" shortcut created on the desktop
)

echo.
echo  ============================================================
echo   INSTALLATION COMPLETE
echo  ============================================================
echo.
echo  To open EasyPAML:
echo    - double-click "EasyPAML" on the desktop
echo    - or double-click EasyPAML.bat in this folder
echo  Command line: easypaml-cli.bat --help
echo.
echo  Example data: examples\
echo.
pause
endlocal
