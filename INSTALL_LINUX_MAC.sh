#!/usr/bin/env bash
# EasyPAML installer for Linux and macOS.
#
# Creates an isolated Python environment in .venv/ inside this folder and
# installs the dependencies there, without --user and without touching the
# system Python (works with PEP 668 on Ubuntu 23.04 and later). Safe to run
# again.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [ -t 1 ]; then
    GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; BOLD='\033[1m'; NC='\033[0m'
else
    GREEN=''; YELLOW=''; RED=''; BOLD=''; NC=''
fi
ok()   { echo -e "${GREEN} OK:${NC} $*"; }
warn() { echo -e "${YELLOW} WARNING:${NC} $*"; }
err()  { echo -e "${RED} ERROR:${NC} $*"; }
cmd()  { echo -e "      ${BOLD}$*${NC}"; }

OS="$(uname -s)"
HAS_APT=0; command -v apt-get >/dev/null 2>&1 && HAS_APT=1
HAS_DNF=0; command -v dnf >/dev/null 2>&1 && HAS_DNF=1

# suggested command to install system packages
pkg_hint() {
    local deb="$1" rpm="$2" brew="$3"
    if [ "$HAS_APT" -eq 1 ]; then cmd "sudo apt update && sudo apt install -y $deb"
    elif [ "$HAS_DNF" -eq 1 ]; then cmd "sudo dnf install -y $rpm"
    elif [ "$OS" = "Darwin" ]; then cmd "brew install $brew"
    else cmd "install: $deb"
    fi
}

echo ""
echo " ============================================================"
echo "  EasyPAML installer (Linux / macOS)"
echo " ============================================================"
echo ""

# ── 1. Python 3.8+ ───────────────────────────────────────────────────────────
echo "[1/5] Checking Python..."
PYTHON=""
for c in python3 python3.13 python3.12 python3.11 python3.10 python3.9 python3.8 python; do
    if command -v "$c" >/dev/null 2>&1 && \
       "$c" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 8) else 1)' 2>/dev/null; then
        PYTHON="$c"; break
    fi
done
if [ -z "$PYTHON" ]; then
    err "Python 3.8 or newer not found. Install it with:"
    pkg_hint "python3" "python3" "python"
    exit 1
fi
ok "$("$PYTHON" --version 2>&1) ($PYTHON)"

# ── 2. venv and tkinter ──────────────────────────────────────────────────────
echo ""
echo "[2/5] Checking venv and tkinter..."
PYVER="$("$PYTHON" -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')"
MISSING_TK=0
if [ ! -x ".venv/bin/python" ]; then
    if ! "$PYTHON" -c 'import venv, ensurepip' >/dev/null 2>&1; then
        err "The Python 'venv' module is not installed (needed to install the dependencies)."
        echo "      Run this command, then run ./INSTALL_LINUX_MAC.sh again:"
        pkg_hint "python3-venv python3-tk" "python3-tkinter" "python-tk@$PYVER"
        exit 1
    fi
fi
ok "venv available"
if ! "$PYTHON" -c 'import tkinter' >/dev/null 2>&1; then
    MISSING_TK=1
    warn "tkinter not installed: the EasyPAML window needs it (the command line works without it)."
    echo "      To install it:"
    pkg_hint "python3-tk" "python3-tkinter" "python-tk@$PYVER"
else
    ok "tkinter available"
fi

# ── 3. .venv and dependencies ────────────────────────────────────────────────
echo ""
echo "[3/5] Installing dependencies in .venv/ (a few minutes the first time)..."
if [ ! -x ".venv/bin/python" ]; then
    if ! "$PYTHON" -m venv .venv; then
        err "Could not create .venv/. Install venv and run again:"
        pkg_hint "python3-venv" "python3" "python"
        exit 1
    fi
fi
VENV_PY="$SCRIPT_DIR/.venv/bin/python"
"$VENV_PY" -m pip install --upgrade pip --quiet --disable-pip-version-check || \
    warn "could not update pip; continuing with the current version"
# the exact versions EasyPAML was tested with; the minimum versions when they have
# no package for this Python
if "$VENV_PY" -m pip install -r tools/requirements-lock.txt --disable-pip-version-check --quiet; then
    ok "Dependencies installed in .venv/ (tested versions, requirements-lock.txt)"
elif "$VENV_PY" -m pip install -r tools/requirements.txt --disable-pip-version-check; then
    warn "the tested versions are not available for this Python; newer ones were installed (requirements.txt)"
    ok "Dependencies installed in .venv/"
else
    err "Installing the dependencies failed (see the message above). Common causes: no internet, a proxy."
    echo "      To try again: ./INSTALL_LINUX_MAC.sh"
    exit 1
fi

# ── 4. CODEML ────────────────────────────────────────────────────────────────
echo ""
echo "[4/5] Checking codeml (PAML)..."
CODEML_OK=0
if [ -x "$SCRIPT_DIR/bin/codeml" ]; then
    ok "codeml in the project: bin/codeml"; CODEML_OK=1
elif command -v codeml >/dev/null 2>&1; then
    ok "system codeml: $(command -v codeml)"; CODEML_OK=1
else
    warn "codeml not found."
    if [ "$HAS_APT" -eq 1 ]; then
        echo "      Installing the 'paml' package (may ask for your administrator password)..."
        if sudo apt-get install -y paml; then CODEML_OK=1; ok "PAML installed (apt)"; fi
    elif [ "$HAS_DNF" -eq 1 ]; then
        if sudo dnf install -y paml; then CODEML_OK=1; ok "PAML installed (dnf)"; fi
    elif [ "$OS" = "Darwin" ] && command -v brew >/dev/null 2>&1; then
        if brew install brewsci/bio/paml; then CODEML_OK=1; ok "PAML installed (brew)"; fi
    fi
    if [ "$CODEML_OK" -eq 0 ] && [ "$OS" = "Linux" ] && [ "$(uname -m)" = "x86_64" ]; then
        echo "      Trying the official PAML binary (GitHub abacus-gene/paml)..."
        PAML_URL="https://github.com/abacus-gene/paml/releases/download/v4.10.10/paml-4.10.10-linux-x86_64.tar.gz"
        TMP_DIR="$(mktemp -d)"
        if (command -v curl >/dev/null 2>&1 && curl -fsSL "$PAML_URL" -o "$TMP_DIR/paml.tgz") || \
           (command -v wget >/dev/null 2>&1 && wget -q "$PAML_URL" -O "$TMP_DIR/paml.tgz"); then
            tar -xzf "$TMP_DIR/paml.tgz" -C "$TMP_DIR" 2>/dev/null || true
            BIN="$(find "$TMP_DIR" -name codeml -type f | head -1)"
            if [ -n "$BIN" ]; then
                mkdir -p "$SCRIPT_DIR/bin" && cp "$BIN" "$SCRIPT_DIR/bin/codeml" && \
                    chmod +x "$SCRIPT_DIR/bin/codeml" && CODEML_OK=1 && ok "codeml copied to bin/codeml"
            fi
        fi
        rm -rf "$TMP_DIR"
    fi
    if [ "$CODEML_OK" -eq 0 ]; then
        warn "codeml was not installed. EasyPAML opens, but needs it to run analyses:"
        pkg_hint "paml" "paml" "brewsci/bio/paml"
    fi
fi

# ── 5. Launchers ─────────────────────────────────────────────────────────────
# written here, so the folder shows only the installer until EasyPAML is installed
echo ""
echo "[5/5] Launchers..."
cat > "$SCRIPT_DIR/EasyPAML.sh" <<'LAUNCHER'
#!/usr/bin/env bash
# Opens the EasyPAML window. Created by INSTALL_LINUX_MAC.sh.
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ ! -x "$DIR/.venv/bin/python" ]; then
    echo "EasyPAML is not installed in this folder. Run ./INSTALL_LINUX_MAC.sh first."
    exit 1
fi
exec "$DIR/.venv/bin/python" "$DIR/src/easypaml_window.py" "$@"
LAUNCHER
cat > "$SCRIPT_DIR/easypaml-cli.sh" <<'LAUNCHER'
#!/usr/bin/env bash
# EasyPAML command line (./easypaml-cli.sh --help). Created by INSTALL_LINUX_MAC.sh.
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "$DIR/.venv/bin/python" "$DIR/src/easypaml_cli.py" "$@"
LAUNCHER
chmod +x "$SCRIPT_DIR/EasyPAML.sh" "$SCRIPT_DIR/easypaml-cli.sh"
ok "EasyPAML.sh (window) and easypaml-cli.sh (command line) created"
if [ "$OS" = "Linux" ] && [ -n "${HOME:-}" ]; then
    DESKTOP_DIR="$HOME/.local/share/applications"
    if mkdir -p "$DESKTOP_DIR" 2>/dev/null; then
        cat > "$DESKTOP_DIR/EasyPAML.desktop" <<EOF
[Desktop Entry]
Version=1.0
Type=Application
Name=EasyPAML
Comment=Positive selection analysis with PAML/codeml
Exec="$SCRIPT_DIR/EasyPAML.sh"
Path=$SCRIPT_DIR
Terminal=false
Categories=Science;Biology;
EOF
        ok "shortcut added to the applications menu (EasyPAML)"
    fi
fi

echo ""
echo " ============================================================"
echo "  INSTALLATION COMPLETE"
echo " ============================================================"
echo ""
echo " To open EasyPAML:"
cmd "./EasyPAML.sh"
echo " Command line:"
cmd "./easypaml-cli.sh --help"
if [ "$MISSING_TK" -eq 1 ]; then
    echo ""
    warn "install tkinter so the window can open:"
    pkg_hint "python3-tk" "python3-tkinter" "python-tk@$PYVER"
fi
if [ "$CODEML_OK" -eq 0 ]; then
    echo ""
    warn "install codeml (PAML) before running analyses:"
    pkg_hint "paml" "paml" "brewsci/bio/paml"
fi
echo ""
