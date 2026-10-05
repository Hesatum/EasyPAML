#!/usr/bin/env bash
# ============================================================
#  EasyPAML - Instalador para Linux / macOS
#
#  Cria um ambiente Python isolado em .venv/ (dentro desta pasta),
#  instala as dependências nele (sem --user, sem mexer no Python do
#  sistema -- funciona no Ubuntu 23.04+/24.04 com PEP 668) e cria o
#  lançador EasyPAML.sh. Pode ser rodado de novo quantas vezes quiser.
# ============================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [ -t 1 ]; then
    GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; BOLD='\033[1m'; NC='\033[0m'
else
    GREEN=''; YELLOW=''; RED=''; BOLD=''; NC=''
fi
ok()   { echo -e "${GREEN} OK:${NC} $*"; }
warn() { echo -e "${YELLOW} AVISO:${NC} $*"; }
err()  { echo -e "${RED} ERRO:${NC} $*"; }
cmd()  { echo -e "      ${BOLD}$*${NC}"; }

OS="$(uname -s)"
HAS_APT=0; command -v apt-get >/dev/null 2>&1 && HAS_APT=1
HAS_DNF=0; command -v dnf >/dev/null 2>&1 && HAS_DNF=1

# comando sugerido para instalar pacotes do sistema
pkg_hint() {
    local deb="$1" rpm="$2" brew="$3"
    if [ "$HAS_APT" -eq 1 ]; then cmd "sudo apt update && sudo apt install -y $deb"
    elif [ "$HAS_DNF" -eq 1 ]; then cmd "sudo dnf install -y $rpm"
    elif [ "$OS" = "Darwin" ]; then cmd "brew install $brew"
    else cmd "instale: $deb"
    fi
}

echo ""
echo " ============================================================"
echo "  EasyPAML - Instalador Linux / macOS"
echo " ============================================================"
echo ""

# ── 1. Python 3.8+ ───────────────────────────────────────────────────────────
echo "[1/5] Verificando Python..."
PYTHON=""
for c in python3 python3.13 python3.12 python3.11 python3.10 python3.9 python3.8 python; do
    if command -v "$c" >/dev/null 2>&1 && \
       "$c" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 8) else 1)' 2>/dev/null; then
        PYTHON="$c"; break
    fi
done
if [ -z "$PYTHON" ]; then
    err "Python 3.8 ou superior não encontrado. Instale com:"
    pkg_hint "python3" "python3" "python"
    exit 1
fi
ok "$("$PYTHON" --version 2>&1) ($PYTHON)"

# ── 2. Pré-requisitos: venv e tkinter ────────────────────────────────────────
echo ""
echo "[2/5] Verificando pré-requisitos (venv, tkinter)..."
PYVER="$("$PYTHON" -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')"
MISSING_TK=0
if [ ! -x ".venv/bin/python" ]; then
    if ! "$PYTHON" -c 'import venv, ensurepip' >/dev/null 2>&1; then
        err "O módulo 'venv' do Python não está instalado (necessário para instalar as dependências)."
        echo "      Rode o comando abaixo e depois rode ./install.sh de novo:"
        pkg_hint "python3-venv python3-tk" "python3-tkinter" "python-tk@$PYVER"
        exit 1
    fi
fi
ok "venv disponível"
if ! "$PYTHON" -c 'import tkinter' >/dev/null 2>&1; then
    MISSING_TK=1
    warn "tkinter não instalado: a janela do EasyPAML não abre sem ele (o modo linha de comando funciona)."
    echo "      Para instalar:"
    pkg_hint "python3-tk" "python3-tkinter" "python-tk@$PYVER"
else
    ok "tkinter disponível"
fi

# ── 3. Ambiente .venv e dependências ─────────────────────────────────────────
echo ""
echo "[3/5] Instalando dependências em .venv/ (pode levar alguns minutos na primeira vez)..."
if [ ! -x ".venv/bin/python" ]; then
    if ! "$PYTHON" -m venv .venv; then
        err "Não foi possível criar .venv/. Instale o venv e rode de novo:"
        pkg_hint "python3-venv" "python3" "python"
        exit 1
    fi
fi
VENV_PY="$SCRIPT_DIR/.venv/bin/python"
"$VENV_PY" -m pip install --upgrade pip --quiet --disable-pip-version-check || \
    warn "não foi possível atualizar o pip; seguindo com a versão atual"
if ! "$VENV_PY" -m pip install -r requirements.txt --disable-pip-version-check; then
    err "Falha ao instalar as dependências (veja a mensagem acima). Causas comuns: sem internet, proxy."
    echo "      Para tentar de novo: ./install.sh"
    exit 1
fi
ok "Dependências instaladas em .venv/"

# ── 4. CODEML ────────────────────────────────────────────────────────────────
echo ""
echo "[4/5] Verificando o CODEML (PAML)..."
CODEML_OK=0
if [ -x "$SCRIPT_DIR/bin/codeml" ]; then
    ok "codeml do projeto: bin/codeml"; CODEML_OK=1
elif command -v codeml >/dev/null 2>&1; then
    ok "codeml do sistema: $(command -v codeml)"; CODEML_OK=1
else
    warn "codeml não encontrado."
    if [ "$HAS_APT" -eq 1 ]; then
        echo "      Instalando o pacote 'paml' (pode pedir a sua senha de administrador)..."
        if sudo apt-get install -y paml; then CODEML_OK=1; ok "PAML instalado (apt)"; fi
    elif [ "$HAS_DNF" -eq 1 ]; then
        if sudo dnf install -y paml; then CODEML_OK=1; ok "PAML instalado (dnf)"; fi
    elif [ "$OS" = "Darwin" ] && command -v brew >/dev/null 2>&1; then
        if brew install brewsci/bio/paml; then CODEML_OK=1; ok "PAML instalado (brew)"; fi
    fi
    if [ "$CODEML_OK" -eq 0 ] && [ "$OS" = "Linux" ] && [ "$(uname -m)" = "x86_64" ]; then
        echo "      Tentando o binário oficial do PAML (GitHub abacus-gene/paml)..."
        PAML_URL="https://github.com/abacus-gene/paml/releases/download/v4.10.10/paml-4.10.10-linux-x86_64.tar.gz"
        TMP_DIR="$(mktemp -d)"
        if (command -v curl >/dev/null 2>&1 && curl -fsSL "$PAML_URL" -o "$TMP_DIR/paml.tgz") || \
           (command -v wget >/dev/null 2>&1 && wget -q "$PAML_URL" -O "$TMP_DIR/paml.tgz"); then
            tar -xzf "$TMP_DIR/paml.tgz" -C "$TMP_DIR" 2>/dev/null || true
            BIN="$(find "$TMP_DIR" -name codeml -type f | head -1)"
            if [ -n "$BIN" ]; then
                mkdir -p "$SCRIPT_DIR/bin" && cp "$BIN" "$SCRIPT_DIR/bin/codeml" && \
                    chmod +x "$SCRIPT_DIR/bin/codeml" && CODEML_OK=1 && ok "codeml copiado para bin/codeml"
            fi
        fi
        rm -rf "$TMP_DIR"
    fi
    if [ "$CODEML_OK" -eq 0 ]; then
        warn "O CODEML não foi instalado. O EasyPAML abre, mas precisa dele para rodar análises:"
        pkg_hint "paml" "paml" "brewsci/bio/paml"
    fi
fi

# ── 5. Lançador ──────────────────────────────────────────────────────────────
echo ""
echo "[5/5] Lançador..."
chmod +x "$SCRIPT_DIR/EasyPAML.sh" 2>/dev/null || true
ok "EasyPAML.sh pronto"
if [ "$OS" = "Linux" ] && [ -n "${HOME:-}" ]; then
    DESKTOP_DIR="$HOME/.local/share/applications"
    if mkdir -p "$DESKTOP_DIR" 2>/dev/null; then
        cat > "$DESKTOP_DIR/EasyPAML.desktop" <<EOF
[Desktop Entry]
Version=1.0
Type=Application
Name=EasyPAML
Comment=Análise de seleção positiva com PAML/codeml
Exec="$SCRIPT_DIR/EasyPAML.sh"
Path=$SCRIPT_DIR
Terminal=false
Categories=Science;Biology;
EOF
        ok "atalho no menu de aplicativos (EasyPAML)"
    fi
fi

echo ""
echo " ============================================================"
echo "  INSTALAÇÃO CONCLUÍDA"
echo " ============================================================"
echo ""
echo " Para abrir o EasyPAML:"
cmd "./EasyPAML.sh"
echo "   (ou: .venv/bin/python EasyPAML.py)"
echo " Modo linha de comando:"
cmd ".venv/bin/python easypaml_cli.py --help"
if [ "$MISSING_TK" -eq 1 ]; then
    echo ""
    warn "lembre de instalar o tkinter para a janela abrir:"
    pkg_hint "python3-tk" "python3-tkinter" "python-tk@$PYVER"
fi
if [ "$CODEML_OK" -eq 0 ]; then
    echo ""
    warn "lembre de instalar o CODEML (PAML) antes de rodar análises:"
    pkg_hint "paml" "paml" "brewsci/bio/paml"
fi
echo ""
