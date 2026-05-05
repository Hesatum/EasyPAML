#!/usr/bin/env bash
# ============================================================
#  EasyPAML - Instalador para Linux / macOS
#  Analise de Selecao Positiva com CODEML/PAML
# ============================================================

set -euo pipefail

# Ir para o diretorio do script, independente de onde foi chamado
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; NC='\033[0m'
ok()   { echo -e "${GREEN} OK:${NC} $*"; }
warn() { echo -e "${YELLOW} AV:${NC} $*"; }
err()  { echo -e "${RED} ERRO:${NC} $*"; }

echo ""
echo " ============================================================"
echo "  EasyPAML - Instalador Linux / macOS"
echo "  Analise de Selecao Positiva com CODEML/PAML"
echo " ============================================================"
echo ""

# ── 1. Verificar Python 3.8+ ─────────────────────────────────────────────────
echo "[1/4] Verificando Python..."

PYTHON=""
for cmd in python3 python python3.12 python3.11 python3.10 python3.9 python3.8; do
    if command -v "$cmd" &>/dev/null; then
        ver=$("$cmd" -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
        major="${ver%%.*}"; minor="${ver##*.}"
        if [ "$major" -ge 3 ] && [ "$minor" -ge 8 ]; then
            PYTHON="$cmd"
            break
        fi
    fi
done

if [ -z "$PYTHON" ]; then
    err "Python 3.8 ou superior nao encontrado!"
    echo ""
    echo " Como instalar:"
    echo "   Ubuntu/Debian : sudo apt-get install python3 python3-pip"
    echo "   Fedora/RHEL   : sudo dnf install python3 python3-pip"
    echo "   macOS (brew)  : brew install python"
    echo "   macOS (site)  : https://www.python.org/downloads/"
    echo ""
    exit 1
fi

PY_VER=$("$PYTHON" --version 2>&1)
ok "$PY_VER encontrado ($PYTHON)"

# ── 2. Instalar dependencias Python ──────────────────────────────────────────
echo ""
echo "[2/4] Instalando dependencias Python..."
echo " (Pode levar alguns minutos na primeira vez)"
echo ""

# Atualizar pip silenciosamente
"$PYTHON" -m pip install --upgrade pip --quiet --user 2>/dev/null || \
    warn "Nao foi possivel atualizar pip, continuando..."

# Instalar requirements
"$PYTHON" -m pip install -r requirements.txt --user
ok "Dependencias instaladas"

# ── 3. Verificar / instalar CODEML ───────────────────────────────────────────
echo ""
echo "[3/4] Verificando CODEML..."

CODEML_FOUND=0
CODEML_BUNDLED="$SCRIPT_DIR/bin/codeml"

if [ -x "$CODEML_BUNDLED" ]; then
    ok "CODEML bundled encontrado em bin/codeml"
    CODEML_FOUND=1
elif command -v codeml &>/dev/null; then
    ok "CODEML encontrado no PATH do sistema: $(which codeml)"
    # Criar symlink/copia local para uso uniforme
    mkdir -p "$SCRIPT_DIR/bin"
    cp "$(which codeml)" "$SCRIPT_DIR/bin/codeml" 2>/dev/null || \
        ln -sf "$(which codeml)" "$SCRIPT_DIR/bin/codeml" 2>/dev/null || true
    CODEML_FOUND=1
else
    warn "CODEML nao encontrado. Tentando instalar automaticamente..."
    OS=$(uname -s)
    ARCH=$(uname -m)

    if [[ "$OS" == "Linux" ]]; then
        # Tentar via gerenciador de pacotes
        if command -v apt-get &>/dev/null; then
            sudo apt-get install -y paml 2>/dev/null && CODEML_FOUND=1 || true
        elif command -v dnf &>/dev/null; then
            sudo dnf install -y paml 2>/dev/null && CODEML_FOUND=1 || true
        elif command -v yum &>/dev/null; then
            sudo yum install -y paml 2>/dev/null && CODEML_FOUND=1 || true
        fi

        # Atualizar link local apos instalacao por pacote
        if [ $CODEML_FOUND -eq 1 ] && command -v codeml &>/dev/null; then
            mkdir -p "$SCRIPT_DIR/bin"
            cp "$(which codeml)" "$SCRIPT_DIR/bin/codeml" 2>/dev/null || \
                ln -sf "$(which codeml)" "$SCRIPT_DIR/bin/codeml" 2>/dev/null || true
        fi

        # Tentar baixar binario pre-compilado se ainda nao encontrado
        if [ $CODEML_FOUND -eq 0 ] && command -v wget &>/dev/null; then
            echo " Tentando baixar binario pre-compilado do PAML..."
            mkdir -p "$SCRIPT_DIR/bin"
            # PAML 4.10.7 Linux x86_64 pre-compiled
            PAML_URL="https://github.com/abacus-gene/paml/releases/download/v4.10.7/paml-4.10.7-linux-x86_64.tar.gz"
            TMP_DIR=$(mktemp -d)
            if wget -q "$PAML_URL" -O "$TMP_DIR/paml.tar.gz"; then
                tar -xzf "$TMP_DIR/paml.tar.gz" -C "$TMP_DIR" 2>/dev/null || true
                CODEML_BIN=$(find "$TMP_DIR" -name "codeml" -type f 2>/dev/null | head -1)
                if [ -n "$CODEML_BIN" ]; then
                    cp "$CODEML_BIN" "$SCRIPT_DIR/bin/codeml"
                    chmod +x "$SCRIPT_DIR/bin/codeml"
                    CODEML_FOUND=1
                    ok "CODEML baixado e instalado em bin/codeml"
                fi
            fi
            rm -rf "$TMP_DIR"
        fi

    elif [[ "$OS" == "Darwin" ]]; then
        if command -v brew &>/dev/null; then
            brew install brewsci/bio/paml 2>/dev/null && CODEML_FOUND=1 || true
            if command -v codeml &>/dev/null; then
                mkdir -p "$SCRIPT_DIR/bin"
                cp "$(which codeml)" "$SCRIPT_DIR/bin/codeml" 2>/dev/null || true
            fi
        fi
    fi

    if [ $CODEML_FOUND -eq 0 ]; then
        warn "CODEML nao foi instalado automaticamente."
        echo ""
        echo " Instale manualmente:"
        echo "   Ubuntu/Debian : sudo apt-get install paml"
        echo "   Fedora        : sudo dnf install paml"
        echo "   macOS         : brew install brewsci/bio/paml"
        echo "   Manual        : http://abacus.gene.ucl.ac.uk/software/paml.html"
        echo ""
        echo " O EasyPAML sera instalado mas precisara do CODEML para rodar analises."
    fi
fi

# ── 4. Criar launcher ─────────────────────────────────────────────────────────
echo ""
echo "[4/4] Criando launcher..."

LAUNCHER="$SCRIPT_DIR/EasyPAML.sh"
cat > "$LAUNCHER" <<EOF
#!/usr/bin/env bash
# EasyPAML launcher — gerado pelo install.sh
cd "\$(dirname "\${BASH_SOURCE[0]}")"
exec $PYTHON EasyPAML.py "\$@"
EOF
chmod +x "$LAUNCHER"
ok "Launcher criado: EasyPAML.sh"

# Atalho .desktop para Linux (GNOME/KDE/XFCE)
if [[ "$(uname -s)" == "Linux" ]]; then
    DESKTOP_DIR="$HOME/.local/share/applications"
    mkdir -p "$DESKTOP_DIR"
    cat > "$DESKTOP_DIR/EasyPAML.desktop" <<EOF
[Desktop Entry]
Version=1.0
Type=Application
Name=EasyPAML
Comment=Analise de Selecao Positiva com CODEML/PAML
Exec=$LAUNCHER
Terminal=false
Categories=Science;Biology;
EOF
    chmod +x "$DESKTOP_DIR/EasyPAML.desktop"
    ok "Atalho de desktop criado (EasyPAML no menu de aplicativos)"
fi

echo ""
echo " ============================================================"
echo "  INSTALACAO CONCLUIDA!"
echo " ============================================================"
echo ""
echo " Para iniciar o EasyPAML:"
echo "   ./EasyPAML.sh"
if [[ "$(uname -s)" == "Linux" ]]; then
    echo "   Ou procure 'EasyPAML' no menu de aplicativos"
fi
echo ""
echo " Dados de exemplo em: exemplos_teste/"
echo ""
