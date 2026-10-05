#!/usr/bin/env bash
# EasyPAML - lançador Linux/macOS. Usa o ambiente .venv/ criado pelo
# install.sh; sem ele, tenta o python3 do sistema e explica o que falta.
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
if [ -x ".venv/bin/python" ]; then
    exec ".venv/bin/python" EasyPAML.py "$@"
fi
if command -v python3 >/dev/null 2>&1 && \
   python3 -c 'import customtkinter, Bio, pandas, scipy, matplotlib' >/dev/null 2>&1; then
    exec python3 EasyPAML.py "$@"
fi
echo "EasyPAML ainda não está instalado nesta pasta."
echo "Rode primeiro:  ./install.sh"
exit 1
