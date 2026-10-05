#!/usr/bin/env python3
"""
EasyPAML - Interface intuitiva para análise de seleção positiva com PAML/CODEML
"""
import sys
import os
from pathlib import Path

# Force UTF-8 output so emoji/Unicode in print() works on Windows (cp1252 consoles).
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
if hasattr(sys.stderr, 'reconfigure'):
    sys.stderr.reconfigure(encoding='utf-8', errors='replace')

# Always run from the directory where this file lives, so relative paths work
# regardless of how the user launched the app (double-click, shortcut, terminal).
_HERE = Path(__file__).resolve().parent
os.chdir(_HERE)
sys.path.insert(0, str(_HERE))

try:
    import tkinter  # noqa: F401
except ImportError:
    print("O EasyPAML precisa do tkinter (janela gráfica do Python), que não está instalado.\n"
          "  Ubuntu/Debian: sudo apt install python3-tk\n"
          "  Fedora:        sudo dnf install python3-tkinter\n"
          "  macOS (brew):  brew install python-tk\n"
          "O modo linha de comando funciona sem ele: python3 easypaml_cli.py --help\n"
          "--- EasyPAML needs tkinter (sudo apt install python3-tk).")
    sys.exit(1)
try:
    from src.gui.main_gui import App
except ImportError as exc:
    print(f"Dependência não encontrada: {exc}\n"
          "Rode o instalador (Linux/macOS: ./install.sh · Windows: install.bat)\n"
          "e abra pelo EasyPAML.sh / EasyPAML.bat (ou .venv/bin/python EasyPAML.py).\n"
          f"--- Missing dependency: {exc}. Run the installer first.")
    sys.exit(1)


def main():
    App.load_language_pref()   # carrega idioma salvo ANTES de construir a janela
    app = App()
    app.mainloop()


if __name__ == "__main__":
    main()
