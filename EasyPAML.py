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

from src.gui.main_gui import App


def main():
    App.load_language_pref()   # carrega idioma salvo ANTES de construir a janela
    app = App()
    app.mainloop()


if __name__ == "__main__":
    main()
