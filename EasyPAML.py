#!/usr/bin/env python3
"""
EasyPAML - Interface intuitiva para análise de seleção positiva com PAML/CODEML
"""
import sys
import os
from pathlib import Path

# Always run from the directory where this file lives, so relative paths work
# regardless of how the user launched the app (double-click, shortcut, terminal).
_HERE = Path(__file__).resolve().parent
os.chdir(_HERE)
sys.path.insert(0, str(_HERE))

from src.gui.main_gui import App


def main():
    app = App()
    app.mainloop()


if __name__ == "__main__":
    main()
