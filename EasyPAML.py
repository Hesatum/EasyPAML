#!/usr/bin/env python3
"""EasyPAML: positive selection analysis with PAML/codeml (window entry point)."""
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
    print("EasyPAML needs tkinter (Python's window toolkit), which is not installed.\n"
          "  Ubuntu/Debian: sudo apt install python3-tk\n"
          "  Fedora:        sudo dnf install python3-tkinter\n"
          "  macOS (brew):  brew install python-tk\n"
          "The command line works without it: python3 easypaml_cli.py --help")
    sys.exit(1)
try:
    from src.gui.main_gui import App
except ImportError as exc:
    print(f"Missing dependency: {exc}\n"
          "Run the installer (Linux/macOS: ./install.sh · Windows: install.bat)\n"
          "and open the program with EasyPAML.sh / EasyPAML.bat (or .venv/bin/python EasyPAML.py).")
    sys.exit(1)


def main():
    App.load_language_pref()   # load the saved language before building the window
    app = App()
    app.mainloop()


if __name__ == "__main__":
    main()
