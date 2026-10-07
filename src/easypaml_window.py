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

# Run from the EasyPAML folder (one above src/), so relative paths work however
# the program was opened (double-click, shortcut, terminal).
_ROOT = Path(__file__).resolve().parents[1]
os.chdir(_ROOT)
sys.path.insert(0, str(_ROOT))

try:
    import tkinter  # noqa: F401
except ImportError:
    print("EasyPAML needs tkinter (Python's window toolkit), which is not installed.\n"
          "  Ubuntu/Debian: sudo apt install python3-tk\n"
          "  Fedora:        sudo dnf install python3-tkinter\n"
          "  macOS (brew):  brew install python-tk\n"
          "The command line works without it: ./easypaml-cli.sh --help")
    sys.exit(1)
def _import_app():
    from src.gui.main_gui import App
    return App


try:
    from src.gui.splash import load_with_splash
    try:
        App = load_with_splash(_import_app)
    except tkinter.TclError:     # no display for a splash window
        App = _import_app()
except ImportError as exc:
    print(f"Missing dependency: {exc}\n"
          "Run the installer (Linux/macOS: ./INSTALL_LINUX_MAC.sh · Windows: INSTALL_WINDOWS.bat),\n"
          "then open the program with the EasyPAML.sh or EasyPAML.bat it creates.")
    sys.exit(1)


def main():
    App.load_language_pref()   # load the saved language before building the window
    app = App()
    app.mainloop()


if __name__ == "__main__":
    main()
