#!/usr/bin/env bash
# EasyPAML launcher for Linux and macOS. Uses the .venv/ created by
# install.sh; without it, tries the system python3.
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
if [ -x ".venv/bin/python" ]; then
    exec ".venv/bin/python" EasyPAML.py "$@"
fi
if command -v python3 >/dev/null 2>&1 && \
   python3 -c 'import customtkinter, Bio, pandas, scipy, matplotlib' >/dev/null 2>&1; then
    exec python3 EasyPAML.py "$@"
fi
echo "EasyPAML is not installed in this folder yet."
echo "Run first:  ./install.sh"
exit 1
