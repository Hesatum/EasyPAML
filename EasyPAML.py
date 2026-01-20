#!/usr/bin/env python3
"""
EasyPAML - Interface intuitiva para análise de seleção positiva com PAML/CODEML
Ponto de entrada único para usuários finais

Uso:
    python main.py
"""

import sys
from pathlib import Path

# Adicionar src ao path para imports relativos
sys.path.insert(0, str(Path(__file__).parent))

from src.gui.main_gui import App


def main():
    """Função principal - inicia a aplicação"""
    app = App()
    app.mainloop()


if __name__ == "__main__":
    main()
