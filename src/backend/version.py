"""Versão do EasyPAML -- fonte única, usada pela GUI (menu Sobre), pelo CLI
(--version) e gravada em run_config.json.

Além do número, o commit exato: do git, quando a pasta é um clone; senão de
_commit.txt, que o GitHub preenche ao gerar o ZIP ou o tar.gz (export-subst
em .gitattributes)."""

import subprocess
from functools import lru_cache
from pathlib import Path
from typing import Optional

__version__ = "0.3.0.dev0"

_ROOT = Path(__file__).resolve().parents[2]
_COMMIT_FILE = Path(__file__).with_name('_commit.txt')


@lru_cache(maxsize=1)
def source_commit() -> Optional[str]:
    """Hash completo do commit do código em uso, com '-dirty' se houver
    mudanças locais não commitadas; None se não der para saber."""
    if (_ROOT / '.git').exists():
        try:
            head = subprocess.run(['git', '-C', str(_ROOT), 'rev-parse', 'HEAD'],
                                  capture_output=True, text=True, timeout=5)
            if head.returncode == 0 and head.stdout.strip():
                dirty = subprocess.run(['git', '-C', str(_ROOT), 'status', '--porcelain',
                                        '--untracked-files=no'],
                                       capture_output=True, text=True, timeout=5)
                return head.stdout.strip() + ('-dirty' if dirty.stdout.strip() else '')
        except (OSError, subprocess.SubprocessError):
            pass
    try:
        text = _COMMIT_FILE.read_text(encoding='utf-8').strip()
    except OSError:
        return None
    return text if text and not text.startswith('$Format') else None


def version_string() -> str:
    """'0.3.0.dev0 (commit e1a7f4a)' -- para citar."""
    c = source_commit()
    return f"{__version__} (commit {c[:7]}{'-dirty' if c and c.endswith('-dirty') else ''})" \
        if c else __version__
