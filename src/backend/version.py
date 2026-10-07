"""EasyPAML version, shown under About, by --version and in run_config.json.

The exact commit comes from git in a clone, or from _commit.txt, which GitHub
fills in when it builds a ZIP or tar.gz (export-subst in .gitattributes)."""

import subprocess
from functools import lru_cache
from pathlib import Path
from typing import Optional

__version__ = "0.4.0"

_ROOT = Path(__file__).resolve().parents[2]
_COMMIT_FILE = Path(__file__).with_name('_commit.txt')


@lru_cache(maxsize=1)
def source_commit() -> Optional[str]:
    """Full hash of the running code's commit, with '-dirty' if there are
    uncommitted changes; None if unknown."""
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


def after_release(describe: str, version: str = __version__) -> str:
    """'0.4.0+3' when `git describe --long` ('v0.4.0-3-gabc1234') says the code is 3
    commits after the release tag of this version; the version alone otherwise."""
    parts = describe.strip().rsplit('-', 2)
    if len(parts) == 3 and parts[0] == f"v{version}" and parts[1].isdigit() and int(parts[1]) > 0:
        return f"{version}+{parts[1]}"
    return version


@lru_cache(maxsize=1)
def _describe() -> str:
    if (_ROOT / '.git').exists():
        try:
            out = subprocess.run(['git', '-C', str(_ROOT), 'describe', '--tags', '--long', '--match', 'v*'],
                                 capture_output=True, text=True, timeout=5)
            if out.returncode == 0:
                return out.stdout.strip()
        except (OSError, subprocess.SubprocessError):
            pass
    return ''


def version_string() -> str:
    """'0.4.0 (commit e1a7f4a)' at the release, '0.4.0+3 (commit …)' three commits
    after it, for citing."""
    c = source_commit()
    ver = after_release(_describe())
    return f"{ver} (commit {c[:7]}{'-dirty' if c and c.endswith('-dirty') else ''})" if c else ver
