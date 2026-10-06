import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'src'))

DATA = ROOT / 'tests' / 'data'


def real_codeml():
    """Real codeml for the integration tests (EASYPAML_TEST_CODEML)."""
    p = os.environ.get('EASYPAML_TEST_CODEML')
    return p if p and Path(p).exists() else None
