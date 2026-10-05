"""Item 6 -- textos PT/EN completos e cores com contraste >= 4,5:1 (sem abrir Tk)."""
import re
import sys
import types
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture(scope='module')
def texts():
    sys.path.insert(0, str(ROOT / 'src'))
    from gui import gui_texts
    return gui_texts


def _strings(value):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for v in value.values():
            yield from _strings(v)
    elif isinstance(value, (list, tuple)):
        for v in value:
            yield from _strings(v)


def test_pt_and_en_have_the_same_keys(texts):
    assert set(texts.TEXTS_PT) == set(texts.TEXTS_EN)


def test_default_language_is_english(texts):
    assert texts._current_lang == 'en'


def test_pt_texts_have_accents(texts):
    """Palavras que apareciam sem acento na interface em PT."""
    unaccented = re.compile(r'\b(Nao|CONFIGURACOES|EXECUCAO|Analise|ANALISE|Parametros|Proposito|'
                            r'Interpretacao|Referencias|arvore|sitios?|Instrucoes)\b')
    for key, value in texts.TEXTS_PT.items():
        for text in _strings(value):
            assert not unaccented.search(text), (key, text)


def test_pt_texts_are_translated(texts):
    """Termos que ficavam em inglês na interface em PT."""
    for key in ('btn_tree_file', 'tab_site_models', 'model_status_default', 'viewer_tab_lrt',
                'viewer_tab_branch'):
        assert texts.TEXTS_PT[key] != texts.TEXTS_EN[key], key


def _lum(h):
    h = h.lstrip('#')
    r, g, b = (int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))
    f = lambda c: c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
    return 0.2126 * f(r) + 0.7152 * f(g) + 0.0722 * f(b)


def contrast(a, b):
    la, lb = sorted((_lum(a), _lum(b)), reverse=True)
    return (la + 0.05) / (lb + 0.05)


CARD_BACKGROUNDS = ('#0c0c0e', '#0d0d11', '#16161a', '#17171f', '#1e1e24', '#20202c', '#1a1a26')


def _palette():
    # ui_helpers importa customtkinter; para o teste basta o dicionário PALETTE
    src = (ROOT / 'src' / 'gui' / 'ui_helpers.py').read_text(encoding='utf-8')
    block = src[src.index('PALETTE = {'):src.index('}', src.index('PALETTE = {')) + 1]
    ns = {}
    exec(block, ns)
    return ns['PALETTE']


@pytest.mark.parametrize('key', ['text_primary', 'text_secondary', 'text_tertiary', 'text_muted',
                                 'accent_text', 'success_text', 'warning_text', 'danger_text', 'info_text'])
def test_text_colors_contrast(key):
    color = _palette()[key]
    for bg in CARD_BACKGROUNDS:
        assert contrast(color, bg) >= 4.5, (key, color, bg, round(contrast(color, bg), 2))


@pytest.mark.parametrize('key', ['accent_fill', 'success_fill', 'warning_fill', 'danger_fill', 'info_fill'])
def test_white_text_on_filled_buttons(key):
    assert contrast('#ffffff', _palette()[key]) >= 4.5


def test_window_colors_in_gui_modules_are_legible():
    """text_tertiary/text_muted das janelas principal e de resultados."""
    for path in ('src/gui/main_gui.py', 'src/gui/results_viewer.py'):
        src = (ROOT / path).read_text(encoding='utf-8')
        for key in ('text_tertiary', 'text_muted'):
            m = re.search(rf"'{key}':\s*'(#[0-9a-fA-F]{{6}})'", src)
            assert m, (path, key)
            assert contrast(m.group(1), '#17171f') >= 4.5, (path, key, m.group(1))


def test_results_window_fits_small_screens():
    src = (ROOT / 'src' / 'gui' / 'results_viewer.py').read_text(encoding='utf-8')
    assert 'fit_to_screen(self, 1400, 900' in src
    assert '<Escape>' in src
