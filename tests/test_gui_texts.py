"""Interface texts in both languages and colour contrast (without opening Tk)."""
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
    """Portuguese words that must keep their accents."""
    unaccented = re.compile(r'\b(Nao|CONFIGURACOES|EXECUCAO|Analise|ANALISE|Parametros|Proposito|'
                            r'Interpretacao|Referencias|arvore|sitios?|Instrucoes)\b')
    for key, value in texts.TEXTS_PT.items():
        for text in _strings(value):
            assert not unaccented.search(text), (key, text)


def test_pt_texts_are_translated(texts):
    """Terms that must be translated in the Portuguese interface."""
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


LAYERS = ('bg_dark', 'bg_card', 'bg_card_hover', 'bg_window', 'bg_panel', 'bg_surface',
          'bg_elevated', 'bg_inset', 'row_alt')


def _palettes():
    """DARK_PALETTE and LIGHT_PALETTE from ui_helpers, without importing customtkinter."""
    src = (ROOT / 'src' / 'gui' / 'ui_helpers.py').read_text(encoding='utf-8')
    out = {}
    for name in ('DARK_PALETTE', 'LIGHT_PALETTE'):
        i = src.index(name + ' = {')
        ns = {}
        exec(src[i:src.index('\n}', i) + 2], ns)
        out[name] = ns[name]
    return out


PALETTES = _palettes()


def _layers(pal):
    return [pal[k] for k in LAYERS]


def test_both_palettes_have_the_same_keys():
    assert set(PALETTES['DARK_PALETTE']) == set(PALETTES['LIGHT_PALETTE'])


@pytest.mark.parametrize('name', sorted(PALETTES))
@pytest.mark.parametrize('key', ['text_primary', 'text_secondary', 'text_tertiary', 'text_muted',
                                 'accent_text', 'success_text', 'warning_text', 'danger_text', 'info_text',
                                 'accent_cyan', 'success_light'])
def test_text_colors_contrast(name, key):
    pal = PALETTES[name]
    for bg in _layers(pal):
        assert contrast(pal[key], bg) >= 4.5, (name, key, pal[key], bg, round(contrast(pal[key], bg), 2))


@pytest.mark.parametrize('name', sorted(PALETTES))
@pytest.mark.parametrize('key', ['accent_fill', 'success_fill', 'warning_fill', 'danger_fill', 'info_fill',
                                 'neutral_fill'])
def test_white_text_on_filled_buttons(name, key):
    assert contrast('#ffffff', PALETTES[name][key]) >= 4.5


@pytest.mark.parametrize('name', sorted(PALETTES))
@pytest.mark.parametrize('kind', ['success', 'warning', 'danger'])
def test_semantic_subtle_pairs(name, kind):
    """Painel de resultados: texto *_fg sobre o fundo *_subtle e sobre as camadas."""
    pal = PALETTES[name]
    fg = pal[f'{kind}_fg']
    for bg in [pal[f'{kind}_subtle']] + _layers(pal):
        assert contrast(fg, bg) >= 4.5, (name, kind, fg, bg, round(contrast(fg, bg), 2))


@pytest.mark.parametrize('name', sorted(PALETTES))
def test_text_on_subtle_backgrounds(name):
    pal = PALETTES[name]
    for key in ('text_primary', 'text_secondary', 'text_tertiary', 'accent_text'):
        for bg in ('success_subtle', 'warning_subtle', 'danger_subtle'):
            assert contrast(pal[key], pal[bg]) >= 4.5, (name, key, bg)


@pytest.mark.parametrize('name', sorted(PALETTES))
def test_plot_and_canvas_text(name):
    pal = PALETTES[name]
    assert contrast(pal['plot_fg'], pal['plot_bg']) >= 4.5
    assert contrast(pal['plot_muted'], pal['plot_bg']) >= 4.5
    assert contrast(pal['canvas_text'], pal['canvas_bg']) >= 4.5


def test_gui_modules_take_text_colors_from_the_palette():
    """Window text colours come from PALETTE, so they work in both themes."""
    for path in ('src/gui/main_gui.py', 'src/gui/results_viewer.py'):
        src = (ROOT / path).read_text(encoding='utf-8')
        for key in ('text_tertiary', 'text_muted'):
            assert re.search(rf"'{key}':\s*PALETTE\['{key}'\]", src), (path, key)


def test_results_window_fits_small_screens():
    src = (ROOT / 'src' / 'gui' / 'results_viewer.py').read_text(encoding='utf-8')
    assert 'fit_to_screen(self, 1400, 900' in src
    assert '<Escape>' in src
