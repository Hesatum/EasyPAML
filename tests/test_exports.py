"""Results panel: Excel and HTML exports (on a fake-codeml run) and the per-gene conclusion."""
import platform
import re

import pytest

from src.gui import results_viewer as rv
from tests.test_runner_failures import _app, fake_codeml  # noqa: F401  (fixture)

pytestmark = pytest.mark.skipif(platform.system() == 'Windows',
                                reason="the fake codeml uses a shebang (Linux/macOS)")


@pytest.fixture
def viewer(tmp_path, fake_codeml, monkeypatch):  # noqa: F811
    app = _app(tmp_path, fake_codeml, models=('M7', 'M8', 'M8a'))
    app.run_batch_analysis()
    w = rv.ResultsViewerWindow.__new__(rv.ResultsViewerWindow)
    w.output_folder = tmp_path / 'out'
    assert w._load_data()
    w._extract_tag_columns()
    w._sort_by_significance()
    messages = []
    monkeypatch.setattr(rv, 'show_message', lambda parent, title, text, kind='info':
                        messages.append((kind, text)))
    w._messages = messages
    return w


def test_excel_export_writes_every_test(viewer, tmp_path, monkeypatch):
    openpyxl = pytest.importorskip('openpyxl')
    out = tmp_path / 'r.xlsx'
    monkeypatch.setattr(rv, 'ask_save_file', lambda *a, **k: str(out))
    viewer._export_excel()
    assert [k for k, _ in viewer._messages] == ['info'], viewer._messages
    wb = openpyxl.load_workbook(out)
    assert 'Summary' in wb.sheetnames and len(wb.sheetnames) >= 3
    head = [c.value for c in wb[wb.sheetnames[0]][1]]
    assert 'lnL (M7)' in head and 'q-value (BH)' in head and 'significant (q < 0.05)' in head


def test_html_export_is_english_and_counts_like_the_table(viewer, tmp_path, monkeypatch):
    out = tmp_path / 'r.html'
    monkeypatch.setattr(rv, 'ask_save_file', lambda *a, **k: str(out))
    viewer._export_html()
    assert [k for k, _ in viewer._messages] == ['info'], viewer._messages
    html = out.read_text(encoding='utf-8')
    assert not re.search(r'\b(analisados|significantes|às)\b', html)
    for meta, table in re.findall(r'<p class="meta">\d+ gene\(s\) tested &nbsp;·&nbsp; (\d+) significant'
                                  r'.*?<tbody>(.*?)</tbody>', html, re.S):
        assert int(meta) == table.count('<tr class="sig">')


def test_m8_vs_m7_alone_is_not_called_positive_selection():
    W = rv.ResultsViewerWindow
    rows = [('M8 vs M7', True, '1e-5', '1e-5', '', 0)]
    text, _ = W._conclusion({('M7', 'M8'): True}, rows)
    assert text.startswith('⚠') and 'M8a' in text
    text, _ = W._conclusion({('M7', 'M8'): True, ('M8a', 'M8'): False},
                            rows + [('M8 vs M8a', False, '0.2', '0.2', '', None)])
    assert text.startswith('⚠')
    text, _ = W._conclusion({('M7', 'M8'): True, ('M8a', 'M8'): True},
                            rows + [('M8 vs M8a', True, '1e-4', '1e-4', '', 3)])
    assert text.startswith('Positive selection')


@pytest.mark.parametrize('n', [3, 40])
def test_chart_figure_is_written_for_few_and_many_genes(tmp_path, n):
    import numpy as np
    from src.gui import charts
    rng = np.random.default_rng(1)
    test = charts.TestData('M8a', 'M8', 1, [charts.Point(f"g{i}", float(rng.chisquare(1) * (1 + 20 * (i % 4 == 0))),
                                                         i % 4 == 0) for i in range(n)])
    whole = [charts.Point(f"g{i}", float(rng.lognormal(-1, 0.5)), False) for i in range(n)]
    pos = [charts.Point(f"g{i}", float(rng.lognormal(1, 0.5)), i % 4 == 0) for i in range(n)]
    out = tmp_path / 'fig.png'
    charts.export_figure(out, [test], whole, pos)
    assert out.stat().st_size > 10_000


def test_significant_without_sites_is_a_weak_signal():
    W = rv.ResultsViewerWindow
    rows = [('M8 vs M7', False, '0.03', '0.1', '', None, 100.6, 0.0115),
            ('M8 vs M8a', True, '0.009', '0.028', '', 0, 100.6, 0.0115)]
    text, _ = W._conclusion({('M7', 'M8'): False, ('M8a', 'M8'): True}, rows)
    assert text.startswith('⚠ Weak signal') and 'ω = 101' in text and '1.1%' in text
