"""Excel and HTML exports of the results panel, on a run with a fake codeml."""
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
