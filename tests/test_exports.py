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
    assert f"with {len(wb.sheetnames)} sheet(s)" in viewer._messages[0][1]
    head = [c.value for c in wb[wb.sheetnames[0]][1]]
    assert 'lnL (M7)' in head and 'q-value (BH)' in head and 'significant (q < 0.05)' in head


def test_html_export_is_english_and_counts_like_the_table(viewer, tmp_path, monkeypatch):
    out = tmp_path / 'r.html'
    monkeypatch.setattr(rv, 'ask_save_file', lambda *a, **k: str(out))
    viewer._export_html()
    assert [k for k, _ in viewer._messages] == ['info'], viewer._messages
    html = out.read_text(encoding='utf-8')
    assert not re.search(r'\b(analisados|significantes|às)\b', html)
    assert '<td>-4100.000</td>' in html and 'np (' not in html
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
    assert text.startswith('Positive selection supported')


@pytest.mark.parametrize('n', [3, 40])
def test_chart_figure_is_written_for_few_and_many_genes(tmp_path, n):
    import numpy as np
    from src.gui import charts
    rng = np.random.default_rng(1)
    test = charts.TestData('M8a', 'M8', 1, [charts.Point(f"g{i}", float(rng.chisquare(1) * (1 + 20 * (i % 4 == 0))),
                                                         i % 4 == 0) for i in range(n)])
    whole = [charts.Point(f"g{i}", float(rng.lognormal(-1, 0.5)), False) for i in range(n)]
    pos = [charts.Point(f"g{i}", float(rng.lognormal(1, 0.5)), i % 4 == 0) for i in range(n)]
    from matplotlib.figure import Figure
    fig = Figure(figsize=(7, 3.2))
    charts.draw_lrt(fig.add_subplot(1, 2, 1), test, charts.LIGHT)
    charts.draw_omega(fig.add_subplot(1, 2, 2), whole, pos, charts.LIGHT)
    out = tmp_path / 'fig.png'
    fig.savefig(out, dpi=300)
    assert out.stat().st_size > 10_000


def test_significant_without_sites_is_a_weak_signal():
    W = rv.ResultsViewerWindow
    rows = [('M8 vs M7', False, '0.03', '0.1', '', None, 100.6, 0.0115),
            ('M8 vs M8a', True, '0.009', '0.028', '', 0, 100.6, 0.0115)]
    text, _ = W._conclusion({('M7', 'M8'): False, ('M8a', 'M8'): True}, rows)
    assert text.startswith('⚠ Weak signal') and 'ω = 101' in text and '1.1%' in text


def test_sites_figure_marks_one_and_two_stars(tmp_path):
    from matplotlib.figure import Figure
    from src.gui import charts
    fig = Figure(figsize=(11, 4))
    charts.draw_sites(fig, [10, 39, 85, 120, 200], [0.986, 0.994, 0.990, 0.7, 0.999],
                      ['*', '**', '*', '', '**'], 300, [150, 151], charts.LIGHT, title="t")
    out = tmp_path / 'sites.png'
    fig.savefig(out, dpi=100)
    assert out.stat().st_size > 5_000
    lollipops = fig.axes[0]
    assert len(lollipops.lines) >= 3     # grey, * and ** markers


def test_hover_box_is_short_and_fast_with_thousands_of_genes():
    import time
    import numpy as np
    from matplotlib.backend_bases import MouseEvent
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from src.gui import charts
    rng = np.random.default_rng(0)
    test = charts.TestData('M8a', 'M8', 1, [charts.Point(f"g{i}", float(rng.chisquare(1)), False)
                                            for i in range(5000)])
    fig = Figure(figsize=(6, 3))
    FigureCanvasAgg(fig)
    ax = fig.add_subplot()
    charts.enable_hover(fig, ax, charts.draw_lrt(ax, test, charts.LIGHT, compact=True), charts.LIGHT)
    fig.canvas.draw()
    x, y = ax.transData.transform((0.5, 0))
    t0 = time.perf_counter()
    fig.canvas.callbacks.process('motion_notify_event', MouseEvent('motion_notify_event', fig.canvas, x, y + 20))
    assert time.perf_counter() - t0 < 0.5
    box = [t for t in ax.texts if t.get_visible() and '…' in t.get_text()]
    assert box and len(box[0].get_text().splitlines()) == 8
    assert box[0].get_text().splitlines()[-1].startswith('… +')


def test_newick_with_codeml_omega_labels_and_branch_marks():
    from src.gui.branch_tab import parse_newick, preorder
    w = parse_newick("((A #0.5 , B #0.5 ) #0.5 , ((C #0.5 , D #0.5 ) #2.1 , E #0.5 ) #0.5 );")
    marks = parse_newick("((A,B),((C,D)#1,E));")
    assert len(preorder(w)) == len(preorder(marks)) == 9
    stem = [n for n in preorder(marks) if n.label == 1]
    assert len(stem) == 1 and [c.name for c in stem[0].children] == ['C', 'D']
    assert [n.w for n in preorder(w)][5] == 2.1
    clade = parse_newick("((A,B),((C,D)$1,E));")
    assert sorted(n.name for n in preorder(clade) if n.label == 1) == ['', 'C', 'D']


def test_summary_answer_says_when_nothing_is_supported(viewer):
    text, supported = viewer._summary_answer_text()
    assert not supported and 'supported in none of the 1 gene(s)' in text


def test_branch_omega_written_once_per_labelled_clade():
    from src.gui.branch_tab import label_tops, parse_newick
    one_clade = parse_newick("(((a#1,b#1),c#1),(d,e));")
    assert [n.children != [] for n in label_tops(one_clade)] == [True]
    apart = parse_newick("(((a#1,b#1),c),((d#1,e#1),f),g#1);")
    tops = label_tops(apart)
    assert len(tops) == 3 and sum(1 for n in tops if not n.children) == 1
