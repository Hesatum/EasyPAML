"""Item 2 -- falhar de forma visível, nunca travar, Parar sem órfãos.

Usa um codeml FALSO (script Python) que imita os modos de falha vistos no
teste de usabilidade; o comportamento é escolhido por FAKE_CODEML_MODE.
"""
import os
import platform
import stat
import sys
import threading
import time
from pathlib import Path

import pytest

from src.backend import messages
from src.backend.codeml_backend import CodemlBatchAnalysis
from tests.conftest import real_codeml

DATA = Path(__file__).resolve().parent / 'data'

pytestmark = pytest.mark.skipif(platform.system() == 'Windows',
                                reason="codeml falso usa shebang (Linux/macOS)")

FAKE = r'''#!{python}
import os, re, sys, time
ctl = open(sys.argv[1]).read()
out = re.search(r'outfile\s*=\s*(\S+)', ctl).group(1)
mode = os.environ.get('FAKE_CODEML_MODE', 'ok')
print("CODONML (in paml version 9.9fake)", flush=True)
if sys.argv[1] == 'p.ctl':     # sonda de versão do EasyPAML
    sys.exit(0)
def write_ok():
    with open(out, 'w') as f:
        f.write("CODONML (in paml version 9.9fake)\nns =  10  ls = 300\n")
        f.write("lnL(ntime: 17  np: 22):  -4100.000000      +0.000000\n")
if mode == 'stdin':
    print("stop codon TGA in seq. #   3 (Gorilla_gorilla), nucleotide site 448")
    print("Press Enter to continue", flush=True)
    sys.stdin.read(1)          # com stdin fechado retorna na hora
    write_ok()
elif mode == 'idle':
    time.sleep(60)             # parado, sem CPU
elif mode == 'rc1':
    print("Error: Number of sequences different in tree and seq files.", flush=True)
    sys.exit(1)
elif mode == 'nolnl':
    open(out, 'w').write("CODONML (in paml version 9.9fake)\nsomething went wrong\n")
elif mode == 'slow':
    t0 = time.time()
    while time.time() - t0 < 60:
        pass                   # usa CPU (não é inativo)
else:
    write_ok()
'''


@pytest.fixture
def fake_codeml(tmp_path):
    path = tmp_path / 'fake_codeml'
    path.write_text(FAKE.format(python=sys.executable))
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return str(path)


def _app(tmp_path, codeml, models=('M7', 'M8'), **extra):
    inp = tmp_path / 'in'
    inp.mkdir(exist_ok=True)
    (inp / 'gene.fasta').write_text((DATA / 'gene_exemplo.fasta').read_text())
    app = CodemlBatchAnalysis()
    app.config = {
        'input_folder': inp, 'tree_file': DATA / 'gene_exemplo.nwk',
        'output_folder': tmp_path / 'out', 'models': list(models), 'n_workers': 1,
        'run_lrt': True, 'codeml_path': codeml, 'timeout': 30, 'idle_timeout': 3,
        'log_callback': lambda level, text: app._test_log.append((level, text)),
    }
    app._test_log = []
    messages.set_language('en')
    return app


def test_success_path(tmp_path, fake_codeml, monkeypatch):
    monkeypatch.setenv('FAKE_CODEML_MODE', 'ok')
    app = _app(tmp_path, fake_codeml)
    summary = app.run_batch_analysis()
    assert summary['ok'] == 1 and summary['failed'] == 0
    assert any('ANALYSIS COMPLETE' in t for _, t in app._test_log)
    # reprodutibilidade: .ctl com caminhos relativos + cópias ao lado
    m8 = tmp_path / 'out' / 'M8'
    ctl = (m8 / 'gene_M8.ctl').read_text()
    assert 'seqfile = gene_M8_seq.fasta' in ctl and 'treefile = gene_M8_tree.nwk' in ctl
    assert (m8 / 'gene_M8_seq.fasta').exists() and (m8 / 'gene_M8_tree.nwk').exists()
    assert (tmp_path / 'out' / 'run_config.json').exists()


def test_codeml_waiting_for_enter_does_not_hang(tmp_path, fake_codeml, monkeypatch):
    monkeypatch.setenv('FAKE_CODEML_MODE', 'stdin')
    app = _app(tmp_path, fake_codeml)
    t0 = time.time()
    summary = app.run_batch_analysis()
    assert time.time() - t0 < 20
    assert summary['ok'] == 1


def test_idle_codeml_is_killed_and_reported(tmp_path, fake_codeml, monkeypatch):
    monkeypatch.setenv('FAKE_CODEML_MODE', 'idle')
    app = _app(tmp_path, fake_codeml, models=('M7',))
    t0 = time.time()
    summary = app.run_batch_analysis()
    assert time.time() - t0 < 30
    assert summary['failed'] == 1
    reason = summary['failures']['gene']
    assert 'no CPU' in reason
    assert any(level == 'error' and 'FINISHED WITH FAILURES' in t for level, t in app._test_log)
    assert not any('ANALYSIS COMPLETE' in t for _, t in app._test_log)


def test_nonzero_exit_reports_last_line(tmp_path, fake_codeml, monkeypatch):
    monkeypatch.setenv('FAKE_CODEML_MODE', 'rc1')
    app = _app(tmp_path, fake_codeml)
    summary = app.run_batch_analysis()
    assert summary['failed'] == 1 and summary['ok'] == 0
    assert 'Number of sequences different' in summary['failures']['gene']
    status = (tmp_path / 'out' / 'genes_status.tsv').read_text()
    assert 'gene\tfailed' in status
    tsv = (tmp_path / 'out' / 'analysis_summary.tsv').read_text().splitlines()
    assert tsv[1].split('\t')[1] == 'failed'


def test_missing_lnl_is_failure(tmp_path, fake_codeml, monkeypatch):
    monkeypatch.setenv('FAKE_CODEML_MODE', 'nolnl')
    summary = _app(tmp_path, fake_codeml, models=('M7',)).run_batch_analysis()
    assert summary['failed'] == 1
    assert 'likelihood' in summary['failures']['gene']


def test_timeout_is_failure(tmp_path, fake_codeml, monkeypatch):
    monkeypatch.setenv('FAKE_CODEML_MODE', 'slow')
    app = _app(tmp_path, fake_codeml, models=('M7',))
    app.config['timeout'] = 2
    summary = app.run_batch_analysis()
    assert summary['failed'] == 1 and 'time limit' in summary['failures']['gene']


def test_stop_kills_codeml_without_orphans(tmp_path, fake_codeml, monkeypatch):
    monkeypatch.setenv('FAKE_CODEML_MODE', 'slow')
    app = _app(tmp_path, fake_codeml, models=('M7', 'M8'))
    app.config['idle_timeout'] = 0
    stop = threading.Event()
    app.config['stop_event'] = stop
    seen = []

    def _stopper():
        for _ in range(100):
            with app._processes_lock:
                procs = list(app._active_processes)
            if procs:
                seen.extend(procs)
                stop.set()
                return
            time.sleep(0.1)
    th = threading.Thread(target=_stopper)
    th.start()
    t0 = time.time()
    summary = app.run_batch_analysis()
    th.join()
    assert time.time() - t0 < 20
    assert seen and all(p.poll() is not None for p in seen)   # recolhidos, sem zumbi
    assert summary['stopped'] is True
    assert not any('ANALYSIS COMPLETE' in t for _, t in app._test_log)


def test_internal_stop_codon_fails_fast_with_position(tmp_path, fake_codeml, monkeypatch):
    monkeypatch.setenv('FAKE_CODEML_MODE', 'ok')
    app = _app(tmp_path, fake_codeml)
    (app.config['input_folder'] / 'gene.fasta').write_text((DATA / 'gene_problematico.fasta').read_text())
    summary = app.run_batch_analysis()
    reason = summary['failures']['gene']
    assert 'Gorilla_gorilla codon 150 (TGA)' in reason


def test_ignore_stop_codons_runs_and_warns_about_excluded_taxon(tmp_path, fake_codeml, monkeypatch):
    monkeypatch.setenv('FAKE_CODEML_MODE', 'ok')
    app = _app(tmp_path, fake_codeml)
    app.config['ignore_stop_codons'] = True
    (app.config['input_folder'] / 'gene.fasta').write_text((DATA / 'gene_problematico.fasta').read_text())
    summary = app.run_batch_analysis()
    assert summary['ok'] == 1
    warns = [t for level, t in app._test_log if level == 'warn']
    assert any('Macaca_mulata' in t and 'EXCLUDED' in t for t in warns)
    assert any('stop codon' in t for t in warns)


@pytest.mark.skipif(real_codeml() is None, reason="defina EASYPAML_TEST_CODEML para rodar com o codeml real")
def test_real_codeml_problematic_data_never_hangs(tmp_path):
    app = _app(tmp_path, real_codeml(), models=('M7',))
    app.config['ignore_stop_codons'] = True
    app.config['timeout'] = 600
    app.config['idle_timeout'] = 120
    (app.config['input_folder'] / 'gene.fasta').write_text((DATA / 'gene_problematico.fasta').read_text())
    t0 = time.time()
    summary = app.run_batch_analysis()
    assert summary['ok'] == 1 and time.time() - t0 < 300
    import json
    sm = json.loads((tmp_path / 'out' / 'M7' / 'gene_M7_sitemap.json').read_text())
    assert sm['verified'] is True and sm['n_codons_used_by_codeml'] == 299
