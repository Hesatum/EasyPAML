"""Failures are visible, codeml never hangs, and Stop leaves no orphan process.

Uses a fake codeml (a Python script) whose behaviour is set by FAKE_CODEML_MODE."""
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
if sys.argv[1] == 'p.ctl':     # EasyPAML's version probe
    sys.exit(0)
def write_ok():
    with open(out, 'w') as f:
        f.write("CODONML (in paml version 9.9fake)\nns =  10  ls = 300\n")
        f.write("lnL(ntime: 17  np: 22):  -4100.000000      +0.000000\n")
fail_model = os.environ.get('FAKE_FAIL_MODEL')
if fail_model and sys.argv[1].endswith('_' + fail_model + '.ctl'):
    write_ok()                 # writes lnL, then stalls (like a codeml stopped during BEB)
    time.sleep(60)
    sys.exit(0)
if mode == 'stdin':
    print("stop codon TGA in seq. #   3 (Gorilla_gorilla), nucleotide site 448")
    print("Press Enter to continue", flush=True)
    sys.stdin.read(1)          # returns at once: stdin is closed
    write_ok()
elif mode == 'idle':
    time.sleep(60)             # idle, no CPU
elif mode == 'rc1':
    print("Error: Number of sequences different in tree and seq files.", flush=True)
    sys.exit(1)
elif mode == 'nolnl':
    open(out, 'w').write("CODONML (in paml version 9.9fake)\nsomething went wrong\n")
elif mode == 'busy':
    t0 = time.time()
    while time.time() - t0 < float(os.environ.get('FAKE_BUSY_S', '5')):
        pass                   # uses CPU, then finishes
    write_ok()
elif mode == 'slow':
    t0 = time.time()
    while time.time() - t0 < 60:
        pass                   # busy, never finishes
else:
    write_ok()
'''


@pytest.fixture
def fake_codeml(tmp_path):
    path = tmp_path / 'fake_codeml'
    path.write_text(FAKE.format(python=sys.executable))
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return str(path)


@pytest.fixture
def wrapped_codeml(tmp_path, fake_codeml):
    """Like /usr/bin/codeml on Debian/Ubuntu: a sh script that runs codeml as a child."""
    path = tmp_path / 'codeml_wrapper'
    path.write_text(f'#!/bin/sh -e\n{fake_codeml} "$@"\n')
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return str(path)


def _app(tmp_path, codeml, models=('M7', 'M8'), **extra):
    inp = tmp_path / 'in'
    inp.mkdir(parents=True, exist_ok=True)
    (inp / 'gene.fasta').write_text((DATA / 'gene_example.fasta').read_text())
    app = CodemlBatchAnalysis()
    app.config = {
        'input_folder': inp, 'tree_file': DATA / 'gene_example.nwk',
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
    # .ctl with relative paths, with the inputs next to it
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


def test_timeout_unset_uses_automatic_limit(tmp_path, fake_codeml, monkeypatch):
    monkeypatch.setenv('FAKE_CODEML_MODE', 'ok')
    app = _app(tmp_path, fake_codeml, models=('M7',))
    app.config['timeout'] = 0
    summary = app.run_batch_analysis()
    assert summary['ok'] == 1
    log = (tmp_path / 'out' / 'batch_analysis_log.txt').read_text()
    assert 'time limit 1800 s' in log and 'automatic' in log


def test_wrapper_script_busy_child_is_not_idle(tmp_path, wrapped_codeml, monkeypatch):
    """A busy codeml under a wrapper script is not taken as idle."""
    monkeypatch.setenv('FAKE_CODEML_MODE', 'busy')
    monkeypatch.setenv('FAKE_BUSY_S', '6')
    app = _app(tmp_path, wrapped_codeml, models=('M7',))
    app.config['idle_timeout'] = 2
    summary = app.run_batch_analysis()
    assert summary['ok'] == 1, summary['failures']


def test_wrapper_script_idle_child_still_detected(tmp_path, wrapped_codeml, monkeypatch):
    monkeypatch.setenv('FAKE_CODEML_MODE', 'idle')
    app = _app(tmp_path, wrapped_codeml, models=('M7',))
    app.config['idle_timeout'] = 2
    t0 = time.time()
    summary = app.run_batch_analysis()
    assert summary['failed'] == 1 and 'CPU' in summary['failures']['gene']
    assert time.time() - t0 < 30


def test_stop_kills_child_of_wrapper_script(tmp_path, wrapped_codeml, monkeypatch):
    monkeypatch.setenv('FAKE_CODEML_MODE', 'slow')
    app = _app(tmp_path, wrapped_codeml, models=('M7',))
    app.config['idle_timeout'] = 0
    app.config['timeout'] = 3
    app.run_batch_analysis()
    time.sleep(0.5)
    import subprocess
    left = subprocess.run(['pgrep', '-f', str(tmp_path / 'fake_codeml')],
                          capture_output=True, text=True).stdout.split()
    assert not left, f"codeml filho ficou rodando: {left}"


def test_failed_run_output_is_not_used_in_lrt(tmp_path, fake_codeml, monkeypatch):
    """Output of a run stopped after writing lnL stays out of the LRT."""
    monkeypatch.setenv('FAKE_CODEML_MODE', 'ok')
    monkeypatch.setenv('FAKE_FAIL_MODEL', 'M8')
    app = _app(tmp_path, fake_codeml, models=('M7', 'M8'))
    app.config['idle_timeout'] = 2
    summary = app.run_batch_analysis()
    assert summary['failed'] == 1
    m8 = tmp_path / 'out' / 'M8'
    assert not list(m8.glob('*_results.txt'))
    assert list(m8.glob('*_results_FAILED.txt'))
    tsv = (tmp_path / 'out' / 'analysis_summary.tsv').read_text().splitlines()
    head, row = tsv[0].split('\t'), tsv[1].split('\t')
    assert row[head.index('status')] == 'failed'
    assert 'lrt_M7_vs_M8' not in head or row[head.index('lrt_M7_vs_M8')] in ('', 'NA', 'nan')
    assert 'Gene: gene' not in (tmp_path / 'out' / 'LRT_results.txt').read_text()


def test_regenerate_skips_genes_marked_failed(tmp_path, fake_codeml, monkeypatch):
    """Older output folders: a failed gene's *_results.txt is ignored."""
    monkeypatch.setenv('FAKE_CODEML_MODE', 'ok')
    app = _app(tmp_path, fake_codeml, models=('M7', 'M8'))
    app.run_batch_analysis()
    out = tmp_path / 'out'
    (out / 'genes_status.tsv').write_text("Gene\tstatus\treason\ngene\tfailed\tM8: inatividade\n")
    CodemlBatchAnalysis.regenerate_summary_files(out)
    tsv = (out / 'analysis_summary.tsv').read_text().splitlines()
    head, row = tsv[0].split('\t'), tsv[1].split('\t')
    assert row[head.index('status')] == 'failed'
    assert 'lrt_M7_vs_M8' not in head
    assert CodemlBatchAnalysis._find_orphaned_analyses(out) == {}


def test_tree_without_branch_lengths_is_not_given_zeros(tmp_path, fake_codeml, monkeypatch):
    """A tree without branch lengths does not get ':0' on every branch."""
    import re
    monkeypatch.setenv('FAKE_CODEML_MODE', 'ok')
    app = _app(tmp_path, fake_codeml, models=('M7',))
    plain = re.sub(r':[0-9.eE+-]+', '', (DATA / 'gene_example.nwk').read_text())
    (tmp_path / 'plain.nwk').write_text(plain)
    app.config['tree_file'] = tmp_path / 'plain.nwk'
    assert app.run_batch_analysis()['ok'] == 1
    tree = next((tmp_path / 'out' / 'M7').glob('*_tree.nwk')).read_text()
    assert ':' not in tree.split('\n', 1)[1]


def test_cli_side_outputs_sites_table_and_m8a_warning(tmp_path, fake_codeml, monkeypatch):
    """The command line writes sites_BEB.tsv and warns about M7 vs M8 without M8a."""
    monkeypatch.setenv('FAKE_CODEML_MODE', 'ok')
    app = _app(tmp_path, fake_codeml, models=('M7', 'M8'))
    app.run_batch_analysis()
    assert (tmp_path / 'out' / 'sites_BEB.tsv').read_text().startswith('gene\tmodel\tmethod\t')
    assert any('M8a' in t for level, t in app._test_log if level == 'warn')


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
    assert seen and all(p.poll() is not None for p in seen)   # reaped, no zombies
    assert summary['stopped'] is True
    assert not any('ANALYSIS COMPLETE' in t for _, t in app._test_log)


def test_stop_button_marks_genes_stopped_not_failed(tmp_path, fake_codeml, monkeypatch):
    """The GUI sets stop_event and kills codeml at once; the killed run (exit -15)
    must count as stopped, not failed."""
    monkeypatch.setenv('FAKE_CODEML_MODE', 'slow')
    app = _app(tmp_path, fake_codeml, models=('M7', 'M8'))
    app.config['idle_timeout'] = 0
    stop = threading.Event()
    app.config['stop_event'] = stop

    def _stopper():
        for _ in range(100):
            with app._processes_lock:
                busy = bool(app._active_processes)
            if busy:
                stop.set()
                app.stop_all_processes()
                return
            time.sleep(0.1)
    th = threading.Thread(target=_stopper)
    th.start()
    summary = app.run_batch_analysis()
    th.join()
    assert summary['stopped'] is True
    assert summary['failed'] == 0
    assert app.gene_status == {'gene': 'stopped'}
    assert app.runs_done == app.runs_total


def test_internal_stop_codon_fails_fast_with_position(tmp_path, fake_codeml, monkeypatch):
    monkeypatch.setenv('FAKE_CODEML_MODE', 'ok')
    app = _app(tmp_path, fake_codeml)
    (app.config['input_folder'] / 'gene.fasta').write_text((DATA / 'gene_problematic.fasta').read_text())
    summary = app.run_batch_analysis()
    reason = summary['failures']['gene']
    assert 'Gorilla_gorilla codon 150 (TGA)' in reason


def test_ignore_stop_codons_runs_and_warns_about_excluded_taxon(tmp_path, fake_codeml, monkeypatch):
    monkeypatch.setenv('FAKE_CODEML_MODE', 'ok')
    app = _app(tmp_path, fake_codeml)
    app.config['ignore_stop_codons'] = True
    (app.config['input_folder'] / 'gene.fasta').write_text((DATA / 'gene_problematic.fasta').read_text())
    summary = app.run_batch_analysis()
    assert summary['ok'] == 1
    warns = [t for level, t in app._test_log if level == 'warn']
    assert any('Macaca_mulata' in t and 'EXCLUDED' in t for t in warns)
    assert any('stop codon' in t for t in warns)
    # the warning is kept after the run
    status = (tmp_path / 'out' / 'genes_status.tsv').read_text().splitlines()
    assert status[0] == 'Gene\tstatus\treason\tnotes'
    gene_line = status[1].split('\t')
    assert gene_line[1] == 'ok' and 'stop codon' in gene_line[3] and 'Macaca_mulata' in gene_line[3]


@pytest.mark.skipif(real_codeml() is None, reason="defina EASYPAML_TEST_CODEML para rodar com o codeml real")
def test_real_codeml_problematic_data_never_hangs(tmp_path):
    app = _app(tmp_path, real_codeml(), models=('M7',))
    app.config['ignore_stop_codons'] = True
    app.config['timeout'] = 600
    app.config['idle_timeout'] = 120
    (app.config['input_folder'] / 'gene.fasta').write_text((DATA / 'gene_problematic.fasta').read_text())
    t0 = time.time()
    summary = app.run_batch_analysis()
    assert summary['ok'] == 1 and time.time() - t0 < 300
    import json
    sm = json.loads((tmp_path / 'out' / 'M7' / 'gene_M7_sitemap.json').read_text())
    assert sm['verified'] is True and sm['n_codons_used_by_codeml'] == 299


def test_per_gene_tree_is_paired_by_file_name(tmp_path, fake_codeml, monkeypatch):
    """GENE.nwk next to the alignment replaces the general tree for that gene."""
    monkeypatch.setenv('FAKE_CODEML_MODE', 'ok')
    app = _app(tmp_path, fake_codeml, models=('M7',))
    inp = app.config['input_folder']
    (inp / 'g2.fasta').write_text((DATA / 'gene_example.fasta').read_text())
    # the gene's own tree has 6 taxa -> 4 sequences excluded
    (inp / 'gene.nwk').write_text(
        "((Homo_sapiens,Pan_troglodytes),Gorilla_gorilla,(Macaca_mulatta,(Papio_anubis,Aotus_nancymaae)));\n")
    app.config['tree_file'] = None
    summary = app.run_batch_analysis()
    # g2 has no tree at all -> clear failure
    assert summary['failures'].get('g2')
    tree_used = (tmp_path / 'out' / 'M7' / 'gene_M7_tree.nwk').read_text()
    assert tree_used.startswith('6  1')
    assert any('EXCLUDED' in t for _, t in app._test_log)

    app2 = _app(tmp_path / 'second', fake_codeml, models=('M7',))
    (app2.config['input_folder'] / 'gene.nwk').write_text(
        "((Homo_sapiens,Pan_troglodytes),Gorilla_gorilla,(Macaca_mulatta,(Papio_anubis,Aotus_nancymaae)));\n")
    (app2.config['input_folder'] / 'g2.fasta').write_text((DATA / 'gene_example.fasta').read_text())
    summary2 = app2.run_batch_analysis()          # with a general tree, g2 uses it
    assert summary2['ok'] == 2
    assert (tmp_path / 'second' / 'out' / 'M7' / 'g2_M7_tree.nwk').read_text().startswith('10  1')
    assert (tmp_path / 'second' / 'out' / 'M7' / 'gene_M7_tree.nwk').read_text().startswith('6  1')


def test_codeml_symlink_is_not_resolved(tmp_path):
    """Debian: /usr/bin/codeml is a link to one script that dispatches by name."""
    from src.backend.codeml_backend import find_codeml
    target = tmp_path / 'baseml'
    target.write_text('#!/bin/sh\n')
    link = tmp_path / 'codeml'
    link.symlink_to(target)
    assert find_codeml(str(link)).endswith('codeml')


def test_branch_site_uses_labeled_tree(tmp_path, fake_codeml, monkeypatch):
    """The #1-labelled tree reaches the branch-site .ctl pruned and without branch
    lengths; the null uses fix_omega = 1."""
    monkeypatch.setenv('FAKE_CODEML_MODE', 'ok')
    app = _app(tmp_path, fake_codeml, models=('Branch-site_null', 'Branch-site'))
    (app.config['input_folder'] / 'gene.fasta').write_text((DATA / 'gene_problematic.fasta').read_text())
    app.config['ignore_stop_codons'] = True
    app.config['labeled_tree_branchsite'] = (
        "((((Homo_sapiens:0.1,Pan_troglodytes:0.1)#1,Gorilla_gorilla),(Pongo_abelii,Hylobates_lar)),"
        "((Macaca_mulatta,Papio_anubis),(Callithrix_jacchus,(Saimiri_boliviensis,Aotus_nancymaae))));")
    summary = app.run_batch_analysis()
    assert summary['ok'] == 1
    out = tmp_path / 'out'
    tree = (out / 'Branch-site' / 'gene_Branch-site_tree.nwk').read_text()
    assert '#1' in tree and ':0.1' not in tree
    assert tree.startswith('9  1')                 # Macaca_mulata missing from the tree -> 9 taxa
    from src.backend.ctl_params import parse_ctl_text
    null = parse_ctl_text((out / 'Branch-site_null' / 'gene_Branch-site_null.ctl').read_text())
    assert null['model'] == '2' and null['NSsites'] == '2' and null['fix_omega'] == '1'
    alt = parse_ctl_text((out / 'Branch-site' / 'gene_Branch-site.ctl').read_text())
    assert alt['fix_omega'] == '0'
    assert 'lrt_Branch-site_null_vs_Branch-site' in (out / 'analysis_summary.tsv').read_text()


def test_pruning_keeps_branch_labels():
    """(Macaca,Papio)#2 without Macaca: the #2 label moves to Papio."""
    from io import StringIO
    from Bio import Phylo
    tree = Phylo.read(StringIO("(((Homo,Pan)#1,Gorilla),((Macaca,Papio)#2,(Aotus,Saimiri)));"), 'newick')
    CodemlBatchAnalysis._prune_keep_labels(tree, {'Homo', 'Pan', 'Gorilla', 'Papio', 'Aotus', 'Saimiri'})
    io = StringIO()
    Phylo.write(tree, io, 'newick')
    text = io.getvalue()
    assert 'Papio#2' in text and '#1' in text and 'Macaca' not in text
    tree = Phylo.read(StringIO("(((Homo,Pan)#1,Gorilla),(Macaca,Papio));"), 'newick')
    CodemlBatchAnalysis._prune_keep_labels(tree, {'Pan', 'Gorilla', 'Macaca', 'Papio'})
    io = StringIO()
    Phylo.write(tree, io, 'newick')
    assert 'Pan#1' in io.getvalue()
