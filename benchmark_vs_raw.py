#!/usr/bin/env python3
"""
benchmark_vs_raw.py -- compara CODEML puro (subprocess.run direto, sem
wrapper) contra o EasyPAML CLI, no mesmo locus/modelo, pra medir overhead
real do wrapper Python e confirmar lnL/np identicos nos dois caminhos.

Reusa a arvore ja podada/desenraizada que o pipeline principal do projeto
deixou em disco (nao reimplementa poda/desenraizamento -- essa logica ja
existe e ja foi validada em codeml_backend.py). So funciona pra loci com
os 24 taxons completos (sem poda necessaria) -- e a maioria.

Uso:
  python3 benchmark_vs_raw.py --loci C10_100132523_g6889,C10_100038951_g6880 --model M1a
"""
import argparse
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent
CODEML_BIN = REPO / 'bin' / 'codeml'
EASYPAML_INPUT = Path('/home/user/Desktop/projetos/matheus_wgs/resultados_analises/13_te_contamination/easypaml_input_5725')
TREE_FILE = Path('/home/user/Desktop/projetos/matheus_wgs/resultados_analises/13_te_contamination/codeml_site_models_5725/tree.nwk')
EASYPAML_TREE = Path('/home/user/Desktop/projetos/matheus_wgs/resultados_analises/13_te_contamination/astral_castles_comparison/astral4_11386.tre')

CTL_TEMPLATE = """\
      seqfile = seq.phy
     treefile = tree.nwk
      outfile = out.txt
        noisy = 1
      verbose = 1
      seqtype = 1
        ndata = 1
        icode = 0
    cleandata = 1
        model = 0
      NSsites = {nssites}
    CodonFreq = 7
      estFreq = 0
        clock = 0
    fix_omega = 0
        omega = 0.5
"""
NSSITES = {'M1a': 1, 'M2a': 2, 'M7': 7, 'M8': 8}
LNL_RE = re.compile(r'lnL\(ntime:\s*(\d+)\s*np:\s*(\d+)\):\s*([+-]?\d+\.\d+)')


def read_fasta(path: Path) -> dict:
    seqs, cur = {}, None
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            cur = line[1:].split()[0]
            seqs[cur] = []
        elif cur:
            seqs[cur].append(line.strip())
    return {k: ''.join(v) for k, v in seqs.items()}


def run_raw_codeml(locus: str, model: str) -> dict:
    seqs = read_fasta(EASYPAML_INPUT / f'{locus}.fasta')
    with tempfile.TemporaryDirectory(prefix='rawcodeml_') as d:
        d = Path(d)
        shutil.copy(TREE_FILE, d / 'tree.nwk')
        with open(d / 'seq.phy', 'w') as fh:
            fh.write(f' {len(seqs)} {len(next(iter(seqs.values())))}\n')
            for name, seq in sorted(seqs.items()):
                fh.write(f'{name}  {seq}\n')
        (d / 'codeml.ctl').write_text(CTL_TEMPLATE.format(nssites=NSSITES[model]))

        t0 = time.time()
        subprocess.run([str(CODEML_BIN), 'codeml.ctl'], cwd=d, input='\n',
                        capture_output=True, text=True, timeout=1600)
        elapsed = time.time() - t0

        out_path = d / 'out.txt'
        m = LNL_RE.search(out_path.read_text(errors='ignore')) if out_path.exists() else None
        return {'elapsed': elapsed, 'lnL': float(m.group(3)) if m else None, 'np': int(m.group(2)) if m else None}


def run_easypaml(locus: str, model: str) -> dict:
    from src.backend.codeml_backend import CodemlBatchAnalysis
    with tempfile.TemporaryDirectory(prefix='easypaml_bench_in_') as ind, \
         tempfile.TemporaryDirectory(prefix='easypaml_bench_out_') as outd:
        shutil.copy(EASYPAML_INPUT / f'{locus}.fasta', ind)
        app = CodemlBatchAnalysis()
        app.config = {
            'input_folder': Path(ind), 'tree_file': EASYPAML_TREE, 'output_folder': Path(outd),
            'models': [model], 'timeout': 1600, 'run_lrt': False, 'n_workers': 1,
            'auto_prune_tree': True,
        }
        t0 = time.time()
        app.run_batch_analysis()
        elapsed = time.time() - t0

        result_file = next(Path(outd, model).glob('*_results.txt'), None)
        m = LNL_RE.search(result_file.read_text(errors='ignore')) if result_file else None
        return {'elapsed': elapsed, 'lnL': float(m.group(3)) if m else None, 'np': int(m.group(2)) if m else None}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--loci', required=True, help='comma-separated locus IDs (must have all 24 taxa)')
    ap.add_argument('--model', default='M1a')
    args = ap.parse_args()

    print(f"{'locus':28s} {'raw_s':>8s} {'ep_s':>8s} {'overhead_s':>10s} {'raw_lnL':>16s} {'ep_lnL':>16s} {'match':>8s}")
    for locus in args.loci.split(','):
        locus = locus.strip()
        raw = run_raw_codeml(locus, args.model)
        ep = run_easypaml(locus, args.model)
        match = raw['lnL'] is not None and ep['lnL'] is not None and abs(raw['lnL'] - ep['lnL']) < 0.5
        print(f"{locus:28s} {raw['elapsed']:8.1f} {ep['elapsed']:8.1f} {ep['elapsed'] - raw['elapsed']:10.1f} "
              f"{raw['lnL']:16.4f} {ep['lnL']:16.4f} {'OK' if match else 'MISMATCH'}")


if __name__ == '__main__':
    main()
