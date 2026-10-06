#!/usr/bin/env python3
"""Compare EasyPAML with plain codeml on the same genes.

For each alignment, runs EasyPAML for one model, then runs the .ctl that
EasyPAML saved directly with codeml, and prints both run times and lnL values.
The lnL values should agree up to codeml's random starting values (see METHODS.md).

    python tools/benchmark_vs_raw.py --input ALIGNMENTS --tree TREE.nwk --model M8 [--codeml PATH]
"""
import argparse
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.backend.codeml_backend import CodemlBatchAnalysis, find_codeml  # noqa: E402
from src.backend.preflight import group_by_gene, list_alignment_files  # noqa: E402

LNL_RE = re.compile(r'lnL\(ntime:\s*\d+\s+np:\s*\d+\):\s*([-\d.]+)')


def lnl(path: Path):
    m = LNL_RE.search(path.read_text(errors='ignore')) if path.exists() else None
    return float(m.group(1)) if m else None


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--input', type=Path, required=True, help='folder with alignments')
    ap.add_argument('--tree', type=Path, required=True)
    ap.add_argument('--model', default='M8')
    ap.add_argument('--codeml', help='codeml executable (default: same search as EasyPAML)')
    args = ap.parse_args()
    codeml = find_codeml(args.codeml)

    genes, _ = group_by_gene(list_alignment_files(args.input))
    print(f"{'gene':30s} {'easypaml_s':>10s} {'codeml_s':>9s} {'easypaml_lnL':>14s} {'codeml_lnL':>14s}")
    for gene, path in sorted(genes.items()):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / 'in').mkdir()
            shutil.copy(path, tmp / 'in')
            app = CodemlBatchAnalysis()
            app.config = {'input_folder': tmp / 'in', 'tree_file': args.tree,
                          'output_folder': tmp / 'out', 'models': [args.model],
                          'n_workers': 1, 'run_lrt': False, 'codeml_path': codeml,
                          'log_callback': lambda *a: None}
            t0 = time.time()
            app.run_batch_analysis()
            t_wrapper = time.time() - t0

            model_dir = tmp / 'out' / args.model
            ctl = next(model_dir.glob('*.ctl'), None)
            ours = next(model_dir.glob('*_results.txt'), None)
            if ctl is None or ours is None:
                print(f"{gene:30s} EasyPAML run failed, see {model_dir}")
                continue
            raw_dir = tmp / 'raw'
            shutil.copytree(model_dir, raw_dir)
            (raw_dir / ours.name).unlink()
            t0 = time.time()
            subprocess.run([codeml, ctl.name], cwd=raw_dir, stdin=subprocess.DEVNULL,
                           capture_output=True, text=True)
            t_raw = time.time() - t0
            print(f"{gene:30s} {t_wrapper:10.1f} {t_raw:9.1f} {lnl(ours) or float('nan'):14.4f} "
                  f"{lnl(raw_dir / ours.name) or float('nan'):14.4f}")


if __name__ == '__main__':
    main()
