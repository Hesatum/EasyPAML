#!/usr/bin/env python3
"""
EasyPAML command-line mode, for servers and many genes.

Examples:

  ./easypaml-cli.sh --input examples/alignments --tree examples/tree.nwk \\
      --output results/ --models M7,M8,M8a --workers 24

  ./easypaml-cli.sh --config my_run.json    # the same options as JSON

(easypaml-cli.bat on Windows; both are created by the installer.)

--skip-beb stops M2a and M8 when codeml starts the BEB step, the slowest part.
lnL, np and omega, which the LRT uses, are already written by then, so the test
stays valid; only the BEB site table is lost (NEB is used instead).

--two-pass runs every gene with --skip-beb, applies the LRT with
Benjamini-Hochberg correction, then reruns with BEB only the genes with
q < --sig-threshold. It needs M2a and/or M8 in --models.

Details: METHODS.md.
"""
import argparse
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.backend import messages
from src.backend.codeml_backend import CodemlBatchAnalysis, codeml_version, find_codeml
from src.backend.ctl_params import CODONFREQ_OPTIONS, DEFAULT_CODONFREQ, codonfreq_label
from src.backend.preflight import run_preflight
from src.backend.version import __version__, version_string

VALID_MODELS = {'M0', 'M1a', 'M2a', 'M7', 'M8', 'M8a', 'Branch', 'Branch-site', 'Branch-site_null'}
_CODONFREQ_HELP = ", ".join(f"{v}={n}" for v, n, _ in CODONFREQ_OPTIONS)


def parse_args():
    ap = argparse.ArgumentParser(
        prog="easypaml-cli",
        description="EasyPAML: batch codeml analyses without the window.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument('--version', action='version', version=f"EasyPAML {version_string()}")
    ap.add_argument('--config', type=Path, help="JSON file with the options below (overrides the flags)")
    ap.add_argument('--input', type=Path, help="folder with .fas/.fasta/.phy/.phylip files, one per gene")
    ap.add_argument('--tree', type=Path, help="Newick tree file (with or without an 'N  1' header); "
                                              "optional if every gene has GENE.nwk next to its alignment")
    ap.add_argument('--tree-folder', type=Path, help="folder with one tree per gene (GENE.nwk/.tree/.tre), "
                                                     "matched by the alignment file name")
    ap.add_argument('--output', type=Path, help="output folder")
    ap.add_argument('--models', default='M1a,M2a,M7,M8,M8a',
                    help="comma-separated models (default: M1a,M2a,M7,M8,M8a)")
    ap.add_argument('--no-m8a', action='store_true',
                    help="do not run M8a (the second null of M8; no M8a vs M8 test)")
    ap.add_argument('--no-auto-nulls', action='store_true',
                    help="run only the listed models; by default the null model of each listed "
                         "alternative is added, as in the window (M8 adds M7 and M8a, M2a adds M1a)")
    ap.add_argument('--codonfreq', type=int, default=DEFAULT_CODONFREQ,
                    help=f"codeml CodonFreq (default: {DEFAULT_CODONFREQ} = F3x4). Options: {_CODONFREQ_HELP}")
    ap.add_argument('--ncatg', type=int, default=10, help="beta categories in M7/M8/M8a (default: 10)")
    ap.add_argument('--kappa', type=float, default=2.0, help="initial kappa (default: 2)")
    ap.add_argument('--omega', type=float, default=0.5, help="initial omega (default: 0.5)")
    ap.add_argument('--cleandata', type=int, choices=(0, 1), default=1,
                    help="1 = drop columns with gaps, ambiguities or stop codons (default); 0 = keep them")
    ap.add_argument('--ignore-stop-codons', action='store_true',
                    help="run genes with an internal stop codon (codeml treats the column as missing data). "
                         "Without it those genes fail and the stop position is reported.")
    ap.add_argument('--codeml', type=Path, help="codeml executable (default: EASYPAML_CODEML, bin/codeml, "
                                                "or codeml on PATH)")
    ap.add_argument('--idle-timeout', type=int, default=300,
                    help="stop codeml after N s without CPU use (default: 300; 0 turns it off)")
    ap.add_argument('--rerun-all', action='store_true',
                    help="run every model again, even those already in OUTPUT with the same "
                         "alignment, tree and .ctl (by default their results are reused)")
    ap.add_argument('--strict', action='store_true',
                    help="run nothing if the data check finds errors or warnings")
    ap.add_argument('--lang', choices=('pt', 'en'), help="language of the messages (default: en)")
    ap.add_argument('--verbose', action='store_true', help="show debug messages")
    ap.add_argument('--workers', type=int, default=4, help="genes run in parallel (default: 4)")
    ap.add_argument('--timeout', type=int, default=0,
                    help="time limit per codeml run, in seconds. 0 (default) = automatic, "
                         "scaled to the model and the gene's taxa and codons (docs/timing_benchmark.md)")
    ap.add_argument('--no-lrt', action='store_true', help="do not compute the LRT at the end")
    ap.add_argument('--skip-beb', action='store_true', help="stop M2a/M8 before BEB (keeps the LRT, loses the BEB site table)")
    ap.add_argument('--no-prune-tree', action='store_true', help="do not prune the tree for each gene (pruning is on by default)")
    ap.add_argument('--two-pass', action='store_true', help="pass 1 without BEB on every gene, pass 2 with BEB only on significant genes")
    ap.add_argument('--sig-threshold', type=float, default=0.05, help="BH q-value cut-off for --two-pass (default: 0.05)")
    ap.add_argument('--warm-start-m0', action='store_true',
                    help="fit M0 first and use its kappa and branch lengths as starting values for the "
                         "site models, with omega started from 0.2, 1.0 and 2.5 (best lnL kept). Faster, "
                         "but not guaranteed to match a fit from scratch. Off by default.")
    args = ap.parse_args()

    if args.config:
        with open(args.config, encoding='utf-8') as fh:
            cfg = json.load(fh)

        missing = [k for k in ('input', 'output') if k not in cfg]
        if missing:
            ap.error(f"--config {args.config}: missing required key(s) {missing}")
        for key in ('input', 'tree', 'output', 'tree_folder'):
            if cfg.get(key):
                cfg[key] = Path(cfg[key])
        cfg.setdefault('tree', None)
        cfg.setdefault('tree_folder', None)

        # same defaults as the flags, so a minimal config (input/tree/output) works
        cfg.setdefault('models', ['M1a', 'M2a', 'M7', 'M8', 'M8a'])
        cfg.setdefault('codonfreq', DEFAULT_CODONFREQ)
        cfg.setdefault('ncatg', 10)
        cfg.setdefault('kappa', 2.0)
        cfg.setdefault('omega', 0.5)
        cfg.setdefault('cleandata', 1)
        cfg.setdefault('ignore_stop_codons', False)
        cfg.setdefault('codeml', None)
        cfg.setdefault('idle_timeout', 300)
        cfg.setdefault('strict', False)
        cfg.setdefault('lang', args.lang)
        cfg.setdefault('verbose', args.verbose)
        cfg.setdefault('no_m8a', args.no_m8a)
        cfg.setdefault('auto_nulls', not args.no_auto_nulls)
        cfg.setdefault('workers', 4)
        cfg.setdefault('timeout', 0)
        cfg.setdefault('run_lrt', True)
        cfg.setdefault('skip_beb', False)
        cfg.setdefault('auto_prune_tree', True)
        cfg.setdefault('two_pass', False)
        cfg.setdefault('sig_threshold', 0.05)
        cfg.setdefault('warm_start_m0', False)
        # models may be a JSON list or a "A,B,C" string as in --models
        if isinstance(cfg['models'], str):
            cfg['models'] = [m.strip() for m in cfg['models'].split(',') if m.strip()]

        return cfg

    if not (args.input and args.output):
        ap.error("--input and --output are required (or use --config)")

    return {
        'input': args.input,
        'tree': args.tree,
        'tree_folder': args.tree_folder,
        'output': args.output,
        'models': [m.strip() for m in args.models.split(',') if m.strip()],
        'workers': args.workers,
        'timeout': args.timeout,
        'run_lrt': not args.no_lrt,
        'skip_beb': args.skip_beb,
        'auto_prune_tree': not args.no_prune_tree,
        'two_pass': args.two_pass,
        'sig_threshold': args.sig_threshold,
        'warm_start_m0': args.warm_start_m0,
        'codonfreq': args.codonfreq,
        'ncatg': args.ncatg,
        'kappa': args.kappa,
        'omega': args.omega,
        'cleandata': args.cleandata,
        'ignore_stop_codons': args.ignore_stop_codons,
        'reuse_results': not args.rerun_all,
        'codeml': str(args.codeml) if args.codeml else None,
        'idle_timeout': args.idle_timeout,
        'strict': args.strict,
        'lang': args.lang,
        'verbose': args.verbose,
        'no_m8a': args.no_m8a,
        'auto_nulls': not args.no_auto_nulls,
    }


def _make_app(cfg, input_folder, output_folder, models, skip_beb):
    app = CodemlBatchAnalysis()
    app.config = {
        'input_folder': input_folder, 'tree_file': cfg['tree'], 'output_folder': output_folder,
        'tree_folder': cfg.get('tree_folder'),
        'models': models, 'timeout': cfg['timeout'], 'run_lrt': cfg['run_lrt'],
        'n_workers': cfg['workers'], 'auto_prune_tree': cfg['auto_prune_tree'], 'skip_beb': skip_beb,
        'warm_start_m0': cfg.get('warm_start_m0', False),
        'CodonFreq': cfg.get('codonfreq', DEFAULT_CODONFREQ),
        'ncatG': cfg.get('ncatg', 10),
        'kappa': cfg.get('kappa', 2.0),
        'omega': cfg.get('omega', 0.5),
        'cleandata': cfg.get('cleandata', 1),
        'ignore_stop_codons': cfg.get('ignore_stop_codons', False),
        'codeml_path': cfg.get('codeml'),
        'idle_timeout': cfg.get('idle_timeout', 300),
        'verbose': cfg.get('verbose', False),
        'interface': 'cli',
    }
    return app


def run_two_pass(cfg):
    """Pass 1 (skip_beb) on every gene, LRT, then pass 2 (full BEB) on the
    significant genes only."""
    beb_models = [m for m in cfg['models'] if m in ('M2a', 'M8')]
    if not beb_models:
        sys.exit("--two-pass needs M2a and/or M8 in --models (the models with BEB).")

    pass1_dir = cfg['output'] / 'pass1_screen'
    pass1_dir.mkdir(parents=True, exist_ok=True)
    print("### Pass 1/2: every gene, without BEB (LRT only) ###\n")
    app1 = _make_app(cfg, cfg['input'], pass1_dir, cfg['models'], skip_beb=True)
    summary1 = app1.run_batch_analysis()

    import pandas as pd
    df = pd.read_csv(pass1_dir / 'analysis_summary.tsv', sep='\t')
    sig_genes = set()
    for null, alt in (('M1a', 'M2a'), ('M7', 'M8'), ('M8a', 'M8')):
        q_col = f'q_{null}_vs_{alt}'
        if null in cfg['models'] and alt in cfg['models'] and q_col in df.columns:
            qvals = pd.to_numeric(df[q_col], errors='coerce')
            sig_genes |= set(df.loc[qvals < cfg['sig_threshold'], 'Gene'])

    print(f"\n### {len(sig_genes)}/{len(df)} genes significant (BH q < {cfg['sig_threshold']}); rerunning them with BEB ###\n")
    if not sig_genes:
        print("No significant gene; pass 2 has nothing to do.")
        return summary1

    pass2_input = cfg['output'] / 'pass2_input'
    pass2_input.mkdir(parents=True, exist_ok=True)
    for gene in sig_genes:
        # the alignment and, if present, the gene's own tree (GENE.nwk)
        for src in cfg['input'].glob(f'{gene}.*'):
            if src.stem == gene and src.is_file():
                shutil.copy(src, pass2_input / src.name)

    pass2_dir = cfg['output'] / 'pass2_beb'
    pass2_dir.mkdir(parents=True, exist_ok=True)
    app2 = _make_app(cfg, pass2_input, pass2_dir, beb_models, skip_beb=False)
    summary2 = app2.run_batch_analysis()

    print(f"\nScreen (every gene, without BEB): {pass1_dir}")
    print(f"BEB for the {len(sig_genes)} significant genes: {pass2_dir}")
    return {'failed': summary1.get('failed', 0) + summary2.get('failed', 0)}


def main():
    cfg = parse_args()
    messages.set_language(cfg.get('lang') or 'en')
    lang = messages.get_language()

    bad_models = set(cfg['models']) - VALID_MODELS
    if bad_models:
        sys.exit(f"Unknown model(s): {sorted(bad_models)}. Valid: {sorted(VALID_MODELS)}")
    if cfg.get('auto_nulls', True):
        listed = list(cfg['models'])
        cfg['models'] = CodemlBatchAnalysis.auto_complete_null_models(
            listed, include_neutral=True, include_m8a=not cfg.get('no_m8a'))
        added = [m for m in cfg['models'] if m not in listed]
        if added:
            if lang == 'pt':
                print(f"Modelos nulos acrescentados para completar os testes: {', '.join(added)}. "
                      "Para rodar só os modelos listados, use --no-auto-nulls.")
            else:
                print(f"Null models added to complete the tests: {', '.join(added)}. "
                      "To run only the models you listed, use --no-auto-nulls.")
    if cfg.get('no_m8a') and 'M8a' in cfg['models']:
        cfg['models'] = [m for m in cfg['models'] if m != 'M8a']
    if not cfg['input'].is_dir():
        sys.exit(f"Input folder not found: {cfg['input']}")
    if cfg.get('tree') and not cfg['tree'].is_file():
        sys.exit(f"Tree file not found: {cfg['tree']}")
    from src.backend.preflight import discover_per_gene_trees
    per_gene = discover_per_gene_trees(cfg['input'], cfg.get('tree_folder'))
    if not cfg.get('tree') and not per_gene:
        sys.exit("Give --tree, or put a GENE.nwk tree for each gene in the folder or in --tree-folder")

    cfg['output'].mkdir(parents=True, exist_ok=True)

    codeml_path = find_codeml(cfg.get('codeml'))
    print("=" * 72)
    print(f"EasyPAML {version_string()} -- CLI")
    print("=" * 72)
    for k, v in cfg.items():
        if k == 'tree' and v is None and per_gene:
            v = (f"{len(per_gene)} árvore(s) por gene (GENE.nwk)" if lang == 'pt'
                 else f"{len(per_gene)} per-gene tree(s) (GENE.nwk)")
        print(f"  {k:18s}: {v}")
    print(f"  {'codeml (resolved)':18s}: {codeml_path} (version {codeml_version(codeml_path)})")
    print(f"  {'CodonFreq':18s}: {codonfreq_label(cfg.get('codonfreq', DEFAULT_CODONFREQ))}")
    print("=" * 72 + "\n")

    # data check: stop codons, names missing from the tree, pruned taxa,
    # duplicate files, length not a multiple of 3
    report = run_preflight(cfg['input'], cfg.get('tree'), auto_prune=cfg['auto_prune_tree'],
                           ignore_stop_codons=cfg.get('ignore_stop_codons', False),
                           per_gene_trees=per_gene, tree_folder=cfg.get('tree_folder'))
    text = report.format_text(lang, include_info=cfg.get('verbose', False))
    if text:
        print("Verificação dos dados:" if lang == 'pt' else "Data check:")
        print(text + "\n")
    if cfg.get('strict') and report.has_problems:
        sys.exit(2)

    try:
        if cfg['two_pass']:
            summary = run_two_pass(cfg) or {}
        else:
            app = _make_app(cfg, cfg['input'], cfg['output'], cfg['models'], cfg['skip_beb'])
            summary = app.run_batch_analysis() or {}
    except KeyboardInterrupt:
        print("\nStopped (Ctrl+C). Partial results are already on disk.")
        sys.exit(130)
    sys.exit(1 if summary.get('failed') else 0)


if __name__ == '__main__':
    main()
