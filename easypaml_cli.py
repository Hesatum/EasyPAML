#!/usr/bin/env python3
"""
easypaml_cli.py -- modo batch headless do EasyPAML, sem GUI.

Pensado para rodar em servidor/HPC via nohup/screen/slurm, com todos os
parametros explicitos na linha de comando (ou num arquivo --config), pra
reprodutibilidade total (o comando usado pode ir direto na secao de
Metodos de um manuscrito).

Exemplos:

  # Direto por flags
  python3 easypaml_cli.py \\
      --input exemplos_teste/amostras --tree exemplos_teste/final-tree.txt \\
      --output resultados/ --models M1a,M2a,M7,M8 --workers 24 --skip-beb

  # Via arquivo de config (equivalente, mais facil de arquivar/citar)
  python3 easypaml_cli.py --config minha_run.json

--skip-beb interrompe M2a/M8 assim que o CODEML imprime "BEBing..." (a
classificacao de sitio por sitio, a etapa mais cara -- o proprio CODEML avisa
que pode levar varios minutos por locus). O lnL/np/omega usados pelo LRT ja
foram escritos no outfile antes disso, entao o teste de selecao positiva
continua valido; so a tabela BEB de sitios fica ausente (o parser de
resultados cai automaticamente para NEB, que roda antes do BEB e fica
preservado). Por padrao BEB roda normalmente (mais confiavel para o proprio
resultado por sitio) -- so pule se o volume de loci tornar isso proibitivo.

--two-pass automatiza a mesma ideia pro dataset inteiro: passada 1 roda
todos os genes com --skip-beb (rapido, so pra ter o LRT); passada 2 reroda
so os genes com LRT significativo (p < --sig-threshold, default 0.05) com
BEB completo. Na maioria dos datasets a maior parte dos genes nao rejeita
o nulo -- essa e a fatia de BEB que fica pulada sem perder nenhum gene de
interesse real. Requer M1a+M2a e/ou M7+M8 em --models (precisa do par pra
calcular LRT).
"""
import argparse
import json
import shutil
import sys
from pathlib import Path

from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from src.backend.codeml_backend import CodemlBatchAnalysis

VALID_MODELS = {'M0', 'M1a', 'M2a', 'M7', 'M8', 'Branch', 'Branch-site', 'Branch-site_null'}
LRT_DF = 2  # M1a vs M2a e M7 vs M8 sempre tem df=2 (ver _CANONICAL_DF em codeml_backend.py)


def parse_args():
    ap = argparse.ArgumentParser(
        description="EasyPAML -- batch de analises CODEML sem interface grafica.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument('--config', type=Path, help="arquivo JSON com todos os parametros abaixo (sobrescreve as flags)")
    ap.add_argument('--input', type=Path, help="pasta com .fas/.fasta/.phy/.phylip (um arquivo por gene)")
    ap.add_argument('--tree', type=Path, help="arquivo de arvore Newick (com ou sem cabecalho 'N  1')")
    ap.add_argument('--output', type=Path, help="pasta de saida")
    ap.add_argument('--models', default='M1a,M2a,M7,M8', help="modelos separados por virgula (default: M1a,M2a,M7,M8)")
    ap.add_argument('--workers', type=int, default=4, help="genes em paralelo (default: 4)")
    ap.add_argument('--timeout', type=int, default=1600, help="timeout por execucao codeml, em segundos (default: 1600)")
    ap.add_argument('--no-lrt', action='store_true', help="nao calcular LRT automaticamente no final")
    ap.add_argument('--skip-beb', action='store_true', help="interrompe M2a/M8 antes do BEB (mantem LRT, perde tabela de sitio BEB -- ver docstring)")
    ap.add_argument('--no-prune-tree', action='store_true', help="desativa poda automatica da arvore por locus (default: poda ativada)")
    ap.add_argument('--two-pass', action='store_true', help="passada 1 sem BEB em todos os genes, passada 2 com BEB so nos LRT-significativos (ver docstring)")
    ap.add_argument('--sig-threshold', type=float, default=0.05, help="p-valor de corte pro --two-pass (default: 0.05)")
    ap.add_argument('--warm-start-m0', action='store_true',
                     help="roda M0 escondido por gene pra usar como ponto de partida nos modelos de sitio, "
                          "com multi-start automatico de omega (3 pontos de partida, warm_start_multistart=True "
                          "por default) pra reduzir risco de otimo local. Medido em 20 loci reais (M1a): "
                          "17.6x mais rapido; 2/20 genes com lnL levemente pior que o from-scratch (maior "
                          "diferenca: 1.47), 2/20 genes com lnL MELHOR (multistart escapou de otimo que o "
                          "from-scratch nao escapou). Sem multi-start (so omega=0.5) o risco era maior: "
                          "~1/3 dos loci de um teste anterior menor. Desligado por padrao mesmo assim -- "
                          "nao e garantia matematica de resultado identico, so estatisticamente muito mais "
                          "raro de divergir. Ligue sabendo do trade-off.")
    args = ap.parse_args()

    if args.config:
        with open(args.config, encoding='utf-8') as fh:
            cfg = json.load(fh)
        for key in ('input', 'tree', 'output'):
            if key in cfg:
                cfg[key] = Path(cfg[key])
        return cfg

    if not (args.input and args.tree and args.output):
        ap.error("--input, --tree e --output sao obrigatorios (ou use --config)")

    return {
        'input': args.input,
        'tree': args.tree,
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
    }


def _make_app(cfg, input_folder, output_folder, models, skip_beb):
    app = CodemlBatchAnalysis()
    app.config = {
        'input_folder': input_folder, 'tree_file': cfg['tree'], 'output_folder': output_folder,
        'models': models, 'timeout': cfg['timeout'], 'run_lrt': cfg['run_lrt'],
        'n_workers': cfg['workers'], 'auto_prune_tree': cfg['auto_prune_tree'], 'skip_beb': skip_beb,
        'warm_start_m0': cfg.get('warm_start_m0', False),
    }
    return app


def run_two_pass(cfg):
    """Passada 1 (skip_beb) em todo mundo -> LRT -> passada 2 (BEB completo)
    so nos genes significativos. Nao reimplementa nada do backend, so chama
    run_batch_analysis() duas vezes com config diferente."""
    beb_models = [m for m in cfg['models'] if m in ('M2a', 'M8')]
    if not beb_models:
        sys.exit("--two-pass so faz sentido com M2a e/ou M8 em --models (sao os unicos com BEB).")

    pass1_dir = cfg['output'] / 'pass1_screen'
    pass1_dir.mkdir(parents=True, exist_ok=True)
    print("### PASSADA 1/2 -- todos os genes, sem BEB (so LRT) ###\n")
    app1 = _make_app(cfg, cfg['input'], pass1_dir, cfg['models'], skip_beb=True)
    app1.run_batch_analysis()

    import pandas as pd
    df = pd.read_csv(pass1_dir / 'analysis_summary.tsv', sep='\t')
    sig_genes = set()
    for null, alt, col in (('M1a', 'M2a', 'lrt_M1a_vs_M2a'), ('M7', 'M8', 'lrt_M7_vs_M8')):
        if null in cfg['models'] and alt in cfg['models'] and col in df.columns:
            stat = df[col].dropna()
            pvals = stats.chi2.sf(stat.clip(lower=0), df=LRT_DF)
            sig_genes |= set(df.loc[stat.index[pvals < cfg['sig_threshold']], 'Gene'])

    print(f"\n### {len(sig_genes)}/{len(df)} genes com LRT significativo (p<{cfg['sig_threshold']}) -- rerodando com BEB ###\n")
    if not sig_genes:
        print("Nenhum gene significativo -- passada 2 nao tem o que fazer.")
        return

    pass2_input = cfg['output'] / 'pass2_input'
    pass2_input.mkdir(parents=True, exist_ok=True)
    for gene in sig_genes:
        src = next(cfg['input'].glob(f'{gene}.*'), None)
        if src:
            shutil.copy(src, pass2_input / src.name)

    pass2_dir = cfg['output'] / 'pass2_beb'
    pass2_dir.mkdir(parents=True, exist_ok=True)
    app2 = _make_app(cfg, pass2_input, pass2_dir, beb_models, skip_beb=False)
    app2.run_batch_analysis()

    print(f"\nScreen completo (todos os genes, sem BEB): {pass1_dir}")
    print(f"BEB detalhado (so os {len(sig_genes)} significativos): {pass2_dir}")


def main():
    cfg = parse_args()

    bad_models = set(cfg['models']) - VALID_MODELS
    if bad_models:
        sys.exit(f"Modelo(s) invalido(s): {sorted(bad_models)}. Validos: {sorted(VALID_MODELS)}")
    if not cfg['input'].is_dir():
        sys.exit(f"Pasta de input nao existe: {cfg['input']}")
    if not cfg['tree'].is_file():
        sys.exit(f"Arquivo de arvore nao existe: {cfg['tree']}")

    cfg['output'].mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("EasyPAML CLI -- configuracao")
    print("=" * 80)
    for k, v in cfg.items():
        print(f"  {k:16s}: {v}")
    print("=" * 80 + "\n")

    # Grava a config efetivamente usada junto com os resultados -- reprodutibilidade
    # (permite citar exatamente esse arquivo nos Metodos, ou refazer o run identico).
    with open(cfg['output'] / 'run_config.json', 'w', encoding='utf-8') as fh:
        json.dump({k: str(v) if isinstance(v, Path) else v for k, v in cfg.items()}, fh, indent=2, ensure_ascii=False)

    try:
        if cfg['two_pass']:
            run_two_pass(cfg)
        else:
            app = _make_app(cfg, cfg['input'], cfg['output'], cfg['models'], cfg['skip_beb'])
            app.run_batch_analysis()
    except KeyboardInterrupt:
        print("\nInterrompido pelo usuario (Ctrl+C). Resultados parciais ja estao em disco.")
        sys.exit(130)


if __name__ == '__main__':
    main()
