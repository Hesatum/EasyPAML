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
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from src.backend.codeml_backend import CodemlBatchAnalysis

VALID_MODELS = {'M0', 'M1a', 'M2a', 'M7', 'M8', 'Branch', 'Branch-site', 'Branch-site_null'}


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
    }


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

    app = CodemlBatchAnalysis()
    app.config = {
        'input_folder': cfg['input'],
        'tree_file': cfg['tree'],
        'output_folder': cfg['output'],
        'models': cfg['models'],
        'timeout': cfg['timeout'],
        'run_lrt': cfg['run_lrt'],
        'n_workers': cfg['workers'],
        'auto_prune_tree': cfg['auto_prune_tree'],
        'skip_beb': cfg['skip_beb'],
    }

    # Grava a config efetivamente usada junto com os resultados -- reprodutibilidade
    # (permite citar exatamente esse arquivo nos Metodos, ou refazer o run identico).
    with open(cfg['output'] / 'run_config.json', 'w', encoding='utf-8') as fh:
        json.dump({k: str(v) if isinstance(v, Path) else v for k, v in cfg.items()}, fh, indent=2, ensure_ascii=False)

    try:
        app.run_batch_analysis()
    except KeyboardInterrupt:
        print("\nInterrompido pelo usuario (Ctrl+C). Resultados parciais ja estao em disco.")
        sys.exit(130)


if __name__ == '__main__':
    main()
