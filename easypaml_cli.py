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
so os genes com LRT significativo apos correcao Benjamini-Hochberg
(q < --sig-threshold, default 0.05) com BEB completo. Na maioria dos
datasets a maior parte dos genes nao rejeita o nulo -- essa e a fatia de
BEB que fica pulada sem perder nenhum gene de interesse real. Requer
M1a+M2a e/ou M7+M8 em --models (precisa do par pra calcular LRT). O
q-valor (nao o p bruto) e o corte porque a passada 1 testa todos os genes
do dataset simultaneamente -- sem correcao de multiplos testes o p bruto
infla falsos positivos nessa escala.
"""
import argparse
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from src.backend import messages
from src.backend.codeml_backend import CodemlBatchAnalysis, codeml_version, find_codeml
from src.backend.ctl_params import CODONFREQ_OPTIONS, DEFAULT_CODONFREQ, codonfreq_label
from src.backend.preflight import run_preflight
from src.backend.version import __version__

VALID_MODELS = {'M0', 'M1a', 'M2a', 'M7', 'M8', 'M8a', 'Branch', 'Branch-site', 'Branch-site_null'}
_CODONFREQ_HELP = ", ".join(f"{v}={n}" for v, n, _ in CODONFREQ_OPTIONS)


def parse_args():
    ap = argparse.ArgumentParser(
        description="EasyPAML -- batch de analises CODEML sem interface grafica.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument('--version', action='version', version=f"EasyPAML {__version__}")
    ap.add_argument('--config', type=Path, help="arquivo JSON com todos os parametros abaixo (sobrescreve as flags)")
    ap.add_argument('--input', type=Path, help="pasta com .fas/.fasta/.phy/.phylip (um arquivo por gene)")
    ap.add_argument('--tree', type=Path, help="arquivo de arvore Newick (com ou sem cabecalho 'N  1'); "
                                              "opcional se cada gene tiver GENE.nwk ao lado do alinhamento")
    ap.add_argument('--tree-folder', type=Path, help="pasta com uma arvore por gene (GENE.nwk/.tree/.tre), "
                                                     "pareada pelo nome do arquivo do alinhamento")
    ap.add_argument('--output', type=Path, help="pasta de saida")
    ap.add_argument('--models', default='M1a,M2a,M7,M8,M8a',
                    help="modelos separados por virgula (default: M1a,M2a,M7,M8,M8a)")
    ap.add_argument('--no-m8a', action='store_true',
                    help="nao roda o M8a (nulo extra do M8; teste M8a vs M8). Por padrao o M8a roda "
                         "junto com M7 e M8")
    ap.add_argument('--codonfreq', type=int, default=DEFAULT_CODONFREQ,
                    help=f"CodonFreq do codeml (default: {DEFAULT_CODONFREQ} = F3x4). Opcoes: {_CODONFREQ_HELP}")
    ap.add_argument('--ncatg', type=int, default=10, help="categorias da beta em M7/M8 (default: 10)")
    ap.add_argument('--kappa', type=float, default=2.0, help="kappa inicial (default: 2)")
    ap.add_argument('--omega', type=float, default=0.5, help="omega inicial (default: 0.5)")
    ap.add_argument('--cleandata', type=int, choices=(0, 1), default=1,
                    help="1 = remove colunas com gap/ambiguidade/stop (default); 0 = mantem")
    ap.add_argument('--ignore-stop-codons', action='store_true',
                    help="roda genes com stop codon interno (o codeml trata a coluna como dado ausente). "
                         "Sem esta opcao esses genes sao marcados como FALHOU, com a posicao do stop.")
    ap.add_argument('--codeml', type=Path, help="caminho do executavel codeml (default: bin/codeml, "
                                                "variavel EASYPAML_CODEML ou codeml no PATH)")
    ap.add_argument('--idle-timeout', type=int, default=300,
                    help="encerra o codeml se ficar N s sem usar CPU (default: 300; 0 desliga)")
    ap.add_argument('--strict', action='store_true',
                    help="nao roda nada se a verificacao inicial encontrar erros ou avisos")
    ap.add_argument('--lang', choices=('pt', 'en'), help="idioma das mensagens (default: idioma do sistema)")
    ap.add_argument('--verbose', action='store_true', help="mostra mensagens de depuracao")
    ap.add_argument('--workers', type=int, default=4, help="genes em paralelo (default: 4)")
    ap.add_argument('--timeout', type=int, default=1600, help="timeout por execucao codeml, em segundos (default: 1600)")
    ap.add_argument('--no-lrt', action='store_true', help="nao calcular LRT automaticamente no final")
    ap.add_argument('--skip-beb', action='store_true', help="interrompe M2a/M8 antes do BEB (mantem LRT, perde tabela de sitio BEB -- ver docstring)")
    ap.add_argument('--no-prune-tree', action='store_true', help="desativa poda automatica da arvore por locus (default: poda ativada)")
    ap.add_argument('--two-pass', action='store_true', help="passada 1 sem BEB em todos os genes, passada 2 com BEB so nos LRT-significativos (ver docstring)")
    ap.add_argument('--sig-threshold', type=float, default=0.05, help="q-valor (BH) de corte pro --two-pass (default: 0.05)")
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

        missing = [k for k in ('input', 'output') if k not in cfg]
        if missing:
            ap.error(f"--config {args.config}: faltando chave(s) obrigatoria(s) {missing}")
        for key in ('input', 'tree', 'output', 'tree_folder'):
            if cfg.get(key):
                cfg[key] = Path(cfg[key])
        cfg.setdefault('tree', None)
        cfg.setdefault('tree_folder', None)

        # Mesmos defaults do caminho via flags -- sem isso, uma config.json
        # minima (so input/tree/output) quebra com KeyError la na frente em
        # vez de rodar com o comportamento padrao esperado.
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
        cfg.setdefault('workers', 4)
        cfg.setdefault('timeout', 1600)
        cfg.setdefault('run_lrt', True)
        cfg.setdefault('skip_beb', False)
        cfg.setdefault('auto_prune_tree', True)
        cfg.setdefault('two_pass', False)
        cfg.setdefault('sig_threshold', 0.05)
        cfg.setdefault('warm_start_m0', False)
        # Aceita tanto lista (forma natural em JSON) quanto string "A,B,C"
        # (pra quem copiar o valor direto de --models)
        if isinstance(cfg['models'], str):
            cfg['models'] = [m.strip() for m in cfg['models'].split(',') if m.strip()]

        return cfg

    if not (args.input and args.output):
        ap.error("--input e --output sao obrigatorios (ou use --config)")

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
        'codeml': str(args.codeml) if args.codeml else None,
        'idle_timeout': args.idle_timeout,
        'strict': args.strict,
        'lang': args.lang,
        'verbose': args.verbose,
        'no_m8a': args.no_m8a,
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
    summary1 = app1.run_batch_analysis()

    import pandas as pd
    df = pd.read_csv(pass1_dir / 'analysis_summary.tsv', sep='\t')
    sig_genes = set()
    for null, alt in (('M1a', 'M2a'), ('M7', 'M8'), ('M8a', 'M8')):
        q_col = f'q_{null}_vs_{alt}'
        if null in cfg['models'] and alt in cfg['models'] and q_col in df.columns:
            qvals = pd.to_numeric(df[q_col], errors='coerce')
            sig_genes |= set(df.loc[qvals < cfg['sig_threshold'], 'Gene'])

    print(f"\n### {len(sig_genes)}/{len(df)} genes com LRT significativo (q BH<{cfg['sig_threshold']}) -- rerodando com BEB ###\n")
    if not sig_genes:
        print("Nenhum gene significativo -- passada 2 nao tem o que fazer.")
        return summary1

    pass2_input = cfg['output'] / 'pass2_input'
    pass2_input.mkdir(parents=True, exist_ok=True)
    for gene in sig_genes:
        # alinhamento e, se houver, a árvore própria do gene (GENE.nwk)
        for src in cfg['input'].glob(f'{gene}.*'):
            if src.stem == gene and src.is_file():
                shutil.copy(src, pass2_input / src.name)

    pass2_dir = cfg['output'] / 'pass2_beb'
    pass2_dir.mkdir(parents=True, exist_ok=True)
    app2 = _make_app(cfg, pass2_input, pass2_dir, beb_models, skip_beb=False)
    summary2 = app2.run_batch_analysis()

    print(f"\nScreen completo (todos os genes, sem BEB): {pass1_dir}")
    print(f"BEB detalhado (so os {len(sig_genes)} significativos): {pass2_dir}")
    return {'failed': summary1.get('failed', 0) + summary2.get('failed', 0)}


def main():
    cfg = parse_args()
    messages.set_language(cfg.get('lang') or messages.system_language())
    lang = messages.get_language()

    if cfg.get('no_m8a') and 'M8a' in cfg['models']:
        cfg['models'] = [m for m in cfg['models'] if m != 'M8a']
    bad_models = set(cfg['models']) - VALID_MODELS
    if bad_models:
        sys.exit(f"Modelo(s) invalido(s): {sorted(bad_models)}. Validos: {sorted(VALID_MODELS)}")
    if not cfg['input'].is_dir():
        sys.exit(f"Pasta de input nao existe: {cfg['input']}")
    if cfg.get('tree') and not cfg['tree'].is_file():
        sys.exit(f"Arquivo de arvore nao existe: {cfg['tree']}")
    from src.backend.preflight import discover_per_gene_trees
    per_gene = discover_per_gene_trees(cfg['input'], cfg.get('tree_folder'))
    if not cfg.get('tree') and not per_gene:
        sys.exit("Informe --tree (ou ponha uma arvore GENE.nwk por gene na pasta / em --tree-folder)")

    cfg['output'].mkdir(parents=True, exist_ok=True)

    codeml_path = find_codeml(cfg.get('codeml'))
    print("=" * 72)
    print(f"EasyPAML {__version__} -- CLI")
    print("=" * 72)
    for k, v in cfg.items():
        print(f"  {k:18s}: {v}")
    print(f"  {'codeml (resolved)':18s}: {codeml_path} (version {codeml_version(codeml_path)})")
    print(f"  {'CodonFreq':18s}: {codonfreq_label(cfg.get('codonfreq', DEFAULT_CODONFREQ))}")
    print("=" * 72 + "\n")

    # Verificacao antes de rodar: stop codons (com posicao), nomes que nao
    # batem com a arvore, taxons podados, duplicados, comprimento % 3.
    report = run_preflight(cfg['input'], cfg.get('tree'), auto_prune=cfg['auto_prune_tree'],
                           ignore_stop_codons=cfg.get('ignore_stop_codons', False),
                           per_gene_trees=per_gene)
    text = report.format_text(lang, include_info=cfg.get('verbose', False))
    if text:
        print("Verificacao dos dados / Data check:" if lang == 'pt' else "Data check:")
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
        print("\nInterrompido pelo usuario (Ctrl+C). Resultados parciais ja estao em disco.")
        sys.exit(130)
    sys.exit(1 if summary.get('failed') else 0)


if __name__ == '__main__':
    main()
