# EasyPAML -- guia para agentes de IA

Wrapper de batch pra CODEML/PAML (testes de seleção positiva dN/dS). Existe
pra eliminar os erros mais comuns de rodar CODEML na mão: árvore enraizada
onde não pode, amostra faltando num locus específico, `.ctl` mal formatado,
LRT com grau de liberdade errado.

## Como rodar (sem GUI)

```bash
python3 easypaml_cli.py --input PASTA --tree ARVORE.nwk --output SAIDA \
    --models M1a,M2a,M7,M8 --workers N [--skip-beb] [--timeout 1600]
```

Ou via `--config config.json` com as mesmas chaves em formato JSON
(`input`, `tree`, `output`, `models`, `workers`, `timeout`, `run_lrt`,
`skip_beb`, `auto_prune_tree`). O comando efetivamente usado sempre é
gravado em `SAIDA/run_config.json` -- se precisar reproduzir ou entender uma
run já feita, comece lendo esse arquivo, não adivinhe pelos flags.

**`--input`**: pasta com um `.fasta`/`.phy` por gene (todas as amostras
daquele gene num arquivo só, alinhadas). **`--tree`**: Newick, enraizada ou
não, com ou sem táxons a mais/a menos que o FASTA -- poda e desenraizamento
pra modelos de sítio são automáticos (`auto_prune_tree`, default ligado).

Pra chamar direto em Python sem o CLI:

```python
from src.backend.codeml_backend import CodemlBatchAnalysis
app = CodemlBatchAnalysis()
app.config = {
    'input_folder': Path(...), 'tree_file': Path(...), 'output_folder': Path(...),
    'models': ['M1a', 'M2a', 'M7', 'M8'], 'n_workers': 8, 'timeout': 1600,
    'run_lrt': True, 'auto_prune_tree': True, 'skip_beb': False,
}
app.run_batch_analysis()
```

## Onde estão os resultados

```
SAIDA/
  run_config.json          -- config efetivamente usada (reprodutibilidade)
  batch_analysis_log.txt   -- log linha a linha, grep por "FINISHED"/"WARN"/"ERROR"
  analysis_summary.tsv     -- uma linha por gene, ja parseavel com pandas
  M1a/GENE_M1a_results.txt -- outfile cru do CODEML, um por gene x modelo
  M2a/ ...
```

`analysis_summary.tsv` já tem tudo que normalmente se quer sem abrir os
outfiles crus: `lnL`, `np`, `ntime`, `omega`, tempo de execução, e colunas
`lrt_*` já calculadas (estatística + p-valor) quando os pares de modelo
nulo/alternativo (M1a/M2a, M7/M8, etc.) foram ambos rodados. Leia com
`pandas.read_csv(path, sep='\t')` -- não precisa de parser novo, não existe
formato JSON separado porque essa tabela já cobre o caso de uso.

Pra ir além do resumo (sítios individuais sob seleção positiva, valores por
ramo), use as classes já prontas em vez de regex no outfile cru:

- `src.backend.sites_parser.SitesParser` -- sítios BEB/NEB, omega robusto,
  filtragem por p-valor (`extract_omega_robust`, `parse_sites_from_file`,
  `filter_sites_by_pvalue`).
- `src.backend.branch_extractor.BranchExtractor` -- omega por ramo, árvore
  anotada em JSON (modelos Branch/Branch-site).

## Coisas que já foram resolvidas, não precisa reimplementar

- **Árvore desenraizada automaticamente** pra M0/M1a/M2a/M7/M8 (CODEML exige
  isso sem clock/marcação de ramo). Ver `_run_single_analysis`.
- **Poda automática por locus**: amostra ausente num gene específico não
  derruba o locus inteiro, só sai da árvore-guia daquele run.
- **`--skip-beb`**: BEB (classificação de sítio por sítio em M2a/M8) é a
  etapa mais cara, pode levar minutos por locus. Desligado por padrão
  (mais confiável). Ligando, o LRT continua valendo -- lnL/np são escritos
  antes do BEB começar -- só a tabela de sítio some (cai pra NEB, que fica
  disponível quando não foi cortado a tempo).
- **Warm-start implícito via M0**: se você pede só modelos de sítio (sem
  M0 na lista), o backend roda M0 escondido uma vez por gene só pra extrair
  κ e branch lengths iniciais, e reusa isso nos modelos pedidos --
  resultado final idêntico (branch lengths são reestimados livremente),
  só converge mais rápido. Não aparece no `analysis_summary.tsv`.

## Antes de "otimizar" ou "consertar" algo aqui

O binário `bin/codeml` (não versionado, cada instalação copia o seu) é a
fonte de qualquer diferença de tempo/resultado -- meça antes de assumir que
o wrapper Python é o gargalo. `benchmark_vs_raw.py` compara CODEML puro
contra o CLI no mesmo locus/modelo pra isolar overhead real.
