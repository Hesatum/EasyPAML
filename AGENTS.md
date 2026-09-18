# EasyPAML -- guia para agentes de IA

Wrapper de batch pra CODEML/PAML (testes de seleção positiva dN/dS). Existe
pra eliminar os erros mais comuns de rodar CODEML na mão: árvore enraizada
onde não pode, amostra faltando num locus específico, `.ctl` mal formatado,
LRT com grau de liberdade errado.

## Como rodar (sem GUI)

```bash
python3 easypaml_cli.py --input PASTA --tree ARVORE.nwk --output SAIDA \
    --models M1a,M2a,M7,M8 --workers N [--skip-beb] [--warm-start-m0] \
    [--two-pass] [--timeout 1600]
```

Ou via `--config config.json` com as mesmas chaves em formato JSON
(`input`, `tree`, `output`, `models`, `workers`, `timeout`, `run_lrt`,
`skip_beb`, `auto_prune_tree`, `warm_start_m0`, `two_pass`,
`sig_threshold`). O comando efetivamente usado sempre é gravado em
`SAIDA/run_config.json` -- se precisar reproduzir ou entender uma run já
feita, comece lendo esse arquivo, não adivinhe pelos flags.

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
- `src.backend.go_enrichment.rank_candidates(summary_tsv, go_annotation_tsv)`
  -- genes com LRT significativo ranqueados por p-valor + termos GO
  enriquecidos (Fisher exato, candidatos vs. todos os genes testados).
  Pensado pra escala de genoma inteiro: reduz "leia N tabelas BEB" pra
  "leia uma lista curta". TSV de anotação esperado: colunas
  `gene_id_full`, `go_biological_process`, `go_cellular_component`,
  `go_molecular_function`, cada uma com `"descrição [GO:XXXXXXX]; ..."`.
  Também acessível na GUI, aba "Interpretação"/"Interpretation" do
  results viewer.

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
- **`--two-pass`**: roda tudo sem BEB primeiro (rápido), calcula LRT, e só
  reroda com BEB completo os genes com p < `--sig-threshold` (default
  0.05). Saída em `SAIDA/pass1_screen/` (todos os genes) e
  `SAIDA/pass2_beb/` (só os significativos). Bom pra escala de genoma
  inteiro, onde a maioria dos genes não rejeita o nulo.
- **`--warm-start-m0`** (opt-in, desligado por padrão): roda M0 escondido
  uma vez por gene pra usar κ/branch-lengths como ponto de partida nos
  modelos de sítio pedidos, com multi-start automático de ômega (3 pontos:
  0.2/1.0/2.5) pra reduzir risco de ótimo local. **Não é garantia
  matemática de resultado idêntico ao from-scratch** -- medido em 20 loci
  reais: 17.6x mais rápido, 2/20 genes com lnL levemente pior (pior caso:
  -1.47), 2/20 com lnL *melhor* (o multistart escapou de ótimo que o
  from-scratch não escapou). Ligue sabendo do trade-off; pra publicação,
  considere rodar uma amostra com e sem pra confirmar que não muda
  conclusões antes de aplicar no dataset inteiro.

## Antes de "otimizar" ou "consertar" algo aqui

O binário `bin/codeml` (não versionado, cada instalação copia o seu) é a
fonte de qualquer diferença de tempo/resultado -- meça antes de assumir que
o wrapper Python é o gargalo. `benchmark_vs_raw.py` compara CODEML puro
contra o CLI no mesmo locus/modelo pra isolar overhead real.
