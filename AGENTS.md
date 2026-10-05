# EasyPAML -- guia para agentes de IA

Wrapper de batch pra CODEML/PAML (testes de seleção positiva dN/dS). O que o
programa faz cientificamente (parâmetros do `.ctl`, LRT, df, BH, numeração de
sítios, warm-start) está em **[METODOS.md](METODOS.md)** -- leia lá, não aqui.
Este arquivo só diz onde as coisas estão no código.

## Como rodar (sem GUI)

```bash
python3 easypaml_cli.py --input PASTA --tree ARVORE.nwk --output SAIDA \
    --models M1a,M2a,M7,M8,M8a --workers N [--codonfreq 2] [--ncatg 10] \
    [--ignore-stop-codons] [--skip-beb] [--two-pass] [--warm-start-m0] [--codeml PATH]
```

`--help` lista tudo; `--config config.json` aceita as mesmas chaves. Código de
saída 1 = algum gene falhou. Comece qualquer investigação de uma run lendo
`SAIDA/run_config.json` (versões + parâmetros do `.ctl` por modelo) e
`SAIDA/genes_status.tsv` (ok/failed + motivo).

Em Python:

```python
from src.backend.codeml_backend import CodemlBatchAnalysis
app = CodemlBatchAnalysis()
app.config = {'input_folder': Path(...), 'tree_file': Path(...), 'output_folder': Path(...),
              'models': ['M7', 'M8', 'M8a'], 'n_workers': 8, 'run_lrt': True}
summary = app.run_batch_analysis()   # {'total','ok','failed','failures',...}
```

## Mapa do código

| Arquivo | O que tem |
|---|---|
| `src/backend/codeml_backend.py` | execução em lote (`run_batch_analysis`, `_process_gene`, `_run_single_analysis`), LRT, TSV |
| `src/backend/ctl_params.py` | todos os parâmetros do `.ctl` e a tabela de CodonFreq |
| `src/backend/lrt_stats.py` | pares de modelos, df, p (chi2.sf), BH -- usado pelo backend E pelo painel |
| `src/backend/alignment_io.py` | FASTA/PHYLIP (estrito/relaxado), stop codons, regra do cleandata |
| `src/backend/preflight.py` | verificação antes de rodar (mensagens PT/EN) |
| `src/backend/site_map.py` | numeração codeml -> alinhamento do usuário |
| `src/backend/sites_parser.py` | BEB/NEB, classe positiva (ω, p₁), ω por ramo |
| `src/backend/messages.py` | mensagens do backend em PT/EN |
| `src/gui/main_gui.py` | janela principal |
| `src/gui/results_viewer.py` | painel de resultados |
| `src/gui/gui_texts.py` | todos os textos da interface (PT/EN) |
| `tests/` | pytest; `EASYPAML_TEST_CODEML=/caminho/codeml` liga o teste com codeml real |

## Coisas que já foram resolvidas, não reimplementar

- Árvore desenraizada automaticamente para M0/M1a/M2a/M7/M8/M8a; poda por locus.
- codeml roda com stdin fechado, timeout e detecção de inatividade (nunca trava).
- Cada gene × modelo grava `.ctl` + alinhamento + árvore com caminhos relativos
  (`cd SAIDA/M8 && codeml GENE_M8.ctl` reproduz).
- `--warm-start-m0` usa `fix_blength = 1` (valores iniciais). **Não** use 2: no
  PAML 2 = comprimentos fixos (era o bug até a 0.2.0). Os números de speedup
  medidos antes da 0.3.0 foram obtidos com o comportamento antigo.

## Antes de "otimizar" ou "consertar" algo aqui

O binário do codeml é a fonte de qualquer diferença de tempo/resultado -- meça
antes de assumir que o wrapper Python é o gargalo. `benchmark_vs_raw.py` compara
CODEML puro contra o CLI no mesmo locus/modelo. Rode `pytest` antes de commitar.
