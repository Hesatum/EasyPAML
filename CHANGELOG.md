# Changelog

Formato: [Keep a Changelog](https://keepachangelog.com/pt-BR/1.1.0/). Versões
seguem a tag do repositório.

## [Não lançada] — 0.3.0.dev0

Correções a partir do teste de usabilidade de 05/10/2026.

### Mudanças que afetam resultados

- **CodonFreq padrão passou de 7 (FMutSel) para 2 (F3x4).** A interface rotulava
  o 7 como "F3×4 (recomendado)"; quem usou o padrão em versões anteriores rodou
  FMutSel. Os lnL mudam e a lista de sítios BEB pode mudar.
- **Todos os parâmetros do `.ctl` agora são escritos explicitamente**, incluindo
  `ncatG = 10`, `kappa = 2`, `fix_kappa = 0`, `fix_blength = 0`, `method = 0`.
  Antes, `ncatG` ficava no default do codeml: 4 categorias quando cada modelo
  roda sozinho (como o EasyPAML faz), contra as 10 da prática comum e do que o
  codeml usa ao rodar `NSsites = 7 8` num mesmo `.ctl`.
- **Warm-start M0 (`--warm-start-m0`, desligado por padrão)**: usava
  `fix_blength = 2`, que no PAML **fixa** os comprimentos de ramo nos valores do M0;
  agora usa `fix_blength = 1` (valores iniciais), como a documentação dizia.
- p-valores calculados com `chi2.sf` (antes `1 − chi2.cdf`, que dava p = 0 para
  estatísticas grandes) e com a mesma regra no arquivo, no painel e nas exportações
  (antes o painel usava a mistura 50:50 no branch-site e o arquivo usava χ²₁).
- Numeração dos sítios: com `cleandata = 1` o codeml numera só as colunas que sobram;
  agora a posição no alinhamento original também é mostrada/exportada.

### Adicionado

- Modelo **M8a** (M8 com ω = 1 fixo) e teste **M8a vs M8** (df = 1, χ²₁; mistura
  50:50 reportada como referência). Ao escolher M8 com modelos nulos automáticos,
  entram M7 e M8a.
- Verificação dos dados antes de rodar (stop codons com posição, nomes fora da árvore
  com sugestão, táxons podados, duplicados, comprimento % 3), em diálogo na janela e
  impressa no CLI (`--strict`).
- Falhas visíveis: gene com código de saída ≠ 0, sem saída, sem lnL, tempo esgotado
  ou codeml sem usar CPU aparece como FALHOU com o motivo; "X de N genes concluídos,
  Y falharam"; `genes_status.tsv`.
- Reprodutibilidade: `.ctl` com caminhos relativos, alinhamento e árvore usados
  copiados ao lado; `run_config.json` também pela janela, com parâmetros do `.ctl`
  por modelo e versões do EasyPAML e do codeml; `--version`; menu Sobre;
  `requirements-lock.txt`.
- Painel de resultados: aba Resumo (uma frase por gene e teste), lnL dos dois
  modelos, q-valor, ω e p₁ da classe positiva, copiar/exportar sítios em TSV, botão
  "Abrir pasta de resultados"; resultados abrem ao terminar.
- CLI: `--codonfreq`, `--ncatg`, `--kappa`, `--omega`, `--cleandata`,
  `--ignore-stop-codons`, `--codeml`, `--idle-timeout`, `--strict`, `--lang`,
  `--verbose`; código de saída 1 se algum gene falhar.
- PHYLIP relaxado (nomes com mais de 10 caracteres), sequencial ou intercalado.
- Testes automatizados (`pytest`) e dados simulados com gabarito (`tests/data`).
- `METODOS.md` (conteúdo técnico antes só no `AGENTS.md`) e este changelog.

### Corrigido

- O codeml podia ficar parado para sempre esperando "Enter" (stop codon): agora roda
  com a entrada padrão fechada, com tempo limite e detecção de inatividade; "Parar"
  encerra e recolhe os processos.
- "Ignorar stop codons" não tinha efeito; agora decide se genes com stop interno
  rodam (coluna tratada como ausente pelo codeml) ou falham com a posição do stop.
- "ANÁLISE CONCLUÍDA" aparecia mesmo com todos os genes falhando.
- Sequências ausentes da árvore eram excluídas sem aviso na janela.
- `.fasta` + `.phy` do mesmo gene contavam como dois genes.
- Nota metodológica do LRT dizia que os modelos NSsites têm ntime = 0 (falso).
- Referência do M7/M8 corrigida para Yang et al. (2000).
- Dados de exemplo: `exemplos_teste/arvore_amostras.nwk` corresponde a `amostras/`
  (a árvore indicada antes, `final-tree.txt`, não tinha nenhum nome em comum).
- Instalação Linux: `install.sh` cria `.venv` (sem `--user`; funciona no Ubuntu
  24.04/PEP 668), checa venv/tkinter com o comando a rodar, cria o lançador mesmo se
  o codeml falhar; `EasyPAML.sh` versionado; README com URL correta e pré-requisitos.
- Windows: `install.bat` cria `.venv` (com volta para `--user`); `EasyPAML.bat` usa
  o `.venv`; `.bat` com CRLF.
- Interface: janela de resultados cabe na tela e fecha com Esc; botões sem texto em
  hover/pressionado; contraste ≥ 4,5:1; roda do mouse não altera mais as CPUs;
  pasta de saída criada pelo diálogo; contagem de alinhamentos; barra de progresso;
  log com uma mensagem por linha e detalhes técnicos escondidos; tradução PT/EN
  completa com acentos; idioma inicial = idioma do sistema (padrão inglês).
- `requirements.txt`: SciPy mínimo 1.11 (exigido pela correção BH) e `psutil`.

## [0.2.0] — 2026-10-05

- Reestrutura em `src/`, modo CLI (`easypaml_cli.py`), correção BH no LRT, aba de
  interpretação GO. (Tag criada sobre o merge da branch `feature/cli-mode-e-skip-beb`.)
