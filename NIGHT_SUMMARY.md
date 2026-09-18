# Sessão noturna 2026-09-17/18 — resumo pra revisão

Branch: `feature/cli-mode-e-skip-beb`. Nada disso tocou o `EasyPAML-main`
que rodou o lote de produção a noite toda (`codeml_easypaml_5725/`) — essa
é uma cópia isolada, com o binário `bin/codeml` também copiado (não
versionado, igual no repo original).

Todo commit desta lista foi testado com dado real do projeto antes de ser
feito — não tem nada aqui que só rodou uma vez na minha cabeça.

## O que decidir

1. **`--warm-start-m0`**: fica ligado por padrão daqui pra frente, ou
   continua opt-in? Números reais (20 loci, M1a): 17,6× mais rápido;
   2/20 genes com lnL levemente pior (pior caso -1,47), 2/20 com lnL
   *melhor*. Não é garantia matemática de resultado idêntico ao
   from-scratch, é estatisticamente raro de divergir depois do
   multi-start. Minha recomendação: vale a pena pra triagem em escala de
   genoma inteiro; pra qualquer gene que vire resultado final de
   publicação, rodar sem warm-start como confirmação antes de citar o
   número.
2. **Aba "Interpretação" (GO)**: lógica testada ponta a ponta com os
   5.725 genes reais do projeto. A parte visual (Tkinter) eu **não
   consegui testar** — sem display gráfico neste ambiente, sem Xvfb, sem
   sudo pra instalar. Abrir e clicar é o próximo passo real de validação.
3. **Recompilar o `codeml` com flags de otimização do compilador**: cogitei
   e não fiz. Trocar o binário científico que gera os números é uma
   categoria de mudança diferente de otimizar o wrapper — não é uma
   decisão pra eu tomar sozinho de madrugada. Fica anotado como ideia, não
   como trabalho pendente.
4. **Worker count**: com 32 cores na máquina, considerar `--workers 28-30`
   pra rodadas futuras (deixa 2-4 livres pro SO). Não testei
   especificamente — é extrapolação, não medição, tratar como tal.

## O que foi feito (ordem cronológica, 1 commit = 1 mudança testada)

| # | O quê | Medido/testado como |
|---|---|---|
| 1 | CLI headless (`easypaml_cli.py`) + `--skip-beb` com truncagem defensiva | Gene real do repo + locus de produção; NEB/BEB nunca aparecem incompletos no outfile |
| 2 | Corte de 313 linhas de código morto/duplicado (audit ponytail) | `run_wgs_analysis`, `interactive_setup`+`main()` sem chamador; `_extract_likelihood/_np/_ntime` viram uma função só |
| 3 | `bin/codeml` tirado do git (não era versionado no original) | — |
| 4 | Warm-start M0 implícito implementado (só documentado antes, não existia) | 3 loci reais × 4 modelos: 2697,7s → 385,97s (7,0×), lnL idêntico |
| 5 | `--two-pass` (screen sem BEB → BEB só nos LRT-significativos) | 6 loci reais: achou os 2 significativos certos, BEB só neles |
| 6 | Aba "Interpretação": candidatos ranqueados + enriquecimento GO | Backend (`go_enrichment.py`) testado com os 5.725 genes reais + saída real de `--two-pass`; GUI não-testável neste ambiente (ver acima) |
| 7 | Warm-start virou opt-in de verdade (`config['warm_start_m0']`, default False) | Confirmado: default agora reproduz exatamente o lnL do from-scratch |
| 8 | Multi-start de ômega (3 pontos: 0.2/1.0/2.5) quando warm-start ativo | 20 loci reais: risco de ótimo local caiu de ~1/3 pra 1/20 grande + 1/20 pequeno, com 2/20 *melhorando* |
| 9 | Auto-limpeza do próprio código de hoje (dead code em `go_enrichment.py`) | Self-check ainda passa depois do corte |
| 10 | Correção Benjamini-Hochberg (FDR) integrada ao `_run_lrt_analysis` (colunas `q_*` no TSV, campo espelhado no `LRT_results.txt`), propagada pro corte do `--two-pass` e do `rank_candidates` (GO) | 1.542-1.573 genes reais já concluídos no lote de produção (leitura, lote nunca tocado): q≥p em 100% dos casos (0 violações); reduz significância de 64,3%→61,5% (M1a/M2a) e 70,3%→68,4% (M7/M8) — confirma que o confundidor dominante é o poder estatístico ligado ao tamanho do alinhamento, não inflação por múltiplos testes |

## Arquivos novos

- `easypaml_cli.py` — entrada headless
- `src/backend/go_enrichment.py` — enriquecimento GO, testável sem GUI (`python3 -m src.backend.go_enrichment` roda o self-check)
- `benchmark_vs_raw.py` — compara CODEML puro vs. wrapper no mesmo locus/modelo
- `AGENTS.md` — guia técnico pra quem (humano ou IA) mexer nisso depois

## Lote de produção (referência, não mexi nele)

Rodando a noite toda com o código de ANTES de todas essas mudanças (não
retroativo). Zero erro, progresso normal — ver `codeml_easypaml_5725/batch_analysis_log.txt`.
