# EasyPAML — Métodos

Este documento descreve exatamente o que o EasyPAML faz, para quem precisa
descrever a análise num artigo ou revisar um manuscrito que usou o programa.
Tudo o que está aqui é registrado, por execução, em `run_config.json`, no
`.ctl` de cada gene × modelo e no `LRT_results.txt`.

Versão descrita: **0.3.0** (ver [CHANGELOG.md](CHANGELOG.md) para diferenças em
relação a versões anteriores).

---

## 1. Entrada

- **Alinhamentos**: um arquivo por gene numa pasta (`.fasta`, `.fas`, `.fa`, `.fna`,
  `.phy`, `.phylip`). PHYLIP é aceito em formato estrito (nomes de 10 colunas) ou
  relaxado (nome separado da sequência por espaço), sequencial ou intercalado.
  Se houver dois arquivos do mesmo gene (`gene.fasta` e `gene.phy`), usa-se o FASTA
  e o outro é ignorado, com aviso.
- **Árvore**: Newick, com ou sem raiz, com ou sem comprimentos de ramo e com ou
  sem cabeçalho `N 1`. Uma árvore por gene é aceita: `GENE.nwk` (`.tree`, `.tre`,
  `.newick`, `.treefile`) na pasta dos alinhamentos ou em `--tree-folder`, pareada
  pelo nome do arquivo; para esse gene ela substitui a árvore geral. Para o modelo Branch/Branch-site, as marcas `#1`, `#2` … são
  feitas na janela "Marcar ramos".
- **Codeml**: o executável é procurado em `--codeml`/configuração, variável
  `EASYPAML_CODEML`, `bin/codeml(.exe)` do projeto e, por fim, no `PATH`. A versão
  é lida do próprio codeml e gravada em `run_config.json`.

## 2. Verificação antes de rodar

Para cada gene, antes de qualquer execução:

| Verificação | Consequência |
|---|---|
| menos de 3 sequências, sequências de comprimentos diferentes, comprimento não múltiplo de 3, nome repetido | o gene não roda (FALHOU, com o motivo) |
| stop codon interno (fora do último códon) | sem "Ignorar stop codons": o gene não roda e o motivo traz a sequência e a posição do códon. Com a opção: o codeml roda e trata a coluna inteira como dado ausente (comportamento do próprio codeml) |
| sequência do alinhamento ausente da árvore | com poda automática (padrão): a sequência é **excluída** da análise do gene (aviso com o nome mais parecido da árvore); sem poda: o codeml falha |
| táxon da árvore ausente do alinhamento | com poda automática: removido da árvore daquele gene |

Na janela, esses itens aparecem num diálogo antes de começar; no CLI, são
impressos antes de começar (`--strict` aborta se houver problemas).

## 3. O que vai para o codeml

O codeml recebe, para cada gene × modelo, uma cópia do alinhamento em FASTA
(só as sequências presentes na árvore) e a árvore podada; para os modelos de
sítio (M0, M1a, M2a, M7, M8, M8a), a árvore é **desenraizada** (raiz bifurcada →
tricotomia), como o PAML exige.

**Todos** os parâmetros relevantes são escritos explicitamente no `.ctl` (nada
fica no default interno do codeml — por exemplo, sem `ncatG` no `.ctl` a beta do
M7/M8 é discretizada em 4 categorias quando o modelo roda sozinho, mas em 10
quando vários modelos rodam num mesmo `.ctl` (`NSsites = 7 8`); conferido no
PAML 4.9j e 4.10.10):

| Parâmetro | Valor padrão | Observação |
|---|---|---|
| `seqtype` | 1 | códons |
| `CodonFreq` | **2 (F3x4)** | 0 Fequal, 1 F1x4, 2 F3x4, 3 F61, 4 F1x4MG, 5 F3x4MG, 6 FMutSel0, 7 FMutSel. Configurável (janela: Configurações / "editar"; CLI: `--codonfreq`) |
| `estFreq` | 0 | frequências observadas |
| `icode` | 0 | código genético universal |
| `cleandata` | 1 | remove colunas com gap, ambiguidade ou stop em qualquer sequência (CLI: `--cleandata`) |
| `kappa` / `fix_kappa` | 2 / 0 | κ inicial 2, estimado |
| `omega` / `fix_omega` | 0,5 / 0 | ω inicial 0,5, estimado (M8a e Branch-site nulo: 1 / 1) |
| `ncatG` | **10** | categorias da beta (M7, M8, M8a) (CLI: `--ncatg`) |
| `fix_alpha` / `alpha` / `Malpha` | 1 / 0 / 0 | sem variação gama |
| `clock` | 0 | sem relógio |
| `fix_blength` | 0 | comprimentos de ramo estimados do zero (1 no warm-start opcional, ver §6) |
| `method` | 0 | otimização simultânea |
| `getSE`, `RateAncestor` | 0, 0 | — |
| `Small_Diff` | 0,5×10⁻⁶ | critério de convergência |

Modelos (`model`, `NSsites`):

| Modelo | model | NSsites | Observação |
|---|---|---|---|
| M0 | 0 | 0 | |
| M1a | 0 | 1 | ω₁ = 1 imposto pelo codeml |
| M2a | 0 | 2 | |
| M7 | 0 | 7 | beta(p, q) |
| M8 | 0 | 8 | beta + classe extra com ω livre |
| M8a | 0 | 8 | `fix_omega = 1`, `omega = 1` (Swanson et al. 2003) |
| Branch | 2 | 0 | árvore com marcas; comprimentos de ramo da árvore marcada removidos |
| Branch-site / nulo | 2 | 2 | nulo com `fix_omega = 1`, `omega = 1` |

O `.ctl`, o alinhamento e a árvore exatamente como o codeml leu ficam em
`SAIDA/MODELO/` com caminhos relativos: `cd SAIDA/M8 && codeml GENE_M8.ctl`
reproduz a execução.

## 4. Execução

- O codeml roda com a entrada padrão **fechada**: quando ele pede "Press Enter"
  (por exemplo, ao encontrar um stop codon), segue na hora em vez de esperar.
- Cada execução tem tempo limite (`--timeout`, padrão 1600 s no CLI e 6 h na
  janela) e é encerrada se ficar `--idle-timeout` segundos sem usar CPU
  (padrão 300 s). Nos dois casos o gene aparece como FALHOU, com o motivo e a
  última linha que o codeml escreveu.
- Um gene **falha** se algum modelo pedido termina com código ≠ 0, sem arquivo de
  saída ou sem lnL. A mensagem final é "X de N genes concluídos, Y falharam";
  `genes_status.tsv` lista o motivo por gene.
- "Parar" encerra e recolhe todos os processos codeml em andamento.

## 5. Teste da razão de verossimilhança (LRT)

- Estatística: 2Δℓ = 2(lnL_alternativo − lnL_nulo). Valores negativos (o
  alternativo não alcançou o nulo na otimização) são truncados em 0 (p = 1) e
  sinalizados no `LRT_results.txt`.
- Graus de liberdade (diferença no número de parâmetros livres; os comprimentos
  de ramo entram igualmente nos dois modelos):

  | Teste | df | Distribuição nula usada |
  |---|---|---|
  | M1a vs M2a | 2 | χ²₂ |
  | M7 vs M8 | 2 | χ²₂ |
  | M8a vs M8 | 1 | χ²₁ (*) |
  | M0 vs M1a | 2 | χ²₂ |
  | M0 vs Branch | nº de grupos foreground | χ²_df |
  | Branch-site nulo vs Branch-site | 1 | χ²₁ (*) |

  (*) Nos testes em que o nulo fica na fronteira do espaço de parâmetros (ω = 1
  fixo), a distribuição assintótica é a mistura 50:50 de χ²₀ e χ²₁ (Self & Liang
  1987). Seguindo a recomendação do manual do PAML para o branch-site, o EasyPAML
  usa χ²₁ puro (mais conservador) para p e q; o p da mistura aparece no
  `LRT_results.txt` só como referência.
- p-valores calculados com a função de sobrevivência (`scipy.stats.chi2.sf`), sem
  arredondar para zero (p ≈ 10⁻²² é mostrado como tal).
- **Correção para múltiplos testes**: Benjamini-Hochberg (FDR,
  `scipy.stats.false_discovery_control`) **dentro de cada par de modelos**; a
  família é o conjunto de genes testados naquele par **naquela execução**. Rodar
  os genes em execuções separadas muda as famílias — informe no artigo como os
  genes foram agrupados. A coluna `q_*` do `analysis_summary.tsv` é o corte
  recomendado quando muitos genes são testados.
- M7 vs M8 pode ser significativo apenas porque parte dos sítios é neutra (ω = 1),
  que a beta do M7 não acomoda; M8a vs M8 não tem esse problema. O painel avisa
  quando M7 vs M8 é significativo e M8a vs M8 não.

## 6. Sítios sob seleção (BEB/NEB)

- Os sítios vêm da seção Bayes Empirical Bayes (BEB; Yang, Wong & Nielsen 2005)
  do M2a/M8 (NEB se o BEB não existir). Pr(ω>1) ≥ 0,95 é marcado com `*` e ≥ 0,99
  com `**`, como no codeml.
- **Numeração**: com `cleandata = 1` o codeml remove colunas e numera os sítios nas
  colunas que sobram. O EasyPAML calcula quais códons foram mantidos (mesma regra
  do codeml: coluna removida se houver, em qualquer sequência, algo diferente de
  A/C/G/T ou um stop codon), confere com o número de sítios que o codeml declara
  (`ls` no arquivo de saída) e grava `GENE_MODELO_sitemap.json`. A interface e as
  exportações mostram a posição no alinhamento do usuário e a do codeml. Se a
  conferência falhar, a numeração original não é mostrada e há aviso.
- O "ω e p₁ da classe positiva" são a última classe da tabela "MLEs of dN/dS (w)
  for site classes" do M2a/M8. O ω médio do gene (Σ pᵢωᵢ) não é usado como critério.

## 7. Opções para muitos genes

- `--skip-beb`: interrompe M2a/M8 quando o codeml começa o BEB. lnL/np/ω (usados
  no LRT) já foram escritos antes do BEB; a tabela de sítios cai para NEB, e
  seções incompletas são cortadas do arquivo.
- `--two-pass`: passada 1 em todos os genes com `--skip-beb`; LRT + BH; passada 2
  com BEB completo só nos genes com q < `--sig-threshold` (padrão 0,05). Saídas em
  `SAIDA/pass1_screen/` e `SAIDA/pass2_beb/`.
- `--warm-start-m0` (desligado por padrão): roda um M0 por gene e usa κ e os
  comprimentos de ramo dele como **valores iniciais** dos modelos de sítio
  (`fix_blength = 1`, "initial" no pamlDOC), com multi-start de ω (0,2 / 1,0 / 2,5)
  e escolha do maior lnL. Não é garantia de resultado idêntico à estimativa do zero.
  **Atenção**: até a versão 0.2.0 esta opção usava `fix_blength = 2`, que no PAML
  significa comprimentos **fixos** nos valores do M0 (não apenas ponto de partida);
  resultados obtidos com `--warm-start-m0` nessas versões devem ser refeitos.

## 8. Saídas

| Arquivo | Conteúdo |
|---|---|
| `run_config.json` | versões do EasyPAML, do codeml e do Python; parâmetros do `.ctl` por modelo; testes LRT; opções |
| `analysis_summary.tsv` | por gene: status; lnL, np, ntime, ω, tempo por modelo; ω e p₁ da classe positiva (M2a/M8); 2Δℓ, p e q por teste |
| `LRT_results.txt` | LRT por gene e por teste, com esta nota metodológica |
| `genes_status.tsv` | ok / failed + motivo |
| `batch_analysis_log.txt` | log completo, inclusive a saída do codeml |
| `MODELO/GENE_MODELO.*` | `.ctl`, alinhamento, árvore, saída bruta (`_results.txt`), `rst`, `sitemap.json` |

## 9. Texto sugerido para Métodos (adapte)

> Positive selection was tested with EasyPAML v0.3.0 (https://github.com/Hesatum/EasyPAML)
> running codeml from PAML vX.Y (Yang 2007). Codon site models M7, M8 and M8a were fitted
> with codon frequencies F3x4 (CodonFreq = 2), the beta distribution discretised into 10
> categories (ncatG = 10), κ and ω estimated (initial values 2 and 0.5), and alignment
> columns with gaps, ambiguities or stop codons removed (cleandata = 1). Input trees were
> pruned to the taxa present in each alignment and unrooted. Models were compared with
> likelihood ratio tests (M8 vs M7, df = 2; M8 vs M8a, df = 1, χ²₁), and p-values were
> corrected for multiple testing with the Benjamini-Hochberg procedure across the N genes
> tested. Sites under positive selection were identified with the Bayes Empirical Bayes
> procedure (posterior probability ≥ 0.95).

## Referências

- Yang Z (2007) PAML 4. *Mol Biol Evol* 24:1586–1591.
- Yang Z, Nielsen R, Goldman N, Pedersen A-MK (2000) Codon-substitution models for
  heterogeneous selection pressure at amino acid sites. *Genetics* 155:431–449. (M7, M8)
- Swanson WJ, Nielsen R, Yang Q (2003) Pervasive adaptive evolution in mammalian
  fertilization proteins. *Mol Biol Evol* 20:18–20. (M8a)
- Wong WSW, Yang Z, Goldman N, Nielsen R (2004) Accuracy and power of statistical
  methods for detecting adaptive evolution in protein coding sequences. *Genetics*
  168:1041–1051. (M1a, M2a, M8a)
- Yang Z, Wong WSW, Nielsen R (2005) Bayes empirical Bayes inference of amino acid
  sites under positive selection. *Mol Biol Evol* 22:1107–1118. (BEB)
- Zhang J, Nielsen R, Yang Z (2005) Evaluation of an improved branch-site likelihood
  method for detecting positive selection at the molecular level. *Mol Biol Evol*
  22:2472–2479. (branch-site)
- Self SG, Liang K-Y (1987) *J Am Stat Assoc* 82:605–610.
- Benjamini Y, Hochberg Y (1995) *J R Stat Soc B* 57:289–300.
