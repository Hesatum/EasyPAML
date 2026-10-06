# EasyPAML

Interface gráfica (e modo linha de comando) para análise de seleção positiva
com **PAML/codeml**: modelos de sítio (M0, M1a, M2a, M7, M8, M8a), modelo de
ramos e branch-site, com teste da razão de verossimilhança (LRT), correção de
Benjamini-Hochberg e tabela de sítios (BEB).

*English summary at the end.*

- Métodos (parâmetros do codeml, LRT, correção BH): [METODOS.md](METODOS.md)
- O que mudou em cada versão: [CHANGELOG.md](CHANGELOG.md)

---

## Instalação

### Baixar

- Com git: `git clone https://github.com/Hesatum/EasyPAML.git`
- Sem git: [baixe o ZIP](https://github.com/Hesatum/EasyPAML/archive/refs/heads/main.zip)
  e extraia, ou pegue uma versão numerada em
  [Releases / Tags](https://github.com/Hesatum/EasyPAML/tags).

### Windows

> Pré-requisito: **Python 3.8+** ([python.org/downloads](https://www.python.org/downloads/)) —
> marque **"Add Python to PATH"** durante a instalação. O `codeml.exe` já vem na pasta `bin/`.

1. Dentro da pasta do EasyPAML, dê **duplo-clique em `install.bat`**
   (cria um ambiente `.venv` com as dependências; 1–3 minutos).
2. Abra pelo atalho **EasyPAML** da Área de Trabalho ou pelo **`EasyPAML.bat`**.

### Linux (Ubuntu/Debian) e macOS

1. Pré-requisitos do sistema (uma vez só):

   ```bash
   sudo apt update && sudo apt install -y git python3-pip python3-venv python3-tk paml
   ```

   (Fedora: `sudo dnf install git python3-pip python3-tkinter paml` · macOS: `brew install python-tk brewsci/bio/paml`)

2. Baixar e instalar:

   ```bash
   git clone https://github.com/Hesatum/EasyPAML.git
   cd EasyPAML
   ./install.sh
   ```

   O `install.sh` cria o ambiente isolado `.venv/` **dentro da pasta do EasyPAML**
   e instala as dependências nele — não usa `pip --user` nem mexe no Python do
   sistema (funciona no Ubuntu 23.04+/24.04, que bloqueia `pip` fora de ambientes
   virtuais). Se faltar algo, ele para e mostra o comando exato a rodar. Pode ser
   rodado de novo quantas vezes quiser.

3. Abrir:

   ```bash
   ./EasyPAML.sh
   ```

   O `EasyPAML.sh` só ativa o `.venv` e abre o programa. Plano B, se ele não
   funcionar: `.venv/bin/python EasyPAML.py` (sempre `python3`, nunca `python`,
   no Ubuntu).

### CODEML (PAML)

- **Windows**: `bin/codeml.exe` (PAML 4.9j) já vem no repositório.
- **Linux**: vem do pacote do sistema, `sudo apt install paml` (Ubuntu 24.04:
  PAML 4.9j). O `install.sh` tenta instalar se faltar; se não houver pacote,
  baixa o binário oficial do PAML 4.10.10 para `bin/codeml`.
- **Versões testadas**: PAML 4.9j e 4.10.10. Todos os parâmetros vão escritos no
  `.ctl` (nada fica no default interno do codeml); ver METODOS.md.
- Outro codeml: defina a variável `EASYPAML_CODEML=/caminho/do/codeml` ou use
  `--codeml` no modo linha de comando.
- A versão do codeml usada aparece em **Sobre** e em `run_config.json`.

---

## Como usar (janela)

1. **Pasta de alinhamentos**: uma pasta com um arquivo por gene (`.fasta`, `.fas`,
   `.phy`, `.phylip`), sequências de códons alinhadas. O programa mostra quantos
   alinhamentos encontrou.
2. **Arquivo de árvore**: árvore Newick (`.nwk`, `.tree`, `.tre`, `.txt`), com ou
   sem raiz; pode ter táxons a mais (são podados em cada gene).
   **Uma árvore por gene**: ponha `GENE.nwk` (ou `.tree`/`.tre`) ao lado de
   `GENE.fasta` na pasta de alinhamentos; para esse gene ela substitui o arquivo de
   árvore (se todos os genes tiverem a sua, o arquivo de árvore é dispensado).
3. **Pasta de resultados**: escolha ou digite o nome de uma pasta nova (ela é criada).
4. Ligue os modelos. Com **Modelos nulos automáticos** ligado, o nulo de cada teste
   entra sozinho (M8 → M7 e M8a; M2a → M1a).
5. Clique em **Iniciar**. Antes de rodar, o EasyPAML verifica os dados e mostra um
   diálogo com o que encontrou: stop codons (sequência e posição do códon), nomes do
   alinhamento que não estão na árvore (com o nome mais parecido), táxons que serão
   podados, arquivos duplicados do mesmo gene e comprimento que não é múltiplo de 3.
   Escolha **Corrigir e voltar** ou **Continuar mesmo assim**.
6. A barra mostra "Gene X de N". Ao terminar, aparece "X de N genes concluídos, Y
   falharam" (genes que falharam ficam em vermelho, com o motivo) e o painel de
   resultados abre sozinho.

**Painel de resultados**

- **Resumo**: uma frase por gene e por teste, por exemplo
  "M8 vs M7: significativo (p = 4,7×10⁻²², q = 4,7×10⁻²²) — 20 sítios com Pr(ω>1) ≥ 0,95 ·
  classe positiva: ω = 3,8, p₁ = 0,11".
- **LRT e p-valores**: lnL de cada modelo, 2Δℓ, p (notação científica), q (BH),
  ω e proporção da classe positiva.
- **Sítios sob seleção**: posição no **seu** alinhamento e posição no arquivo do codeml
  (diferentes quando colunas com gap/stop são removidas), aminoácido, Pr(ω>1),
  `*` (≥ 0,95) / `**` (≥ 0,99), ω médio ± EP; botões para copiar/exportar em TSV.
- **Abrir pasta de resultados** e exportação para Excel, CSV, PNG e HTML.

> O ω médio do gene **não** é critério de seleção positiva: ele fica abaixo de 1
> mesmo quando poucos sítios estão sob seleção forte. Use o LRT (q) e a tabela de sítios.

---

## Modelos

| Modelo | Para que serve |
|---|---|
| **M0** | um ω para o gene inteiro (linha de base; nulo do Branch) |
| **M1a / M2a** | M2a vs M1a: seleção positiva por sítio (df = 2) |
| **M7 / M8** | M8 vs M7: beta + classe extra (df = 2) |
| **M8a** | M8 com a classe extra fixa em ω = 1. **M8 vs M8a** (df = 1) não é enganado por sítios neutros, ao contrário do M8 vs M7 |
| **Branch** | ω por grupo de ramos marcado (vs M0) |
| **Branch-site** | seleção episódica em sítios do ramo foreground (#1) |

Recomendação usual para seleção positiva por sítio: rode **M8** (o M7 e o M8a
entram sozinhos) e confira os dois testes; se só o M8 vs M7 for significativo,
o sinal pode vir de sítios neutros.

---

## Modo linha de comando (servidor / muitos genes)

```bash
.venv/bin/python easypaml_cli.py --input PASTA --tree ARVORE.nwk --output SAIDA \
    --models M1a,M2a,M7,M8,M8a --workers 8
```

Opções principais (`--help` mostra todas):

| Opção | Padrão | Significado |
|---|---|---|
| `--models` | `M1a,M2a,M7,M8,M8a` | modelos a rodar |
| `--no-m8a` | — | não roda o M8a (sem o teste M8a vs M8); na janela: Configurações › "Incluir M8a" |
| `--tree-folder` | — | pasta com uma árvore por gene (`GENE.nwk`), pareada pelo nome; também vale `GENE.nwk` na própria pasta de alinhamentos |
| `--codonfreq` | `2` (F3x4) | CodonFreq do codeml (0 Fequal, 1 F1x4, 2 F3x4, 3 F61, 7 FMutSel …) |
| `--ncatg` | `10` | categorias da beta (M7/M8/M8a) |
| `--cleandata` | `1` | remove colunas com gap/ambiguidade/stop |
| `--ignore-stop-codons` | desligado | sem ela, genes com stop codon interno falham com a posição do stop |
| `--workers` | `4` | genes em paralelo |
| `--timeout` / `--idle-timeout` | automático / 300 s | limite por execução ([como é calculado](docs/benchmark_tempos.md)) / codeml sem usar CPU |
| `--skip-beb`, `--two-pass` | — | para milhares de genes (ver METODOS.md) |
| `--codeml` | — | caminho do codeml |
| `--strict` | — | não roda nada se a verificação inicial achar problemas |
| `--lang pt\|en` | idioma do sistema | idioma das mensagens |
| `--config arquivo.json` | — | as mesmas opções em JSON |
| `--version` | — | versão do EasyPAML |

Código de saída: 0 se todos os genes rodaram, 1 se algum falhou, 2 com `--strict`
e problemas nos dados.

---

## O que fica na pasta de resultados

```
SAIDA/
  run_config.json          versões (EasyPAML, codeml, Python), parâmetros do .ctl por modelo, opções
  batch_analysis_log.txt   log completo (inclui a saída do codeml)
  genes_status.tsv         um gene por linha: ok / failed + motivo
  analysis_summary.tsv     lnL, np, ω por modelo; ω e p₁ da classe positiva; 2Δl, p e q por teste
  LRT_results.txt          LRT por gene, com a nota metodológica
  M8/GENE_M8.ctl           .ctl usado (todos os parâmetros, caminhos relativos)
  M8/GENE_M8_seq.fasta     alinhamento exatamente como o codeml leu
  M8/GENE_M8_tree.nwk      árvore exatamente como o codeml leu (podada/desenraizada)
  M8/GENE_M8_results.txt   saída bruta do codeml
  M8/GENE_M8_sitemap.json  numeração dos sítios: codeml → alinhamento original
```

Para refazer uma execução à mão: `cd SAIDA/M8 && codeml GENE_M8.ctl`.

---

## Dados de exemplo

`exemplos_teste/` tem 25 genes de *Cereus* (cactos), com 8 a 21 sequências cada, e
uma árvore de 21 táxons com exatamente os mesmos nomes dos alinhamentos:

1. **Pasta de alinhamentos**: `exemplos_teste/amostras/`
2. **Arquivo de árvore**: `exemplos_teste/arvore_amostras.nwk`
3. **Pasta de resultados**: uma pasta nova
4. Ligue **M8** e clique em **Iniciar**

Os resultados do exemplo não vêm no repositório; para gerá-los (cerca de 1 h com 12 processos):

```bash
.venv/bin/python easypaml_cli.py --input exemplos_teste/amostras \
    --tree exemplos_teste/arvore_amostras.nwk --output exemplos_teste/resultados --workers 8
```

`tests/data` tem um conjunto simulado com resposta conhecida.

---

## Solução de problemas

**`./install.sh` diz que falta o venv ou o tkinter** → rode o comando que ele mostra
(`sudo apt install python3-venv python3-tk`) e rode `./install.sh` de novo.

**`./EasyPAML.sh: No such file or directory`** → você está fora da pasta do EasyPAML
(`cd EasyPAML`) ou baixou uma versão antiga. Plano B: `.venv/bin/python EasyPAML.py`.

**`ModuleNotFoundError: No module named 'tkinter'`** → `sudo apt install python3-tk`.

**`python: command not found`** → no Ubuntu o comando é `python3`.

**"codeml não encontrado"** → Linux: `sudo apt install paml`. Windows: confira se
`bin/codeml.exe` existe (reinstale com `install.bat`).

**Um gene aparece como FALHOU** → o motivo está na janela, em `genes_status.tsv` e
no `batch_analysis_log.txt` (com a última linha que o codeml escreveu). O EasyPAML
nunca fica "executando" para sempre: codeml parado sem usar CPU por 5 minutos é
encerrado e reportado.

**Windows: a janela abre e fecha rápido** → rode `install.bat` primeiro. Se persistir,
abra o `cmd` na pasta e rode `.venv\Scripts\python.exe EasyPAML.py` para ver a mensagem.

---

## Estrutura

```
EasyPAML/
├── EasyPAML.py           ponto de entrada (janela)
├── easypaml_cli.py       modo linha de comando
├── install.sh / EasyPAML.sh     instalador e lançador Linux/macOS
├── install.bat / EasyPAML.bat   instalador e lançador Windows
├── requirements.txt      dependências (versões mínimas); requirements-lock.txt (exatas testadas)
├── bin/codeml.exe        codeml para Windows (PAML 4.9j)
├── src/                  código
├── tests/                testes (pytest) e dados simulados
├── exemplos_teste/       dados de exemplo
├── METODOS.md            métodos detalhados
└── CHANGELOG.md          histórico de versões
```

---

## English summary

EasyPAML is a GUI and command-line wrapper for PAML/codeml site, branch and
branch-site models, with LRTs (including M8a vs M8), Benjamini-Hochberg
correction and BEB site tables reported in the user's alignment numbering.

- Windows: double-click `install.bat`, then `EasyPAML.bat`.
- Linux: `sudo apt install git python3-pip python3-venv python3-tk paml`, then
  `git clone https://github.com/Hesatum/EasyPAML.git && cd EasyPAML && ./install.sh && ./EasyPAML.sh`.
- CLI: `.venv/bin/python easypaml_cli.py --help`. Methods: [METODOS.md](METODOS.md).
- The interface follows the system language (Portuguese or English; default English).

---

## Licença e citação

MIT. Cite o PAML: Yang Z (2007) *PAML 4: Phylogenetic Analysis by Maximum Likelihood.*
Mol Biol Evol 24:1586–1591. Ao publicar, informe a versão do EasyPAML e do codeml e os
parâmetros (todos estão em `run_config.json`; ver METODOS.md).
