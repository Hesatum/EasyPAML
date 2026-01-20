# 🧬 EasyPAML - Interface Intuitiva para Análise de Seleção Positiva

**EasyPAML** é uma aplicação de fácil uso que permite a análise de seleção positiva em sequências genômicas usando o motor **PAML/CODEML**, sem necessidade de conhecimentos técnicos avançados.

---

## 📋 O que é Seleção Positiva?

Seleção positiva ocorre quando genes evoluem mais rapidamente do que esperado sob neutralidade evolutiva. É medida pela razão **ω (ômega) = dN/dS**, onde:
- **dN** = substituições não-sinônimas (alteram aminoácido)
- **dS** = substituições sinônimas (não alteram aminoácido)

Interpretação:
- **ω < 1**: Seleção purificadora (mantém sequência)
- **ω ≈ 1**: Evolução neutra
- **ω > 1**: **Seleção positiva** (aceita mudanças)

---

## 🚀 Instalação Rápida (2 passos)

### 1️⃣ Baixar Python

Acesse https://www.python.org/downloads/ e baixe **Python 3.8 ou superior** (Windows).

Durante a instalação, **marque a opção "Add Python to PATH"**.

### 2️⃣ Instalar Dependências Python

Abra o terminal (cmd ou PowerShell) e execute:

```bash
pip install -r requirements.txt
```

Este comando instalará automaticamente todas as bibliotecas necessárias.

### 3️⃣ Executar EasyPAML

Duplo-clique em `EasyPAML.py` ou execute no terminal:

```bash
python EasyPAML.py
```

A janela da aplicação abrirá automaticamente.

---

## � Como Citar

### EasyPAML

Se você usar o EasyPAML em seu trabalho acadêmico, cite:

**Formato ABNT:**
```
SILVA, M. V. EasyPAML: Interface gráfica para análise de seleção positiva usando PAML/CODEML. Versão 0.1.0. 2026. Disponível em: https://github.com/Hesatum/EasyPAML. Acesso em: [data].
```

**Formato BibTeX:**
```bibtex
@software{silva2026easypaml,
  author = {Silva, Matheus Vieira da},
  title = {EasyPAML: Interface gráfica intuitiva para análise de seleção positiva usando PAML/CODEML},
  year = {2026},
  version = {0.1.0},
  date = {2026-01-20},
  url = {https://github.com/Hesatum/EasyPAML}
}
```

**Formato APA:**
```
Silva, M. V. (2026). EasyPAML: Interface gráfica intuitiva para análise de seleção positiva usando PAML/CODEML (Versão 0.1.0) [Software]. GitHub. https://github.com/Hesatum/EasyPAML
```

### PAML (obrigatório)

**Você também DEVE citar o PAML:**

```
Yang, Z. (2007). PAML 4: Phylogenetic Analysis by Maximum Likelihood. Molecular Biology and Evolution, 24(8), 1586-1591. https://doi.org/10.1093/molbev/msm088
```

**BibTeX para PAML:**
```bibtex
@article{yang2007paml,
  author = {Yang, Ziheng},
  title = {PAML 4: Phylogenetic Analysis by Maximum Likelihood},
  journal = {Molecular Biology and Evolution},
  volume = {24},
  number = {8},
  pages = {1586--1591},
  year = {2007},
  doi = {10.1093/molbev/msm088}
}
```

> **Nota:** Esta é uma citação provisória para a versão inicial do software.
---

## �📦 Sobre o PAML/CODEML Incluído

**EasyPAML já inclui o programa CODEML (na pasta `bin/`) para sua conveniência!**

- **O que é PAML?** Pacote de programas para análises filogenéticas desenvolvido por Ziheng Yang
- **Licença:** GPL-3.0 (permite redistribuição)
- **Código-fonte:** https://github.com/abacus-gene/paml


**Citação obrigatória para PAML:**
- Yang, Z. (2007). PAML 4: Phylogenetic Analysis by Maximum Likelihood. Molecular Biology and Evolution, 24(8), 1586-1591.

Se preferir usar sua própria versão do PAML, basta substituir o executável na pasta `bin/`.

---

## 💻 Usando a Aplicação

### Passo 1: Selecionar Arquivos de Entrada

1. Clique em **"📂 Selecionar Pasta de Dados"**
2. Navegue até a pasta contendo seus arquivos FASTA

**Formatos esperados:**
- `.fasta` ou `.fas` - Sequências de DNA
- Arquivos podem conter múltiplos genes

### Passo 2: Escolher Modelos de Análise

A aplicação oferece **9 modelos evolutivos diferentes**:

#### Site Models (Seleção de Codon)
- **M0**: Modelo nulo (ω constante)
- **M1a**: 0 ≤ ω₀ < 1, ω₁ = 1
- **M2a**: 0 ≤ ω₀ < 1, ω₁ = 1, ω₂ ≥ 1 (detecta seleção positiva)
- **M7**: Distribuição beta de ω (0,1)
- **M8**: M7 + categoria com ω > 1

#### Branch Models (Seleção em Linhagens)
- **Branch**: ω diferente entre linhagens
- **Branch-Site**: Detecta seleção em ramos específicos
- **Branch-Site null**: Versão nula para comparação

#### 🆕 Nova Feature: Auto-Seleção de Modelos Nulos!

**A partir de agora, você não precisa se preocupar com qual modelo nulo usar!**

Quando você seleciona um modelo alternativo, EasyPAML **automaticamente executa o modelo nulo correspondente** e calcula a **comparação LRT com p-value automático**:

```
Selecionado → Auto-Adiciona → Compara
────────────────────────────────────────
M2a         → M1a           → M1a vs M2a (p-value automático!)
M8          → M7            → M7 vs M8
BranchSite → BranchSite null → null vs BranchSite
Branch      → M0            → M0 vs Branch
```

**Exemplo:** Selecione apenas **M2a** e o sistema automaticamente:
1. ✅ Executa M2a
2. ✅ Executa M1a (nulo correspondente)  
3. ✅ Calcula p-value para comparação
4. ✅ Salva em `LRT_results.txt` com interpretação automática!

**Dica:** Para máxima confiança, selecione **M2a E M8** - o sistema rodará ambos os testes (M1a vs M2a e M7 vs M8) e você terá confirmação dupla!

### Passo 3: Marcar Ramos no Árvore 

Para modelos **Branch** e **BranchSite**:

1. O programa exibirá a **árvore filogenética dos seus genes**
2. **Clique** nos ramos que deseja marcar como "foreground" (sob possível seleção)
3. Os ramos selecionados aparecem coloridos
4. Clique novamente para desmarcar

### Passo 4: Executar Análise

Clique em **"▶️ Iniciar Análise"**

A aplicação:
- ✅ Gerará arquivos de controle CODEML automaticamente
- ✅ Executará análises em **paralelo** para velocidade
- ✅ Analisará resultados estatisticamente
- ✅ Salvará tudo em `resultados/`

### Passo 5: Visualizar Resultados

Clique em **"📊 Ver Resultados"** (ainda em desenvolvimento) para abrir o visualizador com:

**📊 Gráficos**
- Barras de ω por gene
- Heatmap de log-likelihood
- Boxplot de distribuição de ω

**📋 Tabela Interativa**
- Filtre genes por nome
- Veja todos os valores calculados

**🔬 Seleção Positiva**
- Lista de genes com ω > 1
- Automaticamente destacados

**💾 Exportar**
- Excel (.xlsx) com todos os dados
- CSV para análise externa
- Gráficos como PNG

---

## 📁 Estrutura de Arquivos do Projeto

```
EasyPAML/
├── main.py                    ← Clique aqui para iniciar!
├── requirements.txt           ← Dependências (pip install -r)
├── README.md                  ← Este arquivo
│
├── src/
│   ├── backend/
│   │   ├── __init__.py
│   │   └── codeml_backend.py  ← Motor de análise (não editar)
│   │
│   └── gui/
│       ├── __init__.py
│       ├── main_gui.py        ← Interface gráfica (não editar)
│       └── results_viewer.py  ← Visualizador de resultados (não editar)
│
├── bin/
│   └── codeml.exe            ← Executável PAML (incluído)
│
└── resultados/               ← Criado automaticamente
    └── [dados de análises anteriores]
```

---

## 🔧 Troubleshooting (Solução de Problemas)

### ❌ "Python não é reconhecido"

**Solução:** Reinstale Python marcando "Add Python to PATH"

### ❌ "Erro ao instalar dependências"

Tente:
```bash
python -m pip install --upgrade pip
pip install -r requirements.txt --upgrade
```

### ❌ "CODEML não encontrado"

Verifique se a pasta `bin/` existe com `codeml.exe` incluído.

### ❌ "Arquivo analysis_summary.tsv não encontrado"

Isso significa que a análise anterior não completou. Verifique:
- A pasta de dados tem arquivos FASTA válidos?
- Há espaço em disco suficiente?
- Nenhum CODEML ainda está rodando?

### ❌ "Erro de importação (ModuleNotFoundError)"

Execute novamente:
```bash
pip install -r requirements.txt
```

---

## 📊 Interpretação de Resultados

### Valores Retornados

| Campo | Significado |
|-------|------------|
| `lnL` | Log-likelihood - quanto melhor o ajuste do modelo |
| `np` | Número de parâmetros (modelos com mais parâmetros ajustam melhor) |
| `ω` | dN/dS - razão de substituição |
| `p-value` | Significância estatística (< 0.05 é significativo) |

### Teste de Razão de Verossimilhança (LRT)

EasyPAML calcula automaticamente:

```
2 × (ln L_modelo_complexo - ln L_modelo_nulo)
```

Este valor segue distribuição chi-quadrado. Se **p < 0.05**, há evidência de seleção positiva.

### Exemplo de Interpretação

**Gene X com M2a vs M1a:**
- M1a (nulo): lnL = -2000, ω = 0.8
- M2a (seleção): lnL = -1950, ω₂ = 1.3
- LRT: 2 × ((-1950) - (-2000)) = 100
- **p-value = 0.001** ✅ Evidência de seleção positiva!

---

## 🎓 Recursos de Aprendizado

- **Manual PAML**: http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf
- **Papers sobre ω**: Busque "dN/dS positive selection"
- **Tutoriais em vídeo**: YouTube "CODEML tutorial"

---

## 👨‍💻 Para Desenvolvedores

Para expandir ou modificar a aplicação:

1. **Backend** (`src/backend/codeml_backend.py`): Lógica de análise e PAML
2. **Frontend** (`src/gui/main_gui.py`): Interface gráfica CustomTkinter
3. **Visualização** (`src/gui/results_viewer.py`): Gráficos e exportação

Todas as modificações devem preservar as assinaturas de funções públicas.

---

## 📄 Licença

EasyPAML é distribuído como ferramenta educacional.

---

## 🙏 Suporte

Em caso de dúvidas:
1. Verifique este README
2. Procure a seção "Troubleshooting"
3. Verifique se seus arquivos FASTA estão válidos

---

## ✨ Próximas Análises

Após completar uma análise:
- Resultados são salvos em `resultados/`
- Você pode carregar uma **nova pasta de dados** e refazer
- Resultados antigos são preservados automaticamente

**Boa análise!** 🧬✨

