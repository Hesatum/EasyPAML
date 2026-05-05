# EasyPAML

Interface gráfica para análise de seleção positiva com **PAML/CODEML** — sem linha de comando, sem configuração manual.

---

## Instalação rápida

### Windows

> Pré-requisito: **Python 3.8+**  
> Baixe em [python.org/downloads](https://www.python.org/downloads/) e marque **"Add Python to PATH"** durante a instalação.

1. [Baixe o EasyPAML](https://github.com/LAGEVOL/EasyPAML/archive/refs/heads/main.zip) e extraia a pasta
2. Dentro da pasta, dê **duplo-clique em `install.bat`**
3. Aguarde a instalação terminar (1–3 minutos)
4. Use o atalho **EasyPAML** que aparecerá na Área de Trabalho

> Próximas vezes: use o atalho da Área de Trabalho ou dê duplo-clique em `EasyPAML.bat`.

---

### Linux / macOS

```bash
# Clonar o repositório
git clone https://github.com/LAGEVOL/EasyPAML.git
cd EasyPAML

# Tornar o instalador executável e rodar
chmod +x install.sh
./install.sh

# Iniciar
./EasyPAML.sh
```

O instalador tenta instalar o CODEML automaticamente via `apt`, `dnf` ou `brew`.  
Se não funcionar automaticamente: `sudo apt-get install paml`

---

## Como usar

1. **Pasta .fas** → selecione a pasta com os alinhamentos FASTA (`.fas`)
2. **Árvore (.nwk)** → selecione o arquivo de árvore filogenética (Newick)
3. **Pasta Saída** → escolha onde salvar os resultados
4. Marque os **modelos** que deseja rodar
5. Clique em **▶ INICIAR** e aguarde

Os resultados aparecem automaticamente ao final em **Ver Resultados**.

---

## Modelos disponíveis

| Modelo | Para que serve |
|--------|----------------|
| **M0** | dN/dS único (baseline) |
| **M1a** | Modelo neutro (referência para LRT) |
| **M2a** | Detecta seleção positiva global |
| **M7** | Distribuição Beta de ω |
| **M8** | Beta + seleção positiva — **recomendado** |
| **Branch** | ω livre por ramo etiquetado |
| **Branch-site** | Seleção episódica em ramo específico |

---

## Dados de exemplo

A pasta `exemplos_teste/` tem 4 genes e uma árvore prontos para teste:

1. Selecione `exemplos_teste/amostras/` como **Pasta .fas**
2. Selecione `exemplos_teste/final-tree.txt` como **Árvore**
3. Crie uma pasta de saída qualquer
4. Marque **M8** e clique em **▶ INICIAR**

---

## Requisitos

- Python 3.8 ou superior
- Windows 10/11, Ubuntu 20.04+, Fedora 36+, macOS 12+
- CODEML (incluído no Windows; instalado automaticamente no Linux)
- Internet apenas durante a instalação

---

## Solução de problemas

**"Python não encontrado" no install.bat**  
→ Baixe o Python em [python.org/downloads](https://www.python.org/downloads/) e marque **"Add Python to PATH"** durante a instalação.

**A janela abre e fecha rápido**  
→ Abra o `install.bat` primeiro. Se persistir, abra `cmd.exe` na pasta e rode:
```
python EasyPAML.py
```

**"CODEML não encontrado" durante a análise**  
→ Windows: reinstale com `install.bat`. Linux: `sudo apt-get install paml`.

**Erro no pip install**  
→ Tente: `python -m pip install -r requirements.txt --user`

---

## Estrutura

```
EasyPAML/
├── EasyPAML.py       ← ponto de entrada
├── install.bat       ← instalador Windows  (duplo-clique aqui)
├── install.sh        ← instalador Linux/macOS
├── EasyPAML.bat      ← launcher Windows
├── requirements.txt
├── bin/codeml.exe    ← CODEML para Windows (incluído)
├── src/              ← código fonte
└── exemplos_teste/   ← dados de exemplo
```

---

## Licença

MIT. Cite o PAML original:  
Yang Z (2007) *PAML 4: Phylogenetic Analysis by Maximum Likelihood.* Mol Biol Evol 24:1586–1591.
