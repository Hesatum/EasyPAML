# Resultados de exemplo

Gerados pelo EasyPAML 0.3.0.dev0 (modo CLI) com codeml 4.9j (pacote `paml` do
Ubuntu 24.04), padrões novos (CodonFreq = 2/F3x4, ncatG = 10, cleandata = 1):

```bash
python3 easypaml_cli.py --input exemplos_teste/amostras \
    --tree exemplos_teste/arvore_amostras.nwk --output exemplos_teste/resultados \
    --models M1a,M2a,M7,M8,M8a --workers 12
```

Para economizar espaço no repositório, ficaram de fora as cópias do alinhamento
(`MODELO/GENE_MODELO_seq.fasta`, iguais a `amostras/` só com as sequências que
estão na árvore) e os arquivos `rst`. Caminhos absolutos da máquina foram
trocados por caminhos relativos. Configuração completa em `run_config.json`.
