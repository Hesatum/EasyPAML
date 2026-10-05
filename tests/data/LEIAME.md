# Dados de teste (simulados)

Simulados com `evolverNSsites` (PAML 4.10.10), 10 táxons de primatas, árvore
`gene_exemplo.nwk`.

| Arquivo | Conteúdo |
|---|---|
| `gene_exemplo.fasta` | 300 códons; classes p = 0,6/0,3/0,1 com ω = 0,05/1/4 (10% dos sítios sob seleção positiva) |
| `gene_problematico.fasta` | o mesmo alinhamento com um stop TGA no códon 150 de *Gorilla_gorilla* e o nome `Macaca_mulata` (a árvore tem `Macaca_mulatta`) |
| `geneC.fasta` | 250 códons; p = 0,8/0,2 com ω = 0,1/1 -- **sem** seleção positiva (M7 vs M8 dá falso positivo; M8a vs M8 não) |

Sítios simulados com ω = 4 no `gene_exemplo` (numeração do alinhamento):
10 39 41 47 49 50 65 69 77 78 80 81 85 94 109 121 123 131 140 153 156 161 176
178 200 220 221 228 245 249 254 261 281

Referência (codeml 4.10.10, CodonFreq = 2, ncatG = 10, cleandata = 1):
lnL M7 = −4174,2, **lnL M8 = −4122,628318**.
