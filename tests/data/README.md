# Test data (simulated)

Simulated with `evolverNSsites` (PAML 4.10.10), 10 primate taxa, tree
`gene_example.nwk`.

| File | Content |
|---|---|
| `gene_example.fasta` | 300 codons; classes p = 0.6/0.3/0.1 with ω = 0.05/1/4 (10% of sites under positive selection) |
| `gene_problematic.fasta` | the same alignment with a TGA stop at codon 150 of *Gorilla_gorilla* and the name `Macaca_mulata` (the tree has `Macaca_mulatta`) |
| `geneC.fasta` | 250 codons; p = 0.8/0.2 with ω = 0.1/1, **no** positive selection (M7 vs M8 gives a false positive; M8a vs M8 does not) |

Sites simulated with ω = 4 in `gene_example` (alignment numbering):
10 39 41 47 49 50 65 69 77 78 80 81 85 94 109 121 123 131 140 153 156 161 176
178 200 220 221 228 245 249 254 261 281

Reference (codeml 4.10.10, CodonFreq = 2, ncatG = 10, cleandata = 1):
lnL M7 = −4174.2, **lnL M8 = −4122.628318**.
