# Quick example (simulated)

Two genes simulated with `evolverNSsites` (PAML 4.10.10) on the 10-taxon tree
`tree.nwk`; the same files as `tests/data`. "Try the example" in the window loads
this folder with M8 switched on.

| File | Simulated | Expected with M7, M8 and M8a |
|---|---|---|
| `simulated_selection.fasta` | 300 codons, 10% of sites with ω = 4 | M8 vs M8a and M8 vs M7 significant; 21 BEB sites at Pr(ω>1) ≥ 0.95 |
| `simulated_no_selection.fasta` | 250 codons, ω = 0.1 and 1 only | M8 vs M7 significant (false positive), M8 vs M8a not significant |

The sites simulated with ω = 4 are listed in `tests/data/README.md`.
