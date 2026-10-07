# EasyPAML

EasyPAML runs PAML/codeml selection analyses on many genes at once, from a window or
the command line. It fits the site models (M0, M1a, M2a, M7, M8, M8a), the branch and
branch-site models, tests them with likelihood ratio tests corrected for the number
of genes (Benjamini-Hochberg), and lists the sites under selection (BEB).

[METHODS.md](METHODS.md) describes what the program does to the data, for a methods
section or a review. [CHANGELOG.md](CHANGELOG.md) lists the changes in each version.

## Quick start (Linux)

```bash
sudo apt update && sudo apt install -y git python3-pip python3-venv python3-tk paml
git clone https://github.com/Hesatum/EasyPAML.git
cd EasyPAML
./install.sh
./EasyPAML.sh
```

In the window, click "Try the example" and then Run. The example has two simulated
genes with a known answer and takes a few minutes.

The same analysis from the command line:

```bash
.venv/bin/python easypaml_cli.py --input examples/quick \
    --tree examples/quick/tree.nwk --output example_results --models M8
```

## Installation

**Linux.** Run the first line of the quick start once (Fedora: `sudo dnf install git
python3-pip python3-tkinter paml`). `install.sh` puts the tested versions of the
dependencies in a `.venv/` folder and leaves the system Python alone; if something is
missing, it prints the command to run. If `./EasyPAML.sh` fails, run
`.venv/bin/python EasyPAML.py`.

**macOS.** `brew install python-tk brewsci/bio/paml`, then the `git clone`,
`./install.sh` and `./EasyPAML.sh` lines of the quick start.

**Windows.** Install Python 3.8 or newer from [python.org](https://www.python.org/downloads/)
with "Add Python to PATH" ticked. Download the
[ZIP](https://github.com/Hesatum/EasyPAML/archive/refs/heads/main.zip), extract it,
double-click `install.bat` (1 to 3 minutes) and open EasyPAML from the desktop
shortcut. `codeml.exe` (PAML 4.9j) comes in `bin/`.

**codeml.** On Linux it comes from the `paml` package; if there is none, `install.sh`
downloads PAML 4.10.10 to `bin/codeml`. EasyPAML was tested with PAML 4.9j and
4.10.10. Another codeml can be set with `EASYPAML_CODEML=/path/to/codeml` or
`--codeml`.

**A cited version.** To repeat a published analysis, run `git checkout v0.4.0` (or
the commit it cites) in the EasyPAML folder before `./install.sh`. A version shown as
`0.4.0+3` is 3 commits after 0.4.0; cite its commit.

## Using the window

The window opens in English; PT/EN and the theme are at the bottom left.

1. **Alignments folder**: one aligned codon file per gene (`.fasta`, `.fas`, `.phy`,
   `.phylip`).
2. **Tree**: one Newick tree for all genes (extra taxa are pruned for each gene), or
   a folder with one tree per gene, paired with the genes by file name. The names
   IQ-TREE and RAxML give are recognised (`GENE.fasta.treefile`,
   `RAxML_bestTree.GENE`); a near name is never paired, and the program says which
   file to rename.
3. **Output folder**.
4. **Models**. With "Automatic null models" on, M8 adds M7 and M8a and M2a adds M1a.
   Click an added null to leave it out.
5. **Run**. A data check comes first: stop codons, names missing from the tree (with
   the closest name), duplicate files and lengths that are not a multiple of 3.

With 10 sequences a model takes 1 to 5 minutes per gene; with 30, M8 can take over
an hour ([measured times](docs/timing_benchmark.md)). Stop asks first and keeps the
genes already done.

The results panel opens at the end:

- **Summary**: one sentence with the answer, then one test at a time. Each gene has
  its evidence over all tests (Supported, Weak signal, Not confirmed by M8a, Possible,
  Not detected; the rules are in METHODS.md), q, p, 2Δℓ, ω and the number of sites.
  The chart of 2Δℓ or ω per gene can be shown or hidden and exported for a journal
  (PDF, SVG, 600 dpi TIFF). The table, all tests (Excel) and an HTML report can be
  exported.
- **Positive Sites**: where the sites fall along the CDS, and their positions in your
  alignment and in the codeml output, Pr(ω>1) and ω.
- **Branch** (when Branch or Branch-site ran): the tree with the labelled branches and
  their ω.

The mean ω of a gene does not decide anything: it often stays below 1 when a few
sites are under strong selection. Use q and the sites.

## Models

| Model | Used for |
|---|---|
| M0 | one ω for the whole gene; null of Branch |
| M1a / M2a | M2a vs M1a: sites under positive selection (df = 2) |
| M7 / M8 | M8 vs M7 (df = 2); neutral sites alone can make it significant |
| M8a | M8 vs M8a (df = 1): the stricter test of M8 |
| Branch | one ω per group of labelled branches, vs M0 |
| Branch-site | sites under selection on the labelled (#1) branches |

For sites under positive selection, run M8 and read both of its tests.

## Command line

```bash
.venv/bin/python easypaml_cli.py --input FOLDER --tree TREE.nwk --output OUT --models M8 --workers 8
```

| Option | Default | Meaning |
|---|---|---|
| `--models` | `M1a,M2a,M7,M8,M8a` | models to run; the null of each is added unless `--no-auto-nulls` |
| `--tree-folder` | | one tree per gene, paired by file name |
| `--codonfreq` | `2` (F3x4) | codeml CodonFreq |
| `--ncatg` | `10` | beta categories (M7, M8, M8a) |
| `--workers` | `4` | genes run in parallel |
| `--warm-start-m0` | off | start the other models from M0 (see METHODS.md) |
| `--skip-beb`, `--two-pass` | | for thousands of genes (see METHODS.md) |
| `--rerun-all` | | run again models already in OUT with the same input |
| `--lang pt` | `en` | messages in Portuguese |

`--help` lists every option. The exit code is 0 when every gene ran and 1 when any
failed.

## Output folder

```
OUT/
  methods_text.txt       methods paragraph describing this run
  run_config.json        versions, every .ctl parameter, options, per-gene trees
  analysis_summary.tsv   lnL, np and ω per model; 2Δℓ, p and q per test
  LRT_results.txt        the tests per gene
  sites_BEB.tsv          sites with Pr(ω>1) ≥ 0.95 ("none" for a gene without)
  genes_status.tsv       ok or failed per gene, with the reason and warnings
  batch_analysis_log.txt full log, with codeml's output
  M8/GENE_M8.ctl         the .ctl used; cd OUT/M8 && codeml GENE_M8.ctl repeats it
  M8/GENE_M8_results.txt raw codeml output
```

## Example data

`examples/quick/` has the two simulated genes of "Try the example".
`examples/alignments/` has 25 *Cereus* (cactus) genes with the 21-taxon tree
`examples/tree.nwk`; with M8 they take about 1 hour on 12 CPUs.

## Troubleshooting

- `./EasyPAML.sh: No such file or directory`: run it from the EasyPAML folder.
- `No module named 'tkinter'`: `sudo apt install python3-tk`.
- "codeml not found": `sudo apt install paml` on Linux; on Windows, run `install.bat`
  again.
- A gene FAILED: the reason is in the window, in `genes_status.tsv` and in the log. A
  codeml that uses no CPU for 5 minutes is stopped and reported.
- Windows, the window closes at once: run `install.bat`; to see the error, run
  `.venv\Scripts\python.exe EasyPAML.py` in `cmd`.

## License and citation

MIT. Cite PAML: Yang Z (2007) PAML 4: Phylogenetic Analysis by Maximum Likelihood.
*Mol Biol Evol* 24:1586–1591. Report the EasyPAML version and commit (under About or
with `--version`) and the parameters; `methods_text.txt` puts them in one paragraph
for you to check.
