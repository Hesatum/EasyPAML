# EasyPAML

EasyPAML is a desktop program (with a command-line mode) for positive selection
analysis with PAML/codeml. It runs the site models (M0, M1a, M2a, M7, M8, M8a), the
branch model and the branch-site model, computes the likelihood ratio tests (LRT)
with Benjamini-Hochberg correction, and lists the sites under selection (BEB).

- Methods (codeml parameters, LRT, BH correction): [METHODS.md](METHODS.md)
- Changes in each version: [CHANGELOG.md](CHANGELOG.md)

## Installation

### Download

- With git: `git clone https://github.com/Hesatum/EasyPAML.git`
- Without git: [download the ZIP](https://github.com/Hesatum/EasyPAML/archive/refs/heads/main.zip)
  and extract it.
- To repeat a published analysis, use the commit it cites: after cloning, run
  `git checkout COMMIT` inside the EasyPAML folder, then install as below.

### Windows

> Requirement: Python 3.8 or newer ([python.org/downloads](https://www.python.org/downloads/)).
> Tick "Add Python to PATH" during installation. `codeml.exe` is already in the `bin/` folder.

1. In the EasyPAML folder, double-click `install.bat`. It creates a `.venv`
   environment with the dependencies, which takes 1 to 3 minutes.
2. Open EasyPAML from the desktop shortcut or with `EasyPAML.bat`.

### Linux (Ubuntu/Debian) and macOS

1. System requirements (once):

   ```bash
   sudo apt update && sudo apt install -y git python3-pip python3-venv python3-tk paml
   ```

   (Fedora: `sudo dnf install git python3-pip python3-tkinter paml` · macOS: `brew install python-tk brewsci/bio/paml`)

2. Download and install:

   ```bash
   git clone https://github.com/Hesatum/EasyPAML.git
   cd EasyPAML
   ./install.sh
   ```

   `install.sh` installs the dependencies in a `.venv/` folder inside EasyPAML and
   leaves the system Python alone. If something is missing, it prints the command
   to run; then run `./install.sh` again.

3. Open:

   ```bash
   ./EasyPAML.sh
   ```

   If it fails, run `.venv/bin/python EasyPAML.py` instead.

### codeml (PAML)

- Windows: `bin/codeml.exe` (PAML 4.9j) comes with the repository.
- Linux: codeml comes from the system package, `sudo apt install paml` (PAML 4.9j on
  Ubuntu 24.04). `install.sh` tries to install it if it is missing; if there is no
  package, it downloads the official PAML 4.10.10 binary to `bin/codeml`.
- Tested versions: PAML 4.9j and 4.10.10. Every parameter is written to the `.ctl`
  file, so nothing is left to codeml's internal defaults (see METHODS.md).
- To use another codeml, set `EASYPAML_CODEML=/path/to/codeml` or pass `--codeml`
  in command-line mode.
- The codeml version used is shown under About and in `run_config.json`.

## Using the window

The window follows the system language (English or Portuguese) and the system
light or dark theme. The PT/EN buttons and the Theme setting (Auto, Light, Dark)
at the bottom left change them, and the choice is kept for the next session.

1. Alignments folder: a folder with one file per gene (`.fasta`, `.fas`, `.phy`,
   `.phylip`) containing aligned codon sequences. The program shows how many
   alignments it found.
2. Tree file: a Newick tree (`.nwk`, `.tree`, `.tre`, `.txt`), rooted or unrooted.
   It may contain extra taxa, which are pruned for each gene.
   To use one tree per gene, put `GENE.nwk` (or `.tree`/`.tre`) next to `GENE.fasta`
   in the alignments folder. That tree replaces the tree file for that gene, and if
   every gene has its own tree you do not need a tree file at all.
3. Output folder: choose a folder, or create one with "New folder".
4. Switch the models on. With "Automatic null models" on, the null model of each
   test is added for you (M8 adds M7 and M8a; M2a adds M1a). Click the card of an
   added null to leave it out.
5. Click Run. Before running, EasyPAML checks the data and lists what it found:
   stop codons (sequence and codon position), alignment names that are not in the
   tree (with the closest name), taxa that will be pruned, duplicate files for the
   same gene, and lengths that are not a multiple of 3. Choose "Fix and go back" or
   "Continue anyway".
6. A run takes from minutes to hours. With 10 sequences a model takes 1 to 5
   minutes per gene; with 30 sequences M8 takes 15 to 75 minutes, and with 60
   several hours ([measured times](docs/timing_benchmark.md)). The progress line
   shows the model running, the time elapsed and, after the first gene, an estimate
   of the time left. Stop asks first and keeps the genes already completed.
7. At the end the results panel opens. Genes that failed are listed with the reason.

The results panel has these tabs:

- Summary: a chart of the 2Δℓ distribution of a test and one of ω per gene (move the
  mouse along the x axis to see the genes), then a sentence per gene, such as "Positive selection: M8 vs M7, M8 vs M8a
  significant (q < 0.05), 21 site(s) with Pr(ω>1) ≥ 0.95", and below it p, q, ω and
  p₁ of the positive class for each test.
- LRT and p-values: lnL of each model, 2Δℓ, p (in scientific notation), q (BH), and
  ω and proportion of the positive class.
- Positive Sites: the position in your alignment and the position in the codeml
  file (they differ when columns with gaps or stop codons are removed), amino acid,
  Pr(ω>1), `*` (≥ 0.95) or `**` (≥ 0.99), and mean ω ± SE. Buttons copy or export the
  table as TSV.
- Export: Excel, CSV, HTML and a figure (PNG or PDF) with those charts, one panel per
  test. "Open results folder" opens the folder itself.
- Candidates & GO: significant genes, and GO enrichment (Fisher's exact test) if you load a GO annotation TSV.

> The mean ω of a gene is not a criterion for positive selection. It stays below 1
> even when a few sites are under strong selection. Use the LRT (q) and the site table.

## Models

| Model | Used for |
|---|---|
| M0 | one ω for the whole gene (baseline; null of Branch) |
| M1a / M2a | M2a vs M1a: site-wise positive selection (df = 2) |
| M7 / M8 | M8 vs M7: beta plus an extra class (df = 2) |
| M8a | M8 with the extra class fixed at ω = 1 (M8 vs M8a, df = 1) |
| Branch | one ω per group of labelled branches (vs M0) |
| Branch-site | episodic selection at sites on the foreground branch (#1) |

For site-wise positive selection, the usual choice is to run M8 and check both
tests. If only M8 vs M7 is significant, the signal may come from neutral sites
(ω = 1) rather than positive selection.

## Command-line mode (servers, many genes)

```bash
.venv/bin/python easypaml_cli.py --input FOLDER --tree TREE.nwk --output OUT \
    --models M1a,M2a,M7,M8,M8a --workers 8
```

Main options (`--help` lists all of them):

| Option | Default | Meaning |
|---|---|---|
| `--models` | `M1a,M2a,M7,M8,M8a` | models to run |
| `--no-auto-nulls` | | run only the listed models; by default the null of each listed model is added, as in the window (`--models M8` runs M7, M8 and M8a) |
| `--no-m8a` | | do not run M8a (no M8a vs M8 test); in the window, click the M8a card |
| `--tree-folder` | | folder with one tree per gene (`GENE.nwk`), matched by name; `GENE.nwk` in the alignments folder also works |
| `--codonfreq` | `2` (F3x4) | codeml CodonFreq (0 Fequal, 1 F1x4, 2 F3x4, 3 F61, 7 FMutSel …) |
| `--ncatg` | `10` | beta categories (M7/M8/M8a) |
| `--kappa` | `2` | initial κ (estimated) |
| `--cleandata` | `1` | remove columns with gaps, ambiguities or stop codons |
| `--ignore-stop-codons` | off | without it, genes with an internal stop codon fail and the stop position is reported |
| `--workers` | `4` | genes run in parallel |
| `--timeout` / `--idle-timeout` | automatic / 300 s | time limit per run ([how it is calculated](docs/timing_benchmark.md)) / codeml not using CPU |
| `--skip-beb`, `--two-pass` | | for thousands of genes (see METHODS.md) |
| `--codeml` | | path to codeml |
| `--rerun-all` | | run every model again; by default runs already in OUTPUT with the same alignment, tree and `.ctl` are reused |
| `--strict` | | run nothing if the initial check finds problems |
| `--lang pt\|en` | system language | language of the messages |
| `--config file.json` | | the same options as JSON |
| `--version` | | EasyPAML version and commit |

Exit code: 0 if every gene ran, 1 if any failed, 2 with `--strict` and problems in
the data.

## What the output folder contains

```
OUT/
  methods_text.txt         Methods paragraph (English) filled in with what this run did
  run_config.json          EasyPAML version and commit, codeml and Python versions, .ctl parameters, options
  batch_analysis_log.txt   full log (includes the codeml output)
  genes_status.tsv         one gene per line: ok / failed, reason, warnings
  analysis_summary.tsv     lnL, np, ω per model; ω and p₁ of the positive class; 2Δl, p and q per test
  LRT_results.txt          LRT per gene, with the methods note
  sites_BEB.tsv            BEB sites with Pr(ω>1) ≥ 0.95 (M2a, M8, Branch-site), both numberings
  M8/GENE_M8.ctl           the .ctl used (every parameter, relative paths)
  M8/GENE_M8_seq.fasta     the alignment exactly as codeml read it
  M8/GENE_M8_tree.nwk      the tree exactly as codeml read it (pruned/unrooted)
  M8/GENE_M8_results_FAILED.txt  partial output of a run that failed (left out of the LRT)
  M8/GENE_M8_results.txt   raw codeml output
  M8/GENE_M8_sitemap.json  site numbering: codeml → original alignment
```

To repeat a run by hand: `cd OUT/M8 && codeml GENE_M8.ctl`.

## Example data

`examples/` has 25 *Cereus* (cactus) genes, with 8 to 21 sequences each, and a
21-taxon tree whose names match the alignments exactly:

1. Alignments folder: `examples/alignments/`
2. Tree file: `examples/tree.nwk`
3. Output folder: a new folder
4. Switch M8 on and click Run

The 25 genes take about 1 hour with 12 CPUs. For a quick try, copy two or three
alignments to another folder and choose that folder instead. From the command line:

```bash
.venv/bin/python easypaml_cli.py --input examples/alignments \
    --tree examples/tree.nwk --output examples/results --workers 8
```

`tests/data` has a simulated data set with a known answer.

## Troubleshooting

`./EasyPAML.sh: No such file or directory`: you are outside the EasyPAML folder
(`cd EasyPAML`) or you downloaded an old version.

`ModuleNotFoundError: No module named 'tkinter'`: `sudo apt install python3-tk`.

"codeml not found": on Linux, `sudo apt install paml`. On Windows, check that
`bin/codeml.exe` exists (reinstall with `install.bat`).

A gene shows as FAILED: the reason is in the window, in `genes_status.tsv` and in
`batch_analysis_log.txt` (with the last line codeml wrote). EasyPAML never stays
"running" forever: a codeml that uses no CPU for 5 minutes is stopped and reported.

Windows, the window opens and closes at once: run `install.bat` first. If it keeps
happening, open `cmd` in the folder and run `.venv\Scripts\python.exe EasyPAML.py` to
see the message.

## License and citation

MIT. Cite PAML: Yang Z (2007) *PAML 4: Phylogenetic Analysis by Maximum Likelihood.*
Mol Biol Evol 24:1586–1591. When you publish, report the EasyPAML version and commit
(under About or with `--version`), the codeml version and the parameters. The
`methods_text.txt` file of each run puts all of this in one paragraph for you to review
(see METHODS.md).
