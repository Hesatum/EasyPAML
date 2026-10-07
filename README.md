# EasyPAML

EasyPAML runs PAML/codeml selection analyses on many genes at once, from a window or
the command line. It fits the site models (M0, M1a, M2a, M7, M8, M8a), the branch and
branch-site models, tests them with likelihood ratio tests corrected for the number
of genes (Benjamini-Hochberg), and lists the sites under selection (BEB).

[METHODS.md](METHODS.md) describes what the program does to the data, for a methods
section or a review. [CHANGELOG.md](CHANGELOG.md) lists the changes in each version.

## Installation

**Linux (Ubuntu/Debian).**

```bash
sudo apt update && sudo apt install -y git python3-pip python3-venv python3-tk paml
git clone https://github.com/Hesatum/EasyPAML.git
cd EasyPAML
./INSTALL_LINUX_MAC.sh
./EasyPAML.sh
```

On Fedora, the first line is `sudo dnf install git python3-pip python3-tkinter paml`.
The installer keeps everything inside the EasyPAML folder and creates `EasyPAML.sh`
(the window) and `easypaml-cli.sh` (the command line). If something is missing, it
prints the command to run.

**macOS.** Run `brew install python-tk brewsci/bio/paml`, then the last four lines
above.

**Windows.** Install Python 3.8 or newer from [python.org](https://www.python.org/downloads/)
with "Add Python to PATH" ticked. Download the
[ZIP](https://github.com/Hesatum/EasyPAML/archive/refs/heads/main.zip), extract it,
double-click `INSTALL_WINDOWS.bat`, and open EasyPAML from the desktop shortcut or
`EasyPAML.bat`. If Windows shows "Windows protected your PC", click "More info" and
"Run anyway": the installer is a text script you can open and read. codeml comes in
`bin/`.

**A cited version.** To repeat a published analysis, run `git checkout v0.4.0` (or the
commit it cites) in the EasyPAML folder before installing; versions up to 0.4.0 install
with `./install.sh`. A version shown as `0.4.0+3` is 3 commits after 0.4.0, so cite
its commit.

## Using the window

To see how it works, click "Try the example" and then Run: two simulated genes, one
under positive selection and one not.

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
   the closest name), duplicate files and lengths that are not a multiple of 3. Stop
   keeps the genes already done.

The results panel opens at the end:

- **Summary**: one sentence with the answer, then one test at a time. Each gene has
  its evidence over all tests (Supported, Weak signal, Not confirmed by M8a, Possible,
  Not detected; the rules are in METHODS.md), q, p, 2Δℓ, ω and the number of sites.
  The chart can be exported for a journal (PDF, SVG, 600 dpi TIFF), and the table,
  all tests (Excel) and an HTML report can be saved.
- **Positive sites**: where the sites fall along the CDS, with their positions in your
  alignment and in the codeml output.
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
./easypaml-cli.sh --input FOLDER --tree TREE.nwk --output OUT --models M8 --workers 8
```

On Windows, use `easypaml-cli.bat` with the same options.

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

`--help` lists every option.

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

## Troubleshooting

- `./EasyPAML.sh: No such file or directory`: run `./INSTALL_LINUX_MAC.sh` first, in
  the EasyPAML folder; it creates `EasyPAML.sh`.
- `No module named 'tkinter'`: `sudo apt install python3-tk`.
- "codeml not found": `sudo apt install paml` on Linux; on Windows, run
  `INSTALL_WINDOWS.bat` again. Another codeml can be set with
  `EASYPAML_CODEML=/path/to/codeml`.
- A gene FAILED: the reason is in the window, in `genes_status.tsv` and in the log. A
  codeml that uses no CPU for 5 minutes is stopped and reported.
- Windows, the window closes at once: run `INSTALL_WINDOWS.bat` again; to see the
  error, run `EasyPAML.bat` from `cmd`.

## License and citation

MIT. Cite PAML: Yang Z (2007) PAML 4: Phylogenetic Analysis by Maximum Likelihood.
*Mol Biol Evol* 24:1586–1591. Report the EasyPAML version and commit (under About or
with `--version`) and the parameters; `methods_text.txt` puts them in one paragraph
for you to check.
