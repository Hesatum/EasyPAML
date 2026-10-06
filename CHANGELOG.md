# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased] 0.3.0.dev0

### Changes that affect results

- Default CodonFreq is now 2 (F3x4). The interface labelled 7 as "F3×4
  (recommended)" while the generated `.ctl` used FMutSel, so lnL values and BEB
  sites change.
- Every `.ctl` parameter is written explicitly, including `ncatG = 10`. Before,
  codeml used 4 beta categories because each model runs alone.
- M0 vs M1a uses df = 1 (was 2).
- `--warm-start-m0` uses `fix_blength = 1` (starting values). It used 2, which
  fixes branch lengths at the M0 values.
- p-values use `chi2.sf` (`1 − chi2.cdf` rounded large statistics to p = 0), with
  the same rule in the files, the panel and the exports.
- Runs that fail (time limit, idle, error) no longer enter the LRT. Their partial
  output is kept as `*_results_FAILED.txt`.

### Added

- M8a and the M8a vs M8 test. With automatic null models, M8 adds M7 and M8a; the
  M8a can be turned off (`--no-m8a`, "Include M8a").
- Data check before running: stop codons with position, names missing from the
  tree with the closest match, pruned taxa, duplicate files, length not a multiple
  of 3 (`--strict` on the command line).
- Automatic time limit per run, scaled to the model and gene size
  ([docs/timing_benchmark.md](docs/timing_benchmark.md)); `--timeout` and a field in
  Advanced settings override it.
- One tree per gene: `GENE.nwk` next to `GENE.fasta`, or `--tree-folder`.
- Site numbering in the user's alignment next to codeml's numbering.
- `methods_text.txt` (a methods paragraph for each run), `sites_BEB.tsv` (BEB sites of
  every gene), and warnings per gene in `genes_status.tsv`.
- EasyPAML commit in About, `--version`, the log and `run_config.json`, also in ZIP
  downloads.
- Reproducibility: `.ctl` with relative paths, with the alignment and tree codeml
  read; `run_config.json` with versions and parameters per model;
  `requirements-lock.txt`.
- Results panel: Summary tab, q-values, ω and p₁ of the positive class, TSV export of
  sites, failed genes and warnings shown per gene.
- Command line: `--codonfreq`, `--ncatg`, `--kappa`, `--omega`, `--cleandata`,
  `--ignore-stop-codons`, `--codeml`, `--idle-timeout`, `--strict`, `--lang`,
  `--verbose`; exit code 1 if any gene fails.
- Light theme. The window follows the system theme by default; Theme (Auto, Light,
  Dark) changes it and is remembered.
- Relaxed and interleaved PHYLIP.
- Tests (`pytest`) and simulated data with a known answer (`tests/data`).

### Fixed

- codeml could wait forever for "Enter" after a stop codon. It now runs with
  standard input closed, a time limit and an idle check; Stop ends every codeml
  process.
- With the Debian/Ubuntu `paml` package, `/usr/bin/codeml` is a shell script that
  runs codeml as a child process. The idle check measured only the script and
  stopped working runs after 5 minutes; it now counts child processes.
- The results panel counted 0 failed genes and could show a failed gene as
  significant.
- "Ignore stop codons" had no effect.
- "ANALYSIS COMPLETE" appeared when every gene had failed.
- Sequences missing from the tree were dropped without a warning in the window.
- A `.fasta` and a `.phy` of the same gene counted as two genes.
- The Interpretation tab ignored M8a vs M8 and recomputed p with df = 2.
- M7/M8 reference (Yang et al. 2000); Branch and Branch-site references.
- Linux install: `install.sh` creates `.venv` (works with Ubuntu 24.04 and PEP 668)
  and prints the command for anything missing. Windows: `install.bat` creates `.venv`.
- Window: folder picker on Linux, texts no longer cut off, results window fits the
  screen and closes with Esc, contrast of at least 4.5:1, mouse wheel no longer
  changes the CPU count, model help in both languages.
- Example data: one alignment set and one tree with matching names.

## [0.2.0] 2026-10-05

- Code moved to `src/`, command-line mode (`easypaml_cli.py`), BH correction of the
  LRT, GO interpretation tab.
