# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.3.0] 2026-10-07

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
- When M0 runs with other models, each model is fitted from three initial ω
  values and the best lnL is kept; the saved `.ctl` and output were those of the
  last start, so the panel, BEB sites and regenerated files could come from a
  worse fit. The files of the best start are kept now.
- Runs that fail (time limit, idle, error) no longer enter the LRT. Their partial
  output is kept as `*_results_FAILED.txt`.

### Added

- M8a and the M8a vs M8 test. With automatic null models, M8 adds M7 and M8a
  (`--no-m8a`, or a click on the card of an added null, leaves it out).
- Data check before running: stop codons with their position, names missing from
  the tree with the closest match, pruned taxa, duplicate files, lengths that are
  not a multiple of 3.
- Automatic time limit per run ([docs/timing_benchmark.md](docs/timing_benchmark.md)).
- One tree per gene (`GENE.nwk` next to `GENE.fasta`, or `--tree-folder`).
- Site positions in the user's alignment next to codeml's numbering.
- Output files `methods_text.txt`, `sites_BEB.tsv` and `genes_status.tsv`; the
  EasyPAML commit in About, `--version` and `run_config.json`.
- Progress by model with the time elapsed and an estimate of the time left.
- Charts of each test's 2Δℓ distribution and of ω per gene, in the Summary (with the
  genes under the mouse) and in the exported figure.
- Exports include lnL, np and df of each test, without rounding.
- Light theme, following the system by default.
- Command-line options for every `.ctl` setting; exit code 1 if any gene fails.

### Fixed

- codeml could wait forever for input after a stop codon, and Ubuntu's `codeml`
  wrapper made working runs look idle so they were stopped after 5 minutes.
- Stopping a run marked the interrupted genes as failed.
- The results panel counted failed genes as 0 and could show one as significant.
- "Ignore stop codons" had no effect, and sequences missing from the tree were
  dropped without a warning.
- A `.fasta` and a `.phy` of the same gene counted as two genes.
- The Interpretation tab ignored M8a vs M8 and used df = 2 for it.
- The window marked a site with Pr = 0.990 as ** when codeml's value is below 0.99.
- Installers create `.venv` on Linux (Ubuntu 24.04) and Windows.
- References of M7/M8 (Yang et al. 2000), Branch and Branch-site in the model help.

## [0.2.0] 2026-10-05

- Code moved to `src/`, command-line mode (`easypaml_cli.py`), BH correction of the
  LRT, GO interpretation tab.
