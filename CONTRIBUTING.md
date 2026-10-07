# Contributing

## Running from the source

```bash
./INSTALL_LINUX_MAC.sh            # creates .venv and the launchers
.venv/bin/python -m pytest -q     # tests
```

`EASYPAML_TEST_CODEML=/path/to/codeml` also runs the test that uses a real codeml.

The batch runner can be driven from Python:

```python
from pathlib import Path
from src.backend.codeml_backend import CodemlBatchAnalysis

app = CodemlBatchAnalysis()
app.config = {'input_folder': Path('alignments'), 'tree_file': Path('tree.nwk'),
              'output_folder': Path('out'), 'models': ['M7', 'M8', 'M8a'],
              'n_workers': 8, 'run_lrt': True}
summary = app.run_batch_analysis()   # {'total', 'ok', 'failed', 'failures', ...}
```

When investigating a run, start with `OUT/run_config.json` and `OUT/genes_status.tsv`.

## Layout of the repository

```
EasyPAML/
├── INSTALL_LINUX_MAC.sh  installer; writes EasyPAML.sh and easypaml-cli.sh
├── INSTALL_WINDOWS.bat   installer; writes EasyPAML.bat, easypaml-cli.bat and a shortcut
├── bin/codeml.exe        codeml for Windows (PAML 4.9j)
├── src/                  code; entry points easypaml_window.py and easypaml_cli.py
├── tests/                tests (pytest) and simulated data
├── examples/             example data
├── tools/                requirements.txt (minimum versions), requirements-lock.txt
│                         (tested versions) and developer tools
├── docs/                 timing benchmark
├── METHODS.md            detailed methods
├── CONTRIBUTING.md       notes for developers
└── CHANGELOG.md          version history
```

## Where things are

| File | Contents |
|---|---|
| `src/backend/codeml_backend.py` | batch runner (`run_batch_analysis`, `_process_gene`, `_run_single_analysis`), LRT, summary files |
| `src/backend/ctl_params.py` | every `.ctl` parameter and the CodonFreq table |
| `src/backend/lrt_stats.py` | model pairs, df, p-values, BH; shared by the backend and the panel |
| `src/backend/alignment_io.py` | FASTA/PHYLIP reading, stop codons, the cleandata rule |
| `src/backend/preflight.py` | checks before running |
| `src/backend/site_map.py` | codeml site numbering to alignment numbering |
| `src/backend/sites_parser.py` | BEB/NEB tables, positive class, ω per branch |
| `src/backend/timeouts.py` | automatic time limit ([docs/timing_benchmark.md](docs/timing_benchmark.md)) |
| `src/backend/messages.py` | backend messages (English and Portuguese) |
| `src/gui/main_gui.py`, `results_viewer.py` | main window, results panel |
| `src/gui/gui_texts.py` | every interface text (English and Portuguese) |
| `tools/benchmark_vs_raw.py` | compares plain codeml with the EasyPAML command line on the same gene |

## Before changing behaviour

- Scientific behaviour is described in [METHODS.md](METHODS.md). Update it, and the
  changelog, when results can change.
- `--warm-start-m0` must keep `fix_blength = 1`; 2 fixes branch lengths in PAML.
  Its speed-up figures in `--help` were measured with 2 and need to be measured again.
- codeml dominates run time. Measure with `tools/benchmark_vs_raw.py` before
  assuming the Python wrapper is the bottleneck.
- Interface texts go through `gui_texts.py` in both languages.
