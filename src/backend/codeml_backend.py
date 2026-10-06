"""Batch runner for codeml: builds each .ctl, runs codeml per gene and model,
computes the LRTs and writes the summary files."""

import os
import platform
import subprocess
import tempfile
import time
import shutil
import re
import threading
import traceback
import concurrent.futures
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from threading import Thread

import pandas as pd
from scipy import stats

from .sites_parser import SitesParser
from . import lrt_stats
from .alignment_io import (AlignmentError, cleandata_kept_codons, find_stop_codons,
                           read_alignment, to_fasta)
from .ctl_params import (DEFAULT_CODONFREQ, DEFAULT_CTL_PARAMS, build_ctl_text,
                         codonfreq_label)
from .preflight import discover_per_gene_trees, group_by_gene, list_alignment_files
from .site_map import codeml_site_count, write_sitemap
from .timeouts import estimate_seconds, resolve_timeout, user_timeout

# partial output of a failed run (does not match *_results.txt)
FAILED_RESULTS_SUFFIX = '_results_FAILED.txt'
from .version import __version__, source_commit, version_string

_APP_ROOT   = Path(__file__).resolve().parent.parent.parent
_CODEML_BIN = (_APP_ROOT / 'bin' / 'codeml.exe'
               if platform.system() == 'Windows'
               else _APP_ROOT / 'bin' / 'codeml')


def find_codeml(explicit: Optional[str] = None) -> Optional[str]:
    """codeml to use: the argument, EASYPAML_CODEML, bin/codeml(.exe), then PATH."""
    for cand in (explicit, os.environ.get('EASYPAML_CODEML')):
        if cand and Path(cand).exists():
            # absolute(), not resolve(): on Debian/Ubuntu /usr/bin/codeml links to a
            # script that picks the program by the name it was called with
            return str(Path(cand).absolute())
    if _CODEML_BIN.exists():
        return str(_CODEML_BIN)
    return shutil.which('codeml')


_CODEML_VERSION_CACHE: Dict[str, Optional[str]] = {}


def codeml_version(codeml_path: Optional[str]) -> Optional[str]:
    """codeml version ('4.9j', '4.10.10'). codeml prints it only while running, so
    this runs a tiny M0 analysis in a temporary folder."""
    if not codeml_path:
        return None
    if codeml_path in _CODEML_VERSION_CACHE:
        return _CODEML_VERSION_CACHE[codeml_path]
    version = None
    try:
        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            (tdp / 's.fa').write_text(">a\nATGAAACCC\n>b\nATGAAGCCC\n>c\nATGAAACCG\n")
            (tdp / 't.nwk').write_text("(a,b,c);\n")
            (tdp / 'p.ctl').write_text("seqfile = s.fa\ntreefile = t.nwk\noutfile = o.txt\n"
                                       "noisy = 1\nseqtype = 1\nCodonFreq = 0\nmodel = 0\nNSsites = 0\n")
            out = subprocess.run([codeml_path, 'p.ctl'], cwd=td, stdin=subprocess.DEVNULL,
                                 capture_output=True, text=True, timeout=15, errors='replace')
            m = re.search(r'version\s+([0-9][\w.]*)', out.stdout + out.stderr, re.IGNORECASE)
            if m:
                version = m.group(1).rstrip(',')
    except Exception:
        version = None
    _CODEML_VERSION_CACHE[codeml_path] = version
    return version


def _clock(seconds: float) -> str:
    seconds = int(seconds)
    h, rest = divmod(seconds, 3600)
    m, s = divmod(rest, 60)
    return f"{h}:{m:02d}:{s:02d}" if h else f"{m}:{s:02d}"


class CodemlBatchAnalysis:
    """Runs codeml over a folder of genes."""
    
    MODEL_CONFIGS = {
        'M0': {
            'description': 'Homogeneous model - one ω for all sites',
            'model': 0,
            'NSsites': 0,
            'fix_omega': 0,
            'omega': 0.5,
            'CodonFreq': DEFAULT_CODONFREQ
        },
        'M1a': {
            'description': 'Nearly Neutral - ω < 1 or = 1',
            'model': 0,
            'NSsites': 1,
            'fix_omega': 0,
            'omega': 0.5,
            'CodonFreq': DEFAULT_CODONFREQ
        },
        'M2a': {
            'description': 'Positive Selection - adds ω > 1 class',
            'model': 0,
            'NSsites': 2,
            'fix_omega': 0,
            'omega': 0.5,
            'CodonFreq': DEFAULT_CODONFREQ
        },
        'M7': {
            'description': 'Beta distribution - ω < 1',
            'model': 0,
            'NSsites': 7,
            'fix_omega': 0,
            'omega': 0.5,
            'CodonFreq': DEFAULT_CODONFREQ
        },
        'M8': {
            'description': 'Beta + ω - adds ω > 1 class',
            'model': 0,
            'NSsites': 8,
            'fix_omega': 0,
            'omega': 0.5,
            'CodonFreq': DEFAULT_CODONFREQ
        },
        'M8a': {
            'description': 'M8 with the extra class fixed at ω = 1 (null for M8)',
            'model': 0,
            'NSsites': 8,
            'fix_omega': 1,
            'omega': 1.0,
            'CodonFreq': DEFAULT_CODONFREQ
        },
        'Branch': {
            'description': 'Branch model - different ω for foreground',
            'model': 2,
            'NSsites': 0,
            'fix_omega': 0,
            'omega': 0.5,
            'CodonFreq': DEFAULT_CODONFREQ
        },
        'Branch-site': {
            'display_name': 'Branch-site',
            'description': 'Branch-site model - ω varies across sites and branches',
            'model': 2,
            'NSsites': 2,
            'fix_omega': 0,
            'omega': 0.5,
            'CodonFreq': DEFAULT_CODONFREQ
        },
        'Branch-site_null': {
            'description': 'Branch-site null model - fixes ω=1 (null for Branch-site model)',
            'model': 2,
            'NSsites': 2,
            'fix_omega': 1,
            'omega': 1.0,
            'CodonFreq': DEFAULT_CODONFREQ
        }
    }
    
    MODEL_INFO = {
        'M0': {
            'full_name': 'One-Ratio Model',
            'test_type': 'Site Model',
            'parameters': 'ω = dN/dS (constant across all sites)',
            'purpose': 'A single dN/dS ratio for all sites and branches. Null of the Branch model.',
            'interpretation': 'Gives the average ω of the gene. A gene can have ω < 1 here and still '
                              'have a few sites under positive selection.',
            'use_case': 'Baseline; needed as null for the Branch model.',
            'references': 'Goldman & Yang (1994) Mol Biol Evol 11:725-736'
        },
        'M1a': {
            'full_name': 'Nearly Neutral Model',
            'test_type': 'Site Model',
            'parameters': 'Two classes: ω₀ < 1 (purifying), ω₁ = 1 (neutral)',
            'purpose': 'Null hypothesis: allows sites under purifying and neutral selection only.',
            'interpretation': 'If M2a rejects M1a, indicates presence of positive selection (ω > 1).',
            'use_case': 'Compare against M2a to test for positive selection.',
            'references': 'Nielsen & Yang (1998) Genetics 148:929-936; Wong et al. (2004) Genetics 168:1041-1051'
        },
        'M2a': {
            'full_name': 'Positive Selection Model',
            'test_type': 'Site Model',
            'parameters': 'Three classes: ω₀ < 1, ω₁ = 1, ω₂ > 1 (positive selection)',
            'purpose': 'Alternative hypothesis: allows positive selection at specific sites.',
            'interpretation': 'M2a significantly better than M1a (q < 0.05) = evidence for positive selection. '
                              'Which sites: BEB posterior probability Pr(ω > 1) ≥ 0.95.',
            'use_case': 'Compare against M1a to identify sites under positive selection.',
            'references': 'Nielsen & Yang (1998) Genetics 148:929-936; Wong et al. (2004) Genetics 168:1041-1051; '
                          'BEB: Yang, Wong & Nielsen (2005) Mol Biol Evol 22:1107-1118'
        },
        'M7': {
            'full_name': 'Beta Distribution Model',
            'test_type': 'Site Model',
            'parameters': 'ω distributed as beta(p,q), all values ω < 1',
            'purpose': 'Null hypothesis: continuous distribution of selection, ω constrained < 1.',
            'interpretation': 'Provides smooth alternative to discrete M1a for testing positive selection.',
            'use_case': 'Alternative null hypothesis; compare against M8.',
            'references': 'Yang, Nielsen, Goldman & Pedersen (2000) Genetics 155:431-449'
        },
        'M8': {
            'full_name': 'Beta & Positive Selection Model',
            'test_type': 'Site Model',
            'parameters': 'Beta(p,q) for ω < 1, PLUS additional class with ω > 1',
            'purpose': 'Alternative hypothesis: continuous distribution + discrete class for positive selection.',
            'interpretation': 'M8 significantly better than M7 (q < 0.05) suggests positive selection, but M8 '
                              'can also beat M7 because of neutral sites (ω = 1): confirm with M8a vs M8.',
            'use_case': 'Alternative test for positive selection; compare against M7 and M8a.',
            'references': 'Yang, Nielsen, Goldman & Pedersen (2000) Genetics 155:431-449; '
                          'BEB: Yang, Wong & Nielsen (2005) Mol Biol Evol 22:1107-1118'
        },
        'M8a': {
            'full_name': 'Beta & ω = 1 (null for M8)',
            'test_type': 'Site Model',
            'parameters': 'Beta(p,q) for ω < 1, PLUS an extra class with ω fixed at 1',
            'purpose': 'Null hypothesis for M8 that allows neutral sites (ω = 1). '
                       'M7 vs M8 can reject M7 just because some sites are neutral; '
                       'M8a vs M8 only rejects if there are sites with ω > 1.',
            'interpretation': 'Reject M8a in favour of M8 = evidence for positive selection '
                              'that is not explained by neutral sites.',
            'use_case': 'Run together with M8 (added automatically as null).',
            'references': 'Swanson, Nielsen & Yang (2003) Mol Biol Evol 20:18-20; '
                          'Wong et al. (2004) Genetics 168:1041-1051'
        },
        'Branch': {
            'full_name': 'Branch Model',
            'test_type': 'Branch Model',
            'parameters': 'Different ω for designated foreground branch vs. background branches',
            'purpose': 'Tests if one or more branches evolve under different selection pressure.',
            'interpretation': 'Branch significantly better than M0 = the marked branches have a different ω '
                              'from the rest of the tree (not necessarily ω > 1).',
            'use_case': 'Mark the branches with "Label branches" first.',
            'references': 'Yang (1998) Mol Biol Evol 15:568-573'
        },
        'Branch-site': {
            'full_name': 'Branch-site Model',
            'test_type': 'Branch-site Model',
            'parameters': 'ω varies both by site AND by branch (foreground has different classes)',
            'purpose': 'Tests for positive selection affecting specific sites in specific branches.',
            'interpretation': 'Branch-site significantly better than its null = evidence for positive selection '
                              'at some sites on the foreground branch.',
            'use_case': 'Mark the foreground branch first; the null (ω₂ = 1 fixed) is added automatically.',
            'references': 'Yang & Nielsen (2002) Mol Biol Evol 19:908-917; '
                          'Zhang, Nielsen & Yang (2005) Mol Biol Evol 22:2472-2479'
        }
    }
    
    MODEL_INFO_PT = {
        'M0': {
            'full_name': 'Modelo de uma razão',
            'test_type': 'Modelo de sítio',
            'parameters': 'ω = dN/dS (igual em todos os sítios)',
            'purpose': 'Uma única razão dN/dS para todos os sítios e ramos. Nulo do modelo Branch.',
            'interpretation': 'Dá o ω médio do gene. Um gene pode ter ω < 1 aqui e ainda assim ter '
                              'poucos sítios sob seleção positiva.',
            'use_case': 'Referência; necessário como nulo do modelo Branch.',
            'references': 'Goldman & Yang (1994) Mol Biol Evol 11:725-736',
        },
        'M1a': {
            'full_name': 'Quase neutro',
            'test_type': 'Modelo de sítio',
            'parameters': 'Duas classes: ω₀ < 1 (purificadora) e ω₁ = 1 (neutra)',
            'purpose': 'Hipótese nula: só sítios sob seleção purificadora ou neutros.',
            'interpretation': 'Se o M2a for significativamente melhor que o M1a, há indício de seleção positiva (ω > 1).',
            'use_case': 'Comparar com o M2a.',
            'references': 'Nielsen & Yang (1998) Genetics 148:929-936; Wong et al. (2004) Genetics 168:1041-1051',
        },
        'M2a': {
            'full_name': 'Seleção positiva',
            'test_type': 'Modelo de sítio',
            'parameters': 'Três classes: ω₀ < 1, ω₁ = 1 e ω₂ > 1 (seleção positiva)',
            'purpose': 'Hipótese alternativa: permite seleção positiva em alguns sítios.',
            'interpretation': 'M2a significativamente melhor que o M1a (q < 0,05) = indício de seleção positiva. '
                              'Quais sítios: probabilidade posterior do BEB, Pr(ω > 1) ≥ 0,95.',
            'use_case': 'Comparar com o M1a para achar sítios sob seleção positiva.',
            'references': 'Nielsen & Yang (1998) Genetics 148:929-936; Wong et al. (2004) Genetics 168:1041-1051; '
                          'BEB: Yang, Wong & Nielsen (2005) Mol Biol Evol 22:1107-1118',
        },
        'M7': {
            'full_name': 'Distribuição beta',
            'test_type': 'Modelo de sítio',
            'parameters': 'ω segue uma beta(p, q), sempre entre 0 e 1',
            'purpose': 'Hipótese nula: distribuição contínua de ω, sem ω > 1.',
            'interpretation': 'Alternativa contínua ao M1a como nulo para testar seleção positiva.',
            'use_case': 'Comparar com o M8.',
            'references': 'Yang, Nielsen, Goldman & Pedersen (2000) Genetics 155:431-449',
        },
        'M8': {
            'full_name': 'Beta e seleção positiva',
            'test_type': 'Modelo de sítio',
            'parameters': 'Beta(p, q) para ω < 1 MAIS uma classe extra com ω livre (pode ser > 1)',
            'purpose': 'Hipótese alternativa: distribuição contínua + uma classe para seleção positiva.',
            'interpretation': 'M8 significativamente melhor que o M7 (q < 0,05) sugere seleção positiva, mas o M8 '
                              'também pode vencer o M7 por causa de sítios neutros (ω = 1): confirme com M8a vs M8.',
            'use_case': 'Comparar com o M7 e com o M8a.',
            'references': 'Yang, Nielsen, Goldman & Pedersen (2000) Genetics 155:431-449; '
                          'BEB: Yang, Wong & Nielsen (2005) Mol Biol Evol 22:1107-1118',
        },
        'M8a': {
            'full_name': 'Beta e ω = 1 (nulo do M8)',
            'test_type': 'Modelo de sítio',
            'parameters': 'Beta(p, q) para ω < 1 MAIS uma classe extra com ω fixo em 1',
            'purpose': 'Hipótese nula do M8 que admite sítios neutros (ω = 1). O M7 vs M8 pode rejeitar o '
                       'M7 só porque alguns sítios são neutros; o M8a vs M8 só rejeita se houver sítios com ω > 1.',
            'interpretation': 'M8 significativamente melhor que o M8a = indício de seleção positiva que não se '
                              'explica por sítios neutros.',
            'use_case': 'Rodar junto com o M8 (entra automaticamente como nulo).',
            'references': 'Swanson, Nielsen & Yang (2003) Mol Biol Evol 20:18-20; '
                          'Wong et al. (2004) Genetics 168:1041-1051',
        },
        'Branch': {
            'full_name': 'Modelo de ramos',
            'test_type': 'Modelo de ramos',
            'parameters': 'ω diferente nos ramos marcados (foreground) e no resto da árvore',
            'purpose': 'Testa se um ou mais ramos evoluem sob pressão seletiva diferente.',
            'interpretation': 'Branch significativamente melhor que o M0 = os ramos marcados têm ω diferente '
                              'do resto da árvore (não necessariamente ω > 1).',
            'use_case': 'Marque antes os ramos com "Marcar ramos".',
            'references': 'Yang (1998) Mol Biol Evol 15:568-573',
        },
        'Branch-site': {
            'full_name': 'Modelo de ramos e sítios',
            'test_type': 'Modelo de ramos e sítios',
            'parameters': 'ω varia entre sítios E entre ramos (o ramo foreground tem classes próprias)',
            'purpose': 'Testa seleção positiva em alguns sítios de um ramo específico.',
            'interpretation': 'Branch-site significativamente melhor que o nulo = indício de seleção positiva '
                              'em alguns sítios do ramo foreground.',
            'use_case': 'Marque antes o ramo foreground; o nulo (ω₂ = 1 fixo) entra automaticamente.',
            'references': 'Yang & Nielsen (2002) Mol Biol Evol 19:908-917; '
                          'Zhang, Nielsen & Yang (2005) Mol Biol Evol 22:2472-2479',
        },
    }

    LRT_COMPARISONS = {
        'Site Models': [
            ('M0', 'M1a', 'Tests if ω varies among sites'),
            ('M1a', 'M2a', 'Tests for positive selection'),
            ('M7', 'M8', 'Alternative test for positive selection'),
            ('M8a', 'M8', 'Tests for positive selection allowing neutral sites in the null')
        ],
        'Branch Model': [
            ('M0', 'Branch', 'Tests if ω differs in foreground branch')
        ],
        'Branch-Site Models': [
            ('Branch-site_null', 'Branch-site', 'Tests for positive selection in foreground sites')
        ]
    }
    
    NULL_MODEL_PAIRS = {
        'M2a': 'M1a',
        'M8': ['M7', 'M8a'],
        'Branch': 'M0',
        'Branch-site': 'Branch-site_null'
    }
    
    # Nulls that need fix_omega = 1. M1a must not get it: codeml already fixes
    # ω₁ = 1 there, and fix_omega would fix every ω.
    NEUTRAL_MODELS = {
        'M8a': {
            'fix_omega': 1,
            'omega': 1.0,
            'corresponding_alternative': 'M8',
            'reason': 'M8a: extra class with ω = 1 fixed (Swanson et al. 2003)'
        },
        'Branch-site_null': {
            'fix_omega': 1,
            'omega': 1.0,
            'corresponding_alternative': 'Branch-site',
            'reason': 'Branch-site null: ω₂ = 1 fixed (Zhang et al. 2005)'
        }
    }

    # legacy folder names -> current model names
    _LEGACY_MODEL_NAMES: Dict[str, str] = {
        'BranchSite_A':      'Branch-site',
        'BranchSite_A_null': 'Branch-site_null',
    }
    
    @staticmethod
    def available_cores() -> int:
        """Number of logical CPUs."""
        return os.cpu_count() or 1

    @staticmethod
    def _get_fast_tempdir() -> str:
        """Temporary folder for codeml runs: /dev/shm (RAM) on Linux when writable,
        otherwise the system temporary folder."""
        shm = Path('/dev/shm')
        if shm.exists() and shm.is_dir():
            try:
                probe = shm / f'.easypam_probe_{os.getpid()}'
                probe.write_bytes(b'\x00')
                probe.unlink()
                return str(shm)
            except OSError:
                pass
        return tempfile.gettempdir()

    def __init__(self):
        self.results = {}
        self.config = {}
        self.current_stop_count = 0
        self.current_stop_details = []
        self.current_total_genes = 0
        self.current_processed_genes = 0
        self._results_lock = threading.Lock()
        # running codeml processes, for Stop
        self._active_processes: set = set()
        self._processes_lock = threading.Lock()
        self.current_process = None
    
    @staticmethod
    def auto_complete_null_models(selected_models: List[str], include_neutral: bool = True,
                                  include_m8a: bool = True) -> List[str]:
        """Add the null model of each selected alternative (M2a -> M1a, M8 -> M7 and
        M8a, Branch -> M0, Branch-site -> its null). include_m8a=False leaves M8a
        out unless it was selected."""
        completed_models = set(selected_models)
        
        for model in selected_models:
            if model in CodemlBatchAnalysis.NULL_MODEL_PAIRS:
                nulls = CodemlBatchAnalysis.NULL_MODEL_PAIRS[model]
                if isinstance(nulls, str):
                    nulls = [nulls]
                if not include_m8a:   # include_m8a=False keeps an M8a chosen by hand
                    nulls = [m for m in nulls if m != 'M8a']
                completed_models.update(nulls)
        
        if include_neutral:
            if 'M2a' in selected_models and 'M1a' not in completed_models:
                completed_models.add('M1a')
            if 'Branch-site' in selected_models and 'Branch-site_null' not in completed_models:
                completed_models.add('Branch-site_null')
        
        return sorted(list(completed_models), key=lambda x: selected_models.index(x) if x in selected_models else 999)


    def generate_ctl_content(self, seqfile: str, treefile: str, outfile: str,
                             model_config: dict,
                             omega: float = 0.5,
                             cleandata: int = 1,
                             model_name: str = None,
                             kappa: float = None,
                             fix_blength: int = 0) -> str:
        """.ctl text with every parameter written explicitly. Values come from
        ctl_params defaults, then the model, then the global options, then the
        user's edits for this model."""
        params = dict(DEFAULT_CTL_PARAMS)
        cfg = self.config or {}
        skip = ('description', 'display_name')
        for key, val in (model_config or {}).items():
            if key not in skip and val not in (None, ''):
                params[key] = val
        for key in ('CodonFreq', 'ncatG', 'kappa', 'fix_kappa', 'icode', 'method',
                    'Small_Diff', 'getSE', 'estFreq'):
            if cfg.get(key) is not None:
                params[key] = cfg[key]
        custom = cfg.get('custom_model_params') or cfg.get('custom_model_configs') or {}
        for key, val in (custom.get(model_name) or {}).items():
            if key not in skip and val not in (None, ''):
                params[key] = val

        params['seqfile'] = seqfile
        params['treefile'] = treefile
        params['outfile'] = outfile
        params['cleandata'] = int(cleandata)
        params['fix_blength'] = int(fix_blength)

        fix_omega = int(params.get('fix_omega', 0))
        final_omega = omega
        if model_name in self.NEUTRAL_MODELS:
            fix_omega = 1
            final_omega = 1.0
        params['fix_omega'] = fix_omega
        params['omega'] = final_omega

        if kappa is not None and 0.1 <= kappa <= 20:
            params['fix_kappa'] = 0
            params['kappa'] = round(float(kappa), 4)

        return build_ctl_text(params)

    @staticmethod
    def _extract_kappa(output_path: Path) -> Optional[float]:
        """Extract the estimated κ (kappa, ts/tv ratio) from a CODEML output file.

        Used to warm-start subsequent site/branch models with the M0 estimate,
        which significantly reduces the number of optimization iterations needed.
        Returns None if extraction fails or the value is outside a sane range.
        """
        try:
            text = output_path.read_text(encoding='utf-8', errors='ignore')
            m = re.search(r'kappa\s*\(ts/tv\)\s*=\s*([\d.]+)', text, re.IGNORECASE)
            if m:
                v = float(m.group(1))
                if 0.1 <= v <= 20:
                    return v
            m = re.search(r'^\s*kappa\s+([\d.]+)', text, re.MULTILINE)
            if m:
                v = float(m.group(1))
                if 0.1 <= v <= 20:
                    return v
        except Exception:
            pass
        return None
    
    @staticmethod
    def _extract_fitted_tree(output_path: Path) -> Optional[str]:
        """Newick tree with the ML branch lengths from a codeml output file, used as
        starting values by --warm-start-m0 (fix_blength = 1)."""
        try:
            text = output_path.read_text(encoding='utf-8', errors='ignore')
            for line in reversed(text.splitlines()):
                s = line.strip()
                if (s.startswith('(')
                        and s.endswith(';')
                        and ':' in s
                        and re.search(r':\s*\d[\d.]*', s)):
                    return s
        except Exception:
            pass
        return None


    # site models: PAML needs an unrooted tree
    _SITE_MODELS_UNROOT = {'M0', 'M1a', 'M2a', 'M7', 'M8', 'M8a'}
    _SITE_WARMUP = {'M1a', 'M2a', 'M7', 'M8', 'M8a'}

    @staticmethod
    def _t(key: str, **kw) -> str:
        from .messages import t
        return t(key, **kw)

    def _emit(self, level: str, text: str) -> None:
        """One message for the user. level: info, ok, warn, error, debug or header.
        Goes to config['log_callback'] when set (window), else to stdout ('debug'
        only with config['verbose'])."""
        cb = (self.config or {}).get('log_callback')
        if cb is not None:
            try:
                cb(level, text)
                return
            except Exception:
                pass
        if level == 'debug' and not (self.config or {}).get('verbose'):
            return
        print(text, flush=True)

    def _log(self, text: str) -> None:
        """Append a line to batch_analysis_log.txt (thread-safe)."""
        path = getattr(self, '_log_path', None)
        if path is None:
            return
        with self._log_lock:
            with open(path, 'a', encoding='utf-8') as fh:
                fh.write(text.rstrip('\n') + '\n')

    def _progress(self, done: int, total: int, gene: str = '') -> None:
        cb = (self.config or {}).get('progress_callback')
        if cb is not None:
            try:
                cb(done, total, gene)
            except Exception:
                pass


    def effective_ctl_defaults(self) -> Dict[str, object]:
        """Global .ctl parameters used in this run."""
        params = dict(DEFAULT_CTL_PARAMS)
        for key in ('CodonFreq', 'ncatG', 'kappa', 'fix_kappa', 'icode', 'method',
                    'Small_Diff', 'getSE', 'estFreq'):
            if (self.config or {}).get(key) is not None:
                params[key] = self.config[key]
        params['cleandata'] = int((self.config or {}).get('cleandata', 1))
        return params

    def _write_run_config(self, codeml_path: Optional[str], codeml_ver: Optional[str],
                          genes: List[str]) -> None:
        """run_config.json: versions, .ctl parameters per model, LRT settings and
        options, enough to describe and repeat the run."""
        import json
        cfg = self.config
        skip = {'pause_event', 'stop_event', 'manual_continue_event',
                'manual_continue_all_event', 'log_callback', 'progress_callback',
                'labeled_tree_content', 'labeled_tree_branchsite'}
        wrapper = {}
        for k, v in cfg.items():
            if k in skip:
                continue
            if isinstance(v, Path):
                v = str(v)
            elif isinstance(v, dict):
                v = {kk: (str(vv) if isinstance(vv, Path) else vv) for kk, vv in v.items()}
            wrapper[k] = v
        per_model = {}
        for m in cfg['models']:
            mc = dict(self.MODEL_CONFIGS.get(m, {}))
            mc.update((cfg.get('custom_model_params') or {}).get(m, {}))
            ctl = self.effective_ctl_defaults()
            ctl.update({k: v for k, v in mc.items() if k not in ('description', 'display_name')})
            if m in self.NEUTRAL_MODELS:
                ctl['fix_omega'], ctl['omega'] = 1, 1.0
            elif cfg.get('omega') is not None:
                ctl['omega'] = float(cfg['omega'])
            ctl['CodonFreq_name'] = codonfreq_label(ctl.get('CodonFreq'))
            per_model[m] = ctl
        data = {
            'easypaml_version': __version__,
            'easypaml_commit': source_commit(),
            'codeml_path': codeml_path,
            'codeml_version': codeml_ver,
            'python_version': platform.python_version(),
            'platform': platform.platform(),
            'started_at': datetime.now().isoformat(timespec='seconds'),
            'interface': cfg.get('interface', 'cli'),
            'genes': genes,
            'ctl_parameters_by_model': per_model,
            'lrt': {f"{n} vs {a}": {'df': (info['df'] if info['df'] is not None else 'n foreground groups'),
                                    'null_distribution': ('chi2(1) (mixture 50:50 reported as reference)'
                                                          if info['boundary'] else 'chi2(df)'),
                                    'multiple_testing': 'Benjamini-Hochberg within the pair, all genes of this run'}
                    for (n, a), info in lrt_stats.PAIRS.items()
                    if n in cfg['models'] and a in cfg['models']},
            'options': wrapper,
        }
        out = Path(cfg['output_folder']) / 'run_config.json'
        out.write_text(json.dumps(data, indent=2, ensure_ascii=False, default=str), encoding='utf-8')

    def _write_methods_text(self, codeml_ver: Optional[str], n_genes: int) -> None:
        """methods_text.txt: a methods paragraph describing this run."""
        from .methods_text import build_methods_text
        cfg = self.config
        try:
            sizes = {pair: len(q) for pair, q in (getattr(self, '_lrt_qvalues', None) or {}).items()}
            text = build_methods_text(
                version=version_string(), codeml_version=codeml_ver, models=list(cfg['models']),
                ctl=self.effective_ctl_defaults(), omega0=float(cfg.get('omega', 0.5) or 0.5),
                pruned=bool(cfg.get('auto_prune_tree', True)), family_sizes=sizes,
                n_genes=n_genes, beb=not cfg.get('skip_beb'),
                masked_stops=dict(getattr(self, 'masked_stops', {})),
                excluded_taxa=dict(getattr(self, 'excluded_taxa', {})))
            out = Path(cfg['output_folder']) / 'methods_text.txt'
            out.write_text("Suggested Methods text (generated by EasyPAML from this run; "
                           "review before use)\n\n" + text + "\n", encoding='utf-8')
        except Exception as exc:
            self._log(f"[WARN] methods_text.txt: {exc}")

    def run_batch_analysis(self):
        """Run the batch. Returns self.run_summary:
        {'total', 'ok', 'failed', 'stopped', 'failures': {gene: reason}, ...}."""
        if not self.config:
            raise ValueError(
                "self.config is empty: set input_folder, tree_file, output_folder and models "
                "before calling run_batch_analysis() (see easypaml_cli.py or the GUI)."
            )

        cfg = self.config
        cfg['input_folder'] = Path(cfg['input_folder'])
        cfg['output_folder'] = Path(cfg['output_folder'])
        output_folder = cfg['output_folder']
        output_folder.mkdir(parents=True, exist_ok=True)
        log_file = output_folder / "batch_analysis_log.txt"
        self._log_path = log_file
        self._log_lock = threading.Lock()
        self.failures: Dict[str, str] = {}
        self.gene_status: Dict[str, str] = {}
        self.gene_notes: Dict[str, List[str]] = {}
        self.masked_stops: Dict[str, int] = {}
        self.excluded_taxa: Dict[str, List[str]] = {}
        self.results = {}
        self.current_stop_count = 0
        self.current_stop_details = []

        codeml_path = find_codeml(cfg.get('codeml_path'))
        codeml_ver = codeml_version(codeml_path)
        self._codeml_path = codeml_path
        ctl_defaults = self.effective_ctl_defaults()

        files = list_alignment_files(cfg['input_folder'])
        chosen, ignored = group_by_gene(files)
        genes = list(chosen.items())
        # GENE.nwk next to an alignment replaces tree_file for that gene
        if cfg.get('per_gene_trees') is None and cfg.get('auto_per_gene_trees', True):
            cfg['per_gene_trees'] = discover_per_gene_trees(
                cfg['input_folder'], cfg.get('tree_folder'), genes=set(chosen))
        self.current_total_genes = len(genes)
        self.current_processed_genes = 0
        self.runs_total = len(genes) * len(cfg['models'])
        self.runs_done = 0
        self.running: Dict[str, Tuple[str, float]] = {}
        self.run_start_time = time.time()
        self.expected = self._expected_seconds(genes, cfg['models'])
        self.expected_total = sum(sum(m.values()) for m in self.expected.values()) or 1.0
        self.expected_done = 0.0
        self.expected_finished = 0.0     # expected seconds of the runs that ended
        self.actual_finished = 0.0       # their measured seconds
        n_workers = max(1, int(cfg.get('n_workers', 1)))

        with open(log_file, 'w', encoding='utf-8') as log:
            log.write("=" * 80 + "\n")
            log.write("CODEML BATCH ANALYSIS LOG\n")
            log.write("=" * 80 + "\n")
            log.write(f"EasyPAML: {version_string()}\n")
            log.write(f"codeml: {codeml_path} (version {codeml_ver})\n")
            log.write(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            log.write(f"Input folder: {cfg['input_folder']}\n")
            log.write(f"Output folder: {output_folder}\n")
            log.write(f"Tree file: {cfg.get('tree_file')}\n")
            if cfg.get('per_gene_trees'):
                log.write(f"Per-gene trees: {len(cfg['per_gene_trees'])} "
                          f"({', '.join(sorted(cfg['per_gene_trees'])[:10])})\n")
            log.write(f"Models: {', '.join(cfg['models'])}\n")
            log.write(f"Base .ctl parameters: {dict(ctl_defaults)}\n")
            log.write("=" * 80 + "\n\n")

        self._write_run_config(codeml_path, codeml_ver, [g for g, _ in genes])

        self._emit('header', self._t('run_start', n=len(genes), m=len(cfg['models']),
                                     w=n_workers, version=codeml_ver or '?'))
        self._emit('info', self._t('run_models', models=', '.join(cfg['models']),
                                   codonfreq=codonfreq_label(ctl_defaults['CodonFreq']),
                                   ncatg=ctl_defaults['ncatG'], kappa=ctl_defaults['kappa'],
                                   cleandata=ctl_defaults['cleandata']))
        for dup in ignored:
            msg = self._t('duplicate_ignored', ignored=dup.name, used=chosen[dup.stem].name)
            self._emit('warn', msg)
            self._log(f"[WARN] {msg}")

        start_time = time.time()
        self._progress(0, len(genes))

        if not codeml_path:
            self._emit('error', self._t('no_codeml'))
            for gene, _ in genes:
                self._mark_gene_failed(gene, self._t('reason_no_codeml'))
        else:
            indexed = [(i, g, p) for i, (g, p) in enumerate(genes, 1)]
            beat = self._start_heartbeat(float(cfg.get('heartbeat', 60 if cfg.get('interface') == 'cli' else 0)))
            try:
                if n_workers > 1:
                    with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
                        list(executor.map(lambda a: self._process_gene(*a, len(genes)), indexed))
                else:
                    for item in indexed:
                        self._process_gene(*item, len(genes))
            finally:
                if beat is not None:
                    beat.set()

        total_time = time.time() - start_time
        stopped = bool(cfg.get('stop_event') is not None and cfg['stop_event'].is_set())

        # the LRT first: _save_summary() adds its p and q to the TSV
        if cfg.get('run_lrt', True) and len(cfg['models']) > 1:
            self._emit('info', self._t('lrt_start'))
            self._run_lrt_analysis()
        self._save_summary()
        self._write_methods_text(codeml_ver, len(genes))

        n_ok = sum(1 for s in self.gene_status.values() if s == 'ok')
        n_failed = sum(1 for s in self.gene_status.values() if s == 'failed')
        self.run_summary = {
            'total': len(genes), 'ok': n_ok, 'failed': n_failed, 'stopped': stopped,
            'failures': dict(self.failures), 'minutes': total_time / 60,
            'output_folder': str(output_folder), 'codeml_version': codeml_ver,
        }
        self._write_failures_file()
        try:
            self.write_sites_table(output_folder)
        except Exception as exc:
            self._log(f"[WARN] {self.SITES_TABLE}: {exc}")

        if stopped:
            final = self._t('summary_stopped', ok=n_ok, n=len(genes))
            self._emit('warn', final)
        elif n_failed:
            final = self._t('summary_failed', ok=n_ok, n=len(genes), failed=n_failed)
            self._emit('error', final)
            for gene, reason in sorted(self.failures.items()):
                self._emit('error', self._t('summary_failed_item', gene=gene, reason=reason))
        else:
            final = self._t('summary_ok', ok=n_ok, n=len(genes), minutes=total_time / 60)
            self._emit('ok', final)
        self._emit('info', self._t('results_in', path=output_folder))
        self._log("\n" + final)
        self._log(f"Total time: {total_time / 60:.1f} minutes")
        return self.run_summary

    def _add_gene_note(self, gene: str, note: str) -> None:
        with self._results_lock:
            notes = self.gene_notes.setdefault(gene, [])
            if note not in notes:
                notes.append(note)

    def _write_failures_file(self) -> None:
        """genes_status.tsv: ok or the failure reason per gene, plus warnings for genes
        that ran (masked stop codon, excluded sequence)."""
        path = Path(self.config['output_folder']) / 'genes_status.tsv'
        clean = lambda t: t.replace('\t', ' ').replace('\n', ' ')
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write("Gene\tstatus\treason\tnotes\n")
            for gene in sorted(self.gene_status):
                reason = clean(self.failures.get(gene, ''))
                notes = clean(' | '.join(self.gene_notes.get(gene, [])))
                fh.write(f"{gene}\t{self.gene_status[gene]}\t{reason}\t{notes}\n")

    def _mark_gene_failed(self, gene: str, reason: str) -> None:
        with self._results_lock:
            self.failures[gene] = reason
            self.gene_status[gene] = 'failed'
            self.results.setdefault(gene, {})
        self._emit('error', self._t('gene_failed', gene=gene, reason=reason))
        self._log(f"[FAILED] {gene}: {reason}")

    @staticmethod
    def _expected_seconds(genes, models) -> Dict[str, Dict[str, float]]:
        """Expected codeml time of each gene and model from the benchmark fit
        (timeouts.estimate_seconds); used for progress, not for results."""
        out = {}
        for gene, path in genes:
            try:
                aln = read_alignment(path)
                taxa, codons = len(aln.names), aln.length // 3
            except Exception:
                taxa, codons = 10, 300
            out[gene] = {m: estimate_seconds(m, taxa, codons) for m in models}
        return out

    # measured/expected time before any run ends: the benchmark ran 8 jobs at once
    # under nice on a laptop, so a free core is usually faster than the fit
    _SPEED_PRIOR = 0.35

    def progress_snapshot(self) -> Dict:
        """Fraction done, elapsed seconds, seconds left (None before any progress)
        and the runs going on, weighting each run by its expected time. The speed of
        this machine relative to the benchmark is taken from the runs that ended;
        a running model counts with its elapsed time, up to 90% of what it should
        take at that speed."""
        now = time.time()
        with self._results_lock:
            running = dict(self.running)
            done = self.expected_done
            speed = (self.actual_finished / self.expected_finished
                     if getattr(self, 'expected_finished', 0) > 0 else self._SPEED_PRIOR)
        speed = max(speed, 0.01)
        exp = getattr(self, 'expected', {})
        partial = sum(min((now - t0) / speed, 0.9 * exp.get(g, {}).get(m, 60.0))
                      for g, (m, t0) in running.items())
        frac = min(1.0, (done + partial) / getattr(self, 'expected_total', 1.0))
        elapsed = now - self.run_start_time
        left = elapsed / frac * (1 - frac) if frac > 0.002 and elapsed > 5 else None
        return {'frac': frac, 'elapsed': elapsed, 'left': left, 'running': running}

    def _start_heartbeat(self, every: float) -> Optional[threading.Event]:
        """On the command line, print the progress every `every` seconds while
        codeml runs, since nothing else is printed until a model finishes."""
        if every <= 0:
            return None
        done = threading.Event()

        def _beat():
            while not done.wait(every):
                snap = self.progress_snapshot()
                running = ", ".join(f"{m} ({g})" for g, (m, _) in sorted(snap['running'].items())[:3])
                left = self._t('progress_left', left=_clock(snap['left'])) if snap['left'] else ""
                self._emit('info', self._t('progress_line', pct=int(snap['frac'] * 100),
                                           elapsed=_clock(snap['elapsed']), left=left,
                                           running=running or "…"))
        threading.Thread(target=_beat, daemon=True).start()
        return done

    def _process_gene(self, idx: int, gene: str, fas_file: Path, n_total: int):
        cfg = self.config
        pause_event = cfg.get('pause_event')
        stop_event = cfg.get('stop_event')

        runs_counted = [0]
        counted_models = set()
        gene_expected = getattr(self, 'expected', {}).get(gene, {})

        def _run_finished(model, result):
            with self._results_lock:
                runs_counted[0] += 1
                self.runs_done += 1
                counted_models.add(model)
                self.expected_done += gene_expected.get(model, 0.0)
                if result.get('status') == 'success' and not result.get('reused'):
                    t0 = self.running.get(gene, (None, None))[1]
                    if t0 is not None:
                        self.expected_finished += gene_expected.get(model, 0.0)
                        self.actual_finished += time.time() - t0

        def _done():
            with self._results_lock:
                self.runs_done += max(0, len(cfg['models']) - runs_counted[0])
                self.expected_done += sum(v for m, v in gene_expected.items()
                                          if m not in counted_models)
                self.running.pop(gene, None)
                self.current_processed_genes += 1
                done = self.current_processed_genes
            self._progress(done, n_total, gene)

        if stop_event is not None and stop_event.is_set():
            return gene, {}
        if pause_event is not None:
            pause_event.wait()

        self._emit('info', self._t('gene_start', i=idx, n=n_total, gene=gene))

        try:
            aln = read_alignment(fas_file)
        except AlignmentError as exc:
            self._mark_gene_failed(gene, str(exc))
            _done()
            return gene, {}
        from .preflight import Issue, format_issue
        from .messages import get_language
        problem = None
        if len(aln.names) < 3:
            problem = Issue(gene, 'too_few_sequences', 'error', {'n': len(aln.names)})
        elif not aln.is_aligned:
            problem = Issue(gene, 'unaligned', 'error', {'lengths': sorted(set(aln.lengths))})
        elif aln.length % 3 != 0:
            problem = Issue(gene, 'not_multiple_of_3', 'error',
                            {'length': aln.length, 'remainder': aln.length % 3})
        elif aln.duplicate_names:
            problem = Issue(gene, 'duplicate_name', 'error', {'name': aln.duplicate_names[0]})
        if problem is not None:
            self._mark_gene_failed(gene, format_issue(problem, get_language()))
            _done()
            return gene, {}

        stops = find_stop_codons(aln.names, aln.seqs)
        last_codon = aln.length // 3
        internal = [s for s in stops if s[1] != last_codon]
        if stops:
            with self._results_lock:
                self.current_stop_count += len(stops)
                self.current_stop_details.extend(
                    {'gene': gene, 'sequence': n, 'codon_position': p, 'codon': c} for n, p, c in stops)
            details = "; ".join(f"{n} códon {p} ({c})" if get_language() == 'pt'
                                else f"{n} codon {p} ({c})" for n, p, c in stops[:5])
            if len(stops) > 5:
                details += f"; +{len(stops) - 5}"
            if internal and not cfg.get('ignore_stop_codons', False):
                self._mark_gene_failed(gene, self._t('reason_stop_codons', details=details))
                _done()
                return gene, {}
            msg = self._t('warn_stops_masked', gene=gene, count=len(stops), details=details)
            self._emit('warn', msg)
            self._log(f"[WARN] {msg}")
            self._add_gene_note(gene, self._t('note_stops_masked', count=len(stops), details=details))
            with self._results_lock:
                self.masked_stops[gene] = len(stops)

        gene_results: Dict[str, Dict] = {}
        gene_kappa: Optional[float] = None
        gene_fitted_tree: Optional[str] = None

        models_ordered = (['M0'] + [m for m in cfg['models'] if m != 'M0']
                          if 'M0' in cfg['models'] else list(cfg['models']))

        # optional warm start from a hidden M0 (METHODS.md)
        _needs_warmup = bool(set(models_ordered) & self._SITE_WARMUP) and cfg.get('warm_start_m0', False)
        if 'M0' not in models_ordered and _needs_warmup:
            self._emit('debug', f"    {gene} · M0 (implicit warm-start)…")
            _warmup = self._run_single_analysis(fas_file=fas_file, model_name='M0',
                                                log_file=self._log_path, aln=aln,
                                                save_outputs=False)
            if _warmup and _warmup.get('output_file'):
                _wpath = Path(_warmup['output_file'])
                gene_kappa = self._extract_kappa(_wpath)
                gene_fitted_tree = self._extract_fitted_tree(_wpath)
            self._log(f"[M0-implicit] {gene}: warm-start kappa={gene_kappa}")

        for model_name in models_ordered:
            if stop_event is not None and stop_event.is_set():
                break
            if pause_event is not None:
                pause_event.wait()
            kappa_for_this = None if model_name == 'M0' else gene_kappa
            fitted_for_this = None if model_name == 'M0' else gene_fitted_tree
            self._emit('debug', self._t('model_running', gene=gene, model=model_name))
            with self._results_lock:
                self.running[gene] = (model_name, time.time())
            if fitted_for_this is not None and cfg.get('warm_start_multistart', True):
                result = self._run_model_multistart(fas_file, model_name, self._log_path,
                                                    kappa_for_this, fitted_for_this, aln=aln)
            else:
                result = self._run_single_analysis(fas_file=fas_file, model_name=model_name,
                                                   log_file=self._log_path,
                                                   warm_start_kappa=kappa_for_this,
                                                   fitted_tree=fitted_for_this, aln=aln)
            gene_results[model_name] = result
            if result.get('status') != 'stopped':
                _run_finished(model_name, result)
            if result.get('status') == 'success':
                if result.get('reused'):
                    self._emit('ok', self._t('model_reused', gene=gene, model=model_name,
                                             lnl=result['lnL']))
                else:
                    self._emit('ok', self._t('model_ok', gene=gene, model=model_name,
                                             lnl=result['lnL'], t=result.get('execution_time') or 0))
                if model_name == 'M0' and result.get('output_file'):
                    out_path = Path(result['output_file'])
                    k = self._extract_kappa(out_path)
                    if k is not None:
                        gene_kappa = k
                    ft = self._extract_fitted_tree(out_path)
                    if ft:
                        gene_fitted_tree = ft
            elif result.get('status') == 'stopped':
                break
            else:
                self._emit('error', self._t('model_failed', gene=gene, model=model_name,
                                            reason=result.get('fail_reason', '?')))

        stopped = stop_event is not None and stop_event.is_set()
        failed_models = [m for m, r in gene_results.items() if r.get('status') == 'failed']
        with self._results_lock:
            self.results[gene] = gene_results
        if failed_models:
            reasons = "; ".join(f"{m}: {gene_results[m].get('fail_reason', '?')}" for m in failed_models)
            with self._results_lock:
                self.failures[gene] = reasons
                self.gene_status[gene] = 'failed'
            self._log(f"[FAILED] {gene}: {reasons}")
        elif stopped and len(gene_results) < len(models_ordered):
            with self._results_lock:
                self.gene_status[gene] = 'stopped'
        else:
            with self._results_lock:
                self.gene_status[gene] = 'ok'
        _done()
        return gene, gene_results

    @staticmethod
    def _labeled_root_check(nwk_content: str) -> str:
        """Compare the labels of the two branches around the root of a labelled tree:
        'mixed' (one foreground, one background) needs a rooted tree, 'same' does
        not (PAML guide, Figure S1D)."""
        lines = nwk_content.strip().splitlines()
        nwk = ''
        for line in lines:
            stripped = line.strip()
            if stripped and stripped[0].isdigit() and len(stripped.split()) <= 2:
                continue
            nwk += stripped
        nwk = nwk.strip().rstrip(';').strip()
        if not nwk.startswith('('):
            return 'unknown'

        depth = 0
        top_comma = -1
        for i, ch in enumerate(nwk):
            if ch == '(':
                depth += 1
            elif ch == ')':
                depth -= 1
            elif ch == ',' and depth == 1:
                top_comma = i
                break

        if top_comma == -1:
            return 'unknown'

        child1 = nwk[1:top_comma]
        child2_raw = nwk[top_comma + 1:]
        depth = 0
        end_pos = len(child2_raw)
        for i, ch in enumerate(child2_raw):
            if ch == '(':
                depth += 1
            elif ch == ')':
                if depth == 0:
                    end_pos = i
                    break
                depth -= 1
        child2 = child2_raw[:end_pos]

        c1_fg = '#1' in child1
        c2_fg = '#1' in child2

        return 'mixed' if c1_fg != c2_fg else 'same'

    # purifying, neutral and diversifying starting ω
    _WARM_START_OMEGA_TRIALS = (0.2, 1.0, 2.5)

    def _run_model_multistart(self, fas_file: Path, model_name: str, log_file: Path,
                               warm_start_kappa: float, fitted_tree: str,
                               aln=None) -> Dict:
        """Run a model from the warm-start branch lengths and kappa once per initial ω
        in _WARM_START_OMEGA_TRIALS and keep the best lnL. Every start writes the
        same files in MODEL/, so the files of the best start are put back after a
        worse one."""
        out_dir = Path(self.config['output_folder']) / model_name
        safe = re.sub(r'[^\w.\-]+', '_', fas_file.stem)

        def _files():
            names = set(out_dir.glob(f"{safe}_{model_name}*")) | set(out_dir.glob(f"{fas_file.stem}_{model_name}_*"))
            return {p: p.read_bytes() for p in names if p.is_file()}

        best = None
        best_files = {}
        last = None
        for omega0 in self._WARM_START_OMEGA_TRIALS:
            r = self._run_single_analysis(
                fas_file=fas_file, model_name=model_name, log_file=log_file,
                warm_start_kappa=warm_start_kappa, fitted_tree=fitted_tree,
                omega_override=omega0, aln=aln,
            )
            last = r
            if r.get('status') == 'stopped':
                return r
            if r.get('lnL') is not None and (best is None or r['lnL'] > best['lnL']):
                best = r
                best_files = _files()
        if best is not None and best is not last and best_files:
            for path in set(_files()) - set(best_files):
                path.unlink(missing_ok=True)
            for path, data in best_files.items():
                path.write_bytes(data)
            self._log(f"[{model_name}] {fas_file.stem}: kept the start with the best lnL ({best['lnL']})")
        return best if best is not None else last


    @staticmethod
    def _process_cpu_seconds(pid: int) -> Optional[float]:
        """CPU seconds used by a process and its descendants, or None if unknown.

        On Debian/Ubuntu /usr/bin/codeml is a sh script that runs the real codeml
        as a child, so the child must be counted. Uses psutil, else /proc."""
        try:
            import psutil
            proc = psutil.Process(pid)
            total = 0.0
            for p in [proc] + proc.children(recursive=True):
                try:
                    t = p.cpu_times()
                    total += t.user + t.system
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
            return float(total)
        except ImportError:
            pass
        except Exception:
            return None
        try:
            ticks = os.sysconf(os.sysconf_names['SC_CLK_TCK'])
            stats = {}
            for d in os.listdir('/proc'):
                if not d.isdigit():
                    continue
                try:
                    with open(f"/proc/{d}/stat", 'r') as fh:
                        fields = fh.read().rsplit(')', 1)[1].split()
                except OSError:
                    continue
                stats[int(d)] = (int(fields[1]), int(fields[11]) + int(fields[12]))
            if pid not in stats:
                return None
            total, todo = 0, [pid]
            while todo:
                cur = todo.pop()
                total += stats[cur][1]
                todo.extend(c for c, (ppid, _) in stats.items() if ppid == cur)
            return total / ticks
        except Exception:
            return None

    @staticmethod
    def _terminate(process) -> None:
        """End codeml and reap it. Outside Windows codeml runs in its own process
        group, so the signal also reaches a codeml started by a wrapper script."""
        if process.poll() is not None:
            return
        import signal

        def _signal(sig):
            if platform.system() != 'Windows':
                try:
                    os.killpg(process.pid, sig)
                    return
                except (ProcessLookupError, PermissionError, OSError):
                    pass
            (process.terminate if sig == signal.SIGTERM else process.kill)()
        try:
            _signal(signal.SIGTERM)
            process.wait(timeout=3)
        except Exception:
            try:
                _signal(getattr(signal, 'SIGKILL', signal.SIGTERM))
                process.wait(timeout=3)
            except Exception:
                pass

    def stop_all_processes(self) -> int:
        """End every running codeml (Stop button). Returns how many."""
        with self._processes_lock:
            procs = list(self._active_processes)
        for proc in procs:
            self._terminate(proc)
        return len(procs)

    @staticmethod
    def _prune_keep_labels(tree, keep: set) -> None:
        """Prune tips not in `keep` without losing branch labels: when pruning leaves an
        internal node with one child, Bio.Phylo collapses it and its #N label would
        be lost, so the remaining child inherits the label."""
        label_re = re.compile(r'[#$]\d+$')

        def _label(clade):
            m = label_re.search(clade.name or '')
            return m.group(0) if m else None

        def _add_label(clade, lab):
            if lab and not _label(clade):
                clade.name = (clade.name or '') + lab

        def _walk(clade):
            for child in list(clade.clades):
                if child.is_terminal():
                    if label_re.sub('', child.name or '') not in keep:
                        clade.clades.remove(child)
                else:
                    _walk(child)
                    if not child.clades:
                        clade.clades.remove(child)
                    elif len(child.clades) == 1:
                        only = child.clades[0]
                        _add_label(only, _label(child))
                        idx = clade.clades.index(child)
                        clade.clades[idx] = only
        _walk(tree.root)
        while len(tree.root.clades) == 1 and not tree.root.clades[0].is_terminal():
            tree.root = tree.root.clades[0]


    def _failed(self, reason: str, exec_start: float = None, **extra) -> Dict:
        d = {
            'output_file': None, 'results_file': None, 'lnL': None, 'np': None,
            'ntime': None, 'omega': None,
            'execution_time': (time.time() - exec_start) if exec_start else 0.0,
            'status': 'failed', 'fail_reason': reason, 'stop_count': 0, 'beb_skipped': False,
        }
        d.update(extra)
        return d

    def _execute_codeml(self, codeml_path, ctl_filename, temp_dir, model_name, base_name,
                        names, seqs_used):
        """Run codeml in temp_dir with a time limit, an idle check and Stop.
        Returns (returncode, stdout lines, fail reason, stopped, BEB skipped)."""
        cfg = self.config
        pause_event = cfg.get('pause_event')
        stop_event = cfg.get('stop_event')
        cmd = [codeml_path, ctl_filename]
        self._log(f"[{model_name}] {base_name}: Running command: {cmd} in {temp_dir}")
        popen_kw = {}
        if platform.system() == 'Windows':
            popen_kw['creationflags'] = getattr(subprocess, 'CREATE_NO_WINDOW', 0)
        else:
            popen_kw['start_new_session'] = True   # see _terminate
        # stdin closed: after a stop codon codeml waits for "Press Enter";
        # with stdin closed it continues at once
        process = subprocess.Popen(
            cmd, cwd=temp_dir, stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, encoding='utf-8', errors='replace', bufsize=1, **popen_kw)

        stdout_lines: List[str] = []
        skip_beb = bool(cfg.get('skip_beb', False)) and model_name in ('M2a', 'M8')
        beb_was_skipped = [False]

        def read_stream(stream):
            try:
                for line in iter(stream.readline, ''):
                    text_line = line.rstrip()
                    stdout_lines.append(text_line)
                    # skip_beb: lnL, np and ω are already written when BEB starts
                    if skip_beb and 'BEBing' in text_line and not beb_was_skipped[0]:
                        beb_was_skipped[0] = True
                        self._log(f"[{model_name}] {base_name}: skip_beb -- stopping before BEB")
                        try:
                            self._terminate(process)
                        except Exception:
                            pass
            except Exception:
                pass

        reader = Thread(target=read_stream, args=(process.stdout,), daemon=True)
        reader.start()
        with self._processes_lock:
            self._active_processes.add(process)
            self.current_process = process

        n_codons = max((len(sq) for sq in seqs_used.values()), default=0) // 3
        timeout_s = float(resolve_timeout(cfg.get('timeout'), model_name, len(names), n_codons))
        self._log(f"[{model_name}] {base_name}: time limit {int(timeout_s)} s "
                  f"({len(names)} taxa × {n_codons} codons"
                  f"{', user value' if user_timeout(cfg.get('timeout')) else ', automatic'})")
        idle_s = float(cfg.get('idle_timeout', 300) or 0)
        deadline = time.time() + timeout_s
        last_cpu = None
        last_progress = time.time()
        fail_reason = None
        stopped = False
        try:
            while process.poll() is None:
                if stop_event is not None and stop_event.is_set():
                    stopped = True
                    self._terminate(process)
                    break
                paused = pause_event is not None and not pause_event.is_set()
                now = time.time()
                if paused:
                    # paused time counts neither for the limit nor the idle check
                    deadline += 0.25
                    last_progress = now
                elif now > deadline:
                    fail_reason = self._t('reason_timeout', s=int(timeout_s))
                    self._terminate(process)
                    break
                elif idle_s > 0:
                    cpu = self._process_cpu_seconds(process.pid)
                    if cpu is not None:
                        if last_cpu is None or cpu > last_cpu + 0.01:
                            last_cpu = cpu
                            last_progress = now
                        elif now - last_progress > idle_s:
                            last_line = next((l for l in reversed(stdout_lines) if l.strip()), '')
                            fail_reason = self._t('reason_idle', s=int(idle_s), line=last_line[:200])
                            self._terminate(process)
                            break
                time.sleep(0.25)
        finally:
            with self._processes_lock:
                self._active_processes.discard(process)
                if self.current_process is process:
                    self.current_process = None
        reader.join(timeout=5)
        try:
            process.wait(timeout=5)
        except Exception:
            self._terminate(process)
        try:
            if process.stdout:
                process.stdout.close()
        except Exception:
            pass

        rc = process.returncode
        self._log(f"[{model_name}] {base_name}: process returncode={rc}")
        if stdout_lines:
            self._log("  codeml stdout (last 60 lines):\n" + "\n".join(stdout_lines[-60:]))
        return rc, stdout_lines, fail_reason, stopped, beb_was_skipped[0]

    def _reuse_saved_run(self, model_output_dir: Path, temp_dir: Path, ctl_filename: str,
                         seq_filename: str, tree_filename: str, output_filename: str,
                         codeml_outfile: str, prefix: str) -> bool:
        """True when MODEL/ already holds a finished run with the same .ctl,
        alignment and tree; its output is then copied into temp_dir in place of
        running codeml again."""
        saved_out = model_output_dir / output_filename
        try:
            for fname in (ctl_filename, seq_filename, tree_filename):
                saved = model_output_dir / fname
                if not saved.is_file() or saved.read_bytes() != (temp_dir / fname).read_bytes():
                    return False
            if not saved_out.is_file() or self._extract_model_stats(saved_out)['lnL'] is None:
                return False
            shutil.copy2(saved_out, temp_dir / codeml_outfile)
            rst = model_output_dir / f"{prefix}_rst.txt"
            if rst.is_file():
                shutil.copy2(rst, temp_dir / 'rst')
            return True
        except OSError:
            return False

    def _run_single_analysis(self, fas_file: Path, model_name: str,
                             log_file: Path = None,
                             warm_start_kappa: float = None,
                             fitted_tree: str = None,
                             omega_override: float = None,
                             aln=None,
                             save_outputs: bool = True) -> Dict:
        """Run codeml for one gene and model. Returns a dict with status 'success',
        'failed' or 'stopped' (and fail_reason). Inputs and outputs are kept in
        MODEL/ with relative paths, so `cd MODEL && codeml GENE_MODEL.ctl` repeats it."""
        cfg = self.config
        log_file = log_file or getattr(self, '_log_path', None)
        if not hasattr(self, '_log_lock'):
            self._log_lock = threading.Lock()
        if getattr(self, '_log_path', None) is None and log_file is not None:
            self._log_path = Path(log_file)
        base_name = fas_file.stem
        safe = re.sub(r'[^\w.\-]+', '_', base_name)

        model_config = dict(self.MODEL_CONFIGS.get(model_name, {}))
        custom = cfg.get('custom_model_params')
        if custom is None:
            custom = cfg.get('custom_model_configs', {}) or {}
        for k, v in (custom or {}).get(model_name, {}).items():
            model_config[k] = v

        pause_event = cfg.get('pause_event')
        stop_event = cfg.get('stop_event')
        if pause_event is not None:
            pause_event.wait()
        if stop_event is not None and stop_event.is_set():
            return self._failed(self._t('reason_stopped'), status='stopped')

        codeml_path = getattr(self, '_codeml_path', None) or find_codeml(cfg.get('codeml_path'))
        if not codeml_path:
            return self._failed(self._t('reason_no_codeml'))

        model_output_dir = Path(cfg['output_folder']) / model_name
        if save_outputs:
            model_output_dir.mkdir(parents=True, exist_ok=True)

        prefix = f"{safe}_{model_name}"
        output_filename = f"{base_name}_{model_name}_results.txt"
        codeml_outfile = f"{prefix}_results.txt"
        ctl_filename = f"{prefix}.ctl"
        seq_filename = f"{prefix}_seq.fasta"
        tree_filename = f"{prefix}_tree.nwk"

        temp_dir = Path(tempfile.mkdtemp(prefix=f'easypam_{prefix}_', dir=self._get_fast_tempdir()))
        exec_start = time.time()
        try:
            if aln is None:
                aln = read_alignment(fas_file)
            names = list(aln.names)
            excluded: List[str] = []

            from io import StringIO
            from Bio import Phylo
            tree_path = Path(cfg.get('tree_file')) if cfg.get('tree_file') else None
            per_gene = (cfg.get('per_gene_trees') or {}).get(base_name)
            if per_gene:
                tree_path = Path(per_gene)
            if tree_path is None:
                return self._failed(self._t('reason_exception', error='no tree for this gene'), exec_start)
            tree_text = tree_path.read_text(encoding='utf-8', errors='ignore')
            t_lines = tree_text.splitlines()
            if t_lines and t_lines[0].strip() and t_lines[0].strip().split()[0].isdigit() \
                    and not t_lines[0].strip().startswith('('):
                tree_text = '\n'.join(t_lines[1:])

            tree_obj = None
            not_in_fasta: set = set()
            try:
                tree_obj = Phylo.read(StringIO(tree_text), 'newick')
                tree_taxa = {re.sub(r'[#$]\d+$', '', t.name) for t in tree_obj.get_terminals() if t.name}
            except Exception as exc:
                return self._failed(f"tree: {exc}", exec_start)

            if cfg.get('auto_prune_tree', True):
                not_in_tree = [n for n in names if n not in tree_taxa]
                not_in_fasta = tree_taxa - set(names)
                if not_in_tree:
                    excluded = not_in_tree
                    names = [n for n in names if n in tree_taxa]
                    msg = self._t('warn_excluded', gene=base_name, names=', '.join(not_in_tree))
                    with self._results_lock:
                        warned = getattr(self, '_excluded_warned', set())
                        first = base_name not in warned
                        warned.add(base_name)
                        self._excluded_warned = warned
                    if first:   # once per gene, not once per model
                        self._emit('warn', msg)
                        self._add_gene_note(base_name, self._t('note_excluded',
                                                               names=', '.join(not_in_tree)))
                        with self._results_lock:
                            self.excluded_taxa[base_name] = list(not_in_tree)
                    self._log(f"[WARN] {base_name} [{model_name}]: {msg}")
                for tx in not_in_fasta:
                    for term in [t for t in tree_obj.get_terminals()
                                 if t.name and re.sub(r'[#$]\d+$', '', t.name) == tx]:
                        tree_obj.prune(term)
                if not_in_fasta:
                    self._log(f"[tree] {base_name} [{model_name}]: pruned {sorted(not_in_fasta)}")
            if len(names) < 3:
                return self._failed(self._t('reason_exception',
                                            error=f"{len(names)} sequence(s) shared with the tree"),
                                    exec_start)

            seqs_used = {n: aln.seqs[n] for n in names}
            (temp_dir / seq_filename).write_text(to_fasta(names, seqs_used), encoding='utf-8')

            def _newick(tree) -> str:
                io_ = StringIO()
                # a tree without branch lengths stays without them (Bio.Phylo would write 0.00000)
                has_bl = any(c.branch_length for c in tree.find_clades())
                Phylo.write(tree, io_, 'newick', plain=not has_bl)
                txt = io_.getvalue().strip()
                # codeml rejects the root length Bio.Phylo writes (":0.00000;")
                return re.sub(r'\):[0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?;$', ');', txt)

            if model_name in self._SITE_MODELS_UNROOT and len(tree_obj.root.clades) == 2:
                c0, c1 = tree_obj.root.clades
                if c1.clades:
                    tree_obj.root.clades = [c0] + c1.clades
                elif c0.clades:
                    tree_obj.root.clades = c0.clades + [c1]
                tree_note = "unrooted"
            else:
                tree_note = "as given (pruned if needed)"

            labeled_full = cfg.get('labeled_tree_content')
            labeled_bs = cfg.get('labeled_tree_branchsite')
            if model_name.startswith(('BranchSite', 'Branch-site')):
                labeled_content = labeled_bs or (labeled_full if labeled_full and '#1' in labeled_full else None)
            elif model_name == 'Branch':
                labeled_content = labeled_full
            else:
                labeled_content = None

            fix_bl = 0
            if labeled_content:
                lbl = re.sub(r':\s*-?[0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?', '', labeled_content)
                if not_in_fasta or excluded:
                    try:
                        lbl_tree = Phylo.read(StringIO(lbl), 'newick')
                        self._prune_keep_labels(lbl_tree, set(names))
                        lbl = _newick(lbl_tree)
                    except Exception as exc:
                        self._log(f"[WARN] {base_name} [{model_name}]: labeled tree pruning failed: {exc}")
                tree_out = lbl.strip()
                tree_note = "labeled (branch lengths removed)"
                if model_name.startswith(('Branch-site', 'BranchSite')):
                    status = self._labeled_root_check(tree_out)
                    self._log(f"[tree] {base_name} [{model_name}]: root designation {status}")
            elif fitted_tree and not model_name.startswith('Branch') and model_name != 'M0':
                tree_out = fitted_tree.strip()
                # fix_blength = 1: the M0 lengths are starting values (2 would fix them)
                fix_bl = 1
                tree_note = "M0 fitted tree as starting values (fix_blength = 1)"
            else:
                tree_out = _newick(tree_obj)
            n_tips = len(names)
            (temp_dir / tree_filename).write_text(f"{n_tips}  1\n{tree_out}\n", encoding='utf-8')
            self._log(f"[tree] {base_name} [{model_name}]: {tree_filename} ({tree_note})")

            omega_initial = (float(omega_override) if omega_override is not None
                             else float(cfg.get('omega', model_config.get('omega', 0.5) or 0.5)))
            cleandata_val = int(cfg.get('cleandata', 1))
            provided_ctl = (cfg.get('model_ctl_paths', {}) or {}).get(model_name)
            if provided_ctl and Path(provided_ctl).is_file():
                raw = Path(provided_ctl).read_text(encoding='utf-8')

                def _replace_setting(content: str, key: str, newval: str) -> str:
                    return re.sub(rf'(^\s*{re.escape(key)}\s*=).*?$', rf"\1 {newval}",
                                  content, flags=re.MULTILINE)
                ctl_content = _replace_setting(raw, 'seqfile', seq_filename)
                ctl_content = _replace_setting(ctl_content, 'treefile', tree_filename)
                ctl_content = _replace_setting(ctl_content, 'outfile', codeml_outfile)
            else:
                ctl_content = self.generate_ctl_content(
                    seqfile=seq_filename, treefile=tree_filename, outfile=codeml_outfile,
                    model_config=self.MODEL_CONFIGS.get(model_name, {}),
                    omega=omega_initial, cleandata=cleandata_val,
                    model_name=model_name, kappa=warm_start_kappa, fix_blength=fix_bl)
            (temp_dir / ctl_filename).write_text(ctl_content, encoding='utf-8')

            reused = save_outputs and cfg.get('reuse_results', True) and self._reuse_saved_run(
                model_output_dir, temp_dir, ctl_filename, seq_filename, tree_filename,
                output_filename, codeml_outfile, prefix)
            if reused:
                self._log(f"[{model_name}] {base_name}: same .ctl, alignment and tree as the saved "
                          f"run; its output is reused")
                rc, stdout_lines, fail_reason, stopped, beb_skipped = 0, [], None, False, False
            else:
                rc, stdout_lines, fail_reason, stopped, beb_skipped = self._execute_codeml(
                    codeml_path, ctl_filename, temp_dir, model_name, base_name, names, seqs_used)
            beb_was_skipped = [beb_skipped]

            # Stop may end codeml before the loop above sees stop_event
            if not stopped and stop_event is not None and stop_event.is_set() and rc != 0:
                stopped = True
            if stopped:
                return self._failed(self._t('reason_stopped'), exec_start, status='stopped')

            codeml_out = temp_dir / codeml_outfile
            output_path = model_output_dir / output_filename
            last_line = next((l for l in reversed(stdout_lines) if l.strip()), '')

            # skip_beb stops codeml on purpose: cut incomplete sections
            if beb_was_skipped[0] and codeml_out.exists():
                try:
                    text = codeml_out.read_text(encoding='utf-8', errors='ignore')
                    beb_idx = text.find('Bayes Empirical Bayes (BEB) analysis')
                    if beb_idx != -1:
                        text = text[:beb_idx].rstrip() + '\n'
                    neb_marker = 'Naive Empirical Bayes (NEB) analysis'
                    neb_idx = text.find(neb_marker)
                    if neb_idx != -1 and not re.search(r'^\s*\d+\s+[A-Za-z]\s',
                                                       text[neb_idx + len(neb_marker):], re.MULTILINE):
                        text = text[:neb_idx].rstrip() + '\n'
                    codeml_out.write_text(text, encoding='utf-8')
                except Exception as exc:
                    self._log(f"[{model_name}] {base_name}: truncation after skip_beb failed: {exc}")

            if save_outputs:
                (model_output_dir / ctl_filename).write_text(ctl_content, encoding='utf-8')
                for fname in (seq_filename, tree_filename):
                    shutil.copy2(temp_dir / fname, model_output_dir / fname)
                if codeml_out.exists():
                    shutil.move(str(codeml_out), str(output_path))
                for extra in ('rst',):
                    if (temp_dir / extra).exists():
                        try:
                            shutil.move(str(temp_dir / extra), str(model_output_dir / f"{prefix}_{extra}.txt"))
                        except Exception:
                            pass
            else:
                output_path = temp_dir / codeml_outfile
                # hidden M0: copy the result out of the temporary folder
                if output_path.exists():
                    keep = Path(tempfile.mkdtemp(prefix='easypam_m0_')) / codeml_outfile
                    shutil.copy2(output_path, keep)
                    output_path = keep

            lnL = np_params = ntime_params = omega = None
            if output_path.exists():
                st = self._extract_model_stats(output_path)
                lnL, np_params, ntime_params = st['lnL'], st['np'], st['ntime']
                try:
                    omega = SitesParser.extract_omega_robust(output_path)
                except Exception:
                    omega = None

            if save_outputs and output_path.exists():
                try:
                    kept = (cleandata_kept_codons(names, seqs_used) if cleandata_val == 1
                            else list(range(1, aln.length // 3 + 1)))
                    ls = codeml_site_count(output_path)
                    sm = write_sitemap(model_output_dir / f"{base_name}_{model_name}_sitemap.json",
                                       cleandata=cleandata_val, n_codons=aln.length // 3,
                                       kept_codons=kept, sequences=names, codeml_sites=ls)
                    if sm['verified'] is False:
                        msg = self._t('warn_sitemap', gene=base_name, model=model_name,
                                      codeml=ls, expected=len(kept))
                        self._emit('warn', msg)
                        self._log(f"[WARN] {msg}")
                except Exception as exc:
                    self._log(f"[{model_name}] {base_name}: sitemap failed: {exc}")

            execution_time = time.time() - exec_start
            self._log(f"[{model_name}] {base_name}: FINISHED (lnL={lnL}, np={np_params}, "
                      f"ω={omega}, time={execution_time:.1f}s, rc={rc})")

            if fail_reason is None:
                if not output_path.exists():
                    fail_reason = self._t('reason_no_output', line=last_line[:200])
                elif rc not in (0, None) and not beb_was_skipped[0]:
                    fail_reason = self._t('reason_rc', rc=rc, line=last_line[:200])
                elif lnL is None:
                    fail_reason = self._t('reason_no_lnl', line=last_line[:200])

            # keep a failed run's output under another name so no reader of
            # *_results.txt uses its lnL
            if fail_reason and save_outputs and output_path.exists():
                failed_path = output_path.with_name(
                    output_path.name[:-len('_results.txt')] + FAILED_RESULTS_SUFFIX)
                try:
                    output_path.replace(failed_path)
                    output_path = failed_path
                    lnL = np_params = ntime_params = omega = None
                except OSError as exc:
                    self._log(f"[{model_name}] {base_name}: could not rename failed output: {exc}")

            return {
                'output_file': str(output_path) if output_path.exists() else None,
                'results_file': str(output_path) if output_path.exists() else None,
                'lnL': lnL, 'np': np_params, 'ntime': ntime_params, 'omega': omega,
                'execution_time': execution_time,
                'status': 'failed' if fail_reason else 'success',
                'fail_reason': fail_reason,
                'stop_count': 0,
                'beb_skipped': beb_was_skipped[0],
                'excluded_sequences': excluded,
                'reused': bool(reused),
            }

        except Exception as e:
            self._log(f"[{model_name}] {base_name}: EXCEPTION - {e}\n{traceback.format_exc()}")
            return self._failed(self._t('reason_exception', error=e), exec_start)

        finally:
            for attempt in range(8):
                try:
                    shutil.rmtree(temp_dir)
                    break
                except FileNotFoundError:
                    break
                except Exception:
                    time.sleep(0.5)

    # final result line of codeml, standard form first
    _LNL_LINE_PATTERNS = (
        r'lnL[^:]*:\s*([+-]?\d+\.\d+)',
        r'lnL\([^)]*\):\s*([+-]?\d+\.\d+)',
        r'lnL\s*[:=]\s*([+-]?\d+\.\d+)',
    )

    def _extract_model_stats(self, output_file: Path) -> Dict[str, Optional[float]]:
        """lnL, np and ntime from the final "lnL(ntime: X  np: Y): value" line."""
        lnL = np_params = ntime_params = None
        try:
            with open(output_file, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    if 'lnL' not in line:
                        continue
                    for pat in self._LNL_LINE_PATTERNS:
                        match = re.search(pat, line)
                        if not match:
                            continue
                        try:
                            lnL = float(match.group(1))
                        except Exception:
                            continue
                        if 'np:' in line:
                            m_np = re.search(r'np:\s*(\d+)', line)
                            if m_np:
                                np_params = int(m_np.group(1))
                        if 'ntime:' in line:
                            m_nt = re.search(r'ntime:\s*(\d+)', line)
                            if m_nt:
                                ntime_params = int(m_nt.group(1))
                        break
                    if lnL is not None:
                        break
        except Exception:
            pass
        return {'lnL': lnL, 'np': np_params, 'ntime': ntime_params}


    def _lrt_comparisons_for(self, selected_models):
        """LRT pairs available for the selected models (from lrt_stats.PAIRS)."""
        comparisons = [(n, a, lrt_stats.lrt_column(n, a))
                       for n, a in lrt_stats.pairs_for(selected_models)]
        # legacy BranchSite_A folders
        if 'BranchSite_A_null' in selected_models and 'BranchSite_A' in selected_models:
            comparisons.append(('BranchSite_A_null', 'BranchSite_A',
                                'lrt_BranchSite_A_null_vs_BranchSite_A'))
        return comparisons

    def _save_summary(self):
        """Write analysis_summary.tsv: per gene, status, lnL/np/ntime/ω per model, the
        positive class (M2a/M8) and 2Δl, p and q per test."""
        summary_file = self.config['output_folder'] / "analysis_summary.tsv"
        models = self.config['models']
        comparisons = self._lrt_comparisons_for(models)
        qvalues = getattr(self, '_lrt_qvalues', None) or {}
        pvalues = getattr(self, '_lrt_pvalues', None) or {}
        pos_models = [m for m in models if m in ('M2a', 'M8')]

        header = ["Gene", "status"]
        for model in models:
            header.extend([f"{model}_lnL", f"{model}_np", f"{model}_ntime",
                           f"{model}_omega", f"{model}_time", f"{model}_stops"])
        for m in pos_models:
            header.extend([f"{m}_p_pos", f"{m}_w_pos"])
        for null_model, alt_model, col_name in comparisons:
            header.append(col_name)
            if (null_model, alt_model) in pvalues:
                header.append(lrt_stats.p_column(null_model, alt_model))
            if (null_model, alt_model) in qvalues:
                header.append(lrt_stats.q_column(null_model, alt_model))

        def _fmt(v, f="{:.6f}"):
            return f.format(v) if v is not None else 'NA'

        with open(summary_file, 'w', encoding='utf-8') as fh:
            fh.write("\t".join(header) + "\n")
            for gene_name in sorted(self.results.keys()):
                gene_results = self.results[gene_name] or {}
                row = [str(gene_name).replace('\n', '').replace('\r', '').replace('\t', '_'),
                       getattr(self, 'gene_status', {}).get(gene_name, 'ok')]
                for model in models:
                    r = gene_results.get(model)
                    if r and r.get('status') == 'success':
                        row.extend([_fmt(r.get('lnL')), str(r.get('np', 'NA')),
                                    str(r.get('ntime', 'NA')), _fmt(r.get('omega')),
                                    f"{r.get('execution_time', 0):.2f}",
                                    str(r.get('stop_count', 0))])
                    else:
                        row.extend(['NA'] * 6)
                for m in pos_models:
                    r = gene_results.get(m)
                    pc = (SitesParser.extract_positive_class(Path(r['results_file']))
                          if r and r.get('results_file') else None) or {}
                    row.extend([_fmt(pc.get('p')), _fmt(pc.get('omega'))])
                for null_model, alt_model, _ in comparisons:
                    n, a = gene_results.get(null_model), gene_results.get(alt_model)
                    if n and a and n.get('lnL') is not None and a.get('lnL') is not None:
                        row.append(f"{2 * (a['lnL'] - n['lnL']):.6f}")
                    else:
                        row.append('NA')
                    if (null_model, alt_model) in pvalues:
                        pv = pvalues[(null_model, alt_model)].get(gene_name)
                        row.append(f"{pv:.6e}" if pv is not None else 'NA')
                    if (null_model, alt_model) in qvalues:
                        q = qvalues[(null_model, alt_model)].get(gene_name)
                        row.append(f"{q:.6e}" if q is not None else 'NA')
                row = [str(v).replace('\n', '').replace('\r', '') for v in row]
                fh.write("\t".join(row) + "\n")

        self._emit('debug', f"Summary saved: {summary_file}")

    LRT_METHOD_NOTE = (
        "METHODS NOTE:\n"
        "  Statistic: 2*(lnL_alternative - lnL_null); negative values (the\n"
        "  alternative did not improve) are set to 0, giving p = 1.\n"
        "  df = extra free parameters in the alternative (branch lengths count\n"
        "  equally in both models):\n"
        "    M0  vs M1a  : df = 1\n"
        "    M1a vs M2a  : df = 2  (p2 and omega2)\n"
        "    M7  vs M8   : df = 2  (p1 and omega_s)\n"
        "    M8a vs M8   : df = 1  (omega_s free vs fixed at 1)\n"
        "    M0  vs Branch : df = number of foreground groups (#1, #2, ...)\n"
        "  Null distribution: chi2 with that df. For tests whose null lies on the\n"
        "  boundary (M0 vs M1a, M8a vs M8, branch-site) significance uses chi2(1),\n"
        "  as the PAML manual recommends for branch-site; the 0.5*chi2(0) +\n"
        "  0.5*chi2(1) mixture (Self & Liang 1987) is printed for reference.\n"
        "  p-values use the survival function (chi2.sf), never rounded to zero.\n"
        "  q-value = p corrected by Benjamini-Hochberg within each pair of models\n"
        "  (family = every gene tested in that pair in this run).\n"
        "  M7 vs M8 can reject M7 only because some sites are neutral (omega = 1);\n"
        "  M8a vs M8 does not have this problem (Swanson et al. 2003).\n"
    )

    def _run_lrt_analysis(self):
        """Likelihood ratio tests with Benjamini-Hochberg correction. BH needs every
        p-value of a pair first, so each pair is collected, corrected, then written.
        q and p are kept in self._lrt_qvalues and self._lrt_pvalues for
        _save_summary(), which runs after this."""
        lrt_file = self.config['output_folder'] / "LRT_results.txt"
        self._lrt_qvalues = {}
        self._lrt_pvalues = {}

        with open(lrt_file, 'w', encoding='utf-8') as f:
            f.write("="*80 + "\n")
            f.write("LIKELIHOOD RATIO TEST (LRT) RESULTS\n")
            f.write(f"EasyPAML {__version__}\n")
            f.write("="*80 + "\n\n")
            f.write(self.LRT_METHOD_NOTE + "\n")

            selected_models = self.config['models']
            comparisons = lrt_stats.pairs_for(selected_models)
            descriptions = {(n, a): d for grp in self.LRT_COMPARISONS.values() for n, a, d in grp}

            if not comparisons:
                f.write("No valid model comparisons found.\n")
                f.write("For LRT, you need pairs of nested models.\n")
                self._emit('warn', self._t('lrt_none'))
                return

            for null_model, alt_model in comparisons:
                info = lrt_stats.PAIRS[(null_model, alt_model)]
                boundary = info['boundary']
                f.write("\n" + "="*80 + "\n")
                f.write(f"COMPARISON: {null_model} (null) vs {alt_model} (alternative)\n")
                f.write(f"Description: {descriptions.get((null_model, alt_model), '')}\n")
                f.write("="*80 + "\n\n")

                collected = []
                for gene_name in sorted(self.results.keys()):
                    gene_results = self.results[gene_name]
                    null_res = gene_results.get(null_model)
                    alt_res = gene_results.get(alt_model)
                    if not null_res or not alt_res:
                        continue
                    if 'failed' in (null_res.get('status'), alt_res.get('status')):
                        continue
                    lnL_null, lnL_alt = null_res.get('lnL'), alt_res.get('lnL')
                    if lnL_null is None or lnL_alt is None:
                        continue
                    np_null, np_alt = null_res.get('np'), alt_res.get('np')
                    ntime_null, ntime_alt = null_res.get('ntime'), alt_res.get('ntime')

                    df = info['df']
                    if df is None:  # M0 vs Branch: number of foreground groups
                        df = abs(np_alt - np_null) if (np_alt and np_null) else 1
                        if ntime_null is not None and ntime_alt is not None:
                            df = max(1, df - (ntime_alt - ntime_null))

                    raw_stat = 2 * (lnL_alt - lnL_null)
                    lrt_stat = max(0.0, raw_stat)
                    p_value = lrt_stats.p_value(lrt_stat, df, boundary=boundary)
                    p_mix = lrt_stats.p_value_mixture(lrt_stat) if boundary else None
                    df_display = (f"{df} (chi2(1); mixture shown for reference)"
                                  if boundary else str(df))
                    collected.append({
                        'gene': gene_name, 'lnL_null': lnL_null, 'lnL_alt': lnL_alt,
                        'np_null': np_null, 'np_alt': np_alt, 'lrt_stat': lrt_stat,
                        'raw_stat': raw_stat, 'df_display': df_display,
                        'p_value': p_value, 'p_value_mixture': p_mix,
                    })

                for c, q in zip(collected, lrt_stats.bh_qvalues([c['p_value'] for c in collected])):
                    c['q_value'] = q
                self._lrt_qvalues[(null_model, alt_model)] = {c['gene']: c['q_value'] for c in collected}
                self._lrt_pvalues[(null_model, alt_model)] = {c['gene']: c['p_value'] for c in collected}

                sig_count_05 = sig_count_q05 = 0
                for c in sorted(collected, key=lambda c: (c['p_value'], c['gene'])):
                    sig_count_05 += c['p_value'] < 0.05
                    sig_count_q05 += c['q_value'] < 0.05
                    f.write(f"Gene: {c['gene']}\n")
                    f.write(f"  lnL {null_model}: {c['lnL_null']:.6f} (np={c['np_null']})\n")
                    f.write(f"  lnL {alt_model}: {c['lnL_alt']:.6f} (np={c['np_alt']})\n")
                    f.write(f"  2Δl = {c['lrt_stat']:.6f}")
                    if c['raw_stat'] < 0:
                        f.write(f"  (raw {c['raw_stat']:.6f} < 0: the alternative did not reach "
                                f"the null; consider running it again)")
                    f.write("\n")
                    f.write(f"  df = {c['df_display']}\n")
                    f.write(f"  p-value = {c['p_value']:.6e}\n")
                    if c.get('p_value_mixture') is not None:
                        f.write(f"  p-value (50:50 mixture, reference only, not used for q) = "
                                f"{c['p_value_mixture']:.6e}\n")
                    f.write(f"  q-value (BH) = {c['q_value']:.6e}\n")
                    if c['q_value'] < 0.05:
                        f.write(f"  Result: SIGNIFICANT -- {alt_model} better (q < 0.05, BH-corrected)\n")
                    else:
                        f.write(f"  Result: not significant (q >= 0.05, BH-corrected)\n")
                    f.write("\n" + "-"*60 + "\n\n")

                total_valid = len(collected)
                f.write("\nSUMMARY:\n")
                f.write(f"  Genes tested: {total_valid}\n")
                if total_valid > 0:
                    f.write(f"  Significant at raw p < 0.05 (before correction, reference only): {sig_count_05} ({100*sig_count_05/total_valid:.1f}%)\n")
                    f.write(f"  Significant at BH q < 0.05: {sig_count_q05} ({100*sig_count_q05/total_valid:.1f}%)\n")
                f.write("\n")
                self._emit('info', self._t('lrt_pair_done', null=null_model, alt=alt_model,
                                           n=total_valid, sig=sig_count_q05))

        if ('M7', 'M8') in comparisons and ('M8a', 'M8') not in comparisons:
            self._emit('warn', self._t('warn_m7m8_without_m8a'))
        self._emit('debug', f"LRT results saved: {lrt_file}")


    @staticmethod
    def _fasta_to_phylip_block(fas_path: Path) -> Optional[str]:
        """FASTA file as a PHYLIP block, or None if the file is already PHYLIP."""
        text = fas_path.read_text(encoding='utf-8', errors='ignore').strip()
        lines = text.splitlines()
        if not lines:
            return None

        first = lines[0].strip().split()
        if len(first) == 2 and first[0].isdigit() and first[1].isdigit():
            return text + '\n'

        seqs: Dict[str, List[str]] = {}
        order: List[str] = []
        current = None
        for line in lines:
            if line.startswith('>'):
                current = line[1:].split()[0]
                order.append(current)
                seqs[current] = []
            elif current is not None:
                seqs[current].append(line.strip())

        if not seqs:
            return None

        sequences = {k: ''.join(v) for k, v in seqs.items()}
        n_taxa = len(order)
        lengths = {len(s) for s in sequences.values()}
        if len(lengths) != 1:
            print(f"  [WARN] {fas_path.name}: sequences of different lengths, skipped")
            return None
        n_sites = lengths.pop()

        block_lines = [f" {n_taxa} {n_sites}"]
        for name in order:
            padded = name[:10].ljust(10)
            block_lines.append(f"{padded}  {sequences[name]}")
        return '\n'.join(block_lines) + '\n'

    @staticmethod
    def regenerate_summary_files(results_folder: Path) -> Dict[str, str]:
        """Rebuild LRT_results.txt, analysis_summary.tsv, batch_analysis_log.txt and
        sites_BEB.tsv from an existing results folder. Returns {name: path}."""
        results_folder = Path(results_folder)
        
        if not results_folder.exists():
            raise ValueError(f"Results folder not found: {results_folder}")
        
        generated_files = {}

        try:
            # LRT first: the summary needs its q-values
            print("\n[1/3] Generating LRT_results.txt...")
            lrt_file, qvalues = CodemlBatchAnalysis._regenerate_lrt_results(results_folder)
            if lrt_file:
                generated_files['LRT_results'] = str(lrt_file)
                print(f"  OK: {lrt_file.name}")

            print("\n[2/3] Generating analysis_summary.tsv...")
            summary_file = CodemlBatchAnalysis._regenerate_analysis_summary(results_folder, qvalues)
            if summary_file:
                generated_files['analysis_summary'] = str(summary_file)
                print(f"  OK: {summary_file.name}")

            print("\n[3/3] Generating batch_analysis_log.txt...")
            log_file = CodemlBatchAnalysis._regenerate_batch_log(results_folder)
            if log_file:
                generated_files['batch_analysis_log'] = str(log_file)
                print(f"  OK: {log_file.name}")
            sites_file = CodemlBatchAnalysis.write_sites_table(results_folder)
            if sites_file:
                generated_files['sites_BEB'] = str(sites_file)
            
            print(f"\n[SUCCESS] All files regenerated successfully!")
            return generated_files
        
        except Exception as e:
            print(f"[ERROR] Could not regenerate the summary files: {e}")
            traceback.print_exc()
            return {}
    
    SITES_TABLE = 'sites_BEB.tsv'

    @staticmethod
    def write_sites_table(results_folder: Path) -> Optional[Path]:
        """sites_BEB.tsv: BEB sites with Pr(ω>1) ≥ 0.95 (the ones codeml marks * or **)
        for M2a, M8 and Branch-site, with both numberings; NEB only when there is no
        BEB. Failed genes are left out."""
        from .site_map import attach_original_positions
        results_folder = Path(results_folder)
        status = CodemlBatchAnalysis._read_gene_status(results_folder)
        legacy = {v: k for k, v in CodemlBatchAnalysis._LEGACY_MODEL_NAMES.items()}
        frames = []
        for model in ('M2a', 'M8', 'Branch-site'):
            folder = results_folder / model
            if not folder.is_dir() and model in legacy:
                folder = results_folder / legacy[model]
            if not folder.is_dir():
                continue
            for rf in sorted(folder.glob('*_results.txt')):
                gene = rf.name.split(f'_{folder.name}_results')[0]
                if status.get(gene, ('',))[0] == 'failed':
                    continue
                try:
                    method = 'BEB'
                    df = SitesParser.parse_sites_from_file(rf, method='BEB')
                    if df.empty:
                        method, df = 'NEB', SitesParser.parse_sites_from_file(rf, method='NEB')
                    df = df[df['pr_w_gt_1'] >= 0.95]
                    if df.empty:
                        continue
                    df = attach_original_positions(df, rf)
                except Exception:
                    continue
                df = df.assign(gene=gene, model=model, method=method)
                frames.append(df)
        out = results_folder / CodemlBatchAnalysis.SITES_TABLE
        cols = ['gene', 'model', 'method', 'position_original', 'position', 'amino_acid',
                'pr_w_gt_1', 'significance', 'post_mean', 'post_se']
        if frames:
            table = pd.concat(frames, ignore_index=True)
            table = table[[c for c in cols if c in table.columns]].rename(columns={
                'position_original': 'position_alignment', 'position': 'position_codeml'})
        else:
            table = pd.DataFrame(columns=['gene', 'model', 'method', 'position_alignment',
                                          'position_codeml', 'amino_acid', 'pr_w_gt_1',
                                          'significance', 'post_mean', 'post_se'])
        table.to_csv(out, sep='\t', index=False, float_format='%.3f')
        return out

    @staticmethod
    def _regenerate_analysis_summary(results_folder: Path, qvalues: Optional[dict] = None) -> Optional[Path]:
        """Rebuild analysis_summary.tsv. qvalues ({(null, alt): {gene: q}}) come from
        _regenerate_lrt_results()."""
        results_folder = Path(results_folder)
        qvalues = qvalues or {}
        summary_file = results_folder / "analysis_summary.tsv"

        model_name_mapping = CodemlBatchAnalysis._LEGACY_MODEL_NAMES
        
        models = []
        for item in results_folder.iterdir():
            if item.is_dir() and item.name not in ['reports']:
                model_name = model_name_mapping.get(item.name, item.name)
                models.append(model_name)
        
        models = sorted(set(models))
        
        if not models:
            print("  [WARN] No model folders found")
            return None
        
        data = {}
        
        reverse_mapping = {v: k for k, v in model_name_mapping.items()}
        
        for model in models:
            folder_name = reverse_mapping.get(model, model)
            model_folder = results_folder / folder_name
            if not model_folder.exists():
                continue
            
            for results_file in sorted(model_folder.glob("*_results.txt")):
                gene_name = results_file.name.split(f'_{folder_name}_results')[0]
                
                if gene_name not in data:
                    data[gene_name] = {'Gene': gene_name}
                
                try:
                    with open(results_file, 'r', encoding='utf-8', errors='ignore') as f:
                        content = f.read()
                    
                    lnL_match = re.search(r'lnL\(ntime:.*?\):\s+([-\d.]+)', content)
                    lnL = float(lnL_match.group(1)) if lnL_match else None
                    
                    np_match = re.search(r'lnL\(ntime:\s*(\d+)\s+np:\s*(\d+)\)', content)
                    np_val    = int(np_match.group(2)) if np_match else None
                    ntime_val = int(np_match.group(1)) if np_match else None

                    omega = SitesParser.extract_omega_robust(results_file)
                    
                    time_match = re.search(r'Time used:\s+(\d+):(\d+)', content)
                    exec_time = None
                    if time_match:
                        m = int(time_match.group(1))
                        s = int(time_match.group(2))
                        exec_time = m * 60 + s
                    
                    stop_count = content.count('***')
                    
                    data[gene_name][f'{model}_lnL']   = lnL
                    data[gene_name][f'{model}_np']    = np_val
                    data[gene_name][f'{model}_ntime'] = ntime_val
                    data[gene_name][f'{model}_omega'] = omega
                    data[gene_name][f'{model}_time'] = exec_time
                    data[gene_name][f'{model}_stops'] = stop_count

                    # Branch: "w (dN/dS) for branches:" lists background, #1, #2, ...
                    if model == 'Branch':
                        tag_omegas = SitesParser.extract_omega_by_tags(results_file)
                        for tag, tag_omega in tag_omegas.items():
                            data[gene_name][f'{model}_{tag}_omega'] = tag_omega

                    if 'Branch-site' in model:
                        class_data = SitesParser.extract_branchsite_class_data(results_file)
                        if class_data:
                            data[gene_name][f'{model}_class_data'] = class_data
                
                except Exception as e:
                    print(f"  [WARN] Error processing {gene_name} ({model}): {str(e)}")

        
        gene_status = CodemlBatchAnalysis._read_gene_status(results_folder)
        for gene_name in data:
            row = data[gene_name]
            st, reason = gene_status.get(gene_name, ('', ''))
            if st:
                row['status'] = st
            if st == 'failed':
                continue
            for (null_m, alt_m), info in lrt_stats.PAIRS.items():
                a_l, n_l = row.get(f'{alt_m}_lnL'), row.get(f'{null_m}_lnL')
                if a_l is None or n_l is None or pd.isna(a_l) or pd.isna(n_l):
                    continue
                lrt = 2 * (a_l - n_l)
                row[lrt_stats.lrt_column(null_m, alt_m)] = lrt
                df = info['df'] or 1
                if info['df'] is None:
                    na, nn = row.get(f'{alt_m}_np'), row.get(f'{null_m}_np')
                    ta, tn = row.get(f'{alt_m}_ntime'), row.get(f'{null_m}_ntime')
                    if na and nn:
                        df = abs(na - nn)
                        if ta is not None and tn is not None:
                            df = max(1, df - (ta - tn))
                row[lrt_stats.p_column(null_m, alt_m)] = lrt_stats.p_value(
                    max(0.0, lrt), df, boundary=info['boundary'])
            for m in ('M2a', 'M8'):
                rf = results_folder / m / f"{gene_name}_{m}_results.txt"
                if rf.exists():
                    pc = SitesParser.extract_positive_class(rf) or {}
                    row[f'{m}_p_pos'] = pc.get('p')
                    row[f'{m}_w_pos'] = pc.get('omega')

        for (null_model, alt_model), gene_qvals in qvalues.items():
            col = f'q_{null_model}_vs_{alt_model}'
            for gene_name, q in gene_qvals.items():
                if gene_name in data:
                    data[gene_name][col] = q

        for gene_name in data:
            row = data[gene_name]
            
            if 'Branch-site_class_data' in row and row['Branch-site_class_data']:
                class_data = row['Branch-site_class_data']
                
                for cls in ['0', '1', '2a', '2b']:
                    if cls in class_data:
                        row[f'Branch-site_class{cls}_prop'] = class_data[cls].get('prop')
                        row[f'Branch-site_class{cls}_bg_w'] = class_data[cls].get('bg_w')
                        row[f'Branch-site_class{cls}_fg_w'] = class_data[cls].get('fg_w')
            
            if 'Branch-site_class_data' in row:
                del row['Branch-site_class_data']
            if 'Branch-site_null_class_data' in row:
                del row['Branch-site_null_class_data']
        
        # p and q in scientific notation ('%.6f' would print 4e-23 as 0.000000)
        for gene_name, (st, _reason) in gene_status.items():
            if gene_name not in data:   # failed in every model: no *_results.txt
                data[gene_name] = {'Gene': gene_name, 'status': st}
        df = pd.DataFrame(list(data.values()))
        if 'status' in df.columns:   # same column position as in a normal run
            df.insert(1, 'status', df.pop('status').fillna('ok'))
        for col in df.columns:
            if col.startswith(('p_', 'q_')):
                df[col] = [f"{v:.6e}" if isinstance(v, (int, float)) and pd.notna(v) else 'NA'
                           for v in df[col]]
        df.to_csv(summary_file, sep='\t', index=False, float_format='%.6f')

        CodemlBatchAnalysis._find_orphaned_analyses(results_folder)

        return summary_file
    
    @staticmethod
    def _find_orphaned_analyses(results_folder: Path) -> dict:
        """{gene: [models]} whose .ctl exists but whose codeml never finished (crash or
        interrupted session). Failed runs are not counted."""
        results_folder = Path(results_folder)
        _legacy = CodemlBatchAnalysis._LEGACY_MODEL_NAMES
        _reverse = {v: k for k, v in _legacy.items()}

        orphaned: dict = {}
        gene_status = CodemlBatchAnalysis._read_gene_status(results_folder)
        for item in results_folder.iterdir():
            if not item.is_dir() or item.name in {'reports'}:
                continue
            model = _legacy.get(item.name, item.name)
            folder_name = item.name

            for ctl_file in item.glob(f"*_{folder_name}.ctl"):
                gene_name = ctl_file.stem.replace(f"_{folder_name}", "")
                result_file = item / f"{gene_name}_{folder_name}_results.txt"
                failed_file = item / f"{gene_name}_{folder_name}{FAILED_RESULTS_SUFFIX}"
                if not result_file.exists() and not failed_file.exists() \
                        and gene_status.get(gene_name, ('',))[0] != 'failed':
                    orphaned.setdefault(gene_name, []).append(model)

        if orphaned:
            print(f"  [WARN] {len(orphaned)} gene(s) with .ctl but no result "
                  f"(worker crash / interrupted session):")
            for gene, models in sorted(orphaned.items()):
                print(f"    {gene}: {', '.join(sorted(models))}")
        return orphaned

    @staticmethod
    def _regenerate_batch_log(results_folder: Path) -> Optional[Path]:
        """Rebuild batch_analysis_log.txt."""
        results_folder = Path(results_folder)
        log_file = results_folder / "batch_analysis_log.txt"

        model_name_mapping = CodemlBatchAnalysis._LEGACY_MODEL_NAMES
        reverse_mapping = {v: k for k, v in model_name_mapping.items()}

        with open(log_file, 'w', encoding='utf-8') as f:
            f.write("="*80 + "\n")
            f.write("EASYPAML ANALYSIS LOG (REGENERATED)\n")
            f.write("="*80 + "\n")
            f.write(f"Regenerated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Results folder: {results_folder}\n")
            f.write("="*80 + "\n\n")

            f.write("SUMMARY:\n")
            f.write("-"*80 + "\n")

            models = set()
            genes = set()

            for item in results_folder.iterdir():
                if item.is_dir() and item.name not in ['reports']:
                    model_display = model_name_mapping.get(item.name, item.name)
                    models.add(model_display)

                    for results_file in item.glob("*_results.txt"):
                        # split by the folder name, which may differ from the model name
                        gene = results_file.name.split(f'_{item.name}_results')[0]
                        genes.add(gene)

            f.write(f"Models: {', '.join(sorted(models))}\n")
            f.write(f"Genes: {len(genes)}\n")
            f.write(f"  {', '.join(sorted(genes)[:5])}" + ("..." if len(genes) > 5 else "") + "\n")
            f.write("\n")

            f.write("DETAILS:\n")
            f.write("-"*80 + "\n\n")

            for gene in sorted(genes):
                f.write(f"Gene: {gene}\n")
                f.write("-"*40 + "\n")

                for model_display in sorted(models):
                    # paths use the folder name, which may differ from the model name
                    folder_name = reverse_mapping.get(model_display, model_display)
                    model_folder = results_folder / folder_name
                    results_file = model_folder / f"{gene}_{folder_name}_results.txt"

                    if results_file.exists():
                        try:
                            with open(results_file, 'r', encoding='utf-8', errors='ignore') as rf:
                                content = rf.read()

                            lnL_match = re.search(r'lnL\(.*?\):\s+([-\d.]+)', content)
                            np_match = re.search(r'np:\s*(\d+)\)', content)

                            lnL = float(lnL_match.group(1)) if lnL_match else "NA"
                            np_val = np_match.group(1) if np_match else "NA"

                            f.write(f"  {model_display:20s} | lnL = {lnL:>12} | np = {np_val:>2}\n")
                        except Exception:
                            f.write(f"  {model_display:20s} | could not read the file\n")
                    else:
                        f.write(f"  {model_display:20s} | not found\n")

                f.write("\n")

            f.write("="*80 + "\n")
            f.write("END OF LOG\n")
            f.write("="*80 + "\n")

        return log_file
    
    @staticmethod
    def _read_gene_status(results_folder: Path) -> Dict[str, Tuple[str, str]]:
        """{gene: (status, reason)} from genes_status.tsv, or {}."""
        path = Path(results_folder) / 'genes_status.tsv'
        out: Dict[str, Tuple[str, str]] = {}
        try:
            for line in path.read_text(encoding='utf-8').splitlines()[1:]:
                parts = line.split('\t')
                if parts and parts[0]:
                    out[parts[0]] = (parts[1] if len(parts) > 1 else '',
                                     parts[2] if len(parts) > 2 else '')
        except OSError:
            pass
        return out

    @staticmethod
    def _regenerate_lrt_results(results_folder: Path) -> tuple:
        """Rebuild LRT_results.txt with BH correction per pair. Returns (lrt_file,
        {(null, alt): {gene: q}})."""
        results_folder = Path(results_folder)
        lrt_file = results_folder / "LRT_results.txt"
        qvalues: dict = {}
        gene_status = CodemlBatchAnalysis._read_gene_status(results_folder)

        model_name_mapping = CodemlBatchAnalysis._LEGACY_MODEL_NAMES
        reverse_mapping = {v: k for k, v in model_name_mapping.items()}
        
        models = set()
        genes = set()
        
        for item in results_folder.iterdir():
            if item.is_dir() and item.name not in ['reports']:
                model = model_name_mapping.get(item.name, item.name)
                models.add(model)
                folder_name = item.name
                for results_file in item.glob("*_results.txt"):
                    gene = results_file.name.split(f'_{folder_name}_results')[0]
                    genes.add(gene)
        
        _desc = {(n, a): d for grp in CodemlBatchAnalysis.LRT_COMPARISONS.values() for n, a, d in grp}
        comparisons = [(n, a, _desc.get((n, a), ''), info['df'] if info['df'] is not None else 1)
                       for (n, a), info in lrt_stats.PAIRS.items() if n in models and a in models]

        with open(lrt_file, 'w', encoding='utf-8') as f:
            f.write("="*80 + "\n")
            f.write("LIKELIHOOD RATIO TEST (LRT) RESULTS (REGENERATED)\n")
            f.write(f"EasyPAML {__version__}\n")
            f.write("="*80 + "\n\n")
            f.write(CodemlBatchAnalysis.LRT_METHOD_NOTE + "\n")
            
            if not comparisons:
                f.write("No valid comparisons found\n")
                return lrt_file, qvalues

            for null_model, alt_model, description, df in comparisons:
                f.write("\n" + "="*80 + "\n")
                f.write(f"COMPARISON: {null_model} (null) vs {alt_model} (alternative)\n")
                f.write(f"Description: {description}\n")
                f.write("="*80 + "\n\n")

                # collect every gene before the BH correction
                collected = []
                for gene in sorted(genes):
                    if gene_status.get(gene, ('',))[0] == 'failed':
                        continue   # older folders: output of a stopped codeml
                    null_folder = reverse_mapping.get(null_model, null_model)
                    alt_folder = reverse_mapping.get(alt_model, alt_model)
                    
                    null_file = results_folder / null_folder / f"{gene}_{null_folder}_results.txt"
                    alt_file = results_folder / alt_folder / f"{gene}_{alt_folder}_results.txt"
                    
                    if not (null_file.exists() and alt_file.exists()):
                        continue
                    
                    try:
                        with open(null_file, 'r', encoding='utf-8', errors='ignore') as nf:
                            null_content = nf.read()
                        with open(alt_file, 'r', encoding='utf-8', errors='ignore') as af:
                            alt_content = af.read()
                        
                        null_lnL_match = re.search(r'lnL\(.*?\):\s+([-\d.]+)', null_content)
                        alt_lnL_match = re.search(r'lnL\(.*?\):\s+([-\d.]+)', alt_content)
                        
                        if not (null_lnL_match and alt_lnL_match):
                            continue
                        
                        lnL_null = float(null_lnL_match.group(1))
                        lnL_alt = float(alt_lnL_match.group(1))
                        
                        lrt_stat = 2 * (lnL_alt - lnL_null)

                        # M0 vs Branch: df = (np_Branch − np_M0) − (ntime_Branch − ntime_M0), because
                        # M0 uses the unrooted tree (2n-3 branches) and Branch the rooted one (2n-2)
                        gene_df = df
                        if null_model == 'M0' and alt_model == 'Branch':
                            np_null_m    = re.search(r'lnL\(ntime:\s*(\d+)\s+np:\s*(\d+)\)', null_content)
                            np_alt_m     = re.search(r'lnL\(ntime:\s*(\d+)\s+np:\s*(\d+)\)', alt_content)
                            if np_null_m and np_alt_m:
                                raw_df       = abs(int(np_alt_m.group(2)) - int(np_null_m.group(2)))
                                ntime_null_v = int(np_null_m.group(1))
                                ntime_alt_v  = int(np_alt_m.group(1))
                                gene_df      = max(1, raw_df - (ntime_alt_v - ntime_null_v))
                            elif gene_df == 0:
                                gene_df = 1

                        # boundary tests use χ²₁, as the PAML manual recommends for branch-site
                        is_boundary = lrt_stats.PAIRS[(null_model, alt_model)]['boundary']
                        lrt_stat = max(0.0, lrt_stat)
                        p_value = lrt_stats.p_value(lrt_stat, gene_df, boundary=is_boundary)
                        p_value_mixture = lrt_stats.p_value_mixture(lrt_stat) if is_boundary else None

                        collected.append({
                            'gene': gene, 'lnL_null': lnL_null, 'lnL_alt': lnL_alt,
                            'lrt_stat': lrt_stat, 'df_display':
                                f"{gene_df} (chi2(1); mixture shown for reference)" if is_boundary else str(gene_df),
                            'p_value': p_value, 'p_value_mixture': p_value_mixture,
                        })

                    except Exception:
                        continue

                if collected:
                    qvals = stats.false_discovery_control([c['p_value'] for c in collected], method='bh')
                    for c, q in zip(collected, qvals):
                        c['q_value'] = q
                qvalues[(null_model, alt_model)] = {c['gene']: c['q_value'] for c in collected}

                sig_count_05 = sig_count_q05 = 0
                for c in sorted(collected, key=lambda c: (c['p_value'], c['gene'])):
                    if c['p_value'] < 0.05:
                        sig_count_05 += 1
                    if c['q_value'] < 0.05:
                        sig_count_q05 += 1

                    f.write(f"Gene: {c['gene']}\n")
                    f.write(f"  lnL {null_model}: {c['lnL_null']:.6f}\n")
                    f.write(f"  lnL {alt_model}: {c['lnL_alt']:.6f}\n")
                    f.write(f"  2Δl = {c['lrt_stat']:.6f}\n")
                    f.write(f"  df = {c['df_display']}\n")
                    f.write(f"  p-value = {c['p_value']:.6e}\n")
                    if c['p_value_mixture'] is not None:
                        f.write(f"  p-value (50:50 mixture, reference only, not used for q) = "
                                f"{c['p_value_mixture']:.6e}\n")
                    f.write(f"  q-value (BH) = {c['q_value']:.6e}\n")

                    if c['q_value'] < 0.05:
                        f.write(f"  Result: SIGNIFICANT -- {alt_model} better (q < 0.05, BH-corrected)\n")
                    else:
                        f.write(f"  Result: not significant (q >= 0.05, BH-corrected)\n")

                    f.write("\n" + "-"*60 + "\n\n")

                total_valid = len(collected)
                if total_valid > 0:
                    f.write("\nSUMMARY:\n")
                    f.write(f"  Genes tested: {total_valid}\n")
                    f.write(f"  Significant at raw p < 0.05 (before correction, reference only): {sig_count_05} ({100*sig_count_05/total_valid:.1f}%)\n")
                    f.write(f"  Significant at BH q < 0.05: {sig_count_q05} ({100*sig_count_q05/total_valid:.1f}%)\n")
                    f.write("\n")

        return lrt_file, qvalues
