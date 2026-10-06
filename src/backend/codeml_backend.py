"""
CODEML Interactive Batch Analysis System
Sistema interativo para executar análises CODEML em batch
Gera automaticamente os arquivos .ctl necessários
Requer apenas arquivos .fas e .tree
"""

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
from .timeouts import resolve_timeout, user_timeout

# saída parcial de uma execução que falhou (não casa com *_results.txt)
FAILED_RESULTS_SUFFIX = '_results_FAILED.txt'
from .version import __version__, source_commit, version_string

# Absolute path to the bundled codeml binary — works regardless of CWD.
_APP_ROOT   = Path(__file__).resolve().parent.parent.parent
_CODEML_BIN = (_APP_ROOT / 'bin' / 'codeml.exe'
               if platform.system() == 'Windows'
               else _APP_ROOT / 'bin' / 'codeml')


def find_codeml(explicit: Optional[str] = None) -> Optional[str]:
    """Caminho do codeml a usar, em ordem: argumento/config, variável de
    ambiente EASYPAML_CODEML, bin/codeml(.exe) do projeto, codeml no PATH."""
    for cand in (explicit, os.environ.get('EASYPAML_CODEML')):
        if cand and Path(cand).exists():
            # absolute(), NÃO resolve(): no Debian/Ubuntu /usr/bin/codeml é um link
            # para um script único que decide o programa pelo nome com que foi
            # chamado -- resolvendo o link, rodaria o baseml.
            return str(Path(cand).absolute())
    if _CODEML_BIN.exists():
        return str(_CODEML_BIN)
    return shutil.which('codeml')


_CODEML_VERSION_CACHE: Dict[str, Optional[str]] = {}


def codeml_version(codeml_path: Optional[str]) -> Optional[str]:
    """'4.9j' / '4.10.10' -- o codeml só imprime a versão quando roda uma
    análise, então roda uma minúscula (3 sequências, 3 códons, M0) numa
    pasta temporária. Leva milissegundos."""
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


class CodemlBatchAnalysis:
    """Sistema completo para análises CODEML em batch"""
    
    # Templates de configuração para diferentes modelos
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
    
    # Informações detalhadas dos modelos (para exibição no botão "?")
    MODEL_INFO = {
        'M0': {
            'full_name': 'One-Ratio Model',
            'test_type': 'Site Model',
            'parameters': 'ω = dN/dS (constant across all sites)',
            'purpose': 'Null hypothesis: estimates a single dN/dS ratio for all sites. Used as baseline for M1a.',
            'interpretation': 'If M1a rejects M0, suggests variation in selection pressure among codon sites.',
            'use_case': 'Always recommended as baseline comparison.',
            'references': 'Goldman & Yang (1994)'
        },
        'M1a': {
            'full_name': 'Nearly Neutral Model',
            'test_type': 'Site Model',
            'parameters': 'Two classes: ω₀ < 1 (purifying), ω₁ = 1 (neutral)',
            'purpose': 'Null hypothesis: allows sites under purifying and neutral selection only.',
            'interpretation': 'If M2a rejects M1a, indicates presence of positive selection (ω > 1).',
            'use_case': 'Compare against M2a to test for positive selection.',
            'references': 'Wong et al. (2004), Swanson et al. (2003)'
        },
        'M2a': {
            'full_name': 'Positive Selection Model',
            'test_type': 'Site Model',
            'parameters': 'Three classes: ω₀ < 1, ω₁ = 1, ω₂ > 1 (positive selection)',
            'purpose': 'Alternative hypothesis: allows positive selection at specific sites.',
            'interpretation': 'Reject M1a at p < 0.05 = evidence for positive selection. Sites with ω₂ > 1 are under positive selection.',
            'use_case': 'Compare against M1a to identify sites under positive selection.',
            'references': 'Nielsen & Yang (1998)'
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
            'interpretation': 'Reject M7 at p < 0.05 = evidence for positive selection. More flexible than M2a.',
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
            'interpretation': 'Reject M0 at p < 0.05 = foreground branch has different ω than background.',
            'use_case': 'Use with "Marcar Branch" to mark specific branches for comparison.',
            'references': 'Reis et al. (2009)'
        },
        'Branch-site': {
            'full_name': 'Branch-site Model',
            'test_type': 'Branch-site Model',
            'parameters': 'ω varies both by site AND by branch (foreground has different classes)',
            'purpose': 'Tests for positive selection affecting specific sites in specific branches.',
            'interpretation': 'Reject Branch-site_null at p < 0.05 = evidence for positive selection on foreground branch.',
            'use_case': 'Most powerful test when ω varies both spatially (codon sites) and temporally (lineages).',
            'references': 'Zhang et al. (2005), Bielawski & Yang (2004)'
        }
    }
    
    # Comparações LRT comuns
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
    
    # Mapeamento de modelos alternativos -> modelos nulos (auto-seleção)
    NULL_MODEL_PAIRS = {
        'M2a': 'M1a',              # M2a (alternativo) -> M1a (nulo)
        'M8': ['M7', 'M8a'],       # M8 (alternativo) -> M7 e M8a (nulos)
        'Branch': 'M0',            # Branch -> M0 (nulo)
        'Branch-site': 'Branch-site_null'  # Branch-site -> Branch-site_null
    }
    
    # Modelos neutros: apenas Branch-site_null requer fix_omega=1 no .ctl
    # (ω₂=1 fixado conforme Yang et al. 2005, Zhang et al. 2005)
    # M0 e M1a NÃO precisam de fix_omega=1:
    #   - M0 estima ω livremente (null para Branch model)
    #   - M1a (NSsites=1): ω₁=1 é restringido INTERNAMENTE pelo CODEML via NSsites=1;
    #     fix_omega=1 no .ctl fixaria TODOS os ω=1, corrompendo o modelo
    NEUTRAL_MODELS = {
        'M8a': {
            'fix_omega': 1,
            'omega': 1.0,
            'corresponding_alternative': 'M8',
            'reason': 'M8a: classe extra com ω = 1 fixado (Swanson et al. 2003)'
        },
        'Branch-site_null': {
            'fix_omega': 1,
            'omega': 1.0,
            'corresponding_alternative': 'Branch-site',
            'reason': 'Branch-site_null: ω₂=1 fixado (Yang et al. 2005, Zhang et al. 2005)'
        }
    }

    # LRT do Branch-site usa distribuição 50:50 de χ²₀ + χ²₁ (não χ² padrão)
    # Valor crítico a α=0.05: 2.706 (qchisq(0.90, df=1))
    # Valor crítico a α=0.01: 5.412 (qchisq(0.98, df=1))
    BRANCHSITE_MIXTURE_CRITICAL = {0.05: 2.706, 0.01: 5.412}

    # Mapeamento de nomes legados (pastas antigas) → nomes de exibição atuais.
    # Centralizado aqui para evitar repetição em _regenerate_analysis_summary,
    # _regenerate_batch_log e _regenerate_lrt_results.
    _LEGACY_MODEL_NAMES: Dict[str, str] = {
        'BranchSite_A':      'Branch-site',
        'BranchSite_A_null': 'Branch-site_null',
    }
    
    @staticmethod
    def available_cores() -> int:
        """Retorna o número de cores lógicos disponíveis no sistema."""
        return os.cpu_count() or 1

    @staticmethod
    def _get_fast_tempdir() -> str:
        """
        Retorna o diretório temporário mais rápido disponível na plataforma.

        Linux: verifica /dev/shm (tmpfs — filesystem em RAM).  Se existir e for
        gravável, usa-o para que os arquivos intermediários do CODEML (rub, rst,
        2base.t, etc.) nunca toquem o disco, eliminando latência de I/O.

        Windows / Mac / outros: fallback para tempfile.gettempdir(), que em
        instalações modernas geralmente aponta para um SSD NVMe do sistema.

        Nota de segurança: /dev/shm costuma ser limitado a 50 % da RAM, mas os
        arquivos temporários de cada run são pequenos (< 5 MB por gene×modelo) e
        são removidos imediatamente após a execução, então o uso simultâneo máximo
        é de aproximadamente (n_workers × 5 MB) — muitíssimo abaixo do limite.
        """
        shm = Path('/dev/shm')
        if shm.exists() and shm.is_dir():
            try:
                # Verificação de escrita real antes de comprometer
                probe = shm / f'.easypam_probe_{os.getpid()}'
                probe.write_bytes(b'\x00')
                probe.unlink()
                return str(shm)
            except OSError:
                pass          # /dev/shm cheio ou sem permissão → fallback
        return tempfile.gettempdir()

    def __init__(self):
        self.results = {}
        self.config = {}
        # current stop codon count updated during runs (for GUI polling)
        self.current_stop_count = 0
        self.current_stop_details = []
        self.current_total_genes = 0
        self.current_processed_genes = 0
        self._results_lock = threading.Lock()
        # Conjunto thread-safe de processos CODEML ativos; permite stop imediato
        self._active_processes: set = set()
        self._processes_lock = threading.Lock()
        self.current_process = None   # compat. GUI (último processo ativo)
    
    @staticmethod
    def auto_complete_null_models(selected_models: List[str], include_neutral: bool = True,
                                  include_m8a: bool = True) -> List[str]:
        """
        Auto-completa modelos nulos baseado em modelos alternativos selecionados.
        
        Quando um modelo alternativo é selecionado, seu correspondente modelo nulo
        é automaticamente adicionado para permitir comparação LRT. Se include_neutral
        é True, modelos neutros especiais também são incluídos automaticamente.
        
        Args:
            selected_models: Lista de modelos selecionados pelo usuário
            include_neutral: Se True, incluir modelos neutros (M1a, Branch-site_null)
            
        Returns:
            Lista de modelos com os nulos auto-adicionados
        """
        completed_models = set(selected_models)
        
        for model in selected_models:
            if model in CodemlBatchAnalysis.NULL_MODEL_PAIRS:
                nulls = CodemlBatchAnalysis.NULL_MODEL_PAIRS[model]
                if isinstance(nulls, str):
                    nulls = [nulls]
                if not include_m8a:   # opção "não incluir M8a" (o M8a escolhido à mão fica)
                    nulls = [m for m in nulls if m != 'M8a']
                completed_models.update(nulls)
        
        # Se include_neutral está habilitado, adicionar modelos neutros se seus
        # correspondentes alternativos foram selecionados
        if include_neutral:
            # M1a é o neutro para M2a
            if 'M2a' in selected_models and 'M1a' not in completed_models:
                completed_models.add('M1a')
            # Branch-site_null é o neutro para Branch-site
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
        """Gera o .ctl com TODOS os parâmetros relevantes escritos explicitamente.

        Nada fica no default interno do codeml: o mesmo .ctl dá o mesmo
        resultado no PAML 4.9j e no 4.10.x (o ncatG padrão, por exemplo,
        mudou entre versões). Valores base em ctl_params.DEFAULT_CTL_PARAMS;
        model_config (modelo + edições do usuário na janela "cfg") e
        self.config (opções globais: CodonFreq, ncatG, kappa) sobrescrevem.

        - kappa       : κ inicial; se vier do warm-start M0 substitui o padrão
        - fix_blength : 0 = estimar do zero; 1 = usar a árvore como ponto de partida
        Modelos neutros (M8a, Branch-site_null): fix_omega = 1, omega = 1.0 sempre.
        """
        params = dict(DEFAULT_CTL_PARAMS)
        cfg = self.config or {}
        skip = ('description', 'display_name')
        # 1) padrão do modelo (MODEL_CONFIGS)
        for key, val in (model_config or {}).items():
            if key not in skip and val not in (None, ''):
                params[key] = val
        # 2) opções globais (GUI: Configurações; CLI: --codonfreq/--ncatg/--kappa)
        for key in ('CodonFreq', 'ncatG', 'kappa', 'fix_kappa', 'icode', 'method',
                    'Small_Diff', 'getSE', 'estFreq'):
            if cfg.get(key) is not None:
                params[key] = cfg[key]
        # 3) edições do usuário na janela "cfg" deste modelo
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
            # Primary format in M0 output: "  kappa (ts/tv) =  2.54321"
            m = re.search(r'kappa\s*\(ts/tv\)\s*=\s*([\d.]+)', text, re.IGNORECASE)
            if m:
                v = float(m.group(1))
                if 0.1 <= v <= 20:
                    return v
            # Secondary format (parameter table): "  kappa   2.54321"
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
        """
        Extrai a árvore com branch lengths otimizados do arquivo de saída do CODEML.

        O CODEML escreve, perto do final do output, a topologia com os comprimentos
        de ramo estimados por ML em formato Newick.  Essa árvore é usada nos modelos
        de sítio (opção warm_start_m0) como ponto de partida via fix_blength = 1
        ("initial" no pamlDOC: os valores da árvore são só o início da otimização
        e continuam sendo estimados).

        Atenção: fix_blength = 2 é "fixed" -- prende os comprimentos nos valores da
        árvore e muda os resultados. Versões até 0.2.0 usavam 2 aqui por engano.

        Estratégia de extração:
        Varredura reversa das linhas do arquivo (a árvore ajustada aparece após os
        parâmetros ML, próxima ao final).  Critérios:
          – começa com '(' e termina com ';'   (formato Newick)
          – contém ':'                          (branch lengths presentes)
          – contém pelo menos um dígito após ':'(exclui topologias sem comprimentos)

        Returns:
            String Newick com branch lengths, ou None se a extração falhar.
        """
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

    # ══════════════════════════════════════════════════════════════════
    # Log / mensagens
    # ══════════════════════════════════════════════════════════════════

    # Modelos de sítio: árvore desenraizada (o PAML exige sem relógio/marcas)
    _SITE_MODELS_UNROOT = {'M0', 'M1a', 'M2a', 'M7', 'M8', 'M8a'}
    # Modelos que aceitam warm-start do M0
    _SITE_WARMUP = {'M1a', 'M2a', 'M7', 'M8', 'M8a'}

    @staticmethod
    def _t(key: str, **kw) -> str:
        from .messages import t
        return t(key, **kw)

    def _emit(self, level: str, text: str) -> None:
        """Uma mensagem por linha para o usuário.

        level: 'info' | 'ok' | 'warn' | 'error' | 'debug' | 'header'.
        Com config['log_callback'] (GUI) a mensagem vai para lá, já
        classificada; sem callback (CLI) vai para o stdout -- 'debug' só com
        config['verbose']."""
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
        """Linha no batch_analysis_log.txt (thread-safe)."""
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

    # ══════════════════════════════════════════════════════════════════
    # Execução em lote
    # ══════════════════════════════════════════════════════════════════

    def effective_ctl_defaults(self) -> Dict[str, object]:
        """Parâmetros globais do .ctl efetivamente usados nesta execução."""
        params = dict(DEFAULT_CTL_PARAMS)
        for key in ('CodonFreq', 'ncatG', 'kappa', 'fix_kappa', 'icode', 'method',
                    'Small_Diff', 'getSE', 'estFreq'):
            if (self.config or {}).get(key) is not None:
                params[key] = self.config[key]
        params['cleandata'] = int((self.config or {}).get('cleandata', 1))
        return params

    def _write_run_config(self, codeml_path: Optional[str], codeml_ver: Optional[str],
                          genes: List[str]) -> None:
        """run_config.json -- tudo o que é preciso para descrever/reproduzir a
        execução: versões, parâmetros do .ctl por modelo e opções do wrapper.
        Gravado pela GUI e pelo CLI."""
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
        """methods_text.txt: parágrafo de Métodos com o que esta execução fez."""
        from .methods_text import build_methods_text
        cfg = self.config
        try:
            sizes = {pair: len(q) for pair, q in (getattr(self, '_lrt_qvalues', None) or {}).items()}
            text = build_methods_text(
                version=version_string(), codeml_version=codeml_ver, models=list(cfg['models']),
                ctl=self.effective_ctl_defaults(), omega0=float(cfg.get('omega', 0.5) or 0.5),
                pruned=bool(cfg.get('auto_prune_tree', True)), family_sizes=sizes,
                n_genes=n_genes, beb=not cfg.get('skip_beb'))
            out = Path(cfg['output_folder']) / 'methods_text.txt'
            out.write_text("Suggested Methods text (generated by EasyPAML from this run; "
                           "review before use)\n\n" + text + "\n", encoding='utf-8')
        except Exception as exc:
            self._log(f"[WARN] methods_text.txt: {exc}")

    def run_batch_analysis(self):
        """Executa a análise em lote. Retorna self.run_summary:
        {'total', 'ok', 'failed', 'stopped', 'failures': {gene: motivo}, ...}."""
        if not self.config:
            raise ValueError(
                "self.config vazio -- defina input_folder/tree_file/output_folder/models "
                "antes de chamar run_batch_analysis() (ver easypaml_cli.py ou a GUI)."
            )

        cfg = self.config
        cfg['input_folder'] = Path(cfg['input_folder'])
        cfg['output_folder'] = Path(cfg['output_folder'])
        output_folder = cfg['output_folder']
        output_folder.mkdir(parents=True, exist_ok=True)   # cria se não existir
        log_file = output_folder / "batch_analysis_log.txt"
        self._log_path = log_file
        self._log_lock = threading.Lock()
        self.failures: Dict[str, str] = {}
        self.gene_status: Dict[str, str] = {}
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
        # Árvore por gene (GENE.nwk ao lado do alinhamento): substitui tree_file
        if cfg.get('per_gene_trees') is None and cfg.get('auto_per_gene_trees', True):
            cfg['per_gene_trees'] = discover_per_gene_trees(
                cfg['input_folder'], cfg.get('tree_folder'), genes=set(chosen))
        self.current_total_genes = len(genes)
        self.current_processed_genes = 0
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
                                   ncatg=ctl_defaults['ncatG'], cleandata=ctl_defaults['cleandata']))
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
            if n_workers > 1:
                with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
                    list(executor.map(lambda a: self._process_gene(*a, len(genes)), indexed))
            else:
                for item in indexed:
                    self._process_gene(*item, len(genes))

        total_time = time.time() - start_time
        stopped = bool(cfg.get('stop_event') is not None and cfg['stop_event'].is_set())

        # LRT primeiro -- popula q/p, que _save_summary() anexa ao TSV.
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

    def _write_failures_file(self) -> None:
        """genes_status.tsv: uma linha por gene, 'ok' ou o motivo da falha."""
        path = Path(self.config['output_folder']) / 'genes_status.tsv'
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write("Gene\tstatus\treason\n")
            for gene in sorted(self.gene_status):
                reason = self.failures.get(gene, '').replace('\t', ' ').replace('\n', ' ')
                fh.write(f"{gene}\t{self.gene_status[gene]}\t{reason}\n")

    def _mark_gene_failed(self, gene: str, reason: str) -> None:
        with self._results_lock:
            self.failures[gene] = reason
            self.gene_status[gene] = 'failed'
            self.results.setdefault(gene, {})
        self._emit('error', self._t('gene_failed', gene=gene, reason=reason))
        self._log(f"[FAILED] {gene}: {reason}")

    def _process_gene(self, idx: int, gene: str, fas_file: Path, n_total: int):
        cfg = self.config
        pause_event = cfg.get('pause_event')
        stop_event = cfg.get('stop_event')

        def _done():
            with self._results_lock:
                self.current_processed_genes += 1
                done = self.current_processed_genes
            self._progress(done, n_total, gene)

        if stop_event is not None and stop_event.is_set():
            return gene, {}
        if pause_event is not None:
            pause_event.wait()

        self._emit('info', self._t('gene_start', i=idx, n=n_total, gene=gene))

        # ── Validação do alinhamento (mesmas regras do preflight) ──────────
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

        # ── Stop codons: decididos ANTES de rodar ─────────────────────────
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

        gene_results: Dict[str, Dict] = {}
        gene_kappa: Optional[float] = None
        gene_fitted_tree: Optional[str] = None

        models_ordered = (['M0'] + [m for m in cfg['models'] if m != 'M0']
                          if 'M0' in cfg['models'] else list(cfg['models']))

        # Warm-start opcional via M0 implícito (ver AGENTS.md / METODOS.md)
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
            if fitted_for_this is not None and cfg.get('warm_start_multistart', True):
                result = self._run_model_multistart(fas_file, model_name, self._log_path,
                                                    kappa_for_this, fitted_for_this, aln=aln)
            else:
                result = self._run_single_analysis(fas_file=fas_file, model_name=model_name,
                                                   log_file=self._log_path,
                                                   warm_start_kappa=kappa_for_this,
                                                   fitted_tree=fitted_for_this, aln=aln)
            gene_results[model_name] = result
            if result.get('status') == 'success':
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
        """
        Verifica se os dois ramos ao redor da raiz de uma árvore labelada têm a
        mesma designação (ambos foreground #1 ou ambos background) ou designações
        diferentes (um foreground, um background).

        Conforme o guia do PAML (Figura S1D):
        - 'mixed' → ramos com designações diferentes → árvore ENRAIZADA necessária
        - 'same'  → ramos com mesma designação     → árvore não-enraizada pode ser usada
        - 'unknown' → não foi possível determinar (tree com < 2 filhos no root, etc.)
        """
        # Remover cabeçalho PHYLIP ("N  1") se presente
        lines = nwk_content.strip().splitlines()
        nwk = ''
        for line in lines:
            stripped = line.strip()
            if stripped and stripped[0].isdigit() and len(stripped.split()) <= 2:
                continue   # linha de cabeçalho
            nwk += stripped
        nwk = nwk.strip().rstrip(';').strip()
        if not nwk.startswith('('):
            return 'unknown'

        # Encontrar vírgulas de nível 1 (filhos diretos da raiz)
        depth = 0
        top_comma = -1
        for i, ch in enumerate(nwk):
            if ch == '(':
                depth += 1
            elif ch == ')':
                depth -= 1
            elif ch == ',' and depth == 1:
                top_comma = i
                break   # basta a primeira vírgula top-level para separar os dois filhos

        if top_comma == -1:
            return 'unknown'

        child1 = nwk[1:top_comma]          # conteúdo do 1º filho
        child2_raw = nwk[top_comma + 1:]   # restante (2º filho + ")...")
        # Isolar o 2º filho: tudo até a última ')' de nível 0
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

    # Espectro purificadora / neutra / diversificadora -- cobre os regimes
    # onde um otimo local costuma prender a busca de omega. So testado com
    # branch lengths/kappa ja warm-started (o multistart de omega puro e
    # barato; refazer a busca de branch length e que seria caro de repetir).
    _WARM_START_OMEGA_TRIALS = (0.2, 1.0, 2.5)

    def _run_model_multistart(self, fas_file: Path, model_name: str, log_file: Path,
                               warm_start_kappa: float, fitted_tree: str,
                               aln=None) -> Dict:
        """Roda o mesmo modelo com warm-start de branch length/kappa varias
        vezes, cada uma com omega inicial diferente, e fica com o de maior
        lnL. Mitiga o risco medido do warm-start (2026-09-18: em teste com
        loci reais, ~1/3 convergiu pra otimo local pior partindo so de
        omega=0.5) sem pagar o custo total de reotimizar branch length do
        zero em cada tentativa -- so a parte barata (omega) e repetida.
        """
        best = None
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
        return best if best is not None else last

    # ── Monitoramento do processo codeml ─────────────────────────────────

    @staticmethod
    def _process_cpu_seconds(pid: int) -> Optional[float]:
        """Tempo de CPU acumulado do processo E dos descendentes (s).

        No Debian/Ubuntu, /usr/bin/codeml é um script sh que roda
        /usr/lib/paml/bin/codeml como filho (sem exec): medir só o pid do
        script daria 0 s para sempre e a detecção de inatividade mataria um
        codeml que está trabalhando. psutil se houver; /proc no Linux; None
        se não for possível medir (aí só o timeout vale)."""
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
        """Encerra o codeml e recolhe o processo (sem deixar zumbi/órfão).

        Fora do Windows o codeml roda num grupo de processos próprio
        (start_new_session): o sinal vai para o grupo inteiro, assim o codeml
        real também morre quando o executável é um script que o chama."""
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
        """Encerra todos os codeml em execução (botão Parar). Retorna quantos."""
        with self._processes_lock:
            procs = list(self._active_processes)
        for proc in procs:
            self._terminate(proc)
        return len(procs)

    @staticmethod
    def _prune_keep_labels(tree, keep: set) -> None:
        """Poda as folhas que não estão em `keep` sem perder marcas de ramo.

        Quando a poda deixa um nó interno com um só filho, o Bio.Phylo
        colapsa o nó -- e a marca dele (#1, #2...) sumia. Aqui o filho que
        sobra herda a marca (o ramo resultante é a soma dos dois ramos)."""
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

    # ── Uma execução do codeml (gene × modelo) ───────────────────────────

    def _failed(self, reason: str, exec_start: float = None, **extra) -> Dict:
        d = {
            'output_file': None, 'results_file': None, 'lnL': None, 'np': None,
            'ntime': None, 'omega': None,
            'execution_time': (time.time() - exec_start) if exec_start else 0.0,
            'status': 'failed', 'fail_reason': reason, 'stop_count': 0, 'beb_skipped': False,
        }
        d.update(extra)
        return d

    def _run_single_analysis(self, fas_file: Path, model_name: str,
                             log_file: Path = None,
                             warm_start_kappa: float = None,
                             fitted_tree: str = None,
                             omega_override: float = None,
                             aln=None,
                             save_outputs: bool = True) -> Dict:
        """Executa o CODEML para um gene e um modelo. Sempre retorna um dict
        com 'status' = 'success' | 'failed' | 'stopped' (e 'fail_reason').

        Pasta de saída MODELO/ (reprodutível: `cd MODELO && codeml GENE_MODELO.ctl`):
          GENE_MODELO.ctl            todos os parâmetros, caminhos relativos
          GENE_MODELO_seq.fasta      alinhamento exatamente como o codeml leu
          GENE_MODELO_tree.nwk       árvore exatamente como o codeml leu
          GENE_MODELO_results.txt    saída bruta do codeml (mlc)
          GENE_MODELO_sitemap.json   numeração de sítios codeml -> alinhamento

        warm_start_kappa / fitted_tree: ponto de partida vindo do M0
        (fix_blength = 1). omega_override: multi-start de omega.
        """
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
        output_filename = f"{base_name}_{model_name}_results.txt"   # nome esperado pelo painel
        codeml_outfile = f"{prefix}_results.txt"                     # sem espaços para o codeml
        ctl_filename = f"{prefix}.ctl"
        seq_filename = f"{prefix}_seq.fasta"
        tree_filename = f"{prefix}_tree.nwk"

        temp_dir = Path(tempfile.mkdtemp(prefix=f'easypam_{prefix}_', dir=self._get_fast_tempdir()))
        exec_start = time.time()
        try:
            # ── Alinhamento ──────────────────────────────────────────────
            if aln is None:
                aln = read_alignment(fas_file)
            names = list(aln.names)
            excluded: List[str] = []

            # ── Árvore: leitura, poda e desenraizamento ─────────────────
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
                    if first:   # uma vez por gene, não uma por modelo
                        self._emit('warn', msg)
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

            # Cópia FASTA (sempre FASTA, mesmo que o original seja PHYLIP)
            seqs_used = {n: aln.seqs[n] for n in names}
            (temp_dir / seq_filename).write_text(to_fasta(names, seqs_used), encoding='utf-8')

            def _newick(tree) -> str:
                io_ = StringIO()
                # árvore sem comprimentos de ramo: não inventar ':0' em todos
                # os ramos (o Bio.Phylo escreveria 0.00000 no lugar de None)
                has_bl = any(c.branch_length for c in tree.find_clades())
                Phylo.write(tree, io_, 'newick', plain=not has_bl)
                txt = io_.getvalue().strip()
                # Bio.Phylo escreve comprimento no nó raiz (":0.00000;"), que o codeml rejeita
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

            # Árvore marcada (Branch / Branch-site)
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
                # fix_blength = 1 ("initial"): os comprimentos do M0 são só o
                # ponto de partida e continuam sendo estimados. (2 = "fixed"
                # os prenderia nos valores do M0 -- ver pamlDOC, fix_blength.)
                fix_bl = 1
                tree_note = "M0 fitted tree as starting values (fix_blength = 1)"
            else:
                tree_out = _newick(tree_obj)
            n_tips = len(names)
            (temp_dir / tree_filename).write_text(f"{n_tips}  1\n{tree_out}\n", encoding='utf-8')
            self._log(f"[tree] {base_name} [{model_name}]: {tree_filename} ({tree_note})")

            # ── .ctl ─────────────────────────────────────────────────────
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

            # ── Executar o codeml ───────────────────────────────────────
            cmd = [codeml_path, ctl_filename]
            self._log(f"[{model_name}] {base_name}: Running command: {cmd} in {temp_dir}")
            popen_kw = {}
            if platform.system() == 'Windows':
                popen_kw['creationflags'] = getattr(subprocess, 'CREATE_NO_WINDOW', 0)
            else:
                popen_kw['start_new_session'] = True   # ver _terminate
            # stdin FECHADO: quando o codeml encontra um stop codon ele imprime
            # "Press Enter to continue" e chama getchar(); com a entrada padrão
            # fechada ele segue na hora (tratando a coluna como dado ausente)
            # em vez de esperar para sempre por um Enter que nunca chega.
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
                        # skip_beb: o lnL/np/omega usados no LRT já foram escritos
                        # no outfile antes do BEB começar (ver AGENTS.md / METODOS.md)
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
                        # tempo pausado não conta para timeout nem inatividade
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

            if stopped:
                return self._failed(self._t('reason_stopped'), exec_start, status='stopped')

            codeml_out = temp_dir / codeml_outfile
            output_path = model_output_dir / output_filename
            last_line = next((l for l in reversed(stdout_lines) if l.strip()), '')

            # skip_beb mata o processo de propósito: truncar seções incompletas
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

            # ── Guardar arquivos (reprodutibilidade) ─────────────────────
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
                # resultado temporário (M0 implícito): copiar para fora do sandbox
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

            # ── Mapa de sítios (numeração original) ──────────────────────
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

            # Execução que falhou (inatividade, tempo limite, código de erro):
            # a saída parcial fica para diagnóstico, mas com outro nome, para
            # nenhum leitor de *_results.txt (LRT, painel, Atualizar
            # Resultados) usar um lnL de um codeml interrompido.
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

    # Padroes que identificam a linha de resultado final do CODEML, em ordem
    # de preferencia -- cobre as variantes de formato conhecidas ("lnL(ntime:
    # X np: Y): valor" eh a forma padrao; as outras sao formatos mais antigos/
    # alternativos). So essa linha (nunca a "lnL0 = ..." pre-otimizacao) tem
    # lnL, ntime e np juntos.
    _LNL_LINE_PATTERNS = (
        r'lnL[^:]*:\s*([+-]?\d+\.\d+)',
        r'lnL\([^)]*\):\s*([+-]?\d+\.\d+)',
        r'lnL\s*[:=]\s*([+-]?\d+\.\d+)',
    )

    def _extract_model_stats(self, output_file: Path) -> Dict[str, Optional[float]]:
        """Le o outfile do CODEML uma unica vez e extrai lnL, np e ntime da
        mesma linha de resultado final ("lnL(ntime: X  np: Y): valor").

        Substitui os antigos _extract_likelihood/_extract_np/_extract_ntime,
        que abriam e varriam o arquivo tres vezes separadas pra ler tres
        valores que sempre estao na mesma linha.
        """
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
        """Pares (null, alt, nome_da_coluna) validos dado o conjunto de modelos
        selecionados -- fonte unica em lrt_stats.PAIRS (usada tambem pelo
        painel de resultados)."""
        comparisons = [(n, a, lrt_stats.lrt_column(n, a))
                       for n, a in lrt_stats.pairs_for(selected_models)]
        # Compatibilidade com pastas antigas (BranchSite_A)
        if 'BranchSite_A_null' in selected_models and 'BranchSite_A' in selected_models:
            comparisons.append(('BranchSite_A_null', 'BranchSite_A',
                                'lrt_BranchSite_A_null_vs_BranchSite_A'))
        return comparisons

    def _save_summary(self):
        """Salva analysis_summary.tsv: uma linha por gene com status, lnL/np/
        ntime/omega por modelo, ω e p₁ da classe positiva (M2a/M8), e para
        cada par de LRT a estatística 2Δl, o p e o q (BH)."""
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
        "NOTA METODOLOGICA / METHODS NOTE:\n"
        "  Estatistica: 2*(lnL_alternativo - lnL_nulo); valores negativos (o\n"
        "  alternativo nao melhorou) sao truncados em 0 e dao p = 1.\n"
        "  Graus de liberdade = parametros livres a mais no alternativo\n"
        "  (os comprimentos de ramo entram igualmente nos dois modelos):\n"
        "    M0  vs M1a  : df = 2  (p0 e omega0)\n"
        "    M1a vs M2a  : df = 2  (p2 e omega2)\n"
        "    M7  vs M8   : df = 2  (p1 e omega_s)\n"
        "    M8a vs M8   : df = 1  (omega_s livre vs fixo em 1)\n"
        "    M0  vs Branch : df = numero de grupos foreground (#1, #2, ...)\n"
        "  Distribuicao nula: chi2 com esses df. Para M8a vs M8 e Branch-site\n"
        "  (nulo na fronteira, omega = 1 fixo) a significancia usa chi2(1) puro,\n"
        "  como o manual do PAML recomenda para o branch-site; a mistura\n"
        "  0.5*chi2(0) + 0.5*chi2(1) (Self & Liang 1987) e reportada so como\n"
        "  referencia. p-valores calculados com a funcao de sobrevivencia\n"
        "  (chi2.sf), sem arredondar para zero.\n"
        "  q-value = p corrigido por Benjamini-Hochberg (FDR) DENTRO de cada\n"
        "  comparacao (familia = todos os genes testados nesse par de modelos\n"
        "  nesta execucao).\n"
        "  M7 vs M8 pode rejeitar M7 so porque ha sitios neutros (omega = 1);\n"
        "  M8a vs M8 nao tem esse problema (Swanson et al. 2003).\n"
    )

    def _run_lrt_analysis(self):
        """Executa Likelihood Ratio Tests + correcao Benjamini-Hochberg (FDR).

        BH precisa da familia COMPLETA de p-valores de uma comparacao antes de
        corrigir -- por isso o metodo e em duas fases por comparacao. Os
        q-valores ficam em self._lrt_qvalues[(null_model, alt_model)][gene] e
        os p em self._lrt_pvalues, que _save_summary() anexa como colunas
        q_*/p_* no TSV (por isso este metodo roda ANTES de _save_summary()).
        """
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

                # Fase 1: coletar todos os genes validos ANTES de corrigir por BH
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
                    if df is None:  # M0 vs Branch: grupos foreground
                        df = abs(np_alt - np_null) if (np_alt and np_null) else 1
                        if ntime_null is not None and ntime_alt is not None:
                            df = max(1, df - (ntime_alt - ntime_null))

                    raw_stat = 2 * (lnL_alt - lnL_null)
                    lrt_stat = max(0.0, raw_stat)
                    p_value = lrt_stats.p_value(lrt_stat, df, boundary=boundary)
                    p_mix = lrt_stats.p_value_mixture(lrt_stat) if boundary else None
                    df_display = (f"{df} (chi2(1) puro / pure; mistura so referencia)"
                                  if boundary else str(df))
                    collected.append({
                        'gene': gene_name, 'lnL_null': lnL_null, 'lnL_alt': lnL_alt,
                        'np_null': np_null, 'np_alt': np_alt, 'lrt_stat': lrt_stat,
                        'raw_stat': raw_stat, 'df_display': df_display,
                        'p_value': p_value, 'p_value_mixture': p_mix,
                    })

                # Fase 2: BH na familia completa desta comparacao
                for c, q in zip(collected, lrt_stats.bh_qvalues([c['p_value'] for c in collected])):
                    c['q_value'] = q
                self._lrt_qvalues[(null_model, alt_model)] = {c['gene']: c['q_value'] for c in collected}
                self._lrt_pvalues[(null_model, alt_model)] = {c['gene']: c['p_value'] for c in collected}

                # Fase 3: escrever
                sig_count_05 = sig_count_01 = sig_count_q05 = 0
                for c in collected:
                    sig_count_05 += c['p_value'] < 0.05
                    sig_count_01 += c['p_value'] < 0.01
                    sig_count_q05 += c['q_value'] < 0.05
                    f.write(f"Gene: {c['gene']}\n")
                    f.write(f"  lnL {null_model}: {c['lnL_null']:.6f} (np={c['np_null']})\n")
                    f.write(f"  lnL {alt_model}: {c['lnL_alt']:.6f} (np={c['np_alt']})\n")
                    f.write(f"  2Δl = {c['lrt_stat']:.6f}")
                    if c['raw_stat'] < 0:
                        f.write(f"  (bruto {c['raw_stat']:.6f} < 0: otimizacao do alternativo nao "
                                f"alcancou o nulo; considere rodar de novo)")
                    f.write("\n")
                    f.write(f"  df = {c['df_display']}\n")
                    f.write(f"  p-value = {c['p_value']:.6e}\n")
                    if c.get('p_value_mixture') is not None:
                        f.write(f"  p-value (mistura 50:50, referencia -- NAO usado pro q-valor) = "
                                f"{c['p_value_mixture']:.6e}\n")
                    f.write(f"  q-value (BH) = {c['q_value']:.6e}\n")
                    if c['q_value'] < 0.01:
                        f.write(f"  Result: SIGNIFICANT -- {alt_model} better (q < 0.01, BH-corrected)\n")
                    elif c['q_value'] < 0.05:
                        f.write(f"  Result: SIGNIFICANT -- {alt_model} better (q < 0.05, BH-corrected)\n")
                    else:
                        f.write(f"  Result: not significant (q >= 0.05, BH-corrected)\n")
                    f.write("\n" + "-"*60 + "\n\n")

                total_valid = len(collected)
                f.write("\nRESUMO:\n")
                f.write(f"  Total de genes analisados: {total_valid}\n")
                if total_valid > 0:
                    f.write(f"  Significativo em p bruto < 0.05: {sig_count_05} ({100*sig_count_05/total_valid:.1f}%)\n")
                    f.write(f"  Significativo em p bruto < 0.01: {sig_count_01} ({100*sig_count_01/total_valid:.1f}%)\n")
                    f.write(f"  Significativo em q (BH) < 0.05: {sig_count_q05} ({100*sig_count_q05/total_valid:.1f}%)\n")
                f.write("\n")
                self._emit('info', self._t('lrt_pair_done', null=null_model, alt=alt_model,
                                           n=total_valid, sig=sig_count_q05))

        self._emit('debug', f"LRT results saved: {lrt_file}")

    # ══════════════════════════════════════════════════════════════════
    # WGS / ndata MODE  (genome-scale multi-gene analysis)
    # ══════════════════════════════════════════════════════════════════

    @staticmethod
    def _fasta_to_phylip_block(fas_path: Path) -> Optional[str]:
        """Converte um arquivo FASTA para um bloco no formato PHYLIP do CODEML.

        Retorna None se o arquivo já estiver em formato PHYLIP (primeira linha
        com '<ntaxa> <nsite>').
        """
        text = fas_path.read_text(encoding='utf-8', errors='ignore').strip()
        lines = text.splitlines()
        if not lines:
            return None

        # Detectar se já é PHYLIP (primeira linha = dois inteiros)
        first = lines[0].strip().split()
        if len(first) == 2 and first[0].isdigit() and first[1].isdigit():
            return text + '\n'

        # Parsear FASTA
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
            print(f"  [WARN] {fas_path.name}: sequências com tamanhos diferentes — pulando")
            return None
        n_sites = lengths.pop()

        # Montar bloco PHYLIP
        block_lines = [f" {n_taxa} {n_sites}"]
        for name in order:
            # PHYLIP: nome com 10 chars (padded/truncated)
            padded = name[:10].ljust(10)
            block_lines.append(f"{padded}  {sequences[name]}")
        return '\n'.join(block_lines) + '\n'

    @staticmethod
    def regenerate_summary_files(results_folder: Path) -> Dict[str, str]:
        """
        Atualiza os 3 arquivos de síntese a partir de resultados já existentes
        
        Detecta automaticamente quais modelos estão presentes na pasta e regenera:
        - analysis_summary.tsv: Tabela com lnL, np, ω, e LRTs
        - batch_analysis_log.txt: Log consolidado de todas as análises
        - LRT_results.txt: Resultados detalhados dos testes LRT
        
        Considera modelos neutros:
        - M1a é neutro de M2a
        - M7 é neutro de M8
        - BranchSite_A_null é neutro de BranchSite_A
        - M0 é neutro de Branch
        
        Args:
            results_folder: Pasta contendo os subdirectórios de modelos (M0, M1a, etc.)
        
        Returns:
            Dict com paths dos arquivos gerados: {'analysis_summary', 'batch_analysis_log', 'LRT_results'}
        """
        results_folder = Path(results_folder)
        
        if not results_folder.exists():
            raise ValueError(f"Results folder not found: {results_folder}")
        
        generated_files = {}

        try:
            # ═══ 1. REGENERAR LRT_results.txt (+ q-valores BH) ═══
            # Roda primeiro porque analysis_summary.tsv precisa dos q-valores
            # pra anexar as colunas q_* -- mesma dependência que uma run ao
            # vivo tem (_run_lrt_analysis antes de _save_summary).
            print("\n[1/3] Generating LRT_results.txt...")
            lrt_file, qvalues = CodemlBatchAnalysis._regenerate_lrt_results(results_folder)
            if lrt_file:
                generated_files['LRT_results'] = str(lrt_file)
                print(f"  OK: {lrt_file.name}")

            # ═══ 2. REGENERAR analysis_summary.tsv ═══
            print("\n[2/3] Generating analysis_summary.tsv...")
            summary_file = CodemlBatchAnalysis._regenerate_analysis_summary(results_folder, qvalues)
            if summary_file:
                generated_files['analysis_summary'] = str(summary_file)
                print(f"  OK: {summary_file.name}")

            # ═══ 3. REGENERAR batch_analysis_log.txt ═══
            print("\n[3/3] Generating batch_analysis_log.txt...")
            log_file = CodemlBatchAnalysis._regenerate_batch_log(results_folder)
            if log_file:
                generated_files['batch_analysis_log'] = str(log_file)
                print(f"  OK: {log_file.name}")
            
            print(f"\n[SUCCESS] All files regenerated successfully!")
            return generated_files
        
        except Exception as e:
            print(f"[ERRO] Falha ao regenerar arquivos: {str(e)}")
            traceback.print_exc()
            return {}
    
    @staticmethod
    def _regenerate_analysis_summary(results_folder: Path, qvalues: Optional[dict] = None) -> Optional[Path]:
        """Regenera analysis_summary.tsv.

        qvalues (opcional): {(null_model, alt_model): {gene: q_value}}, vindo
        de _regenerate_lrt_results() -- anexa colunas q_{null}_vs_{alt} do
        mesmo jeito que uma run normal faz via _save_summary(). Sem isso o
        TSV regenerado ficaria sem correção de múltiplos testes, divergindo
        do schema de uma run ao vivo."""
        results_folder = Path(results_folder)
        qvalues = qvalues or {}
        summary_file = results_folder / "analysis_summary.tsv"

        # Mapeamento de nomes de pasta (legados) para nomes de modelo (atuais)
        model_name_mapping = CodemlBatchAnalysis._LEGACY_MODEL_NAMES
        
        # Descobrir quais modelos estão presentes
        models = []
        for item in results_folder.iterdir():
            if item.is_dir() and item.name not in ['reports']:
                # Mapear nomes antigos para novos
                model_name = model_name_mapping.get(item.name, item.name)
                models.append(model_name)
        
        models = sorted(set(models))  # Remove duplicatas e ordena
        
        if not models:
            print("  [WARN] No model folders found")
            return None
        
        # Coletar dados de todos os genes
        data = {}
        
        # Mapa reverso: nome do modelo novo -> nome da pasta antiga
        reverse_mapping = {v: k for k, v in model_name_mapping.items()}
        
        for model in models:
            # Usar o nome da pasta original (se existir) para encontrar os arquivos
            folder_name = reverse_mapping.get(model, model)
            model_folder = results_folder / folder_name
            if not model_folder.exists():
                continue
            
            for results_file in sorted(model_folder.glob("*_results.txt")):
                # Extrair nome do gene - precisa usar o nome da pasta original nos arquivos
                gene_name = results_file.name.split(f'_{folder_name}_results')[0]
                
                if gene_name not in data:
                    data[gene_name] = {'Gene': gene_name}
                
                # Extrair valores
                try:
                    with open(results_file, 'r', encoding='utf-8', errors='ignore') as f:
                        content = f.read()
                    
                    # Extrair lnL
                    lnL_match = re.search(r'lnL\(ntime:.*?\):\s+([-\d.]+)', content)
                    lnL = float(lnL_match.group(1)) if lnL_match else None
                    
                    # Extrair np e ntime
                    np_match = re.search(r'lnL\(ntime:\s*(\d+)\s+np:\s*(\d+)\)', content)
                    np_val    = int(np_match.group(2)) if np_match else None
                    ntime_val = int(np_match.group(1)) if np_match else None

                    # Extrair omega
                    omega = SitesParser.extract_omega_robust(results_file)
                    
                    # Extrair tempo de execução
                    time_match = re.search(r'Time used:\s+(\d+):(\d+)', content)
                    exec_time = None
                    if time_match:
                        m = int(time_match.group(1))
                        s = int(time_match.group(2))
                        exec_time = m * 60 + s
                    
                    # Contar STOPs
                    stop_count = content.count('***')
                    
                    # Guardar dados
                    data[gene_name][f'{model}_lnL']   = lnL
                    data[gene_name][f'{model}_np']    = np_val
                    data[gene_name][f'{model}_ntime'] = ntime_val
                    data[gene_name][f'{model}_omega'] = omega
                    data[gene_name][f'{model}_time'] = exec_time
                    data[gene_name][f'{model}_stops'] = stop_count

                    # Para Branch model: guardar omegas por tag (#1, #2, ... e background)
                    # A linha "w (dN/dS) for branches:" lista os grupos na ordem:
                    #   [0]=background, [1]=#1, [2]=#2, ... conforme definido no PAML
                    if model == 'Branch':
                        tag_omegas = SitesParser.extract_omega_by_tags(results_file)
                        for tag, tag_omega in tag_omegas.items():
                            data[gene_name][f'{model}_{tag}_omega'] = tag_omega

                    # Se é Branch-site ou Branch-site_null, extrair dados de classes de sítios
                    if 'Branch-site' in model:
                        class_data = SitesParser.extract_branchsite_class_data(results_file)
                        if class_data:
                            # Armazenar os dados de classe para exibição estruturada
                            data[gene_name][f'{model}_class_data'] = class_data
                
                except Exception as e:
                    print(f"  [WARN] Error processing {gene_name} ({model}): {str(e)}")

        
        # Calcular LRTs (pares em lrt_stats.PAIRS) e p-valores
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

        # Anexar q-valores (BH), se fornecidos por _regenerate_lrt_results()
        for (null_model, alt_model), gene_qvals in qvalues.items():
            col = f'q_{null_model}_vs_{alt_model}'
            for gene_name, q in gene_qvals.items():
                if gene_name in data:
                    data[gene_name][col] = q

        # ═══ PÓS-PROCESSAMENTO: Expandir dados de classes Branch-site ═══
        # Adicionar colunas de foreground omega para cada classe
        for gene_name in data:
            row = data[gene_name]
            
            # Se existe dados de classe do Branch-site, extrair e adicionar colunas
            if 'Branch-site_class_data' in row and row['Branch-site_class_data']:
                class_data = row['Branch-site_class_data']
                
                # Classes em ordem: 0, 1, 2a, 2b
                for cls in ['0', '1', '2a', '2b']:
                    if cls in class_data:
                        # Adicionar colunas com proporção, background omega e foreground omega
                        row[f'Branch-site_class{cls}_prop'] = class_data[cls].get('prop')
                        row[f'Branch-site_class{cls}_bg_w'] = class_data[cls].get('bg_w')
                        row[f'Branch-site_class{cls}_fg_w'] = class_data[cls].get('fg_w')
            
            # Remover a coluna temporária class_data (não salvar no TSV)
            if 'Branch-site_class_data' in row:
                del row['Branch-site_class_data']
            if 'Branch-site_null_class_data' in row:
                del row['Branch-site_null_class_data']
        
        # Converter para DataFrame e salvar. p/q em notação científica: com
        # '%.6f' um p de 4e-23 virava "0.000000" (o "p = 0" de volta).
        for gene_name, (st, _reason) in gene_status.items():
            if gene_name not in data:   # falhou em todos os modelos: sem *_results.txt
                data[gene_name] = {'Gene': gene_name, 'status': st}
        df = pd.DataFrame(list(data.values()))
        if 'status' in df.columns:   # mesma posição que numa execução normal
            df.insert(1, 'status', df.pop('status').fillna('ok'))
        for col in df.columns:
            if col.startswith(('p_', 'q_')):
                df[col] = [f"{v:.6e}" if isinstance(v, (int, float)) and pd.notna(v) else 'NA'
                           for v in df[col]]
        df.to_csv(summary_file, sep='\t', index=False, float_format='%.6f')

        # Alertar sobre genes com .ctl mas sem resultado (worker crash / sessão interrompida)
        CodemlBatchAnalysis._find_orphaned_analyses(results_folder)

        return summary_file
    
    @staticmethod
    def _find_orphaned_analyses(results_folder: Path) -> dict:
        """Detecta genes com arquivo .ctl mas sem resultado (_results.txt).

        Retorna dict  { gene_name: [model, ...] }  listando, para cada gene,
        os modelos cujo CODEML foi iniciado (ctl gravado) mas nunca concluiu.
        Causas típicas: worker paralelo morreu por pressão de memória ou
        crash numérico do CODEML, sessão interrompida pelo usuário.

        Use _regenerate_analysis_summary() depois de corrigir os órfãos para
        atualizar o TSV.
        """
        results_folder = Path(results_folder)
        _legacy = CodemlBatchAnalysis._LEGACY_MODEL_NAMES
        _reverse = {v: k for k, v in _legacy.items()}

        orphaned: dict = {}
        gene_status = CodemlBatchAnalysis._read_gene_status(results_folder)
        for item in results_folder.iterdir():
            if not item.is_dir() or item.name in {'reports'}:
                continue
            model = _legacy.get(item.name, item.name)
            folder_name = item.name  # nome real da pasta

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
        """Regenera batch_analysis_log.txt"""
        results_folder = Path(results_folder)
        log_file = results_folder / "batch_analysis_log.txt"

        # Mapeamento de nomes legados → atuais (centralizado na constante de classe)
        model_name_mapping = CodemlBatchAnalysis._LEGACY_MODEL_NAMES
        # Mapeamento inverso: nome de exibição → nome da pasta no disco
        reverse_mapping = {v: k for k, v in model_name_mapping.items()}

        with open(log_file, 'w', encoding='utf-8') as f:
            f.write("="*80 + "\n")
            f.write("LOG DE ANÁLISE CODEML (REGENERADO)\n")
            f.write("="*80 + "\n")
            f.write(f"Regenerado em: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Pasta de resultados: {results_folder}\n")
            f.write("="*80 + "\n\n")

            f.write("RESUMO DA ANÁLISE:\n")
            f.write("-"*80 + "\n")

            # Descobrir modelos e genes — usa item.name (pasta real) para o split do gene
            models = set()
            genes = set()

            for item in results_folder.iterdir():
                if item.is_dir() and item.name not in ['reports']:
                    model_display = model_name_mapping.get(item.name, item.name)
                    models.add(model_display)

                    for results_file in item.glob("*_results.txt"):
                        # CORREÇÃO: split pelo nome real da pasta (item.name), não pelo
                        # nome de exibição (model_display), que pode ser diferente.
                        gene = results_file.name.split(f'_{item.name}_results')[0]
                        genes.add(gene)

            f.write(f"Modelos encontrados: {', '.join(sorted(models))}\n")
            f.write(f"Genes encontrados: {len(genes)} genes\n")
            f.write(f"  {', '.join(sorted(genes)[:5])}" + ("..." if len(genes) > 5 else "") + "\n")
            f.write("\n")

            # Detalhes de cada gene/modelo
            f.write("RESULTADOS DETALHADOS:\n")
            f.write("-"*80 + "\n\n")

            for gene in sorted(genes):
                f.write(f"Gene: {gene}\n")
                f.write("-"*40 + "\n")

                for model_display in sorted(models):
                    # CORREÇÃO: usar o nome real da pasta (folder_name) para construir
                    # os caminhos — o nome de exibição pode não corresponder ao nome no disco.
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
                            f.write(f"  {model_display:20s} | Erro ao ler arquivo\n")
                    else:
                        f.write(f"  {model_display:20s} | Não encontrado\n")

                f.write("\n")

            f.write("="*80 + "\n")
            f.write("FIM DO LOG\n")
            f.write("="*80 + "\n")

        return log_file
    
    @staticmethod
    def _read_gene_status(results_folder: Path) -> Dict[str, Tuple[str, str]]:
        """{gene: (status, motivo)} de genes_status.tsv ({} se não existir)."""
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
        """Regenera LRT_results.txt, com correcao Benjamini-Hochberg (FDR) por
        comparacao -- mesma logica de duas fases que _run_lrt_analysis (coleta
        todos os p-valores da familia, corrige, so depois escreve).

        Retorna (lrt_file, qvalues) onde qvalues e
        {(null_model, alt_model): {gene: q_value}} -- usado por
        regenerate_summary_files() pra anexar colunas q_* no TSV, do mesmo
        jeito que uma run normal faz via _save_summary()."""
        results_folder = Path(results_folder)
        lrt_file = results_folder / "LRT_results.txt"
        qvalues: dict = {}
        gene_status = CodemlBatchAnalysis._read_gene_status(results_folder)

        # Mapeamento de nomes legados → atuais (centralizado na constante de classe)
        model_name_mapping = CodemlBatchAnalysis._LEGACY_MODEL_NAMES
        reverse_mapping = {v: k for k, v in model_name_mapping.items()}
        
        # Descobrir quais modelos estão presentes
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
        
        # Pares e df: fonte unica em lrt_stats.PAIRS
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

                # Fase 1: coletar todos os genes validos ANTES de corrigir por BH
                # (q-valor de um gene depende do rank do seu p-valor entre todos
                # os outros da mesma comparacao -- nao da pra escrever linha a
                # linha como o p bruto).
                collected = []
                for gene in sorted(genes):
                    if gene_status.get(gene, ('',))[0] == 'failed':
                        continue   # pasta antiga: saída de um codeml interrompido
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
                        
                        # Extrair lnL
                        null_lnL_match = re.search(r'lnL\(.*?\):\s+([-\d.]+)', null_content)
                        alt_lnL_match = re.search(r'lnL\(.*?\):\s+([-\d.]+)', alt_content)
                        
                        if not (null_lnL_match and alt_lnL_match):
                            continue
                        
                        lnL_null = float(null_lnL_match.group(1))
                        lnL_alt = float(alt_lnL_match.group(1))
                        
                        # Calcular LRT
                        lrt_stat = 2 * (lnL_alt - lnL_null)

                        # Para M0 vs Branch: df = (np_Branch − np_M0) − (ntime_Branch − ntime_M0)
                        # M0 usa árvore não-enraizada (ntime = 2n-3); Branch usa a árvore
                        # rotulada/enraizada (ntime = 2n-2).  Sem a correção, df = k+1.
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
                                gene_df = 1  # fallback seguro

                        # Branch-site: a mistura 50:50 de χ²(0)+χ²(1) e a nula
                        # assintotica correta (Self & Liang 1987), mas o manual do
                        # PAML recomenda explicitamente usar χ²(1) puro em vez dela
                        # ("guard against violations of model assumptions") -- ver
                        # nota completa em _run_lrt_analysis. Consistente com la.
                        is_boundary = lrt_stats.PAIRS[(null_model, alt_model)]['boundary']
                        lrt_stat = max(0.0, lrt_stat)
                        p_value = lrt_stats.p_value(lrt_stat, gene_df, boundary=is_boundary)
                        p_value_mixture = lrt_stats.p_value_mixture(lrt_stat) if is_boundary else None

                        collected.append({
                            'gene': gene, 'lnL_null': lnL_null, 'lnL_alt': lnL_alt,
                            'lrt_stat': lrt_stat, 'df_display':
                                f"{gene_df} (chi2(1) puro; mistura so referencia)" if is_boundary else str(gene_df),
                            'p_value': p_value, 'p_value_mixture': p_value_mixture,
                        })

                    except Exception:
                        continue

                # Fase 2: corrigir por BH usando a familia completa desta comparacao
                if collected:
                    qvals = stats.false_discovery_control([c['p_value'] for c in collected], method='bh')
                    for c, q in zip(collected, qvals):
                        c['q_value'] = q
                qvalues[(null_model, alt_model)] = {c['gene']: c['q_value'] for c in collected}

                # Fase 3: escrever (p bruto e q-valor BH lado a lado)
                sig_count_05 = sig_count_01 = sig_count_q05 = 0
                for c in collected:
                    if c['p_value'] < 0.05:
                        sig_count_05 += 1
                    if c['p_value'] < 0.01:
                        sig_count_01 += 1
                    if c['q_value'] < 0.05:
                        sig_count_q05 += 1

                    f.write(f"Gene: {c['gene']}\n")
                    f.write(f"  lnL {null_model}: {c['lnL_null']:.6f}\n")
                    f.write(f"  lnL {alt_model}: {c['lnL_alt']:.6f}\n")
                    f.write(f"  2Δl = {c['lrt_stat']:.6f}\n")
                    f.write(f"  df = {c['df_display']}\n")
                    f.write(f"  p-value = {c['p_value']:.6e}\n")
                    if c['p_value_mixture'] is not None:
                        f.write(f"  p-value (mistura 50:50, referencia -- NAO usado pro q-valor) = "
                                f"{c['p_value_mixture']:.6e}\n")
                    f.write(f"  q-value (BH) = {c['q_value']:.6e}\n")

                    if c['q_value'] < 0.01:
                        f.write(f"  Result: [OK][OK] {alt_model} significantly better (q < 0.01, BH-corrected)\n")
                    elif c['q_value'] < 0.05:
                        f.write(f"  Result: [OK] {alt_model} significantly better (q < 0.05, BH-corrected)\n")
                    else:
                        f.write(f"  Result: [ERROR] No significant difference (q >= 0.05, BH-corrected)\n")

                    f.write("\n" + "-"*60 + "\n\n")

                # Resumo
                total_valid = len(collected)
                if total_valid > 0:
                    f.write("\nRESUMO:\n")
                    f.write(f"  Total de genes analisados: {total_valid}\n")
                    f.write(f"  Significativo em p bruto < 0.05: {sig_count_05} ({100*sig_count_05/total_valid:.1f}%)\n")
                    f.write(f"  Significativo em p bruto < 0.01: {sig_count_01} ({100*sig_count_01/total_valid:.1f}%)\n")
                    f.write(f"  Significativo em q (BH) < 0.05: {sig_count_q05} ({100*sig_count_q05/total_valid:.1f}%)\n")
                    f.write("\n")

        return lrt_file, qvalues
