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

# Absolute path to the bundled codeml binary — works regardless of CWD.
_APP_ROOT   = Path(__file__).resolve().parent.parent.parent
_CODEML_BIN = (_APP_ROOT / 'bin' / 'codeml.exe'
               if platform.system() == 'Windows'
               else _APP_ROOT / 'bin' / 'codeml')


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
            'CodonFreq': 7
        },
        'M1a': {
            'description': 'Nearly Neutral - ω < 1 or = 1',
            'model': 0,
            'NSsites': 1,
            'fix_omega': 0,
            'omega': 0.5,
            'CodonFreq': 7
        },
        'M2a': {
            'description': 'Positive Selection - adds ω > 1 class',
            'model': 0,
            'NSsites': 2,
            'fix_omega': 0,
            'omega': 0.5,
            'CodonFreq': 7
        },
        'M7': {
            'description': 'Beta distribution - ω < 1',
            'model': 0,
            'NSsites': 7,
            'fix_omega': 0,
            'omega': 0.5,
            'CodonFreq': 7
        },
        'M8': {
            'description': 'Beta + ω - adds ω > 1 class',
            'model': 0,
            'NSsites': 8,
            'fix_omega': 0,
            'omega': 0.5,
            'CodonFreq': 7
        },
        'Branch': {
            'description': 'Branch model - different ω for foreground',
            'model': 2,
            'NSsites': 0,
            'fix_omega': 0,
            'omega': 0.5,
            'CodonFreq': 7
        },
        'Branch-site': {
            'display_name': 'Branch-site',
            'description': 'Branch-site model - ω varies across sites and branches',
            'model': 2,
            'NSsites': 2,
            'fix_omega': 0,
            'omega': 0.5,
            'CodonFreq': 7
        },
        'Branch-site_null': {
            'description': 'Branch-site null model - fixes ω=1 (null for Branch-site model)',
            'model': 2,
            'NSsites': 2,
            'fix_omega': 1,
            'omega': 1.0,
            'CodonFreq': 7
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
            'references': 'Yang et al. (2005)'
        },
        'M8': {
            'full_name': 'Beta & Positive Selection Model',
            'test_type': 'Site Model',
            'parameters': 'Beta(p,q) for ω < 1, PLUS additional class with ω > 1',
            'purpose': 'Alternative hypothesis: continuous distribution + discrete class for positive selection.',
            'interpretation': 'Reject M7 at p < 0.05 = evidence for positive selection. More flexible than M2a.',
            'use_case': 'Alternative test for positive selection; compare against M7.',
            'references': 'Yang et al. (2005)'
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
            ('M7', 'M8', 'Alternative test for positive selection')
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
        'M8': 'M7',                # M8 (alternativo) -> M7 (nulo)
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
    def auto_complete_null_models(selected_models: List[str], include_neutral: bool = True) -> List[str]:
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
                null_model = CodemlBatchAnalysis.NULL_MODEL_PAIRS[model]
                completed_models.add(null_model)
        
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
                             fix_kappa_heuristic: bool = False,
                             fix_blength: int = 0) -> str:
        """Gera conteúdo do arquivo .ctl baseado no modelo.

        Parameters:
        - seqfile, treefile, outfile : caminhos/nomes a inserir no .ctl
        - model_config               : dict com parâmetros do modelo
        - omega                      : valor inicial de ω
        - cleandata                  : 0/1 — remover sítios ambíguos
        - model_name                 : nome do modelo (para verificar modelos neutros)
        - kappa                      : κ estimado pelo M0 (warm-start ou fixado)
        - fix_kappa_heuristic        : se True, insere fix_kappa=1 (modo heurístico —
                                       κ fixado no valor do M0; acelera otimização mas
                                       é uma aproximação.  LRT ainda válido se ambos os
                                       modelos do par usarem o mesmo κ fixado.)
        - fix_blength                : 0=estimar do zero  2=warm-start do M0 (safe —
                                       branch lengths re-estimados livremente a partir
                                       de valores iniciais melhores; sem impacto nos
                                       resultados finais)

        Branch-site_null: fix_omega=1, omega=1.0 (ω₂=1 fixado, Yang et al. 2005)
        """
        # Para Branch-site_null, garantir fix_omega=1 mesmo se o usuário tiver editado
        fix_omega = model_config['fix_omega']
        final_omega = omega

        if model_name in self.NEUTRAL_MODELS:
            fix_omega = 1
            final_omega = 1.0
        elif fix_omega == 0:
            final_omega = omega

        # ── Bloco de kappa ────────────────────────────────────────────────────
        # fix_kappa_heuristic=True  → fix_kappa=1, kappa=<valor M0>  (modo rápido)
        # fix_kappa_heuristic=False → fix_kappa=0, kappa=<warm-start> (padrão)
        if kappa is not None and 0.1 <= kappa <= 20:
            if fix_kappa_heuristic:
                kappa_block = (
                    f"\n    fix_kappa = 1              * κ fixado no valor M0 (modo heurístico)"
                    f"\n        kappa = {kappa:.4f}     * ts/tv ratio estimado pelo M0"
                )
            else:
                kappa_block = (
                    f"\n    fix_kappa = 0              * κ livre (re-estimado)"
                    f"\n        kappa = {kappa:.4f}     * ts/tv warm-start do M0 (ponto de partida)"
                )
        else:
            kappa_block = ""

        # ── Linha fix_blength ─────────────────────────────────────────────────
        # fix_blength=2: warm-start dos branch lengths do M0; re-estimados livremente.
        # Matematicamente equivalente a fix_blength=0, porém converge mais rápido.
        blength_line = f"\n  fix_blength = {fix_blength}              * 0=estimar do zero  2=warm-start dos branch lengths" if fix_blength != 0 else ""

        ctl_template = f"""      seqfile = {seqfile}
     treefile = {treefile}
      outfile = {outfile}

        noisy = 1              * 0-9: output detail (1=minimal stdout, model fit written to outfile)
      verbose = 1              * More or less detailed report in outfile
      seqtype = 1              * Data type
        ndata = 1              * Number of data sets or loci
        icode = 0              * Genetic code
    cleandata = {cleandata}              * Remove sites with ambiguity data?

        model = {model_config['model']}         * Models for omega varying across lineages
      NSsites = {model_config['NSsites']}          * Models for omega varying across sites
    CodonFreq = {model_config['CodonFreq']}        * Codon frequencies
      estFreq = 0              * Use observed freqs or estimate freqs by ML
        clock = 0              * Clock model
    fix_omega = {fix_omega}         * Estimate or fix omega
        omega = {final_omega}        * Initial or fixed omega{kappa_block}{blength_line}
"""
        return ctl_template

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
        complexos (M1a, M2a, M7, M8, Branch-site) como ponto de partida via
        fix_blength = 2 ("inicializar a partir dos valores fornecidos na árvore mas
        ainda re-estimar livremente").

        fix_blength = 2  ≠  fix_blength = 1
          • fix_blength = 1: branch lengths FIXADOS (aproximação, altera resultados).
          • fix_blength = 2: branch lengths usados só como WARM-START; o otimizador
            ainda os re-estima livremente.  Os resultados finais são matematicamente
            idênticos a qualquer outro ponto de partida — apenas convergem mais rápido.

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

    def interactive_setup(self):
        """Configuração interativa via input do usuário"""
        
        print("\n" + "="*80)
        print("CODEML INTERACTIVE BATCH ANALYSIS")
        print("="*80 + "\n")
        
        # 1. Pasta com arquivos .fas / .fasta / .phy / .phylip
        while True:
            input_folder = input("Enter path to folder with .fas/.fasta/.phy/.phylip files: ").strip().strip('"')
            input_path = Path(input_folder)
            if input_path.exists() and input_path.is_dir():
                fas_files = (
                    list(input_path.glob("*.fas"))
                    + list(input_path.glob("*.fasta"))
                    + list(input_path.glob("*.phy"))
                    + list(input_path.glob("*.phylip"))
                )
                if fas_files:
                    print(f"   [OK] Found {len(fas_files)} sequence files (.fas/.fasta/.phy/.phylip)")
                    self.config['input_folder'] = input_path
                    break
                else:
                    print("   [ERROR] No .fas/.fasta/.phy/.phylip files found in this folder. Try again.")
            else:
                print("   [ERROR] Folder not found. Try again.")
        
        # 2. Arquivo de árvore
        while True:
            tree_file = input("\nEnter path to tree file (.tree or .txt): ").strip().strip('"')
            tree_path = Path(tree_file)
            if tree_path.exists() and tree_path.is_file():
                print(f"   [OK] Tree file loaded: {tree_path.name}")
                self.config['tree_file'] = tree_path
                break
            else:
                print("   [ERROR] Tree file not found. Try again.")
        
        # 3. Pasta de saída
        output_folder = input("\nEnter path for output folder: ").strip().strip('"')
        output_path = Path(output_folder)
        output_path.mkdir(parents=True, exist_ok=True)
        print(f"   [OK] Output folder: {output_path}")
        self.config['output_folder'] = output_path
        
        # 4. Selecionar modelos
        print("\n" + "="*80)
        print("AVAILABLE MODELS")
        print("="*80 + "\n")
        
        print("SITE MODELS (variation among sites):")
        for i, (code, info) in enumerate([
            ('M0', self.MODEL_CONFIGS['M0']),
            ('M1a', self.MODEL_CONFIGS['M1a']),
            ('M2a', self.MODEL_CONFIGS['M2a']),
            ('M7', self.MODEL_CONFIGS['M7']),
            ('M8', self.MODEL_CONFIGS['M8'])
        ], 1):
            print(f"  {i}. {code:5s} - {info['description']}")
        
        print("\nBRANCH MODEL (variation among branches):")
        print(f"  7. Branch - {self.MODEL_CONFIGS['Branch']['description']}")
        
        print("\nBRANCH-SITE MODELS (variation in both):")
        print(f"  8. Branch-site      - {self.MODEL_CONFIGS['Branch-site']['description']}")
        print(f"  9. Branch-site_null - {self.MODEL_CONFIGS['Branch-site_null']['description']}")
        
        print("\nEnter model numbers separated by spaces (e.g., '1 4 5' for M0, M7, M8)")
        print("Or enter 'all' for all site models (recommended for testing positive selection)")
        
        model_input = input("\nSelect models: ").strip().lower()
        
        model_map = {
            '1': 'M0', '2': 'M1a', '3': 'M2a', '4': 'M7',
            '5': 'M8', '6': 'Branch',
            '7': 'Branch-site', '8': 'Branch-site_null'
        }
        
        if model_input == 'all':
            selected_models = ['M0', 'M1a', 'M2a', 'M7', 'M8']
            print("   [OK] Selected all site models")
        else:
            numbers = model_input.split()
            selected_models = [model_map[n] for n in numbers if n in model_map]
            if not selected_models:
                print("   [WARN] No valid models selected. Using M0 and M8 as default.")
                selected_models = ['M0', 'M8']
        
        self.config['models'] = selected_models
        print(f"   [OK] Models to run: {', '.join(selected_models)}")
        
        # 5. Timeout
        print("\nSet timeout per analysis (in seconds)")
        print("   Recommended: 1600 (≈27 minutes)")
        timeout_input = input("   Timeout [1600]: ").strip()
        self.config['timeout'] = int(timeout_input) if timeout_input else 1600
        print(f"   [OK] Timeout set to {self.config['timeout']} seconds")
        
        # 6. LRT
        print("\nPerform Likelihood Ratio Tests (LRT)?")
        lrt_input = input("   Run LRT? [Y/n]: ").strip().lower()
        self.config['run_lrt'] = lrt_input != 'n'
        print(f"   [OK] LRT: {'Yes' if self.config['run_lrt'] else 'No'}")
        
        # Resumo da configuração
        print("\n" + "="*80)
        print("CONFIGURATION SUMMARY")
        print("="*80)
        print(f"Input folder:  {self.config['input_folder']}")
        print(f"Tree file:     {self.config['tree_file']}")
        print(f"Output folder: {self.config['output_folder']}")
        print(f"Models:        {', '.join(self.config['models'])}")
        print(f"Timeout:       {self.config['timeout']}s")
        print(f"Run LRT:       {self.config['run_lrt']}")
        print(f"Total genes:   {len(fas_files)}")
        print(f"Total runs:    {len(fas_files) * len(self.config['models'])}")
        print("="*80 + "\n")
        
        confirm = input("Proceed with analysis? [Y/n]: ").strip().lower()
        if confirm == 'n':
            print("Analysis cancelled.")
            return False
        
        return True
    
    def run_batch_analysis(self):
        """Executa análise em batch"""
        
        if not self.config:
            if not self.interactive_setup():
                return
        
        output_folder = self.config['output_folder']
        log_file = output_folder / "batch_analysis_log.txt"
        
        # Criar log inicial
        with open(log_file, 'w', encoding='utf-8') as log:
            log.write("="*80 + "\n")
            log.write("CODEML BATCH ANALYSIS LOG\n")
            log.write("="*80 + "\n")
            log.write(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            log.write(f"Input folder: {self.config['input_folder']}\n")
            log.write(f"Output folder: {output_folder}\n")
            log.write(f"Tree file: {self.config['tree_file']}\n")
            log.write(f"Models: {', '.join(self.config['models'])}\n")
            log.write("="*80 + "\n\n")
        
        # Obter arquivos de sequência (.fas, .fasta, .phy, .phylip)
        _input = self.config['input_folder']
        fas_files = sorted(
            list(_input.glob("*.fas"))
            + list(_input.glob("*.fasta"))
            + list(_input.glob("*.phy"))
            + list(_input.glob("*.phylip")),
            key=lambda p: p.name.lower()
        )
        self.current_total_genes = len(fas_files)
        self.current_processed_genes = 0

        n_workers = max(1, int(self.config.get('n_workers', 1)))

        print("\n" + "="*80)
        print("STARTING BATCH ANALYSIS")
        print("="*80)
        print(f"Processing {len(fas_files)} genes × {len(self.config['models'])} models  |  workers: {n_workers}\n")

        start_time = time.time()

        def _process_gene(args):
            idx, fas_file = args
            pause_event = self.config.get('pause_event')
            stop_event = self.config.get('stop_event')

            if stop_event is not None and stop_event.is_set():
                return fas_file.stem, {}

            if pause_event is not None:
                pause_event.wait()

            # Resetar flag manual-all para cada gene
            manual_all_ev = self.config.get('manual_continue_all_event')
            if manual_all_ev is not None:
                try:
                    manual_all_ev.clear()
                except Exception:
                    pass  # evento já limpo ou inválido — não é crítico

            print(f"\n{'='*60}")
            print(f"[{idx}/{len(fas_files)}] Gene: {fas_file.stem}")
            print(f"{'='*60}")

            # ── Validar arquivo de sequência antes de passar ao CODEML ───────────
            # CODEML rejeita silenciosamente arquivos com sequências de tamanhos
            # diferentes (cria outfile vazio e sai com código -1).  Detectar aqui
            # evita arquivos de resultado vazios e dá ao usuário uma mensagem clara.
            # Suporta FASTA (>..) e PHYLIP sequential/interleaved (N  L na 1ª linha).
            try:
                _raw_text  = fas_file.read_text(encoding='utf-8', errors='ignore')
                _raw_lines = [l for l in _raw_text.splitlines() if l.strip()]
                _seqs: dict[str, str] = {}

                _first_clean = _raw_lines[0].strip() if _raw_lines else ''
                _is_phylip   = (
                    bool(_first_clean)
                    and _first_clean.split()[0].lstrip('-').isdigit()
                    and not _first_clean.startswith('>')
                )

                if _is_phylip:
                    # PHYLIP sequential: "N  L\nname10+seq\n..."
                    _ph_parts = _first_clean.split()
                    _ph_ns    = int(_ph_parts[0])
                    _ph_ls    = int(_ph_parts[1]) if len(_ph_parts) > 1 else 0
                    # Cada sequência ocupa uma ou mais linhas; nome = primeiros 10 chars
                    _seq_lines = _raw_lines[1:]
                    _cur_name: str | None = None
                    _cur_seq:  list[str]  = []
                    for _sln in _seq_lines:
                        # Nova sequência: linha com nome no início (não é espaço ou continuação)
                        if not _sln.startswith(' ') and len(_seqs) < _ph_ns:
                            if _cur_name is not None:
                                _seqs[_cur_name] = ''.join(_cur_seq)
                            _cur_name = _sln[:10].strip() or f'seq{len(_seqs)+1}'
                            _cur_seq  = [re.sub(r'\s', '', _sln[10:])]
                        elif _cur_name is not None:
                            _cur_seq.append(re.sub(r'\s', '', _sln))
                    if _cur_name is not None:
                        _seqs[_cur_name] = ''.join(_cur_seq)
                else:
                    # FASTA: > header lines
                    _cur: str | None = None
                    _parts: list[str] = []
                    for _ln in _raw_lines:
                        if _ln.startswith('>'):
                            if _cur is not None:
                                _seqs[_cur] = ''.join(_parts)
                            _cur   = _ln[1:].split()[0]
                            _parts = []
                        elif _cur is not None:
                            _parts.append(_ln.strip())
                    if _cur is not None:
                        _seqs[_cur] = ''.join(_parts)

                _fmt_label = "PHYLIP" if _is_phylip else "FASTA"

                if len(_seqs) < 2:
                    print(f"[SKIP] {fas_file.name}: arquivo {_fmt_label} com menos de 2 sequencias — pulando")
                    with open(log_file, 'a', encoding='utf-8') as _log:
                        _log.write(f"[SKIP] {fas_file.stem}: menos de 2 sequencias no {_fmt_label}\n")
                    with self._results_lock:
                        self.results[fas_file.stem] = {}
                        self.current_processed_genes += 1
                    return fas_file.stem, {}

                _lengths = {len(s) for s in _seqs.values()}
                if len(_lengths) != 1:
                    _sorted = sorted(_lengths)
                    print(
                        f"[SKIP] {fas_file.name}: sequencias nao alinhadas "
                        f"(tamanhos: {_sorted}) — pulando"
                    )
                    with open(log_file, 'a', encoding='utf-8') as _log:
                        _log.write(
                            f"[SKIP] {fas_file.stem}: sequencias nao alinhadas "
                            f"(tamanhos distintos: {_sorted})\n"
                        )
                    with self._results_lock:
                        self.results[fas_file.stem] = {}
                        self.current_processed_genes += 1
                    return fas_file.stem, {}

                _seq_len = _lengths.pop()
                if _seq_len < 6:
                    print(f"[SKIP] {fas_file.name}: sequencias muito curtas ({_seq_len} bp) — pulando")
                    with open(log_file, 'a', encoding='utf-8') as _log:
                        _log.write(f"[SKIP] {fas_file.stem}: sequencias com {_seq_len} bp (minimo 6 bp)\n")
                    with self._results_lock:
                        self.results[fas_file.stem] = {}
                        self.current_processed_genes += 1
                    return fas_file.stem, {}

            except Exception as _val_err:
                # Erro ao ler o arquivo — deixar o CODEML tentar e lidar com a falha
                with open(log_file, 'a', encoding='utf-8') as _log:
                    _log.write(f"[WARN] {fas_file.stem}: nao foi possivel validar arquivo de sequencia: {_val_err}\n")

            gene_results = {}
            gene_kappa:       Optional[float] = None   # κ estimado pelo M0 → warm-start
            gene_fitted_tree: Optional[str]   = None   # árvore ajustada M0 → warm-start branch lengths

            # Modo heurístico (opcional, ativado pela GUI):
            # fix_kappa=1 fixa κ no valor do M0 em vez de apenas usá-lo como ponto de partida.
            # Economiza ~20-30 % de iterações por modelo mas é uma aproximação.
            heuristic_mode = bool(self.config.get('heuristic_mode', False))

            # Sempre rodar M0 primeiro (se selecionado):
            #  • κ estimado pelo M0 é usado como warm-start ou fixado nos modelos seguintes
            #  • árvore ajustada pelo M0 (branch lengths ML) é usada como warm-start via fix_blength=2
            models_ordered = (
                ['M0'] + [m for m in self.config['models'] if m != 'M0']
                if 'M0' in self.config['models']
                else list(self.config['models'])
            )

            for model_name in models_ordered:
                if stop_event is not None and stop_event.is_set():
                    break
                if pause_event is not None:
                    pause_event.wait()

                # M0 não usa warm-start (ele É a fonte)
                kappa_for_this       = None if model_name == 'M0' else gene_kappa
                fitted_tree_for_this = None if model_name == 'M0' else gene_fitted_tree
                fix_kappa_for_this   = heuristic_mode and (model_name != 'M0') and (kappa_for_this is not None)

                print(f"  - Running {model_name}...", end=" ", flush=True)
                result = self._run_single_analysis(
                    fas_file=fas_file,
                    model_name=model_name,
                    log_file=log_file,
                    warm_start_kappa=kappa_for_this,
                    fitted_tree=fitted_tree_for_this,
                    fix_kappa_heuristic=fix_kappa_for_this,
                )
                if result:
                    gene_results[model_name] = result
                    lnL = result.get('lnL')
                    t   = result.get('execution_time')

                    # Extrair κ e árvore ajustada do M0 para warm-start dos modelos seguintes
                    if model_name == 'M0' and result.get('output_file'):
                        out_path = Path(result['output_file'])

                        extracted_k = self._extract_kappa(out_path)
                        if extracted_k is not None:
                            gene_kappa = extracted_k

                        extracted_tree = self._extract_fitted_tree(out_path)
                        if extracted_tree:
                            gene_fitted_tree = extracted_tree

                        # Montar sufixo de status para o log
                        ws_parts = []
                        if gene_kappa        is not None: ws_parts.append(f"k={gene_kappa:.3f}")
                        if gene_fitted_tree  is not None: ws_parts.append("bl=ok")
                        ws_tag = "  →warm-start[" + ", ".join(ws_parts) + "]" if ws_parts else ""
                        print(f"[OK]  lnL={lnL:.2f}  t={t:.1f}s{ws_tag}" if lnL is not None and t is not None else f"[OK]{ws_tag}")
                    else:
                        tags = []
                        if kappa_for_this is not None:
                            tags.append(f"k0={'fix' if fix_kappa_for_this else 'warm'}={kappa_for_this:.3f}")
                        if fitted_tree_for_this is not None:
                            tags.append("bl=warm")
                        ws = ("  (" + ", ".join(tags) + ")") if tags else ""
                        print(f"[OK]  lnL={lnL:.2f}  t={t:.1f}s{ws}" if lnL is not None and t is not None else "[OK]")
                else:
                    print("[ERROR]")

            with self._results_lock:
                self.results[fas_file.stem] = gene_results
                self.current_processed_genes += 1

            return fas_file.stem, gene_results

        indexed = list(enumerate(fas_files, 1))

        if n_workers > 1:
            with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
                list(executor.map(_process_gene, indexed))
        else:
            for item in indexed:
                _process_gene(item)
        
        total_time = time.time() - start_time
        
        # Salvar sumário
        print(f"\n{'='*80}")
        print("SAVING RESULTS")
        print(f"{'='*80}")
        self._save_summary()
        
        # Executar LRT
        if self.config['run_lrt'] and len(self.config['models']) > 1:
            print(f"\n{'='*80}")
            print("PERFORMING LIKELIHOOD RATIO TESTS")
            print(f"{'='*80}")
            self._run_lrt_analysis()
        
        # Sumário final
        print(f"\n{'='*80}")
        print("ANALYSIS COMPLETE!")
        print(f"{'='*80}")
        print(f"Total time: {total_time/60:.1f} minutes")
        print(f"Results saved to: {output_folder}")
        print(f"Log file: {log_file}")
        print(f"{'='*80}\n")
    
    def _run_single_analysis(self, fas_file: Path, model_name: str,
                            log_file: Path,
                            warm_start_kappa: float = None,
                            fitted_tree: str = None,
                            fix_kappa_heuristic: bool = False) -> Optional[Dict]:
        """Executa análise CODEML para um arquivo e modelo.

        Parâmetros de otimização de velocidade (sem impacto nos resultados):
          warm_start_kappa  : κ estimado pelo M0 — usado como ponto de partida para
                              a busca de κ em modelos subsequentes (fix_kappa=0).
          fitted_tree       : árvore Newick com branch lengths otimizados pelo M0 —
                              escrita no sandbox e referenciada com fix_blength=2,
                              de forma que o otimizador parte de valores já próximos
                              do ótimo.  Os branch lengths são re-estimados livremente;
                              os resultados finais são matematicamente idênticos.
          fix_kappa_heuristic: se True, insere fix_kappa=1 no .ctl (modo heurístico —
                              κ fixado no valor M0; acelera ~20-30 % por modelo mas
                              é uma aproximação.  Ativado pelo toggle na GUI.)
        """
        base_name = fas_file.stem
        # start from default config and allow GUI-provided custom overrides
        model_config = dict(self.MODEL_CONFIGS.get(model_name, {}))
        try:
            # GUI uses 'custom_model_params'; keep backward-compatible key 'custom_model_configs'
            custom_configs = self.config.get('custom_model_params', None)
            if custom_configs is None:
                custom_configs = self.config.get('custom_model_configs', {}) or {}
            else:
                custom_configs = custom_configs or {}

            if model_name in custom_configs:
                for k, v in custom_configs[model_name].items():
                    model_config[k] = v
        except Exception:
            pass  # parâmetros customizados inválidos — usa configuração padrão
        
        # Pause support: if provided, wait before creating output dir / starting work
        pause_event = self.config.get('pause_event')
        if pause_event is not None:
            pause_event.wait()

        # Criar diretório para o modelo
        model_output_dir = self.config['output_folder'] / model_name
        model_output_dir.mkdir(exist_ok=True)
        
        # Nomes dos arquivos
        output_filename = f"{base_name}_{model_name}_results.txt"
        ctl_filename = f"{base_name}_{model_name}.ctl"
        ctl_path = model_output_dir / ctl_filename
        
        # ── Preparar sandbox de execução ────────────────────────────────────────
        # Cria diretório temporário isolado usando tempfile.mkdtemp().
        # • No Linux, usa /dev/shm (tmpfs em RAM) quando disponível → zero I/O de disco.
        # • No Windows/Mac, usa o temp dir padrão do sistema (normalmente SSD NVMe).
        # • Nome único gerado pelo tempfile → sem colisões, sem necessidade de retry.
        temp_dir = Path(tempfile.mkdtemp(
            prefix=f'easypam_{model_name}_{base_name}_',
            dir=self._get_fast_tempdir()
        ))

        try:
            # ── Sanitizar headers / preparar cópia do arquivo de sequência ─────────
            # Para FASTA: CODEML 4.9j tem um limite interno de ~90 chars por linha de
            # header.  Headers mais longos corrompem o parser e causam:
            #   "Error in sequence data file: O at 10 seq 1."
            # Solução: cópia no sandbox com headers truncados ao nome da espécie.
            # Para PHYLIP: arquivo já está no formato correto; cópia direta no sandbox.
            sanitized_fas = temp_dir / fas_file.name
            _seqfile_ref  = str(fas_file.absolute())   # fallback: arquivo original
            try:
                _raw_seq = fas_file.read_text(encoding='utf-8', errors='replace')
                _raw_seq_lines = _raw_seq.splitlines()
                _first_sq = (_raw_seq_lines[0].strip() if _raw_seq_lines else '')
                _seq_is_phylip = (
                    bool(_first_sq)
                    and _first_sq.split()[0].lstrip('-').isdigit()
                    and not _first_sq.startswith('>')
                )
                if _seq_is_phylip:
                    # PHYLIP: copiar sem modificar (formato já adequado para CODEML)
                    sanitized_fas.write_text(_raw_seq, encoding='utf-8')
                else:
                    # FASTA: truncar headers ao primeiro token (nome da espécie)
                    _lines_out: list[str] = []
                    for _fline in _raw_seq_lines:
                        if _fline.startswith('>'):
                            _spname = _fline[1:].split()[0] if _fline[1:].strip() else 'seq'
                            _lines_out.append(f'>{_spname}')
                        else:
                            _lines_out.append(_fline)
                    sanitized_fas.write_text('\n'.join(_lines_out) + '\n', encoding='utf-8')
                _seqfile_ref = str(sanitized_fas)
            except Exception as _san_err:
                with open(log_file, 'a', encoding='utf-8') as _log:
                    _log.write(
                        f"[WARN] {base_name}: preparação do arquivo de sequência falhou: "
                        f"{_san_err}; usando arquivo original\n"
                    )

            # ── Podar árvore para corresponder ao FASTA ──────────────────────────
            # CODEML exige que o número de sequências no FASTA seja igual ao número
            # de taxons declarado no cabeçalho do arquivo de árvore ("N  1").
            # Controlado pelo toggle 'auto_prune_tree' (padrão: True).
            #   1. Identificar taxons da árvore ausentes no FASTA → podar
            #   2. Identificar sequências do FASTA ausentes na árvore → excluir
            #   3. Atualizar o contador na primeira linha do arquivo de árvore
            _pruned_tree_path: Optional[Path] = None
            if not self.config.get('auto_prune_tree', True):
                pass  # Poda desativada pelo usuário — usar árvore original
            else:
                try:
                    from io import StringIO as _SIO
                    from Bio import Phylo as _Phylo

                    # Taxons presentes no arquivo de sequência sanitizado
                    # (suporta FASTA com '>' e PHYLIP sequential com nome nos primeiros 10 chars)
                    _fasta_taxa: set = set()
                    _san_lines = sanitized_fas.read_text(encoding='utf-8', errors='ignore').splitlines()
                    _san_first = _san_lines[0].strip() if _san_lines else ''
                    _san_is_phy = (
                        bool(_san_first)
                        and _san_first.split()[0].lstrip('-').isdigit()
                        and not _san_first.startswith('>')
                    )
                    if _san_is_phy:
                        # PHYLIP: extrair nomes (primeiros 10 chars não-espaço de cada linha de sequência)
                        _phy_ns = int(_san_first.split()[0])
                        for _pln in _san_lines[1:]:
                            if _pln and not _pln.startswith(' ') and len(_fasta_taxa) < _phy_ns:
                                _tx = _pln[:10].strip()
                                if _tx:
                                    _fasta_taxa.add(_tx)
                    else:
                        # FASTA: extrair nomes dos headers '>'
                        for _fln in _san_lines:
                            if _fln.startswith('>'):
                                _tx = _fln[1:].split()[0] if _fln[1:].strip() else ''
                                if _tx:
                                    _fasta_taxa.add(_tx)

                    # Ler e parsear a árvore original
                    _orig_tree_file = Path(self.config['tree_file'])
                    _orig_raw = _orig_tree_file.read_text(encoding='utf-8', errors='ignore')
                    _orig_lines = _orig_raw.splitlines()
                    # Cabeçalho PHYLIP opcional ("N  k") na primeira linha
                    _has_header = (
                        _orig_lines
                        and _orig_lines[0].strip()
                        and _orig_lines[0].strip().split()[0].isdigit()
                    )
                    _nwk_str = '\n'.join(_orig_lines[1:]) if _has_header else _orig_raw

                    _bio_tree = _Phylo.read(_SIO(_nwk_str), 'newick')
                    _tree_taxa: set = {t.name for t in _bio_tree.get_terminals() if t.name}

                    _not_in_tree  = _fasta_taxa - _tree_taxa   # sequências FASTA sem match na árvore
                    _not_in_fasta = _tree_taxa - _fasta_taxa   # taxons da árvore ausentes no FASTA

                    if _not_in_tree:
                        with open(log_file, 'a', encoding='utf-8') as _log:
                            _log.write(
                                f"[WARN] {base_name} [{model_name}]: {len(_not_in_tree)} "
                                f"sequencia(s) do FASTA nao encontrada(s) na arvore (serao "
                                f"excluidas da analise): {sorted(_not_in_tree)}\n"
                            )

                    # Podar árvore: remover taxons ausentes no FASTA
                    for _tx in _not_in_fasta:
                        _bio_tree.prune(_tx)

                    # Filtrar FASTA: manter apenas taxons presentes na árvore
                    _matching = _fasta_taxa & _tree_taxa
                    if _not_in_tree:
                        _raw_san = sanitized_fas.read_text(encoding='utf-8', errors='ignore')
                        _kept_lines: list = []
                        _include = False
                        for _fln in _raw_san.splitlines():
                            if _fln.startswith('>'):
                                _tx = _fln[1:].split()[0] if _fln[1:].strip() else ''
                                _include = _tx in _matching
                            if _include:
                                _kept_lines.append(_fln)
                        sanitized_fas.write_text('\n'.join(_kept_lines) + '\n', encoding='utf-8')

                    # Escrever árvore podada no sandbox
                    _pnwk_io = _SIO()
                    _Phylo.write(_bio_tree, _pnwk_io, 'newick')
                    _pnwk = _pnwk_io.getvalue().strip()
                    # Bio.Phylo adiciona branch length no nó raiz (ex: "...):0.00000;")
                    # que CODEML não aceita — remover esse artefato
                    _pnwk = re.sub(r'\):[0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?;$', ');', _pnwk)
                    _n_match = _bio_tree.count_terminals()
                    _pruned_content = f"{_n_match}  1\n{_pnwk}\n"
                    _pruned_tree_path = temp_dir / 'pruned_tree.nwk'
                    _pruned_tree_path.write_text(_pruned_content, encoding='utf-8')

                except ImportError:
                    pass  # Biopython indisponível; CODEML pode falhar com contagem diferente
                except Exception as _prune_err:
                    with open(log_file, 'a', encoding='utf-8') as _log:
                        _log.write(
                            f"[WARN] {base_name} [{model_name}]: poda automatica da arvore "
                            f"falhou: {_prune_err}\n"
                        )

            # ── Desenraizar árvore para site models ──────────────────────────────
            # O guia do PAML é explícito: site models (M0, M1a, M2a, M7, M8)
            # exigem árvore NÃO-enraizada. Uma árvore enraizada tem o nó raiz
            # com exatamente 2 filhos (bifurcação). Desraizamos colapsando um
            # dos filhos do root para criar uma tricotomia no root.
            # Branch/Branch-site são excluídos pois exigem árvore enraizada com
            # marcação de ramo.
            _SITE_MODELS_UNROOT = {'M0', 'M1a', 'M2a', 'M7', 'M8'}
            if model_name in _SITE_MODELS_UNROOT:
                try:
                    from io import StringIO as _SIO_u
                    from Bio import Phylo as _Phylo_u

                    # Ler a árvore que está sendo usada (podada ou original)
                    _usrc = _pruned_tree_path if _pruned_tree_path is not None \
                            else Path(self.config['tree_file'])
                    _uraw  = _usrc.read_text(encoding='utf-8', errors='ignore')
                    _ulines = _uraw.splitlines()
                    _uhdr   = (
                        _ulines
                        and _ulines[0].strip()
                        and _ulines[0].strip().split()[0].isdigit()
                    )
                    _unwk = '\n'.join(_ulines[1:]) if _uhdr else _uraw

                    _utree = _Phylo_u.read(_SIO_u(_unwk), 'newick')

                    # Enraizada ↔ root com exatamente 2 filhos diretos
                    if len(_utree.root.clades) == 2:
                        _uc0, _uc1 = _utree.root.clades
                        if _uc1.clades:
                            # _uc1 é nó interno → elevar seus filhos ao root
                            _utree.root.clades = [_uc0] + _uc1.clades
                        elif _uc0.clades:
                            # _uc0 é nó interno → elevar seus filhos ao root
                            _utree.root.clades = _uc0.clades + [_uc1]
                        # (se ambos forem folhas não há como desraizar — ignorar)

                        _uio = _SIO_u()
                        _Phylo_u.write(_utree, _uio, 'newick')
                        _upnwk = _uio.getvalue().strip()
                        # Remover artefato de branch length no root gerado pelo Bio.Phylo
                        _upnwk = re.sub(
                            r'\):[0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?;$', ');', _upnwk
                        )
                        _un = _utree.count_terminals()
                        _unrooted_path = temp_dir / 'unrooted_tree.nwk'
                        _unrooted_path.write_text(
                            f"{_un}  1\n{_upnwk}\n", encoding='utf-8'
                        )
                        _pruned_tree_path = _unrooted_path   # substituir referência
                except ImportError:
                    pass  # Bio.Phylo indisponível; prosseguir com árvore original
                except Exception as _unroot_err:
                    with open(log_file, 'a', encoding='utf-8') as _log:
                        _log.write(
                            f"[WARN] {base_name} [{model_name}]: desenraizamento "
                            f"automatico falhou: {_unroot_err}\n"
                        )

            # ── Determinar conteúdo da árvore para este modelo ────────────────
            labeled_full      = self.config.get('labeled_tree_content')
            labeled_branchsite = self.config.get('labeled_tree_branchsite')

            # Suporte ao nome antigo (BranchSite*) e novo (Branch-site*)
            if model_name.startswith('BranchSite') or model_name.startswith('Branch-site'):
                labeled_content = labeled_branchsite or (
                    labeled_full if labeled_full and '#1' in labeled_full else None)
            elif model_name == 'Branch':
                labeled_content = labeled_full
            else:
                labeled_content = None

            # ── Decidir treefile no .ctl e o que escrever no sandbox ──────────
            # Prioridade:
            #   1. Árvore marcada pelo usuário (labeled_content) — obrigatória para Branch/Branch-site
            #   2. Árvore ajustada do M0 (fitted_tree) — warm-start via fix_blength=2
            #      (apenas para modelos de sítios; Branch exige marcação específica)
            #   3. Árvore original do usuário — caminho absoluto; nenhuma cópia necessária
            if labeled_content:
                # Escrever árvore marcada no sandbox; .ctl referencia pelo nome relativo
                (temp_dir / 'labeled.nwk').write_text(labeled_content, encoding='utf-8')
                tree_ref  = 'labeled.nwk'   # relativo ao CWD (temp_dir)
                fix_bl    = 0                # branch lengths estimados normalmente
            elif (fitted_tree
                  and not model_name.startswith('Branch')
                  and model_name != 'M0'):
                # Warm-start: escrever árvore M0 no sandbox; usar fix_blength=2
                # Os branch lengths serão RE-ESTIMADOS livremente — sem impacto nos resultados.
                (temp_dir / 'warm_tree.nwk').write_text(fitted_tree, encoding='utf-8')
                tree_ref  = 'warm_tree.nwk' # relativo ao CWD (temp_dir)
                fix_bl    = 2                # inicializar a partir dos valores do M0
            elif _pruned_tree_path is not None:
                # Árvore podada/desenraizada — referenciada pelo nome do arquivo no sandbox
                tree_ref  = _pruned_tree_path.name   # ex: 'pruned_tree.nwk' ou 'unrooted_tree.nwk'
                fix_bl    = 0
            else:
                # Árvore original referenciada por caminho absoluto → sem cópia
                tree_ref  = str(Path(self.config['tree_file']).absolute())
                fix_bl    = 0

            # Nota: seqfile agora aponta para a cópia sanitizada no sandbox.

            # ── Gerar conteúdo do .ctl ────────────────────────────────────────
            custom_paths  = self.config.get('model_ctl_paths', {}) or {}
            provided_ctl  = custom_paths.get(model_name)
            omega_initial = float(self.config.get('omega', model_config.get('omega', 0.5) or 0.5))
            cleandata_val = int(self.config.get('cleandata', 1))

            if provided_ctl:
                provided_path = Path(provided_ctl)
                if provided_path.exists() and provided_path.is_file():
                    raw = provided_path.read_text(encoding='utf-8')

                    def _replace_setting(content: str, key: str, newval: str) -> str:
                        pat  = rf'(^\s*{re.escape(key)}\s*=).*?$'
                        repl = rf"\1 {newval}"
                        return re.sub(pat, repl, content, flags=re.MULTILINE)

                    ctl_content = _replace_setting(raw,         'seqfile', _seqfile_ref)
                    ctl_content = _replace_setting(ctl_content, 'treefile', tree_ref)
                    ctl_content = _replace_setting(ctl_content, 'outfile',  output_filename)
                else:
                    with open(log_file, 'a', encoding='utf-8') as log:
                        log.write(f"Warning: provided .ctl for {model_name} not found: "
                                  f"{provided_ctl}; generating default .ctl\n")
                    ctl_content = self.generate_ctl_content(
                        seqfile=_seqfile_ref,
                        treefile=tree_ref,
                        outfile=output_filename,
                        model_config=model_config,
                        omega=omega_initial,
                        cleandata=cleandata_val,
                        model_name=model_name,
                        kappa=warm_start_kappa,
                        fix_kappa_heuristic=fix_kappa_heuristic,
                        fix_blength=fix_bl,
                    )
            else:
                ctl_content = self.generate_ctl_content(
                    seqfile=_seqfile_ref,
                    treefile=tree_ref,
                    outfile=output_filename,
                    model_config=model_config,
                    omega=omega_initial,
                    cleandata=cleandata_val,
                    model_name=model_name,
                    kappa=warm_start_kappa,
                    fix_kappa_heuristic=fix_kappa_heuristic,
                    fix_blength=fix_bl,
                )

            # Escrever .ctl no sandbox (para o CODEML) e em model_output_dir (para referência)
            ctl_in_sandbox = temp_dir / ctl_filename
            ctl_in_sandbox.write_text(ctl_content, encoding='utf-8')
            try:
                ctl_path.write_text(ctl_content, encoding='utf-8')
            except Exception:
                pass  # falha ao salvar referência não impede a execução

            exec_start = time.time()

            # Executar CODEML — use absolute bundled binary, fall back to system PATH
            if _CODEML_BIN.exists():
                cmd = [str(_CODEML_BIN), ctl_filename]
            else:
                cmd = ["codeml", ctl_filename]

            with open(log_file, 'a', encoding='utf-8') as log:
                log.write(f"[{model_name}] {base_name}: Running command: {cmd} in {temp_dir}\n")

            process = subprocess.Popen(
                cmd,
                cwd=temp_dir,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                encoding='utf-8',
                bufsize=1
            )

            # Capturar output
            stdout_lines = []
            stderr_lines = []
            # stop codon tracking
            stop_count = [0]
            stop_details = []

            def read_stream(stream, storage):
                try:
                    for line in iter(stream.readline, ''):
                        if line is None:
                            continue
                        text_line = line.rstrip()
                        storage.append(text_line)

                        # Detect possible stop-codon prompt or pause requiring Enter
                        lower = text_line.lower()
                        is_stop_codon = 'stop' in lower and 'codon' in lower
                        is_press_enter = 'press' in lower and 'enter' in lower

                        if is_stop_codon or is_press_enter:
                            with open(log_file, 'a', encoding='utf-8') as log:
                                log.write(f"[{model_name}] {base_name}: Detected prompt line: {text_line}\n")

                            # increment stop counter and try to parse details
                            stop_count[0] += 1
                            # update object-level counter for GUI polling
                            try:
                                self.current_stop_count = stop_count[0]
                            except Exception:
                                pass
                            # try regex: stop codon TAG in seq. #   1 (...), nucleotide site 214
                            m = re.search(r"stop codon\s+(\w+)\s+in seq\.\s*#\s*(\d+).*?site\s*(\d+)", text_line, re.IGNORECASE)
                            if m:
                                codon = m.group(1)
                                seqnum = m.group(2)
                                site = m.group(3)
                                stop_details.append({'codon': codon, 'seqnum': int(seqnum), 'site': int(site), 'line': text_line})
                            else:
                                stop_details.append({'line': text_line})
                            # publish details to object for GUI polling
                            try:
                                self.current_stop_details = list(stop_details)
                            except Exception:
                                pass
                        

                            # Notify user (prints are redirected to GUI when used from GUI)
                            try:
                                print(f"[{model_name}] {base_name}: Detected stop codon (count={stop_count[0]}).")
                            except Exception:
                                pass

                            auto_continue = bool(self.config.get('auto_continue_stop_codons', False))
                            n_workers = int(self.config.get('n_workers', 1))

                            # In parallel mode, multiple workers share the same
                            # manual_continue_event.  Blocking on a shared event
                            # causes all-but-one worker to deadlock → genes appear
                            # skipped.  Force auto-continue whenever n_workers > 1.
                            if n_workers > 1 and not auto_continue:
                                auto_continue = True
                                with open(log_file, 'a', encoding='utf-8') as log:
                                    log.write(
                                        f"[{model_name}] {base_name}: Parallel mode — "
                                        f"forced auto-continue for stop codon.\n"
                                    )

                            manual_all_ev = self.config.get('manual_continue_all_event')

                            def _send_enter(reason: str) -> None:
                                """Write a newline to the subprocess stdin reliably on Windows."""
                                try:
                                    if process.stdin:
                                        # Use the binary buffer when available so the
                                        # write bypasses TextIOWrapper internal buffering.
                                        if hasattr(process.stdin, 'buffer'):
                                            process.stdin.buffer.write(b'\n')
                                            process.stdin.buffer.flush()
                                        else:
                                            process.stdin.write('\n')
                                            process.stdin.flush()
                                    with open(log_file, 'a', encoding='utf-8') as log:
                                        log.write(f"[{model_name}] {base_name}: {reason}\n")
                                except Exception as exc:
                                    with open(log_file, 'a', encoding='utf-8') as log:
                                        log.write(
                                            f"[{model_name}] {base_name}: "
                                            f"Failed to send Enter ({reason}): {exc}\n"
                                        )

                            # If user requested auto-continue, or a manual-all event is set, send Enter
                            if auto_continue:
                                _send_enter("Auto-sent Enter to subprocess (stop codon).")
                            elif manual_all_ev is not None and getattr(manual_all_ev, 'is_set') and manual_all_ev.is_set():
                                _send_enter("manual-all event set; sent Enter.")
                            else:
                                # Wait for GUI/user to signal continuation via event
                                event = self.config.get('manual_continue_event')
                                if event is None:
                                    # no event provided -> fallback to auto
                                    _send_enter("No manual event provided; auto-sent Enter.")
                                else:
                                    with open(log_file, 'a', encoding='utf-8') as log:
                                        log.write(f"[{model_name}] {base_name}: Waiting for manual continue event...\n")
                                    # Wait until GUI sets the event
                                    event.wait()
                                    # clear event for next prompt
                                    try:
                                        event.clear()
                                    except Exception:
                                        pass
                                    _send_enter("Manual continue event received; sent Enter.")
                except Exception:
                    pass

            stdout_thread = Thread(target=read_stream, args=(process.stdout, stdout_lines))
            stderr_thread = Thread(target=read_stream, args=(process.stderr, stderr_lines))

            stdout_thread.start()
            stderr_thread.start()

            # Registrar processo como ativo (stop imediato e monitoramento paralelo)
            with self._processes_lock:
                self._active_processes.add(process)
                self.current_process = process

            stop_event  = self.config.get('stop_event')
            pause_event = self.config.get('pause_event')
            timeout_s   = self.config.get('timeout', 1600)
            deadline    = time.time() + timeout_s
            timed_out   = False
            stopped     = False

            try:
                while process.poll() is None:
                    # Verificar stop imediato
                    if stop_event is not None and stop_event.is_set():
                        stopped = True
                        try:
                            process.terminate()
                            process.wait(timeout=3)
                        except Exception:
                            try:
                                process.kill()
                            except Exception:
                                pass
                        break
                    # Verificar timeout
                    if time.time() > deadline:
                        timed_out = True
                        try:
                            process.kill()
                        except Exception:
                            pass
                        break
                    # Aguardar; checar a cada 0.5 s para responsividade
                    time.sleep(0.5)
            finally:
                with self._processes_lock:
                    self._active_processes.discard(process)
                    if self.current_process is process:
                        self.current_process = None

            stdout_thread.join(timeout=1)
            stderr_thread.join(timeout=1)

            if timed_out:
                with open(log_file, 'a', encoding='utf-8') as log:
                    log.write(f"[{model_name}] {base_name}: TIMEOUT after {timeout_s}s\n")
                    log.write(f"  Captured stdout (last 200 lines):\n")
                    for L in stdout_lines[-200:]:
                        log.write(L + "\n")
                    log.write(f"  Captured stderr (last 200 lines):\n")
                    for L in stderr_lines[-200:]:
                        log.write(L + "\n")
                return None

            if stopped:
                with open(log_file, 'a', encoding='utf-8') as log:
                    log.write(f"[{model_name}] {base_name}: STOPPED by user\n")
                return None

            # Wait for reader threads to finish
            stdout_thread.join(timeout=1)
            stderr_thread.join(timeout=1)

            # Close streams to release file handles on Windows
            try:
                if process.stdout:
                    process.stdout.close()
            except Exception:
                pass
            try:
                if process.stderr:
                    process.stderr.close()
            except Exception:
                pass
            try:
                if process.stdin:
                    process.stdin.close()
            except Exception:
                pass

            rc = process.returncode
            with open(log_file, 'a', encoding='utf-8') as log:
                log.write(f"[{model_name}] {base_name}: process returncode={rc}\n")
                if stdout_lines:
                    log.write(f"  stdout (last 200 lines):\n")
                    for L in stdout_lines[-200:]:
                        log.write(L + "\n")
                if stderr_lines:
                    log.write(f"  stderr (last 200 lines):\n")
                    for L in stderr_lines[-200:]:
                        log.write(L + "\n")

            if rc != 0:
                with open(log_file, 'a', encoding='utf-8') as log:
                    log.write(f"[{model_name}] {base_name}: Non-zero return code {rc}\n")
                # continue to attempt to find outputs

            # Mover arquivos de saída
            for src_file in temp_dir.glob("*"):
                # move expected outputs (output file, .rst, .txt), leave input files
                if src_file.name == ctl_filename:
                    continue
                if src_file.name == output_filename or src_file.suffix in ['.rst', '.txt']:
                    dest = model_output_dir / src_file.name
                    try:
                        shutil.move(src_file, dest)
                    except Exception as e:
                        with open(log_file, 'a', encoding='utf-8') as log:
                            log.write(f"Failed to move {src_file} -> {dest}: {e}\n")

            output_path = model_output_dir / output_filename

            # Extrair informações
            lnL = None
            np_params = None
            omega = None

            if output_path.exists():
                lnL = self._extract_likelihood(output_path)
                np_params = self._extract_np(output_path)
                omega = self._extract_omega(output_path)
            else:
                with open(log_file, 'a', encoding='utf-8') as log:
                    log.write(f"[{model_name}] {base_name}: expected output file not found: {output_path}\n")

            execution_time = time.time() - exec_start

            # Log sucesso (ou parcial) including stop count
            with open(log_file, 'a', encoding='utf-8') as log:
                log.write(f"[{model_name}] {base_name}: FINISHED (lnL={lnL}, np={np_params}, ω={omega}, time={execution_time:.1f}s, stop_count={stop_count[0]})\n")
                if stop_details:
                    log.write(f"  Stop details:\n")
                    for d in stop_details:
                        log.write(f"    {d}\n")

            return {
                'output_file': str(output_path) if output_path.exists() else None,
                'results_file': str(output_path) if output_path.exists() else None,  # Alias para compatibilidade
                'lnL': lnL,
                'np': np_params,
                'omega': omega,
                'execution_time': execution_time,
                'status': 'success' if rc == 0 and output_path.exists() else 'partial',
                'stop_count': stop_count[0]
            }

        except Exception as e:
            tb = traceback.format_exc()
            with open(log_file, 'a', encoding='utf-8') as log:
                log.write(f"[{model_name}] {base_name}: EXCEPTION - {e}\n")
                log.write(tb + "\n")
            return None

        finally:
            # Final cleanup: try to remove temp_dir with retries; if fails, log and continue
            if temp_dir.exists():
                for attempt in range(8):
                    try:
                        shutil.rmtree(temp_dir)
                        break
                    except PermissionError as pe:
                        with open(log_file, 'a', encoding='utf-8') as log:
                            log.write(f"Cleanup: could not remove {temp_dir} (attempt {attempt+1}/8): {pe}\n")
                        time.sleep(0.5)
                    except Exception as e:
                        with open(log_file, 'a', encoding='utf-8') as log:
                            log.write(f"Cleanup: unexpected error removing {temp_dir}: {e}\n")
                        break
                else:
                    with open(log_file, 'a', encoding='utf-8') as log:
                        log.write(f"Cleanup: failed to remove {temp_dir} after retries; leaving it in place.\n")
    
    def _extract_likelihood(self, output_file: Path) -> Optional[float]:
        """Extrai log-likelihood"""
        try:
            with open(output_file, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    if 'lnL' in line:
                        # Try a few regex patterns to capture common CODEML formats
                        patterns = [
                            r'lnL[^:]*:\s*([+-]?\d+\.\d+)',
                            r'lnL\([^)]*\):\s*([+-]?\d+\.\d+)',
                            r'lnL\s*[:=]\s*([+-]?\d+\.\d+)'
                        ]
                        for pat in patterns:
                            match = re.search(pat, line)
                            if match:
                                try:
                                    return float(match.group(1))
                                except Exception:
                                    continue
        except Exception:
            pass
        return None
    
    def _extract_np(self, output_file: Path) -> Optional[int]:
        """Extrai número de parâmetros"""
        try:
            with open(output_file, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    if 'lnL' in line and 'np:' in line:
                        match = re.search(r'np:\s*(\d+)', line)
                        if match:
                            return int(match.group(1))
        except Exception:
            pass
        return None

    def _extract_omega(self, output_file: Path) -> Optional[float]:
        """Extrai omega usando SitesParser (suporta Branch/Branch-Site/Site models)"""
        try:
            omega = SitesParser.extract_omega_robust(output_file)
            return omega
        except Exception:
            return None
    
    def _save_summary(self):
        """Salva sumário em TSV com colunas de LRT e omegas extraídos robustamente"""
        summary_file = self.config['output_folder'] / "analysis_summary.tsv"
        
        with open(summary_file, 'w', encoding='utf-8') as f:
            # Header
            header = ["Gene"]
            for model in self.config['models']:
                header.extend([f"{model}_lnL", f"{model}_np", f"{model}_omega", f"{model}_time", f"{model}_stops"])
            
            # Adicionar colunas de LRT
            selected_models = self.config['models']
            if 'M0' in selected_models and 'M1a' in selected_models:
                header.append("lrt_M0_vs_M1a")
            if 'M1a' in selected_models and 'M2a' in selected_models:
                header.append("lrt_M1a_vs_M2a")
            if 'M7' in selected_models and 'M8' in selected_models:
                header.append("lrt_M7_vs_M8")
            if 'M0' in selected_models and 'Branch' in selected_models:
                header.append("lrt_M0_vs_Branch")
            # Support both old name (BranchSite_A) and new name (Branch-site)
            if ('BranchSite_A_null' in selected_models and 'BranchSite_A' in selected_models) or \
               ('Branch-site_null' in selected_models and 'Branch-site' in selected_models):
                if 'Branch-site_null' in selected_models and 'Branch-site' in selected_models:
                    header.append("lrt_Branch-site_null_vs_Branch-site")
                else:
                    header.append("lrt_BranchSite_A_null_vs_BranchSite_A")
            
            f.write("\t".join(header) + "\n")
            
            # Data
            for gene_name in sorted(self.results.keys()):
                gene_results = self.results[gene_name]
                # Remover caracteres que corrompem o formato TSV
                row = [str(gene_name).replace('\n', '').replace('\r', '').replace('\t', '_')]
                
                for model in self.config['models']:
                    if model in gene_results and gene_results[model]:
                        result = gene_results[model]
                        
                        # Extrair omega robustamente do arquivo de resultados
                        omega_value = result.get('omega')
                        if omega_value is None or omega_value == 'NA':
                            # Tentar extrair do arquivo de resultados
                            results_file = result.get('results_file')
                            if results_file:
                                try:
                                    from pathlib import Path
                                    omega_value = SitesParser.extract_omega_robust(Path(results_file))
                                except Exception:
                                    omega_value = None
                        
                        row.extend([
                            f"{result.get('lnL', 'NA'):.6f}" if result.get('lnL') else 'NA',
                            str(result.get('np', 'NA')),
                            f"{omega_value:.6f}" if omega_value is not None and omega_value != 'NA' else 'NA',
                            f"{result.get('execution_time', 0):.2f}",
                            str(result.get('stop_count', 0))
                        ])
                    else:
                        row.extend(['NA', 'NA', 'NA', 'NA', '0'])
                
                # Calcular LRTs para este gene
                lrt_comparisons = []
                if 'M0' in selected_models and 'M1a' in selected_models:
                    lrt_comparisons.append(('M0', 'M1a', 'lrt_M0_vs_M1a'))
                if 'M1a' in selected_models and 'M2a' in selected_models:
                    lrt_comparisons.append(('M1a', 'M2a', 'lrt_M1a_vs_M2a'))
                if 'M7' in selected_models and 'M8' in selected_models:
                    lrt_comparisons.append(('M7', 'M8', 'lrt_M7_vs_M8'))
                if 'M0' in selected_models and 'Branch' in selected_models:
                    lrt_comparisons.append(('M0', 'Branch', 'lrt_M0_vs_Branch'))
                # Support both old name (BranchSite_A) and new name (Branch-site)
                if ('BranchSite_A_null' in selected_models and 'BranchSite_A' in selected_models) or \
                   ('Branch-site_null' in selected_models and 'Branch-site' in selected_models):
                    if 'Branch-site_null' in selected_models and 'Branch-site' in selected_models:
                        lrt_comparisons.append(('Branch-site_null', 'Branch-site', 'lrt_Branch-site_null_vs_Branch-site'))
                    else:
                        lrt_comparisons.append(('BranchSite_A_null', 'BranchSite_A', 'lrt_BranchSite_A_null_vs_BranchSite_A'))
                
                for null_model, alt_model, _ in lrt_comparisons:
                    if (null_model in gene_results and gene_results[null_model] and 
                        alt_model in gene_results and gene_results[alt_model]):
                        null_lnL = gene_results[null_model].get('lnL')
                        alt_lnL = gene_results[alt_model].get('lnL')
                        if null_lnL is not None and alt_lnL is not None:
                            lrt_stat = 2 * (alt_lnL - null_lnL)
                            row.append(f"{lrt_stat:.6f}")
                        else:
                            row.append('NA')
                    else:
                        row.append('NA')
                
                # Limpar newlines de todos os valores antes de escrever
                row = [str(v).replace('\n', '').replace('\r', '') for v in row]
                f.write("\t".join(row) + "\n")
        
        print(f"  [OK] Summary saved: {summary_file}")
    
    def _run_lrt_analysis(self):
        """Executa Likelihood Ratio Tests"""
        
        lrt_file = self.config['output_folder'] / "LRT_results.txt"
        
        with open(lrt_file, 'w', encoding='utf-8') as f:
            f.write("="*80 + "\n")
            f.write("LIKELIHOOD RATIO TEST (LRT) RESULTS\n")
            f.write("="*80 + "\n\n")
            
            # Determinar comparações relevantes
            comparisons = []
            selected_models = self.config['models']
            
            # Site models comparisons
            if 'M0' in selected_models and 'M1a' in selected_models:
                comparisons.append(('M0', 'M1a', 'Tests if ω varies among sites'))
            if 'M1a' in selected_models and 'M2a' in selected_models:
                comparisons.append(('M1a', 'M2a', 'Tests for positive selection'))
            if 'M7' in selected_models and 'M8' in selected_models:
                comparisons.append(('M7', 'M8', 'Alternative test for positive selection'))
            
            # Branch models
            if 'M0' in selected_models and 'Branch' in selected_models:
                comparisons.append(('M0', 'Branch', 'Tests if ω differs in foreground'))
            
            # Branch-site models
            if 'Branch-site_null' in selected_models and 'Branch-site' in selected_models:
                comparisons.append(('Branch-site_null', 'Branch-site',
                                  'Tests for positive selection in foreground sites (50:50 mixture χ²)'))
            
            if not comparisons:
                f.write("No valid model comparisons found.\n")
                f.write("For LRT, you need pairs of nested models.\n")
                print("  [WARN] No valid LRT comparisons found")
                return
            
            print(f"\n  Running {len(comparisons)} LRT comparison(s):\n")
            
            # Realizar cada comparação
            for null_model, alt_model, description in comparisons:
                print(f"    • {null_model} vs {alt_model}")
                
                f.write("\n" + "="*80 + "\n")
                f.write(f"COMPARISON: {null_model} (null) vs {alt_model} (alternative)\n")
                f.write(f"Description: {description}\n")
                f.write("="*80 + "\n\n")
                
                sig_count_05 = 0
                sig_count_01 = 0
                total_valid = 0
                
                # Comparar cada gene
                for gene_name in sorted(self.results.keys()):
                    gene_results = self.results[gene_name]
                    
                    if null_model not in gene_results or alt_model not in gene_results:
                        continue
                    
                    null_res = gene_results[null_model]
                    alt_res = gene_results[alt_model]
                    
                    if not null_res or not alt_res:
                        continue
                    
                    lnL_null = null_res.get('lnL')
                    lnL_alt = alt_res.get('lnL')
                    np_null = null_res.get('np')
                    np_alt = alt_res.get('np')
                    
                    if lnL_null is None or lnL_alt is None:
                        continue
                    
                    # Calcular LRT
                    lrt_stat = 2 * (lnL_alt - lnL_null)
                    df = abs(np_alt - np_null)

                    if df == 0:
                        continue
                    if lrt_stat < 0:
                        lrt_stat = 0.0

                    # Calcular p-value
                    is_branchsite = (null_model == 'Branch-site_null' and alt_model == 'Branch-site')

                    if is_branchsite:
                        # Distribuição nula: mistura 50:50 de χ²(0) e χ²(1)
                        # P(2Δl > x) = 0.5 * P(χ²(1) > x)  para x > 0
                        # Valor crítico α=0.05: 2.706  (qchisq(0.90, df=1))
                        # Valor crítico α=0.01: 5.412  (qchisq(0.98, df=1))
                        if lrt_stat <= 0:
                            p_value = 1.0
                        else:
                            p_value = 0.5 * stats.chi2.sf(lrt_stat, df=1)
                        df_display = "mixture(0,1)"
                    else:
                        p_value = 1 - stats.chi2.cdf(lrt_stat, df)
                        df_display = str(df)
                    
                    total_valid += 1
                    
                    if p_value < 0.05:
                        sig_count_05 += 1
                    if p_value < 0.01:
                        sig_count_01 += 1
                    
                    # Escrever resultado
                    f.write(f"Gene: {gene_name}\n")
                    f.write(f"  lnL {null_model}: {lnL_null:.6f} (np={np_null})\n")
                    f.write(f"  lnL {alt_model}: {lnL_alt:.6f} (np={np_alt})\n")
                    f.write(f"  2Δl = {lrt_stat:.6f}\n")
                    f.write(f"  df = {df_display}\n")
                    f.write(f"  p-value = {p_value:.6e}\n")
                    
                    if p_value < 0.01:
                        f.write(f"  Result: [OK][OK] {alt_model} significantly better (p < 0.01)\n")
                    elif p_value < 0.05:
                        f.write(f"  Result: [OK] {alt_model} significantly better (p < 0.05)\n")
                    else:
                        f.write(f"  Result: [ERROR] No significant difference\n")
                    
                    f.write("\n" + "-"*60 + "\n\n")
                
                # Sumário da comparação
                f.write("\nRESUMO:\n")
                f.write(f"  Total de genes analisados: {total_valid}\n")
                if total_valid > 0:
                    f.write(f"  Significativo em p < 0.05: {sig_count_05} ({100*sig_count_05/total_valid:.1f}%)\n")
                    f.write(f"  Significativo em p < 0.01: {sig_count_01} ({100*sig_count_01/total_valid:.1f}%)\n")
                else:
                    f.write(f"  Significativo em p < 0.05: {sig_count_05}\n")
                    f.write(f"  Significativo em p < 0.01: {sig_count_01}\n")
                f.write("\n")
        
        print(f"\n  [OK] LRT results saved: {lrt_file}")

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

    def run_wgs_analysis(self) -> None:
        """Modo WGS: combina todos os .fas em um único arquivo PHYLIP e executa
        CODEML com  ndata = N maintree 1  para analisar todos os genes de uma vez.

        Suporta apenas modelos de sítios (site models) sem marcação de ramos.
        Para modelos de ramo use o modo batch padrão.

        Configuração esperada em self.config (além das chaves padrão):
          - 'wgs_nsites': lista/string de NSsites, ex: [0, 1, 2, 7, 8]
          - 'wgs_model' : valor de model=, default 0 (site models)
        """
        output_folder = Path(self.config['output_folder'])
        output_folder.mkdir(parents=True, exist_ok=True)
        log_file = output_folder / "wgs_analysis_log.txt"

        _wgs_input = Path(self.config['input_folder'])
        fas_files = sorted(
            list(_wgs_input.glob("*.fas"))
            + list(_wgs_input.glob("*.fasta"))
            + list(_wgs_input.glob("*.phy"))
            + list(_wgs_input.glob("*.phylip")),
            key=lambda p: p.name.lower()
        )
        if not fas_files:
            print("[WGS] Nenhum arquivo .fas / .fasta / .phy / .phylip encontrado.")
            return

        print(f"\n[WGS] Convertendo {len(fas_files)} genes para PHYLIP combinado...")

        # Combinar todos em um único arquivo
        combined_phy = output_folder / "wgs_combined.phy"
        gene_names = []
        n_valid = 0
        with open(combined_phy, 'w', encoding='utf-8') as out:
            for fas in fas_files:
                block = self._fasta_to_phylip_block(fas)
                if block is None:
                    print(f"  [SKIP] {fas.name}")
                    continue
                out.write(block)
                out.write('\n')
                gene_names.append(fas.stem)
                n_valid += 1

        if n_valid == 0:
            print("[WGS] Nenhum arquivo válido para converter.")
            return
        print(f"[WGS] {n_valid} genes combinados → {combined_phy.name}")

        # Salvar lista de genes na ordem
        (output_folder / "wgs_gene_order.txt").write_text('\n'.join(gene_names), encoding='utf-8')

        # Construir NSsites
        nsites_raw = self.config.get('wgs_nsites', [0, 1, 2, 7, 8])
        if isinstance(nsites_raw, (list, tuple)):
            nsites_str = ' '.join(str(n) for n in nsites_raw)
        else:
            nsites_str = str(nsites_raw)

        model_val = int(self.config.get('wgs_model', 0))
        omega_val = float(self.config.get('omega', 0.5))
        cleandata_val = int(self.config.get('cleandata', 0))
        codon_freq = int(self.config.get('wgs_codonfreq', 7))

        tree_path = Path(self.config['tree_file'])
        outfile_name = "wgs_results.txt"

        ctl_content = (
            f"      seqfile = {combined_phy.name}\n"
            f"     treefile = {tree_path.name}\n"
            f"      outfile = {outfile_name}\n\n"
            f"        noisy = 1\n"
            f"      verbose = 1\n"
            f"      seqtype = 1\n"
            f"        ndata = {n_valid} maintree 1\n"
            f"        icode = 0\n"
            f"    cleandata = {cleandata_val}\n\n"
            f"        model = {model_val}\n"
            f"      NSsites = {nsites_str}\n"
            f"    CodonFreq = {codon_freq}\n"
            f"      estFreq = 0\n"
            f"        clock = 0\n"
            f"    fix_omega = 0\n"
            f"        omega = {omega_val}\n"
        )

        ctl_path = output_folder / "wgs_analysis.ctl"
        ctl_path.write_text(ctl_content, encoding='utf-8')
        print(f"[WGS] .ctl gerado: {ctl_path.name}")
        print(f"[WGS] ndata = {n_valid} maintree 1  |  NSsites = {nsites_str}")

        # Copiar arquivos para temp_dir e executar
        temp_dir = output_folder / "wgs_temp"
        if temp_dir.exists():
            shutil.rmtree(temp_dir, ignore_errors=True)
        temp_dir.mkdir()

        shutil.copy(combined_phy, temp_dir)
        shutil.copy(tree_path, temp_dir)
        shutil.copy(ctl_path, temp_dir)

        cmd = [str(_CODEML_BIN), ctl_path.name] if _CODEML_BIN.exists() else ['codeml', ctl_path.name]

        print(f"\n[WGS] Executando: {' '.join(cmd)}")
        print(f"[WGS] Isso pode demorar muito para grandes datasets WGS...\n")

        stop_event = self.config.get('stop_event')
        with open(log_file, 'w', encoding='utf-8') as log:
            log.write(f"WGS Analysis started: {datetime.now()}\n")
            log.write(f"ndata = {n_valid}  NSsites = {nsites_str}\n")
            log.write(f"Genes: {', '.join(gene_names)}\n\n")

        try:
            process = subprocess.Popen(
                cmd, cwd=temp_dir,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True, encoding='utf-8', bufsize=1
            )
            self.current_process = process

            with open(log_file, 'a', encoding='utf-8') as log:
                for line in iter(process.stdout.readline, ''):
                    if stop_event and stop_event.is_set():
                        process.terminate()
                        break
                    stripped = line.rstrip()
                    log.write(stripped + '\n')
                    if stripped:
                        print(f"  {stripped}")
                    # Auto-responder stop codons
                    if 'stop' in stripped.lower() and 'codon' in stripped.lower():
                        try:
                            process.stdin.write('\n')
                            process.stdin.flush()
                        except Exception:
                            pass

            process.wait(timeout=7200)

        except subprocess.TimeoutExpired:
            process.kill()
            print("[WGS] TIMEOUT após 2h — processo encerrado.")
        except Exception as e:
            print(f"[WGS] ERRO: {e}")
        finally:
            self.current_process = None

        # Mover resultados
        results_file = temp_dir / outfile_name
        if results_file.exists():
            dest = output_folder / outfile_name
            shutil.copy(results_file, dest)
            print(f"\n[WGS] Resultado salvo: {dest}")
        else:
            print("[WGS] Arquivo de resultado não encontrado.")

        # Limpeza
        try:
            shutil.rmtree(temp_dir, ignore_errors=True)
        except Exception:
            pass

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
            # ═══ 1. REGENERAR analysis_summary.tsv ═══
            print("\n[1/3] Generating analysis_summary.tsv...")
            summary_file = CodemlBatchAnalysis._regenerate_analysis_summary(results_folder)
            if summary_file:
                generated_files['analysis_summary'] = str(summary_file)
                print(f"  OK: {summary_file.name}")
            
            # ═══ 2. REGENERAR batch_analysis_log.txt ═══
            print("\n[2/3] Generating batch_analysis_log.txt...")
            log_file = CodemlBatchAnalysis._regenerate_batch_log(results_folder)
            if log_file:
                generated_files['batch_analysis_log'] = str(log_file)
                print(f"  OK: {log_file.name}")
            
            # ═══ 3. REGENERAR LRT_results.txt ═══
            print("\n[3/3] Generating LRT_results.txt...")
            lrt_file = CodemlBatchAnalysis._regenerate_lrt_results(results_folder)
            if lrt_file:
                generated_files['LRT_results'] = str(lrt_file)
                print(f"  OK: {lrt_file.name}")
            
            print(f"\n[SUCCESS] All files regenerated successfully!")
            return generated_files
        
        except Exception as e:
            print(f"[ERRO] Falha ao regenerar arquivos: {str(e)}")
            traceback.print_exc()
            return {}
    
    @staticmethod
    def _regenerate_analysis_summary(results_folder: Path) -> Optional[Path]:
        """Regenera analysis_summary.tsv"""
        results_folder = Path(results_folder)
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
                    
                    # Extrair np
                    np_match = re.search(r'lnL\(ntime:\s*(\d+)\s+np:\s*(\d+)\)', content)
                    np_val = int(np_match.group(2)) if np_match else None
                    
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
                    data[gene_name][f'{model}_lnL'] = lnL
                    data[gene_name][f'{model}_np'] = np_val
                    data[gene_name][f'{model}_omega'] = omega
                    data[gene_name][f'{model}_time'] = exec_time
                    data[gene_name][f'{model}_stops'] = stop_count
                    
                    # Se é Branch-site ou Branch-site_null, extrair dados de classes de sítios
                    if 'Branch-site' in model:
                        class_data = SitesParser.extract_branchsite_class_data(results_file)
                        if class_data:
                            # Armazenar os dados de classe para exibição estruturada
                            data[gene_name][f'{model}_class_data'] = class_data
                
                except Exception as e:
                    print(f"  [WARN] Error processing {gene_name} ({model}): {str(e)}")

        
        # Calcular LRTs
        for gene_name in data:
            row = data[gene_name]
            
            # M0 vs M1a
            if f'M0_lnL' in row and f'M1a_lnL' in row and row[f'M0_lnL'] and row[f'M1a_lnL']:
                lrt = 2 * (row[f'M1a_lnL'] - row[f'M0_lnL'])
                row['lrt_M0_vs_M1a'] = lrt
            
            # M1a vs M2a
            if f'M1a_lnL' in row and f'M2a_lnL' in row and row[f'M1a_lnL'] and row[f'M2a_lnL']:
                lrt = 2 * (row[f'M2a_lnL'] - row[f'M1a_lnL'])
                row['lrt_M1a_vs_M2a'] = lrt
            
            # M7 vs M8
            if f'M7_lnL' in row and f'M8_lnL' in row and row[f'M7_lnL'] and row[f'M8_lnL']:
                lrt = 2 * (row[f'M8_lnL'] - row[f'M7_lnL'])
                row['lrt_M7_vs_M8'] = lrt
            
            # M0 vs Branch
            if f'M0_lnL' in row and f'Branch_lnL' in row and row[f'M0_lnL'] and row[f'Branch_lnL']:
                lrt = 2 * (row[f'Branch_lnL'] - row[f'M0_lnL'])
                row['lrt_M0_vs_Branch'] = lrt
            
            # Branch-site_null vs Branch-site
            if f'Branch-site_null_lnL' in row and f'Branch-site_lnL' in row and row[f'Branch-site_null_lnL'] and row[f'Branch-site_lnL']:
                lrt = 2 * (row[f'Branch-site_lnL'] - row[f'Branch-site_null_lnL'])
                row['lrt_Branch-site_null_vs_Branch-site'] = lrt
        
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
        
        # Converter para DataFrame e salvar
        df = pd.DataFrame(list(data.values()))
        df.to_csv(summary_file, sep='\t', index=False, float_format='%.6f')
        
        return summary_file
    
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
    def _regenerate_lrt_results(results_folder: Path) -> Optional[Path]:
        """Regenera LRT_results.txt"""
        results_folder = Path(results_folder)
        lrt_file = results_folder / "LRT_results.txt"

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
        
        # Definir comparações possíveis
        comparisons = []
        if 'M0' in models and 'M1a' in models:
            comparisons.append(('M0', 'M1a', 'Tests if ω varies among sites', 1))
        if 'M1a' in models and 'M2a' in models:
            comparisons.append(('M1a', 'M2a', 'Tests for positive selection', 1))
        if 'M7' in models and 'M8' in models:
            comparisons.append(('M7', 'M8', 'Tests for positive selection (alternative)', 1))
        if 'M0' in models and 'Branch' in models:
            comparisons.append(('M0', 'Branch', 'Tests branch model (independent evolution rates)', 1))
        if 'Branch-site_null' in models and 'Branch-site' in models:
            # df=1 é usado no teste de mistura 50:50 χ²(0)+χ²(1); ver cálculo abaixo.
            comparisons.append(('Branch-site_null', 'Branch-site', 'Testa seleção positiva no ramo foreground (mistura 50:50 χ²)', 1))
        
        with open(lrt_file, 'w', encoding='utf-8') as f:
            f.write("="*80 + "\n")
            f.write("LIKELIHOOD RATIO TEST (LRT) RESULTS (REGENERATED)\n")
            f.write("="*80 + "\n\n")
            
            if not comparisons:
                f.write("No valid comparisons found\n")
                return lrt_file
            
            for null_model, alt_model, description, df in comparisons:
                f.write("\n" + "="*80 + "\n")
                f.write(f"COMPARISON: {null_model} (null) vs {alt_model} (alternative)\n")
                f.write(f"Description: {description}\n")
                f.write("="*80 + "\n\n")
                
                sig_count_05 = 0
                sig_count_01 = 0
                total_valid = 0
                
                # Comparar cada gene
                for gene in sorted(genes):
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

                        # Branch-site usa distribuição nula 50:50 de χ²(0)+χ²(1),
                        # não χ² padrão — consistente com _run_lrt_analysis.
                        # Referência: Yang et al. (2005), Zhang et al. (2005).
                        is_branchsite = (null_model == 'Branch-site_null'
                                         and alt_model == 'Branch-site')
                        if is_branchsite:
                            p_value = 0.5 * stats.chi2.sf(lrt_stat, df=1) if lrt_stat > 0 else 1.0
                        else:
                            p_value = 1 - stats.chi2.cdf(lrt_stat, df)

                        total_valid += 1
                        
                        if p_value < 0.05:
                            sig_count_05 += 1
                        if p_value < 0.01:
                            sig_count_01 += 1
                        
                        # Escrever resultado
                        f.write(f"Gene: {gene}\n")
                        f.write(f"  lnL {null_model}: {lnL_null:.6f}\n")
                        f.write(f"  lnL {alt_model}: {lnL_alt:.6f}\n")
                        f.write(f"  2Δl = {lrt_stat:.6f}\n")
                        f.write(f"  df = {df}\n")
                        f.write(f"  p-value = {p_value:.6e}\n")
                        
                        if p_value < 0.01:
                            f.write(f"  Result: [OK][OK] {alt_model} significantly better (p < 0.01)\n")
                        elif p_value < 0.05:
                            f.write(f"  Result: [OK] {alt_model} significantly better (p < 0.05)\n")
                        else:
                            f.write(f"  Result: [ERROR] No significant difference\n")
                        
                        f.write("\n" + "-"*60 + "\n\n")
                    
                    except Exception:
                        continue

                # Resumo
                if total_valid > 0:
                    f.write("\nRESUMO:\n")
                    f.write(f"  Total de genes analisados: {total_valid}\n")
                    f.write(f"  Significativo em p < 0.05: {sig_count_05} ({100*sig_count_05/total_valid:.1f}%)\n")
                    f.write(f"  Significativo em p < 0.01: {sig_count_01} ({100*sig_count_01/total_valid:.1f}%)\n")
                    f.write("\n")
        
        return lrt_file


def main():
    """Função principal"""
    try:
        analysis = CodemlBatchAnalysis()
        analysis.run_batch_analysis()
    except KeyboardInterrupt:
        print("\n\n[WARN]  Analysis interrupted by user")
    except Exception as e:
        print(f"\n\n[ERRO] {str(e)}")
        traceback.print_exc()


if __name__ == "__main__":
    main()
