"""
CODEML Interactive Batch Analysis System
Sistema interativo para executar análises CODEML em batch
Gera automaticamente os arquivos .ctl necessários
Requer apenas arquivos .fas e .tree
"""

import os
import subprocess
import time
import shutil
import re
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from threading import Thread
from .sites_parser import SitesParser


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
            'CodonFreq': 2
        },
        'M1a': {
            'description': 'Nearly Neutral - ω < 1 or = 1',
            'model': 0,
            'NSsites': 1,
            'fix_omega': 0,
            'omega': 0.5,
            'CodonFreq': 2
        },
        'M2a': {
            'description': 'Positive Selection - adds ω > 1 class',
            'model': 0,
            'NSsites': 2,
            'fix_omega': 0,
            'omega': 0.5,
            'CodonFreq': 2
        },
        'M7': {
            'description': 'Beta distribution - ω < 1',
            'model': 0,
            'NSsites': 7,
            'fix_omega': 0,
            'omega': 0.5,
            'CodonFreq': 2
        },
        'M8': {
            'description': 'Beta + ω - adds ω > 1 class',
            'model': 0,
            'NSsites': 8,
            'fix_omega': 0,
            'omega': 0.5,
            'CodonFreq': 2
        },
        'Branch': {
            'description': 'Branch model - different ω for foreground',
            'model': 2,
            'NSsites': 0,
            'fix_omega': 0,
            'omega': 0.5,
            'CodonFreq': 2
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
    
    # Modelos neutros (ω=1 FIXADO conforme PAML oficial)
    # Wong et al. 2004, Swanson et al. 2003, Yang et al. 2005
    NEUTRAL_MODELS = {
        'M0': {
            'fix_omega': 1,
            'omega': 1.0,
            'corresponding_alternative': 'Branch',
            'reason': 'M0: ω=1 fixado (modelo nulo para Branch site models)'
        },
        'M1a': {
            'fix_omega': 1,
            'omega': 1.0,
            'corresponding_alternative': 'M2a',
            'reason': 'M1a: ω₁=1 fixado (Wong et al. 2004)'
        },
        'BranchSite_A_null': {
            'fix_omega': 1,
            'omega': 1.0,
            'corresponding_alternative': 'BranchSite_A',
            'reason': 'BranchSite_A_null: ω₂=1 fixado (Yang et al. 2005)'
        }
    }
    
    def __init__(self):
        self.results = {}
        self.config = {}
        # current stop codon count updated during runs (for GUI polling)
        self.current_stop_count = 0
        self.current_stop_details = []
    
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
                             model_name: str = None) -> str:
        """Gera conteúdo do arquivo .ctl baseado no modelo

        Parameters:
        - seqfile, treefile, outfile: paths/names to write into the ctl
        - model_config: dict with model parameters
        - omega: initial omega value (float)
        - cleandata: 0/1 whether to remove ambiguous sites
        - model_name: name of the model (to check if it's a neutral model)
        
        For neutral models (M1a, Branch-site_null), omega is forced to 1.0
        with fix_omega=1, per PAML specifications (Wong et al. 2004, Swanson et al. 2003, Yang et al. 2005)
        """
        # Enforce omega=1 for neutral models
        fix_omega = model_config['fix_omega']
        final_omega = omega
        
        if model_name in self.NEUTRAL_MODELS:
            fix_omega = 1
            final_omega = 1.0
        
        ctl_template = f"""      seqfile = {seqfile}
     treefile = {treefile}
      outfile = {outfile}
   
        noisy = 3              * How much rubbish on the screen
      verbose = 1              * More or less detailed report
      seqtype = 1              * Data type
        ndata = 1              * Number of data sets or loci
        icode = 0              * Genetic code 
    cleandata = {cleandata}              * Remove sites with ambiguity data?
		
        model = {model_config['model']}         * Models for ω varying across lineages
	  NSsites = {model_config['NSsites']}          * Models for ω varying across sites
    CodonFreq = {model_config['CodonFreq']}        * Codon frequencies
	  estFreq = 0              * Use observed freqs or estimate freqs by ML
        clock = 0              * Clock model
    fix_omega = {fix_omega}         * Estimate or fix omega
        omega = {final_omega}        * Initial or fixed omega
"""
        return ctl_template
    
    def interactive_setup(self):
        """Configuração interativa via input do usuário"""
        
        print("\n" + "="*80)
        print("CODEML INTERACTIVE BATCH ANALYSIS")
        print("="*80 + "\n")
        
        # 1. Pasta com arquivos .fas
        while True:
            input_folder = input("📁 Enter path to folder with .fas files: ").strip().strip('"')
            input_path = Path(input_folder)
            if input_path.exists() and input_path.is_dir():
                fas_files = list(input_path.glob("*.fas"))
                if fas_files:
                    print(f"   [OK] Found {len(fas_files)} .fas files")
                    self.config['input_folder'] = input_path
                    break
                else:
                    print("   [ERROR] No .fas files found in this folder. Try again.")
            else:
                print("   [ERROR] Folder not found. Try again.")
        
        # 2. Arquivo de árvore
        while True:
            tree_file = input("\n🌳 Enter path to tree file (.tree or .txt): ").strip().strip('"')
            tree_path = Path(tree_file)
            if tree_path.exists() and tree_path.is_file():
                print(f"   [OK] Tree file loaded: {tree_path.name}")
                self.config['tree_file'] = tree_path
                break
            else:
                print("   [ERROR] Tree file not found. Try again.")
        
        # 3. Pasta de saída
        output_folder = input("\n💾 Enter path for output folder: ").strip().strip('"')
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
        print(f"  9. BranchSite_A_null - {self.MODEL_CONFIGS['BranchSite_A_null']['description']}")
        
        print("\nEnter model numbers separated by spaces (e.g., '1 4 5' for M0, M7, M8)")
        print("Or enter 'all' for all site models (recommended for testing positive selection)")
        
        model_input = input("\n🔢 Select models: ").strip().lower()
        
        model_map = {
            '1': 'M0', '2': 'M1a', '3': 'M2a', '4': 'M7', 
            '5': 'M8', '6': 'Branch',
            '7': 'BranchSite_A', '8': 'BranchSite_A_null'
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
        print("\n⏱️  Set timeout per analysis (in seconds)")
        print("   Recommended: 1600 (≈27 minutes)")
        timeout_input = input("   Timeout [1600]: ").strip()
        self.config['timeout'] = int(timeout_input) if timeout_input else 1600
        print(f"   [OK] Timeout set to {self.config['timeout']} seconds")
        
        # 6. LRT
        print("\n📊 Perform Likelihood Ratio Tests (LRT)?")
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
        
        # Obter arquivos .fas
        fas_files = sorted(self.config['input_folder'].glob("*.fas"))
        # expose total/processed counters for GUI
        try:
            self.current_total_genes = len(fas_files)
            self.current_processed_genes = 0
        except Exception:
            self.current_total_genes = 0
            self.current_processed_genes = 0
        
        print("\n" + "="*80)
        print("STARTING BATCH ANALYSIS")
        print("="*80)
        print(f"Processing {len(fas_files)} genes with {len(self.config['models'])} models each\n")
        
        start_time = time.time()
        
        # Processar cada arquivo
        for idx, fas_file in enumerate(fas_files, 1):
            # Pause support: wait here if pause_event is provided and cleared
            pause_event = self.config.get('pause_event')
            if pause_event is not None:
                pause_event.wait()

            # Reset the 'manual all' flag at the start of each new gene (exon)
            try:
                manual_all_ev = self.config.get('manual_continue_all_event')
                if manual_all_ev is not None:
                    try:
                        manual_all_ev.clear()
                    except Exception:
                        pass
            except Exception:
                pass

            print(f"\n{'='*80}")
            print(f"[{idx}/{len(fas_files)}] Gene: {fas_file.stem}")
            print(f"{'='*80}")
            
            gene_results = {}
            
            # Executar cada modelo
            for model_name in self.config['models']:
                # Pause support before starting each model
                pause_event = self.config.get('pause_event')
                if pause_event is not None:
                    pause_event.wait()

                print(f"\n  - Running {model_name}...", end=" ", flush=True)
                
                result = self._run_single_analysis(
                    fas_file=fas_file,
                    model_name=model_name,
                    log_file=log_file
                )
                
                if result:
                    gene_results[model_name] = result
                    lnL = result.get('lnL')
                    exec_time = result.get('execution_time')
                    lnL_str = f"{lnL:.2f}" if lnL is not None else 'NA'
                    time_str = f"{exec_time:.1f}s" if exec_time is not None else 'NA'
                    print(f"[OK] [lnL: {lnL_str}, time: {time_str}]")
                else:
                    print("[ERROR] Failed")
            
            self.results[fas_file.stem] = gene_results
            # update processed genes counter for GUI
            try:
                self.current_processed_genes = idx
            except Exception:
                pass
        
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
                            log_file: Path) -> Optional[Dict]:
        """Executa análise CODEML para um arquivo e modelo"""
        
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
                # override keys present in custom config
                for k, v in custom_configs[model_name].items():
                    model_config[k] = v
        except Exception:
            pass
        
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
        
        # Gerar arquivo .ctl ou usar .ctl custom do usuário
        try:
            # determine the tree filename that will be present in the temp_dir
            labeled_full = self.config.get('labeled_tree_content')
            labeled_branchsite = self.config.get('labeled_tree_branchsite')

            # Choose which labeled content to use depending on model
            # Support both old name (BranchSite*) and new name (Branch-site*)
            if model_name.startswith('BranchSite') or model_name.startswith('Branch-site'):
                labeled_content = labeled_branchsite or (labeled_full if labeled_full and '#1' in labeled_full else None)
            elif model_name == 'Branch':
                labeled_content = labeled_full
            else:
                labeled_content = None

            if labeled_content:
                tree_filename = 'labeled.nwk'
            else:
                tree_filename = Path(self.config['tree_file']).name

            custom_paths = self.config.get('model_ctl_paths', {}) or {}
            provided_ctl = custom_paths.get(model_name)

            if provided_ctl:
                # If GUI passed a Path object or string, normalize
                provided_path = Path(provided_ctl)
                if provided_path.exists() and provided_path.is_file():
                    # Read user ctl and replace seqfile/treefile/outfile only
                    raw = provided_path.read_text(encoding='utf-8')

                    # Patterns: replace the right-hand side of seqfile/treefile/outfile
                    def _replace_setting(content: str, key: str, newval: str) -> str:
                        pat = rf'(^\s*{re.escape(key)}\s*=).*?$'
                        repl = rf"\1 {newval}"
                        return re.sub(pat, repl, content, flags=re.MULTILINE)

                    raw2 = _replace_setting(raw, 'seqfile', str(fas_file.absolute()))
                    raw2 = _replace_setting(raw2, 'treefile', tree_filename)
                    raw2 = _replace_setting(raw2, 'outfile', output_filename)

                    with open(ctl_path, 'w', encoding='utf-8') as f:
                        f.write(raw2)
                else:
                    # fall back to generated ctl if provided path invalid
                    with open(log_file, 'a', encoding='utf-8') as log:
                        log.write(f"Warning: provided .ctl for {model_name} not found: {provided_ctl}; generating default .ctl\n")
                    ctl_content = self.generate_ctl_content(
                        seqfile=str(fas_file.absolute()),
                        treefile=tree_filename,
                        outfile=output_filename,
                        model_config=model_config,
                        omega=float(self.config.get('omega', model_config.get('omega', 0.5) or 0.5)),
                        cleandata=int(self.config.get('cleandata', 1)),
                        model_name=model_name
                    )
                    with open(ctl_path, 'w', encoding='utf-8') as f:
                        f.write(ctl_content)
            else:
                ctl_content = self.generate_ctl_content(
                    seqfile=str(fas_file.absolute()),
                    treefile=tree_filename,
                    outfile=output_filename,
                    model_config=model_config,
                    omega=float(self.config.get('omega', model_config.get('omega', 0.5) or 0.5)),
                    cleandata=int(self.config.get('cleandata', 1)),
                    model_name=model_name
                )
                with open(ctl_path, 'w', encoding='utf-8') as f:
                    f.write(ctl_content)

        except Exception as e:
            with open(log_file, 'a', encoding='utf-8') as log:
                log.write(f"ERROR creating .ctl [{model_name}] {base_name}: {str(e)}\n")
            return None
        
        # Executar CODEML
        temp_dir = model_output_dir / f"temp_{base_name}"
        
        try:
            # Preparar diretório temporário (remover com retries no Windows se necessário)
            if temp_dir.exists():
                removed = False
                for attempt in range(5):
                    try:
                        shutil.rmtree(temp_dir)
                        removed = True
                        break
                    except PermissionError as pe:
                        with open(log_file, 'a', encoding='utf-8') as log:
                            log.write(f"Warning: could not remove {temp_dir} (attempt {attempt+1}/5): {pe}\n")
                        time.sleep(0.5)
                if not removed:
                    # If we couldn't remove the directory (likely locked files), avoid reusing it.
                    # Create a new unique temp dir to avoid conflicts with locked files.
                    timestamp = int(time.time())
                    new_temp_dir = model_output_dir / f"temp_{base_name}_{timestamp}"
                    with open(log_file, 'a', encoding='utf-8') as log:
                        log.write(f"Warning: failed to remove {temp_dir}; using new temp dir {new_temp_dir}\n")
                    temp_dir = new_temp_dir
            temp_dir.mkdir(exist_ok=True)

            # Copiar arquivos necessários
            shutil.copy(ctl_path, temp_dir)
            shutil.copy(fas_file, temp_dir)
            # If a labeled tree was selected for this model, write it into temp_dir
            if labeled_content:
                try:
                    with open(temp_dir / 'labeled.nwk', 'w', encoding='utf-8') as tf:
                        tf.write(labeled_content)
                except Exception:
                    try:
                        shutil.copy(self.config['tree_file'], temp_dir)
                    except Exception:
                        pass
            else:
                shutil.copy(self.config['tree_file'], temp_dir)

            exec_start = time.time()

            # Executar CODEML
            # Prefer explicit bin/codeml.exe if present, otherwise try 'codeml' on PATH
            codeml_bin = Path('bin') / 'codeml.exe'
            if codeml_bin.exists():
                cmd = [str(codeml_bin), ctl_filename]
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

                            auto_continue = bool(self.config.get('auto_continue_stop_codons', True))
                            manual_all_ev = self.config.get('manual_continue_all_event')

                            # If user requested auto-continue, or a manual-all event is set, send Enter
                            if auto_continue:
                                try:
                                    if process.stdin:
                                        process.stdin.write('\n')
                                        process.stdin.flush()
                                    with open(log_file, 'a', encoding='utf-8') as log:
                                        log.write(f"[{model_name}] {base_name}: Auto-sent Enter to subprocess (stop codon).\n")
                                except Exception as e:
                                    with open(log_file, 'a', encoding='utf-8') as log:
                                        log.write(f"[{model_name}] {base_name}: Failed to auto-send Enter: {e}\n")
                            elif manual_all_ev is not None and getattr(manual_all_ev, 'is_set') and manual_all_ev.is_set():
                                try:
                                    if process.stdin:
                                        process.stdin.write('\n')
                                        process.stdin.flush()
                                    with open(log_file, 'a', encoding='utf-8') as log:
                                        log.write(f"[{model_name}] {base_name}: manual-all event set; sent Enter.\n")
                                except Exception as e:
                                    with open(log_file, 'a', encoding='utf-8') as log:
                                        log.write(f"[{model_name}] {base_name}: Failed to send Enter for manual-all: {e}\n")
                            else:
                                # Wait for GUI/user to signal continuation via event
                                event = self.config.get('manual_continue_event')
                                if event is None:
                                    # no event provided -> fallback to auto
                                    try:
                                        if process.stdin:
                                            process.stdin.write('\n')
                                            process.stdin.flush()
                                        with open(log_file, 'a', encoding='utf-8') as log:
                                            log.write(f"[{model_name}] {base_name}: No manual event provided; auto-sent Enter.\n")
                                    except Exception as e:
                                        with open(log_file, 'a', encoding='utf-8') as log:
                                            log.write(f"[{model_name}] {base_name}: Failed fallback auto-send Enter: {e}\n")
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
                                    try:
                                        if process.stdin:
                                            process.stdin.write('\n')
                                            process.stdin.flush()
                                        with open(log_file, 'a', encoding='utf-8') as log:
                                            log.write(f"[{model_name}] {base_name}: Manual continue event received; sent Enter.\n")
                                    except Exception as e:
                                        with open(log_file, 'a', encoding='utf-8') as log:
                                            log.write(f"[{model_name}] {base_name}: Failed to send Enter after manual event: {e}\n")
                except Exception:
                    pass

            stdout_thread = Thread(target=read_stream, args=(process.stdout, stdout_lines))
            stderr_thread = Thread(target=read_stream, args=(process.stderr, stderr_lines))

            stdout_thread.start()
            stderr_thread.start()

            try:
                process.wait(timeout=self.config['timeout'])
            except subprocess.TimeoutExpired:
                # Timeout: try to capture what we have and kill process
                try:
                    if process.poll() is None:
                        process.kill()
                except Exception:
                    pass
                stdout_thread.join(timeout=1)
                stderr_thread.join(timeout=1)
                with open(log_file, 'a', encoding='utf-8') as log:
                    log.write(f"[{model_name}] {base_name}: TIMEOUT after {self.config['timeout']}s\n")
                    log.write(f"  Captured stdout (last 200 lines):\n")
                    for L in stdout_lines[-200:]:
                        log.write(L + "\n")
                    log.write(f"  Captured stderr (last 200 lines):\n")
                    for L in stderr_lines[-200:]:
                        log.write(L + "\n")
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
            import traceback
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
        except:
            pass
        return None
    
    def _extract_omega(self, output_file: Path) -> Optional[float]:
        """Extrai omega usando SitesParser (suporta Branch/Branch-Site/Site models)"""
        try:
            # Usar função robusta que tenta múltiplas estratégias
            omega = SitesParser.extract_omega_robust(output_file)
            return omega
        except:
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
                row = [str(gene_name).replace('\n', '').replace('\r', '')]  # Remover newlines
                
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
                                except:
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
            if 'BranchSite_A_null' in selected_models and 'BranchSite_A' in selected_models:
                comparisons.append(('BranchSite_A_null', 'BranchSite_A', 
                                  'Tests for positive selection in foreground sites'))
            
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
                    
                    if df == 0 or lrt_stat < 0:
                        continue
                    
                    # Calcular p-value
                    from scipy import stats
                    
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
                f.write("\nSUMMARY:\n")
                f.write(f"  Total genes analyzed: {total_valid}\n")
                f.write(f"  Significant at p < 0.05: {sig_count_05} ({100*sig_count_05/total_valid:.1f}%)\n")
                f.write(f"  Significant at p < 0.01: {sig_count_01} ({100*sig_count_01/total_valid:.1f}%)\n")
                f.write("\n")
        
        print(f"\n  [OK] LRT results saved: {lrt_file}")
    
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
            print(f"[ERROR] Error regenerating files: {str(e)}")
            import traceback
            traceback.print_exc()
            return {}
    
    @staticmethod
    def _regenerate_analysis_summary(results_folder: Path) -> Optional[Path]:
        """Regenera analysis_summary.tsv"""
        import pandas as pd
        
        results_folder = Path(results_folder)
        summary_file = results_folder / "analysis_summary.tsv"
        
        # Mapeamento de nomes de pasta (antigos) para nomes de modelo (novos)
        model_name_mapping = {
            'BranchSite_A': 'Branch-site',
            'BranchSite_A_null': 'Branch-site_null'
        }
        
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
        
        # Mapeamento de nomes
        model_name_mapping = {
            'BranchSite_A': 'Branch-site',
            'BranchSite_A_null': 'Branch-site_null'
        }
        
        with open(log_file, 'w', encoding='utf-8') as f:
            f.write("="*80 + "\n")
            f.write("CODEML BATCH ANALYSIS LOG (REGENERATED)\n")
            f.write("="*80 + "\n")
            f.write(f"Regenerated at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Results folder: {results_folder}\n")
            f.write("="*80 + "\n\n")
            
            f.write("ANALYSIS SUMMARY:\n")
            f.write("-"*80 + "\n")
            
            # Descobrir modelos e genes
            models = set()
            genes = set()
            
            for item in results_folder.iterdir():
                if item.is_dir() and item.name not in ['reports']:
                    model = model_name_mapping.get(item.name, item.name)
                    models.add(model)
                    
                    for results_file in item.glob("*_results.txt"):
                        gene = results_file.name.split(f'_{model}_results')[0]
                        genes.add(gene)
            
            f.write(f"Models found: {', '.join(sorted(models))}\n")
            f.write(f"Genes found: {len(genes)} genes\n")
            f.write(f"  {', '.join(sorted(genes)[:5])}" + ("..." if len(genes) > 5 else "") + "\n")
            f.write("\n")
            
            # Detalhes de cada gene/modelo
            f.write("DETAILED RESULTS:\n")
            f.write("-"*80 + "\n\n")
            
            for gene in sorted(genes):
                f.write(f"Gene: {gene}\n")
                f.write("-"*40 + "\n")
                
                for model in sorted(models):
                    model_folder = results_folder / model
                    results_file = model_folder / f"{gene}_{model}_results.txt"
                    
                    if results_file.exists():
                        try:
                            with open(results_file, 'r', encoding='utf-8', errors='ignore') as rf:
                                content = rf.read()
                            
                            lnL_match = re.search(r'lnL\(.*?\):\s+([-\d.]+)', content)
                            np_match = re.search(r'np:\s*(\d+)\)', content)
                            
                            lnL = float(lnL_match.group(1)) if lnL_match else "NA"
                            np_val = np_match.group(1) if np_match else "NA"
                            
                            f.write(f"  {model:20s} | lnL = {lnL:>12} | np = {np_val:>2}\n")
                        except:
                            f.write(f"  {model:20s} | Error reading file\n")
                    else:
                        f.write(f"  {model:20s} | Not found\n")
                
                f.write("\n")
            
            f.write("="*80 + "\n")
            f.write("END OF LOG\n")
            f.write("="*80 + "\n")
        
        return log_file
    
    @staticmethod
    def _regenerate_lrt_results(results_folder: Path) -> Optional[Path]:
        """Regenera LRT_results.txt"""
        from scipy import stats
        
        results_folder = Path(results_folder)
        lrt_file = results_folder / "LRT_results.txt"
        
        # Mapeamento de nomes
        model_name_mapping = {
            'BranchSite_A': 'Branch-site',
            'BranchSite_A_null': 'Branch-site_null'
        }
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
            comparisons.append(('Branch-site_null', 'Branch-site', 'Tests branch-site model', 2))
        
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
                    
                    except Exception as e:
                        continue
                
                # Sumário
                if total_valid > 0:
                    f.write("\nSUMMARY:\n")
                    f.write(f"  Total genes analyzed: {total_valid}\n")
                    f.write(f"  Significant at p < 0.05: {sig_count_05} ({100*sig_count_05/total_valid:.1f}%)\n")
                    f.write(f"  Significant at p < 0.01: {sig_count_01} ({100*sig_count_01/total_valid:.1f}%)\n")
                    f.write("\n")
        
        return lrt_file


def main():
    """Função principal"""
    try:
        analysis = CodemlBatchAnalysis()
        analysis.run_batch_analysis()
    except KeyboardInterrupt:
        print("\n\n[WARN]️  Analysis interrupted by user")
    except Exception as e:
        print(f"\n\n[ERROR] Error: {str(e)}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
