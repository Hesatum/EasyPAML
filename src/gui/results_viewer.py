"""
Visualizador de Resultados CODEML - EasyPAML
=============================================
Versão 2.0 - CORRIGIDA com busca recursiva em subpastas
- Fix na busca de arquivos em subpastas (M0/, M1a/, M2a/, M7/, M8/, Branch/)
- Extração correta de ω de cada tipo de modelo
- Detecção adequada de seleção positiva
"""

import customtkinter as ctk
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy import stats
from tkinter import filedialog, messagebox
import re
import sys
from src.backend.branch_extractor import BranchExtractor
from .gui_texts import TEXTS


class ResultsViewerWindow(ctk.CTkToplevel):
    """Janela de visualização profissional de resultados"""
    
    COLORS = {
        # Backgrounds — near-black, like Linear / Discord
        'bg_dark':        '#0c0c0e',
        'bg_card':        '#16161a',
        'bg_card_hover':  '#1e1e24',
        'bg_feed':        '#111115',
        'bg_sidebar':     '#111115',
        'bg_input':       '#1e1e24',

        # Text hierarchy
        'text_primary':   '#ededef',
        'text_secondary': '#9898a6',
        'text_tertiary':  '#5e5e6e',
        'text_muted':     '#3a3a48',

        # Accent — indigo (Linear-inspired)
        'accent_blue':        '#6366f1',
        'accent_blue_hover':  '#4f46e5',
        'accent_blue_light':  '#818cf8',

        # Secondary accents
        'accent_cyan':    '#22d3ee',
        'accent_cyan_hover': '#06b6d4',
        'accent_purple':  '#a78bfa',
        'accent_purple_hover': '#7c3aed',

        # Status
        'success':        '#22c55e',
        'success_hover':  '#16a34a',
        'success_light':  '#86efac',
        'warning':        '#f59e0b',
        'warning_hover':  '#d97706',
        'danger':         '#f87171',
        'danger_hover':   '#ef4444',
        'info':           '#22d3ee',

        # Borders
        'border':         '#222228',
        'border_hover':   '#32323e',
    }
    
    def __init__(self, parent, output_folder: Path):
        super().__init__(parent)
        self.title(TEXTS["viewer_window_title"])
        w, h = 1400, 900
        self.update_idletasks()
        sw = self.winfo_screenwidth()
        sh = self.winfo_screenheight()
        x  = max(0, (sw - w) // 2)
        y  = max(0, (sh - h) // 2)
        self.geometry(f"{w}x{h}+{x}+{y}")
        
        self.configure(fg_color=self.COLORS['bg_dark'])
        self.after(100, self.lift)
        self.after(150, self.focus_force)
        self.output_folder = output_folder
        self.df = None
        self.tag_columns = {}
        
        if not self._load_data():
            self._show_error(TEXTS["viewer_error_no_tsv"])
            return
        
        self._extract_tag_columns()
        self.setup_ui()
    
    def _load_data(self) -> bool:
        """Carrega dados do TSV com tratamento robusto"""
        tsv_file = self.output_folder / "analysis_summary.tsv"
        if not tsv_file.exists():
            return False

        try:
            # ── 1. Detect model subfolders not yet in the TSV ──────────────────
            _legacy = {'BranchSite_A': 'Branch-site', 'BranchSite_A_null': 'Branch-site_null'}
            _exclude = {'reports'}
            found_models = set()
            for _item in self.output_folder.iterdir():
                if _item.is_dir() and _item.name not in _exclude:
                    found_models.add(_legacy.get(_item.name, _item.name))

            # Models already described in the TSV (read header only)
            import csv
            with open(tsv_file, newline='', encoding='utf-8', errors='ignore') as _f:
                _header = next(csv.reader(_f, delimiter='\t'), [])
            known_models = {col.rsplit('_', 1)[0] for col in _header
                            if col.endswith('_lnL')}

            missing_models = found_models - known_models
            if missing_models:
                print(f"[INFO] Models in subfolders missing from TSV: {missing_models}")
                print("[INFO] Regenerating analysis_summary.tsv to include all models...")
                from src.backend.codeml_backend import CodemlBatchAnalysis
                CodemlBatchAnalysis._regenerate_analysis_summary(self.output_folder)
                print("[OK] analysis_summary.tsv regenerated successfully.")

            # ── 2. Load (possibly updated) TSV ────────────────────────────────
            self.df = pd.read_csv(tsv_file, sep='\t')
            self.df = self.df.replace(['NA', 'nan', '', 'None'], np.nan)

            numeric_cols = [col for col in self.df.columns if col != 'Gene']
            for col in numeric_cols:
                self.df[col] = pd.to_numeric(self.df[col], errors='coerce')

            # Tentar recuperar omegas faltantes dos arquivos de resultados
            self._recover_missing_omegas()

            print(f"[OK] Dados carregados: {len(self.df)} genes")
            print(f"[OK] Colunas: {list(self.df.columns)}")
            return True
        except Exception as e:
            print(f"[ERR] Erro ao carregar dados: {e}")
            return False
    
    def _recover_missing_omegas(self):
        """Recupera omegas faltantes diretamente dos arquivos de resultados"""
        from src.backend.sites_parser import SitesParser
        from pathlib import Path
        
        # Identificar colunas de omega que estão vazias
        omega_cols = [col for col in self.df.columns if '_omega' in col]
        
        for omega_col in omega_cols:
            # Extrair nome do modelo (ex: M2a_omega -> M2a)
            model_name = omega_col.replace('_omega', '')
            
            # Procurar por arquivos faltantes
            missing_rows = self.df[self.df[omega_col].isna()].index
            
            if len(missing_rows) == 0:
                continue
            
            print(f"\n[INFO] Recuperando omegas faltantes para {model_name}...")
            
            for idx in missing_rows:
                gene_name = self.df.loc[idx, 'Gene']
                
                # Procurar arquivo de resultados com suporte a múltiplas variações
                results_file = self._find_results_file(gene_name, model_name)
                
                if results_file:
                    try:
                        omega = SitesParser.extract_omega_robust(results_file)
                        if omega is not None:
                            self.df.loc[idx, omega_col] = omega
                            print(f"  [OK] {gene_name} ({model_name}): w = {omega:.4f}")
                        else:
                            print(f"  [SKIP] {gene_name} ({model_name}): nao foi possivel extrair")
                    except Exception as e:
                        print(f"  [ERR] {gene_name} ({model_name}): erro - {e}")
                else:
                    print(f"  [INFO] {gene_name} ({model_name}): arquivo nao encontrado")
    
    def _find_results_file(self, gene_name: str, model_name: str):
        """
        Busca arquivo de resultados para um gene e modelo específico.
        
        A estrutura esperada é:
        - results_folder/M8/gene_name_M8_results.txt
        - results_folder/M2a/gene_name_M2a_results.txt
        - results_folder/Branch-site/gene_name_Branch-site_results.txt
        - results_folder/BranchSite_A/gene_name_BranchSite_A_results.txt (compatibilidade com versão antiga)
        """
        # Procurar no padrão padrão
        results_file = self.output_folder / model_name / f"{gene_name}_{model_name}_results.txt"
        if results_file.exists():
            return results_file
        
        # Para Branch-site, tentar também o nome antigo (BranchSite_A)
        if model_name == 'Branch-site':
            results_file = self.output_folder / 'BranchSite_A' / f"{gene_name}_BranchSite_A_results.txt"
            if results_file.exists():
                return results_file
        
        return None
    
    def _extract_tag_columns(self):
        """Extrai colunas dinâmicas de tags"""
        tag_pattern = r'(.+?)_([^_]+)_(omega|lnL)$'
        
        self.tag_columns = {}
        
        for col in self.df.columns:
            match = re.match(tag_pattern, col)
            if match:
                model, tag, metric = match.groups()
                
                if model not in self.tag_columns:
                    self.tag_columns[model] = {'tags': set(), 'omega': {}, 'lnL': {}}
                
                self.tag_columns[model]['tags'].add(tag)
                if metric == 'omega':
                    self.tag_columns[model]['omega'][tag] = col
                elif metric == 'lnL':
                    self.tag_columns[model]['lnL'][tag] = col
        
        print(f"Tags detected: {self.tag_columns}")
    
    def _format_branchsite_class_data(self, gene_idx: int) -> str:
        """
        Formata dados de classes de Branch-site para exibição estruturada
        Retorna string com formatação multilinea
        """
        row = self.df.iloc[gene_idx]
        
        # Procurar colunas de classe
        class_cols = [col for col in self.df.columns if 'Branch-site_class' in col and 'null' not in col]
        
        if not class_cols:
            return "N/A"
        
        # Agrupar por classe
        classes = {}
        for col in class_cols:
            # Parse: Branch-site_class0_fg_w -> classe='0', métrica='fg_w'
            parts = col.replace('Branch-site_class', '').split('_', 1)
            if len(parts) == 2:
                cls, metric = parts
                if cls not in classes:
                    classes[cls] = {}
                classes[cls][metric] = row[col]
        
        # Construir string formatada
        lines = ["Branch-site Classes:"]
        
        # Ordenar classes: 0, 1, 2a, 2b
        for cls in ['0', '1', '2a', '2b']:
            if cls in classes:
                data = classes[cls]
                fg_w = data.get('fg_w', 'N/A')
                prop = data.get('prop', 'N/A')
                
                if isinstance(fg_w, float):
                    fg_w_str = f"{fg_w:.5f}"
                else:
                    fg_w_str = str(fg_w)
                
                if isinstance(prop, float):
                    prop_str = f"{prop:.5f}"
                else:
                    prop_str = str(prop)
                
                lines.append(f"  Class {cls}: prop={prop_str}, fg_w={fg_w_str}")
        
        return "\n".join(lines)
    
    @staticmethod
    def _fmt_pval(p: float) -> str:
        """Format p-value as human-readable decimal (no scientific notation)."""
        if p <= 0:        return "0.00000000"
        if p < 0.000001:  return f"{p:.8f}"
        if p < 0.0001:    return f"{p:.7f}"
        if p < 0.001:     return f"{p:.6f}"
        if p < 0.01:      return f"{p:.5f}"
        return f"{p:.4f}"

    def _show_error(self, message: str):
        """Exibe tela de erro"""
        error_frame = ctk.CTkFrame(self, fg_color=self.COLORS['bg_dark'])
        error_frame.pack(fill='both', expand=True, padx=20, pady=20)
        
        ctk.CTkLabel(error_frame, text="[X]", font=("Roboto", 48)).pack(pady=20)
        ctk.CTkLabel(error_frame, text=message,
                    font=("Roboto", 14, "bold"),
                    text_color=self.COLORS['danger']).pack(pady=10)
        ctk.CTkLabel(error_frame, text=TEXTS["viewer_error_run_analysis"],
                    font=("Roboto", 11),
                    text_color=self.COLORS['text_tertiary']).pack()
    
    def setup_ui(self):
        """Setup da interface premium"""

        # ── HEADER ──────────────────────────────────────────────────
        header = ctk.CTkFrame(self, fg_color=self.COLORS['bg_sidebar'], corner_radius=0, height=64)
        header.pack(fill='x', padx=0, pady=0)
        header.pack_propagate(False)

        left = ctk.CTkFrame(header, fg_color='transparent')
        left.pack(side="left", padx=24, pady=0, fill='y')

        ctk.CTkLabel(left, text="EP", font=("Roboto", 24)).pack(side="left", padx=(0, 12))
        title_block = ctk.CTkFrame(left, fg_color='transparent')
        title_block.pack(side="left", fill='y', pady=14)
        ctk.CTkLabel(title_block, text=TEXTS["viewer_header_title"],
                     font=("Roboto", 16, "bold"),
                     text_color=self.COLORS['text_primary']).pack(anchor="w")
        ctk.CTkLabel(title_block, text=TEXTS["viewer_header_subtitle"],
                     font=("Roboto", 9),
                     text_color=self.COLORS['text_muted']).pack(anchor="w")

        right = ctk.CTkFrame(header, fg_color='transparent')
        right.pack(side="right", padx=24, pady=16, fill='y')
        ctk.CTkLabel(right, text=TEXTS["viewer_genes_loaded"].format(n=len(self.df)),
                     font=("Roboto", 10), text_color=self.COLORS['accent_blue']).pack()
        
        # PAINEL DE ESTATÍSTICAS
        stats_frame = ctk.CTkFrame(self, fg_color='transparent')
        stats_frame.pack(fill='x', padx=20, pady=(16, 4))
        self._create_stats_panel(stats_frame)
        
        # ABAS PRINCIPAIS
        tabs = ctk.CTkTabview(self, fg_color=self.COLORS['bg_card'],
                              segmented_button_fg_color=self.COLORS['bg_sidebar'],
                              segmented_button_selected_color=self.COLORS['accent_blue'],
                              segmented_button_unselected_color=self.COLORS['bg_sidebar'],
                              text_color=self.COLORS['text_tertiary'],
                              segmented_button_selected_hover_color=self.COLORS['accent_blue_hover'],
                              corner_radius=12,
                              border_width=1,
                              border_color=self.COLORS['border'])
        tabs.pack(fill='both', expand=True, padx=20, pady=(8, 20))
        
        tabs.add(TEXTS["viewer_tab_lrt"])
        tabs.add(TEXTS["viewer_tab_omega"])
        tabs.add(TEXTS["viewer_tab_sites"])

        # Verificar se há dados de Branch-site para adicionar aba especial
        branchsite_cols = [col for col in self.df.columns if 'Branch-site_class' in col]
        if branchsite_cols:
            tabs.add(TEXTS["viewer_tab_branchsite_classes"])

        tabs.add(TEXTS["viewer_tab_branch"])
        tabs.add(TEXTS["viewer_tab_export"])

        self._create_lrt_stats_tab(tabs.tab(TEXTS["viewer_tab_lrt"]))
        self._create_positive_selection_tab(tabs.tab(TEXTS["viewer_tab_omega"]))
        self._create_sites_tab(tabs.tab(TEXTS["viewer_tab_sites"]))

        if branchsite_cols:
            self._create_branchsite_class_tab(tabs.tab(TEXTS["viewer_tab_branchsite_classes"]))

        self._create_tree_tab(tabs.tab(TEXTS["viewer_tab_branch"]))
        self._create_export_tab(tabs.tab(TEXTS["viewer_tab_export"]))
    
    def _create_stats_panel(self, parent):
        """Painel com estatísticas gerais — cards premium"""
        positive_genes = self._detect_positive_selection()

        stats_data = [
            (TEXTS["stats_total_genes"],        str(len(self.df)),               self.COLORS['accent_blue'], ""),
            (TEXTS["stats_models_run"],         self._count_models(),            self.COLORS['accent_cyan'], ""),
            (TEXTS["stats_positive_selection"], str(len(positive_genes)),        self.COLORS['success'],     ""),
            (TEXTS["stats_avg_omega"],          f"{self._calc_avg_omega():.3f}", self.COLORS['warning'],     ""),
        ]

        for label, value, color, icon in stats_data:
            card = ctk.CTkFrame(parent, fg_color=self.COLORS['bg_card'],
                                corner_radius=10, border_width=1,
                                border_color=self.COLORS['border'])
            card.pack(side="left", fill="both", expand=True, padx=5)

            inner = ctk.CTkFrame(card, fg_color='transparent')
            inner.pack(fill='both', expand=True, padx=16, pady=14)

            ctk.CTkLabel(inner, text=f"{label}" if not icon else f"{icon}  {label}",
                         font=("Roboto", 10),
                         text_color=self.COLORS['text_tertiary'],
                         anchor='w').pack(anchor='w')
            ctk.CTkLabel(inner, text=value,
                         font=("Roboto", 28, "bold"),
                         text_color=color,
                         anchor='w').pack(anchor='w', pady=(4, 0))
    
    def _create_lrt_stats_tab(self, parent):
        """Aba de Tabela LRT com p-valores"""
        comparisons, descriptions = self._get_available_lrt_columns()
        if not comparisons:
            ctk.CTkLabel(parent, text=TEXTS["lrt_no_comparisons"],
                        font=("Roboto", 12),
                        text_color=self.COLORS['warning']).pack(pady=50)
            return

        # ── selector card ──────────────────────────────────────────────
        ctrl_frame = ctk.CTkFrame(parent, fg_color=self.COLORS['bg_card'], corner_radius=8)
        ctrl_frame.pack(fill='x', padx=10, pady=(10, 4))

        row1 = ctk.CTkFrame(ctrl_frame, fg_color='transparent')
        row1.pack(fill='x', padx=14, pady=(12, 4))

        ctk.CTkLabel(row1, text=TEXTS["lrt_label_model"],
                     font=("Roboto", 11, "bold"),
                     text_color=self.COLORS['text_secondary']).pack(side='left', padx=(0, 10))

        comp_combo = ctk.CTkComboBox(row1, values=list(comparisons.keys()), width=380,
                                     font=("Roboto", 11))
        comp_combo.pack(side='left')
        comp_combo.set(list(comparisons.keys())[0])

        # dynamic null-hypothesis description
        desc_lbl = ctk.CTkLabel(ctrl_frame, text="",
                                font=("Roboto", 9),
                                text_color=self.COLORS['text_tertiary'],
                                anchor='w', justify='left')
        desc_lbl.pack(fill='x', padx=16, pady=(0, 10))

        table_frame = ctk.CTkScrollableFrame(parent, fg_color=self.COLORS['bg_feed'],
                                            corner_radius=8)
        table_frame.pack(fill='both', expand=True, padx=10, pady=(4, 10))

        def update_lrt_table(*args):
            for widget in table_frame.winfo_children():
                widget.destroy()
            selected_comp = comp_combo.get()
            col_name = comparisons[selected_comp]
            desc = descriptions.get(col_name, '')
            desc_lbl.configure(text=desc)
            self._render_lrt_table(table_frame, col_name, selected_comp)

        comp_combo.configure(command=update_lrt_table)
        update_lrt_table()
    
    def _create_positive_selection_tab(self, parent):
        """Genes com ω global > 1 (seleção positiva ao nível do gene inteiro)"""
        # ── Info banner ────────────────────────────────────────────────
        info = ctk.CTkFrame(parent, fg_color='#1a1a26', corner_radius=8)
        info.pack(fill='x', padx=10, pady=(10, 2))

        ctk.CTkLabel(info,
                     text=TEXTS["pos_sel_tab_title"],
                     font=("Roboto", 11, "bold"),
                     text_color=self.COLORS['accent_blue']).pack(side="left", padx=14, pady=(10, 2))

        ctk.CTkLabel(info,
                     text=TEXTS["pos_sel_tab_criterion"],
                     font=("Roboto", 9),
                     text_color=self.COLORS['text_tertiary'],
                     wraplength=900,
                     justify='left').pack(anchor='w', padx=14, pady=(0, 10))

        positive_data = self._detect_positive_selection()

        if not positive_data:
            empty = ctk.CTkFrame(parent, fg_color='transparent')
            empty.pack(expand=True)
            ctk.CTkLabel(empty, text="—", font=("Roboto", 32),
                         text_color=self.COLORS['text_muted']).pack(pady=(60, 8))
            ctk.CTkLabel(empty, text=TEXTS["pos_sel_none_found"],
                         font=("Roboto", 13, "bold"),
                         text_color=self.COLORS['text_tertiary']).pack()
            ctk.CTkLabel(empty, text=TEXTS["pos_sel_criterion_short"],
                         font=("Roboto", 10),
                         text_color=self.COLORS['text_muted']).pack(pady=(4, 0))
            return

        scroll_frame = ctk.CTkScrollableFrame(parent, fg_color='transparent',
                                              corner_radius=8)
        scroll_frame.pack(fill='both', expand=True, padx=10, pady=10)

        for gene_name, signals in positive_data.items():
            card = ctk.CTkFrame(scroll_frame, fg_color='#0b2016',
                                corner_radius=12, border_width=1,
                                border_color='#10b981')
            card.pack(fill='x', pady=6, padx=4)

            # Left accent
            ctk.CTkFrame(card, fg_color='#10b981', width=5,
                         corner_radius=2).pack(side="left", fill="y", padx=(6, 0), pady=8)

            content = ctk.CTkFrame(card, fg_color='transparent')
            content.pack(side="left", fill="both", expand=True, padx=14, pady=12)

            ctk.CTkLabel(content, text=f"{gene_name}",
                         font=("Roboto", 13, "bold"),
                         text_color='#6ee7b7').pack(anchor='w')

            for signal_type, signal_data in signals.items():
                signal_text = (
                    f"  {signal_type}  ·  "
                    f"ω = {signal_data['omega']:.4f}  ·  "
                    f"2Δℓ = {signal_data['lrt']:.3f}  ·  "
                    f"p = {self._fmt_pval(signal_data['p_value'])}"
                )
                ctk.CTkLabel(content, text=signal_text,
                             font=("Roboto", 10),
                             text_color='#a7f3d0').pack(anchor='w', pady=(3, 0))

            badge = ctk.CTkFrame(card, fg_color='#10b981',
                                 corner_radius=8, width=70, height=36)
            badge.pack(side="right", padx=14)
            badge.pack_propagate(False)
            ctk.CTkLabel(badge, text=TEXTS["pos_sel_badge"],
                         font=("Roboto", 9, "bold"),
                         text_color="white").pack(expand=True)
    
    def _create_sites_tab(self, parent):
        """Aba de visualização de sítios sob seleção positiva"""
        ctrl_frame = ctk.CTkFrame(parent, fg_color=self.COLORS['bg_card'],
                                 corner_radius=8, height=120)
        ctrl_frame.pack(fill='x', padx=10, pady=10)
        ctrl_frame.pack_propagate(False)
        
        line1 = ctk.CTkFrame(ctrl_frame, fg_color='transparent')
        line1.pack(fill='x', padx=15, pady=(12, 8))
        
        # 1º: Selecionar Modelo primeiro
        ctk.CTkLabel(line1, text=TEXTS["sites_label_model"], font=("Roboto", 10, "bold")).pack(side='left', padx=(0, 8))

        model_combo = ctk.CTkComboBox(line1, values=['M8', 'M2a', 'Branch-site'], width=150)
        model_combo.pack(side='left', padx=(0, 30))
        model_combo.set('M8')

        # 2º: Gene (será atualizado quando modelo mudar)
        ctk.CTkLabel(line1, text=TEXTS["sites_label_gene"], font=("Roboto", 10, "bold")).pack(side='left', padx=(0, 8))

        gene_combo = ctk.CTkComboBox(line1, values=[], width=200)
        gene_combo.pack(side='left', padx=(0, 30))

        # 3º: Análise
        ctk.CTkLabel(line1, text=TEXTS["sites_label_analysis"], font=("Roboto", 10, "bold")).pack(side='left', padx=(0, 8))
        
        method_combo = ctk.CTkComboBox(line1, values=['BEB', 'NEB'], width=100)
        method_combo.pack(side='left')
        method_combo.set('BEB')
        
        # Função para atualizar genes quando modelo mudar
        def update_gene_list(*args):
            """Atualiza lista de genes baseado no modelo selecionado"""
            selected_model = model_combo.get()
            
            # Tentar encontrar a pasta correta
            model_folder = self.output_folder / selected_model
            
            # Compatibilidade para Branch-site (testar nome novo e antigo)
            if selected_model == 'Branch-site' and not model_folder.exists():
                fallback = self.output_folder / 'BranchSite_A'
                if fallback.exists():
                    model_folder = fallback
            
            # Extrair genes daquela pasta específica
            genes = []
            if model_folder.exists():
                # Regex atualizado para suportar hifens no nome do modelo (ex: Branch-site)
                genes = sorted(set(
                    re.match(r'(.+?)_[A-Za-z0-9\-]+_results\.txt', f.name).group(1)
                    for f in model_folder.glob('*_results.txt')
                    if re.match(r'(.+?)_[A-Za-z0-9\-]+_results\.txt', f.name)
                ))
            
            gene_combo.configure(values=genes)
            if genes:
                gene_combo.set(genes[0])
            else:
                gene_combo.set('')
        
        # Conectar callback para quando modelo mudar
        model_combo.configure(command=update_gene_list)
        
        # Popular inicial de genes com M8
        update_gene_list()
        
        line2 = ctk.CTkFrame(ctrl_frame, fg_color='transparent')
        line2.pack(fill='x', padx=15, pady=(0, 12))
        
        ctk.CTkLabel(line2, text=TEXTS["sites_label_filter"], font=("Roboto", 10, "bold")).pack(side='left', padx=(0, 8))
        
        p_filter = ctk.CTkEntry(line2, placeholder_text="0.95", width=80)
        p_filter.pack(side='left', padx=(0, 15))
        p_filter.insert(0, "0.95")
        
        def on_p_filter_change(*args):
            try:
                float(p_filter.get())
            except ValueError:
                p_filter.delete(0, 'end')
                p_filter.insert(0, '0.95')
        
        p_filter.bind('<KeyRelease>', on_p_filter_change)
        
        def update_sites_table(*args):
            for widget in table_frame.winfo_children():
                widget.destroy()
            
            try:
                p_threshold = float(p_filter.get())
            except ValueError:
                p_threshold = 0.95
            
            gene_name = gene_combo.get()
            model_name = model_combo.get()
            method = method_combo.get()
            
            self._render_sites_table(table_frame, gene_name, model_name, method, p_threshold)
        
        btn_update = ctk.CTkButton(line2, text=TEXTS["sites_btn_update"], width=120,
                                  fg_color=self.COLORS['accent_blue'],
                                  hover_color=self.COLORS['accent_blue_hover'],
                                  command=update_sites_table,
                                  font=("Roboto", 10, "bold"))
        btn_update.pack(side='left')
        
        table_frame = ctk.CTkScrollableFrame(parent, fg_color=self.COLORS['bg_feed'],
                                            corner_radius=8)
        table_frame.pack(fill='both', expand=True, padx=10, pady=10)
        
        gene_combo.configure(command=update_sites_table)
        method_combo.configure(command=update_sites_table)
        
        update_sites_table()
    
    def _render_sites_table(self, parent, gene_name: str, model_name: str, method: str, p_threshold: float = 0.95):
        """Renderiza tabela de sítios sob seleção"""
        
        # Para Branch-site, procurar ambas as variações (novo nome e nome antigo)
        search_patterns = [f"*{gene_name}*{model_name}*results.txt"]
        if model_name == 'Branch-site':
            search_patterns.append(f"*{gene_name}*BranchSite_A*results.txt")
        elif model_name == 'BranchSite_A':
            search_patterns.append(f"*{gene_name}*Branch-site*results.txt")
        
        # Buscar arquivo recursivamente
        results_file = []
        for pattern in search_patterns:
            results_file = list(self.output_folder.rglob(pattern))
            if results_file:
                break
        
        if not results_file:
            display_name = f"{gene_name}_{model_name}_results.txt"
            ctk.CTkLabel(parent,
                        text=TEXTS["sites_file_not_found"].format(filename=display_name),
                        font=("Roboto", 11),
                        text_color=self.COLORS['warning']).pack(pady=50)
            return
        
        # Parse com a classe SitesParser (se disponível)
        try:
            from src.backend.sites_parser import SitesParser
            df_sites = SitesParser.parse_sites_from_file(results_file[0], method=method)
            df_filtered = SitesParser.filter_sites_by_pvalue(df_sites, p_threshold)
            # Usar função robusta que tenta múltiplas estratégias
            omega_global = SitesParser.extract_omega_robust(results_file[0])
        except ImportError:
            # Fallback: parser manual básico
            df_sites, omega_global = self._parse_sites_manual(results_file[0], method)
            df_filtered = df_sites[df_sites['pr_w_gt_1'] >= p_threshold] if not df_sites.empty else pd.DataFrame()
        except Exception as e:
            ctk.CTkLabel(parent,
                        text=TEXTS["sites_parse_error"].format(error=str(e)),
                        font=("Roboto", 11),
                        text_color=self.COLORS['danger']).pack(pady=50)
            return

        if df_filtered.empty:
            ctk.CTkLabel(parent,
                        text=TEXTS["sites_no_sites"].format(threshold=p_threshold),
                        font=("Roboto", 11),
                        text_color=self.COLORS['text_tertiary']).pack(pady=50)
            return
        
        omega_text = f"ω = {omega_global:.4f}" if omega_global else "ω = N/A"

        # ── Info banner ─────────────────────────────────────────────
        banner = ctk.CTkFrame(parent, fg_color='#131326', corner_radius=8,
                              border_width=1, border_color=self.COLORS['accent_blue'])
        banner.pack(fill='x', padx=8, pady=(8, 10))

        banner_left = ctk.CTkFrame(banner, fg_color='transparent')
        banner_left.pack(side="left", padx=14, pady=10)

        ctk.CTkLabel(banner_left,
                     text=f"{gene_name}",
                     font=("Roboto", 13, "bold"),
                     text_color=self.COLORS['text_primary']).pack(anchor="w")
        ctk.CTkLabel(banner_left,
                     text=TEXTS["viewer_sites_subtitle"].format(model=model_name, method=method, omega=omega_text),
                     font=("Roboto", 9),
                     text_color=self.COLORS['text_tertiary']).pack(anchor="w")

        badge = ctk.CTkFrame(banner, fg_color=self.COLORS['accent_blue'],
                             corner_radius=8, width=80)
        badge.pack(side="right", padx=14, pady=10)
        badge.pack_propagate(False)
        ctk.CTkLabel(badge, text=TEXTS["viewer_sites_count"].format(n=len(df_filtered)),
                     font=("Roboto", 11, "bold"),
                     text_color="white").pack(expand=True)

        # ── Table header ─────────────────────────────────────────────
        th = ctk.CTkFrame(parent, fg_color='#1a1a26', corner_radius=6)
        th.pack(fill='x', padx=8, pady=(0, 2))

        widths  = [70, 55, 100, 55, 105, 105, 105]
        headers = list(zip(TEXTS["sites_table_headers"], widths))

        for h_text, width in headers:
            ctk.CTkLabel(th, text=h_text,
                         font=("Roboto", 11, "bold"),
                         text_color=self.COLORS['accent_blue'],
                         width=width).pack(side='left', padx=4, pady=8)

        # ── Rows ─────────────────────────────────────────────────────
        for i, (_, row) in enumerate(df_filtered.iterrows()):
            is_99 = row.get('is_significant_99', False)
            is_95 = row.get('is_significant_95', False)

            if is_99:
                bg_color    = '#0b2016'
                text_color  = '#6ee7b7'
                border_color = '#10b981'
                sig_text    = "**"
                sig_color   = '#10b981'
            elif is_95:
                bg_color    = '#0d1a10'
                text_color  = '#a7f3d0'
                border_color = '#059669'
                sig_text    = "*"
                sig_color   = '#6ee7b7'
            else:
                bg_color    = self.COLORS['bg_card'] if i % 2 == 0 else self.COLORS['bg_sidebar']
                text_color  = self.COLORS['text_secondary']
                border_color = self.COLORS['border']
                sig_text    = "—"
                sig_color   = self.COLORS['text_muted']

            row_frame = ctk.CTkFrame(parent, fg_color=bg_color, corner_radius=4,
                                     border_width=1, border_color=border_color)
            row_frame.pack(fill='x', padx=8, pady=1)

            cells = [
                (str(int(row['position'])), 70,  text_color),
                (row.get('amino_acid', 'X'), 55, text_color),
                (f"{row['pr_w_gt_1']:.4f}",  100, text_color),
                (sig_text,                   55,  sig_color),
                (f"{row['post_mean']:.3f}",  105, text_color),
                (f"{row.get('omega_lower',0):.3f}", 105, self.COLORS['text_tertiary']),
                (f"{row.get('omega_upper',0):.3f}", 105, self.COLORS['text_tertiary']),
            ]

            for cell_text, width, clr in cells:
                ctk.CTkLabel(row_frame, text=cell_text,
                             font=("Roboto", 11),
                             text_color=clr, width=width).pack(side='left', padx=4, pady=8)

        # ── Footer stats ──────────────────────────────────────────────
        footer = ctk.CTkFrame(parent, fg_color=self.COLORS['bg_sidebar'], corner_radius=6)
        footer.pack(fill='x', padx=8, pady=(10, 8))

        sig99 = int(df_filtered.get('is_significant_99', False).sum()) if 'is_significant_99' in df_filtered else 0
        sig95 = int(df_filtered.get('is_significant_95', False).sum()) if 'is_significant_95' in df_filtered else 0

        stats_parts = [
            f"ω médio: {df_filtered['post_mean'].mean():.3f}",
            f"ω max: {df_filtered['post_mean'].max():.3f}",
            f"Pr(ω>1) médio: {df_filtered['pr_w_gt_1'].mean():.3f}",
            f"** (p>=0.99): {sig99}",
            f"* (p>=0.95): {sig95}",
        ]
        ctk.CTkLabel(footer, text="  ·  ".join(stats_parts),
                     font=("Roboto", 10),
                     text_color=self.COLORS['text_tertiary']).pack(pady=8, padx=12)
    
    def _parse_sites_manual(self, filepath: Path, method: str):
        """Parser manual básico caso SitesParser não esteja disponível"""
        df_sites = pd.DataFrame()
        omega_global = None
        
        try:
            with open(filepath, 'r') as f:
                content = f.read()
            
            # Extrair ω global - usar função robusta de extração
            try:
                from src.backend.sites_parser import SitesParser
                omega_global = SitesParser.extract_omega_robust(str(filepath))
            except Exception:
                # Fallback: padrão direto (para compatibilidade)
                omega_match = re.search(r'omega \(dN/dS\)\s*=\s*([\d.]+)', content)
                if omega_match:
                    omega_global = float(omega_match.group(1))
            
            # Extrair sites (pattern simplificado)
            if method == 'BEB':
                pattern = r'(\d+)\s+([A-Z])\s+([\d.]+)\*{0,2}\s+([\d.]+)\+?-\s+([\d.]+)'
            else:
                pattern = r'(\d+)\s+([A-Z])\s+([\d.]+)'
            
            sites_data = []
            for match in re.finditer(pattern, content):
                if method == 'BEB':
                    pos, aa, prob, mean, se = match.groups()
                    sites_data.append({
                        'position': int(pos),
                        'amino_acid': aa,
                        'pr_w_gt_1': float(prob),
                        'post_mean': float(mean),
                        'omega_lower': float(mean) - float(se),
                        'omega_upper': float(mean) + float(se),
                        'is_significant_95': float(prob) >= 0.95,
                        'is_significant_99': float(prob) >= 0.99
                    })
                else:
                    pos, aa, omega = match.groups()
                    sites_data.append({
                        'position': int(pos),
                        'amino_acid': aa,
                        'pr_w_gt_1': 1.0 if float(omega) > 1 else 0.0,
                        'post_mean': float(omega),
                        'is_significant_95': False,
                        'is_significant_99': False
                    })
            
            df_sites = pd.DataFrame(sites_data)
        except Exception as e:
            print(f"[ERR] Erro no parser manual: {e}")
        
        return df_sites, omega_global
    
    def _create_branchsite_class_tab(self, parent):
        """Aba de visualizacao estruturada de classes de Branch-site"""
        ctrl_frame = ctk.CTkFrame(parent, fg_color=self.COLORS['bg_card'],
                                 corner_radius=8, height=60)
        ctrl_frame.pack(fill='x', padx=10, pady=10)
        ctrl_frame.pack_propagate(False)
        
        ctk.CTkLabel(ctrl_frame, text=TEXTS["branchsite_classes_gene_label"], font=("Roboto", 11, "bold")).pack(side='left', padx=15, pady=10)
        
        genes = self.df['Gene'].tolist()
        gene_combo = ctk.CTkComboBox(ctrl_frame, values=genes, width=300)
        gene_combo.pack(side='left', padx=(0, 20))
        gene_combo.set(genes[0] if genes else "")
        
        table_frame = ctk.CTkScrollableFrame(parent, fg_color=self.COLORS['bg_feed'],
                                            corner_radius=8)
        table_frame.pack(fill='both', expand=True, padx=10, pady=10)
        
        def update_branchsite_table(*args):
            for widget in table_frame.winfo_children():
                widget.destroy()
            
            selected_gene = gene_combo.get()
            gene_row = self.df[self.df['Gene'] == selected_gene]
            
            if gene_row.empty:
                ctk.CTkLabel(table_frame, text=TEXTS["branchsite_classes_not_found"],
                           font=("Roboto", 11),
                           text_color=self.COLORS['warning']).pack(pady=50)
                return
            
            idx = gene_row.index[0]
            self._render_branchsite_class_table(table_frame, idx)
        
        gene_combo.configure(command=update_branchsite_table)
        update_branchsite_table()
    
    def _render_branchsite_class_table(self, parent, gene_idx: int):
        """Renderiza tabela estruturada de classes de Branch-site"""
        row = self.df.iloc[gene_idx]
        gene = row['Gene']
        
        # Header
        header_frame = ctk.CTkFrame(parent, fg_color=self.COLORS['bg_card_hover'],
                                   corner_radius=8)
        header_frame.pack(fill='x', padx=8, pady=(8, 12))
        
        ctk.CTkLabel(header_frame, text=TEXTS["branchsite_classes_header"].format(gene=gene),
                    font=("Roboto", 12, "bold"),
                    text_color=self.COLORS['accent_cyan']).pack(pady=8)
        
        # Table header
        table_header_frame = ctk.CTkFrame(parent, fg_color=self.COLORS['bg_card_hover'],
                                         corner_radius=6)
        table_header_frame.pack(fill='x', padx=8, pady=(0, 4))
        
        bs_widths = [120, 150, 150, 150]
        headers = list(zip(TEXTS["branchsite_classes_table_headers"], bs_widths))

        for h_text, width in headers:
            ctk.CTkLabel(table_header_frame, text=h_text,
                        font=("Roboto", 11, "bold"),
                        text_color=self.COLORS['accent_blue'],
                        width=width).pack(side='left', padx=8, pady=8)
        
        # Dados
        for cls in ['0', '1', '2a', '2b']:
            class_col = f'Branch-site_class{cls}_fg_w'
            prop_col = f'Branch-site_class{cls}_prop'
            
            if class_col not in self.df.columns or prop_col not in self.df.columns:
                continue
            
            fg_w = row[class_col]
            prop = row[prop_col]
            
            # Also get background w from omega values
            # For Branch-site we can parse from the file, but for now use omega
            # Note: This is simplified, could be improved with full parsing
            
            row_frame = ctk.CTkFrame(parent, fg_color=self.COLORS['bg_card'],
                                    corner_radius=4, border_width=1,
                                    border_color=self.COLORS['bg_card_hover'])
            row_frame.pack(fill='x', padx=8, pady=2)
            
            # Format values
            if isinstance(prop, float):
                prop_str = f"{prop:.5f}"
            else:
                prop_str = "N/A"
            
            if isinstance(fg_w, float):
                fg_w_str = f"{fg_w:.5f}"
            else:
                fg_w_str = "N/A"
            
            # Background w (would need full parsing - simplified here)
            bg_w_str = TEXTS["branchsite_classes_bg_w_placeholder"]
            
            cells = [
                (f"Class {cls}", 120),
                (prop_str, 150),
                (bg_w_str, 150),
                (fg_w_str, 150)
            ]
            
            for cell_text, width in cells:
                ctk.CTkLabel(row_frame, text=cell_text, font=("Roboto", 11),
                           text_color=self.COLORS['text_secondary'], width=width).pack(side='left', padx=8, pady=8)
        
        # Footer with interpretation
        footer_frame = ctk.CTkFrame(parent, fg_color=self.COLORS['bg_card_hover'],
                                   corner_radius=6)
        footer_frame.pack(fill='x', padx=8, pady=(12, 8))
        
        ctk.CTkLabel(footer_frame, text=TEXTS["branchsite_classes_footer"],
                    font=("Roboto", 9),
                    text_color=self.COLORS['text_tertiary'],
                    wraplength=400).pack(pady=8, padx=8)
    
    def _create_tree_tab(self, parent):
        """Branch Analysis — Cladograma vertical com ramos coloridos por dN/dS."""
        has_branch_omega = ('Branch_omega' in self.df.columns and
                            self.df['Branch_omega'].notna().any())
        has_branchsite   = ('Branch-site_omega' in self.df.columns and
                            self.df['Branch-site_omega'].notna().any())
        has_branch_model = has_branch_omega or bool(
            self.tag_columns.get('Branch', {}).get('omega'))

        # ── Info banner ────────────────────────────────────────────────
        info = ctk.CTkFrame(parent, fg_color='#1a1a26', corner_radius=8)
        info.pack(fill='x', padx=10, pady=(10, 4))
        ctk.CTkLabel(info,
                     text=TEXTS["branch_tab_title"],
                     font=("Roboto", 11, "bold"),
                     text_color=self.COLORS['accent_blue']).pack(side="left", padx=14, pady=(10, 2))
        ctk.CTkLabel(info,
                     text=TEXTS["branch_tab_legend"],
                     font=("Roboto", 9),
                     text_color=self.COLORS['text_tertiary']).pack(side="left", padx=(0, 14), pady=(10, 2))

        if not has_branch_model and not has_branchsite:
            empty = ctk.CTkFrame(parent, fg_color='transparent')
            empty.pack(expand=True)
            ctk.CTkLabel(empty, text=TEXTS["branch_no_data_title"],
                         font=("Roboto", 14, "bold"),
                         text_color=self.COLORS['text_tertiary']).pack(pady=(60, 6))
            ctk.CTkLabel(empty,
                         text=TEXTS["branch_no_data_hint"],
                         font=("Roboto", 10),
                         text_color=self.COLORS['text_muted']).pack()
            return

        # ── Gene + Outgroup selector ───────────────────────────────────
        ctrl = ctk.CTkFrame(parent, fg_color=self.COLORS['bg_card'], corner_radius=8)
        ctrl.pack(fill='x', padx=10, pady=(0, 4))

        ctk.CTkLabel(ctrl, text=TEXTS["branch_label_gene"], font=("Roboto", 11, "bold")).pack(
            side='left', padx=(15, 4), pady=10)

        if has_branch_omega:
            valid_genes = self.df[self.df['Branch_omega'].notna()]['Gene'].tolist()
        elif has_branchsite:
            valid_genes = self.df[self.df['Branch-site_omega'].notna()]['Gene'].tolist()
        else:
            valid_genes = self.df['Gene'].tolist()
        if not valid_genes:
            valid_genes = self.df['Gene'].tolist()

        gene_combo = ctk.CTkComboBox(ctrl, values=valid_genes, width=260)
        gene_combo.pack(side='left', padx=(0, 18), pady=10)
        if valid_genes:
            gene_combo.set(valid_genes[0])

        ctk.CTkLabel(ctrl, text=TEXTS["branch_label_outgroup"], font=("Roboto", 11, "bold")).pack(
            side='left', padx=(0, 4), pady=10)
        outgroup_combo = ctk.CTkComboBox(ctrl, values=[TEXTS["branch_outgroup_none"]], width=210)
        outgroup_combo.pack(side='left', padx=(0, 16), pady=10)
        outgroup_combo.set(TEXTS["branch_outgroup_none"])

        lrt_col  = 'lrt_M0_vs_Branch'
        has_lrt  = lrt_col in self.df.columns
        info_lbl = ctk.CTkLabel(ctrl, text="", font=("Roboto", 9),
                                text_color=self.COLORS['text_tertiary'])
        info_lbl.pack(side='left', padx=10, pady=10)

        # ── PNG export button ──────────────────────────────────────────
        current_fig = [None]

        def _export_png():
            if current_fig[0] is None:
                messagebox.showwarning(TEXTS["msg_warning"], TEXTS["msg_no_figure"])
                return
            fp = filedialog.asksaveasfilename(
                title=TEXTS["dialog_export_cladogram"],
                defaultextension=".png",
                filetypes=[("PNG", "*.png"), ("Todos os arquivos", "*.*")]
            )
            if fp:
                try:
                    current_fig[0].savefig(fp, dpi=200, bbox_inches='tight',
                                           facecolor='#111115')
                    messagebox.showinfo(TEXTS["msg_success"], TEXTS["msg_exported_to"].format(path=fp))
                except Exception as e:
                    messagebox.showerror(TEXTS["msg_error"], TEXTS["msg_export_err"].format(error=e))

        btn_bar = ctk.CTkFrame(parent, fg_color='transparent')
        btn_bar.pack(fill='x', padx=10, pady=(0, 4))
        ctk.CTkButton(btn_bar, text=TEXTS["branch_btn_export_png"],
                      font=("Roboto", 10),
                      fg_color=self.COLORS['accent_blue'],
                      hover_color=self.COLORS['accent_blue_hover'],
                      width=130, height=28,
                      command=_export_png).pack(side='right')

        # ── Chart frame ───────────────────────────────────────────────
        chart_frame = ctk.CTkFrame(parent, fg_color=self.COLORS['bg_feed'], corner_radius=8)
        chart_frame.pack(fill='both', expand=True, padx=10, pady=(0, 10))

        # ── Shared mutable state (gene, outgroup, rotated nodes) ───────
        state = {
            'gene':          valid_genes[0] if valid_genes else None,
            'outgroup':      TEXTS["branch_outgroup_none"],
            'rotated':       set(),   # set of node_ids whose children are reversed
            'node_pos':      {},      # node_id → (data_x, data_y) for click detection
            'internals':     set(),   # set of internal node_ids
        }

        # ── Helper: species map ────────────────────────────────────────
        def _parse_species_map(filepath):
            mapping = {}
            try:
                with open(filepath, 'r', encoding='utf-8', errors='ignore') as fh:
                    content = fh.read()
                m = re.search(r'^\s*(\d+)\s+\d+\s*$', content, re.MULTILINE)
                if not m:
                    return mapping
                ntaxa = int(m.group(1))
                rest  = content[m.end():].lstrip('\n')
                idx   = 1
                for line in rest.split('\n'):
                    line = line.strip()
                    if not line:
                        continue
                    first = line.split()[0]
                    if first in ('Printing', 'CODONML', 'BASEML', 'lnL',
                                 'tree', 'TREE', 'Sequence'):
                        break
                    if re.match(r'^[A-Za-z][A-Za-z0-9._\-]*$', first):
                        mapping[idx] = first
                        idx += 1
                    if idx > ntaxa:
                        break
            except Exception:
                pass
            return mapping

        # ── Helper: re-root ────────────────────────────────────────────
        def _reroot(orig_children_map, outgroup_leaf):
            """Re-root at parent of outgroup_leaf.
            Preserves bifurcating structure: if the new root would end up with >2
            children (because its original siblings get merged with the path-up
            node), we move those siblings under the 'path-up' node so the root
            always has exactly 2 children: [outgroup, rest_of_tree].
            """
            adj = {}
            for p, kids in orig_children_map.items():
                for c, t, w in kids:
                    adj.setdefault(p, []).append((c, t, w))
                    adj.setdefault(c, []).append((p, t, w))

            # Find the direct parent of outgroup_leaf in the original tree
            new_root = None
            for p, kids in orig_children_map.items():
                for c, t, w in kids:
                    if c == outgroup_leaf:
                        new_root = p
                        break
                if new_root is not None:
                    break
            if new_root is None:
                return orig_children_map, next(iter(orig_children_map))

            # Standard BFS re-root (reverses all edges away from new_root)
            new_ch  = {}
            visited = {new_root}
            queue   = [new_root]
            while queue:
                node = queue.pop(0)
                new_ch.setdefault(node, [])
                for nbr, t, w in adj.get(node, []):
                    if nbr not in visited:
                        visited.add(nbr)
                        new_ch[node].append((nbr, t, w))
                        queue.append(nbr)

            # Fix potential polytomy at new_root:
            # BFS gives new_root children = [outgroup, original_siblings..., original_parent].
            # For a bifurcating tree we want new_root = [outgroup, original_parent],
            # and original_parent absorbs the original_siblings.
            root_kids    = new_ch.get(new_root, [])
            og_kid       = [(c, t, w) for c, t, w in root_kids if c == outgroup_leaf]
            non_og       = [(c, t, w) for c, t, w in root_kids if c != outgroup_leaf]

            if len(non_og) > 1:
                # Identify which child was the original parent (path-up)
                orig_ch_ids = {c for c, _, _ in orig_children_map.get(new_root, [])}
                path_up  = [(c, t, w) for c, t, w in non_og if c not in orig_ch_ids]
                siblings = [(c, t, w) for c, t, w in non_og if c in orig_ch_ids]
                if path_up and siblings:
                    path_node = path_up[0][0]
                    # Move original siblings under the path-up node
                    new_ch[path_node].extend(siblings)
                    # New root now has exactly 2 children: outgroup + path-up
                    new_ch[new_root] = og_kid + path_up

            return new_ch, new_root

        # ── Render (horizontal cladogram: root at left, tips at right) ──
        def render_tree():
            gene_name     = state['gene']
            outgroup_name = state['outgroup']
            rotated_nodes = state['rotated']

            old = current_fig[0]
            for w in chart_frame.winfo_children():
                w.destroy()
            if old is not None:
                plt.close(old)
            current_fig[0] = None

            if not gene_name:
                return

            row_data = self.df[self.df['Gene'] == gene_name]
            if row_data.empty:
                ctk.CTkLabel(chart_frame, text=TEXTS["viewer_gene_not_found"],
                             text_color=self.COLORS['text_tertiary']).pack(expand=True)
                return
            row = row_data.iloc[0]

            if has_lrt:
                lrt_val = row.get(lrt_col, np.nan)
                if pd.notna(lrt_val):
                    # df = (np_Branch − np_M0) − (ntime_Branch − ntime_M0)
                    # Correção necessária porque M0 usa árvore não-enraizada (ntime=2n-3)
                    # enquanto Branch usa a árvore rotulada/enraizada (ntime=2n-2).
                    m0_np_val        = row.get('M0_np', np.nan)
                    branch_np_val    = row.get('Branch_np', np.nan)
                    m0_ntime_val     = row.get('M0_ntime', np.nan)
                    branch_ntime_val = row.get('Branch_ntime', np.nan)
                    if pd.notna(m0_np_val) and pd.notna(branch_np_val):
                        raw_df = abs(int(branch_np_val) - int(m0_np_val))
                        if pd.notna(m0_ntime_val) and pd.notna(branch_ntime_val):
                            df_branch = max(1, raw_df - (int(branch_ntime_val) - int(m0_ntime_val)))
                        else:
                            df_branch = max(1, raw_df)
                    else:
                        df_branch = 1
                    p   = 1 - stats.chi2.cdf(lrt_val, df=df_branch) if lrt_val > 0 else 1.0
                    sig = "  * p < 0.05" if p < 0.05 else ""
                    info_lbl.configure(
                        text=f"LRT (M0 vs Branch): 2Df = {lrt_val:.3f}  ·  df = {df_branch}  ·  p = {self._fmt_pval(p)}{sig}",
                        text_color=self.COLORS['success'] if p < 0.05 else self.COLORS['text_tertiary']
                    )

            results_file = self._find_results_file(gene_name, 'Branch')
            if not results_file:
                ctk.CTkLabel(chart_frame,
                             text=TEXTS["viewer_branch_no_file"],
                             text_color=self.COLORS['text_tertiary']).pack(expand=True)
                return
            try:
                df_br = BranchExtractor.extract_branch_table(results_file)
            except Exception as exc:
                ctk.CTkLabel(chart_frame, text=TEXTS["viewer_branch_read_err"].format(error=exc),
                             text_color=self.COLORS['text_tertiary']).pack(expand=True)
                return
            if df_br.empty:
                ctk.CTkLabel(chart_frame,
                             text=TEXTS["viewer_branch_no_table"],
                             text_color=self.COLORS['text_tertiary']).pack(expand=True)
                return

            species_map = _parse_species_map(results_file)

            # Refresh outgroup combo values
            sp_names = [TEXTS["branch_outgroup_none"]] + sorted(species_map.values())
            outgroup_combo.configure(values=sp_names)
            if outgroup_name not in sp_names:
                outgroup_combo.set(TEXTS["branch_outgroup_none"])
                state['outgroup'] = TEXTS["branch_outgroup_none"]
                outgroup_name = TEXTS["branch_outgroup_none"]

            # Undirected omega lookup (works after any re-root)
            branch_omega_map = {}
            for _, r in df_br.iterrows():
                try:
                    ps, cs = str(r['branch']).split('..')
                    a, b   = int(ps), int(cs)
                    w      = float(r['dN_dS'])
                    branch_omega_map[(a, b)] = w
                    branch_omega_map[(b, a)] = w
                except Exception:
                    continue

            # Build directed topology
            children_map = {}
            children_set = set()
            all_nodes    = set()
            for _, r in df_br.iterrows():
                try:
                    ps, cs = str(r['branch']).split('..')
                    p_node, c_node = int(ps), int(cs)
                except Exception:
                    continue
                all_nodes.update([p_node, c_node])
                children_map.setdefault(p_node, []).append(
                    (c_node, float(r['t']), float(r['dN_dS'])))
                children_set.add(c_node)

            if not all_nodes:
                ctk.CTkLabel(chart_frame, text=TEXTS["viewer_branch_invalid"],
                             text_color=self.COLORS['text_tertiary']).pack(expand=True)
                return

            root = next(n for n in all_nodes if n not in children_set)

            # Re-root if outgroup selected; place outgroup at BOTTOM (first in DFS)
            if outgroup_name != TEXTS["branch_outgroup_none"]:
                og_leaf = next((k for k, v in species_map.items()
                                if v == outgroup_name), None)
                if og_leaf is not None and og_leaf in all_nodes:
                    children_map, root = _reroot(children_map, og_leaf)
                    all_nodes    = set()
                    children_set = set()
                    for p, kids in children_map.items():
                        all_nodes.add(p)
                        for c, t, w in kids:
                            all_nodes.add(c)
                            children_set.add(c)
                    # Move outgroup child first → DFS puts it at bottom (Y=0)
                    root_kids = children_map.get(root, [])
                    og_entry  = next((e for e in root_kids if e[0] == og_leaf), None)
                    if og_entry:
                        root_kids.remove(og_entry)
                        root_kids.insert(0, og_entry)
                        children_map[root] = root_kids

            # Apply rotations: flip children order at toggled internal nodes
            for nid in rotated_nodes:
                if nid in children_map and children_map[nid]:
                    children_map[nid] = list(reversed(children_map[nid]))

            # DFS leaf ordering
            tip_order = []
            def _dfs(node):
                kids = children_map.get(node, [])
                if not kids:
                    tip_order.append(node)
                else:
                    for child, _, _ in kids:
                        _dfs(child)
            _dfs(root)

            n_leaves = len(tip_order)

            # Y: evenly spaced leaves; internal = midpoint of children Y
            leaf_y = {leaf: float(i) for i, leaf in enumerate(tip_order)}
            node_y = {}
            def _cy(node):
                if node in leaf_y:
                    node_y[node] = leaf_y[node]
                    return leaf_y[node]
                ys = [_cy(c) for c, _, _ in children_map.get(node, [])]
                node_y[node] = (min(ys) + max(ys)) / 2.0
                return node_y[node]
            _cy(root)
            for n in all_nodes:
                if n not in node_y:
                    node_y[n] = 0.0

            # X: depth from root (root=0 on left, tips=max_depth on right)
            # Leaves are forced to max_depth so all tips align flush right
            # (cladogram style — like ITOL "Ignore branch lengths")
            depth_fr = {}
            def _dfr(node, d=0):
                depth_fr[node] = d
                for c, t, _ in children_map.get(node, []):
                    _dfr(c, d + 1)
            _dfr(root)
            max_depth = max(depth_fr.values()) if depth_fr else 1
            is_leaf = {n for n in all_nodes
                       if not children_map.get(n)}
            node_x = {n: float(max_depth) if n in is_leaf
                         else float(depth_fr.get(n, 0))
                      for n in all_nodes}

            # Store positions for click detection
            state['node_pos']  = {n: (node_x.get(n, 0), node_y.get(n, 0))
                                  for n in all_nodes}
            state['internals'] = {n for n, kids in children_map.items() if kids}

            # Color map
            all_omegas = list(branch_omega_map.values())
            vmax = max(max(all_omegas, default=2.0), 2.0)
            cmap = LinearSegmentedColormap.from_list(
                'omega_ramp', ['#ef4444', '#fbbf24', '#3b82f6'])
            norm = TwoSlopeNorm(vmin=0.0, vcenter=1.0, vmax=vmax)

            # Figure size: height by #leaves, width fixed
            fig_h = max(4.5, n_leaves * 0.28)
            fig_w = max(8.0, max_depth * 1.2 + 5.5)
            fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor='#111115')
            ax.set_facecolor('#111115')
            LW      = 2.0
            done_vc = set()

            # Build parent map for tip ω lookup
            parent_of = {}
            for p_node, kids in children_map.items():
                for c_node, t, _ in kids:
                    parent_of[c_node] = p_node

            # Draw horizontal branches + vertical connectors
            for p_node, kids in children_map.items():
                if not kids:
                    continue
                px = node_x[p_node]
                for c_node, t, _ in kids:
                    cx    = node_x[c_node]
                    cy    = node_y[c_node]
                    omega = branch_omega_map.get((p_node, c_node), 0.5)
                    color = cmap(norm(omega))
                    ax.plot([px, cx], [cy, cy],       # horizontal branch
                            color=color, linewidth=LW,
                            solid_capstyle='round', zorder=2)
                if p_node not in done_vc:
                    child_ys = [node_y[c] for c, _, _ in kids]
                    # Color connector by the incoming branch (parent → p_node)
                    gp = parent_of.get(p_node)
                    if gp is not None:
                        conn_omega = branch_omega_map.get((gp, p_node), 0.5)
                        conn_color = cmap(norm(conn_omega))
                    else:
                        conn_color = '#5a5a6a'  # root has no incoming branch
                    ax.plot([px, px],
                            [min(child_ys), max(child_ys)],   # vertical connector
                            color=conn_color, linewidth=LW,
                            solid_capstyle='round', zorder=1)
                    done_vc.add(p_node)

            # Internal node markers (click targets; square = rotated, circle = normal)
            for nid in state['internals']:
                nx, ny = node_x.get(nid, 0), node_y.get(nid, 0)
                mk = 's' if nid in rotated_nodes else 'o'
                ax.scatter([nx], [ny], color='#2a2a36', s=48, zorder=5,
                           edgecolors='#5a5a6a', linewidths=0.8, marker=mk)

            # Tip dots + labels on the right (ω value + species name)
            for leaf in tip_order:
                lx = node_x[leaf]   # = max_depth
                ly = node_y[leaf]

                p     = parent_of.get(leaf)
                omega = branch_omega_map.get((p, leaf)) if p is not None else None

                if omega is not None:
                    dot_color = cmap(norm(omega))
                    ax.scatter([lx], [ly], color=dot_color, s=55, zorder=6,
                               edgecolors='#222228', linewidths=0.6)
                    omega_str = f'{omega:.3f}  '
                else:
                    omega_str = ''

                name = species_map.get(leaf, str(leaf))
                if len(name) > 30:
                    name = name[:27] + '...'
                label = omega_str + name
                ax.text(max_depth + 0.15, ly, label,
                        ha='left', va='center',
                        color='#c8c8d4', fontsize=7.5, fontfamily='monospace')

            # Colorbar (horizontal, bottom-left)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax, orientation='horizontal',
                                fraction=0.04, pad=0.06, shrink=0.30,
                                anchor=(0.0, 1.0))
            cbar.set_label('w (dN/dS)', color='#9898a6', fontsize=8)
            cbar.ax.tick_params(colors='#9898a6', labelsize=7)
            # Smart tick generator: always include 0.0, 1.0, and vmax;
            # spread intermediate ticks proportionally to the actual range.
            if vmax <= 3.0:
                tick_vals = sorted({0.0, 0.5, 1.0, round(vmax * 0.75, 2), vmax})
            elif vmax <= 10.0:
                step = vmax / 4.0
                tick_vals = sorted({0.0, 1.0,
                                    round(step, 1), round(step * 2, 1),
                                    round(step * 3, 1), round(vmax, 1)})
            else:
                # Large range: pick ~4 round intermediates between 1 and vmax
                import math as _math
                magnitude = 10 ** _math.floor(_math.log10(vmax))
                step = magnitude / 2 if vmax / magnitude < 3 else magnitude
                intermediates = []
                v = step
                while v < vmax:
                    if v > 1.0:
                        intermediates.append(round(v, 1) if step < 1 else int(round(v)))
                    v += step
                tick_vals = sorted({0.0, 1.0, *intermediates, round(vmax, 1)})
                # Cap at 6 ticks to avoid crowding: keep 0, 1, last-3, vmax
                if len(tick_vals) > 6:
                    keep = sorted({tick_vals[0], tick_vals[1],
                                   *tick_vals[-4:]})
                    tick_vals = keep
            tick_fmt = [f'{v:.0f}' if v == int(v) else f'{v:.1f}'
                        for v in tick_vals]
            cbar.set_ticks(tick_vals)
            cbar.set_ticklabels(tick_fmt)

            # Axis limits and style
            ax.set_xlim(-0.3, max_depth + 4.5)
            ax.set_ylim(-0.5, n_leaves - 0.5)
            ax.set_yticks([])
            ax.set_xticks([])
            ax.set_title(f'Cladograma  —  {gene_name}',
                         color='#ededef', fontsize=11, fontweight='bold', pad=8)
            for spine in ax.spines.values():
                spine.set_visible(False)

            plt.tight_layout(pad=0.5)

            # Click handler: detect click near internal node → toggle rotation
            def on_click(event):
                if event.inaxes != ax or event.xdata is None:
                    return
                ex, ey = event.xdata, event.ydata
                best, best_d = None, float('inf')
                for nid in state['internals']:
                    nx, ny = state['node_pos'].get(nid, (0, 0))
                    if abs(ex - nx) > 0.8 or abs(ey - ny) > 1.2:
                        continue
                    d = ((ex - nx) ** 2 + (ey - ny) ** 2) ** 0.5
                    if d < best_d:
                        best_d, best = d, nid
                if best is not None:
                    if best in state['rotated']:
                        state['rotated'].discard(best)
                    else:
                        state['rotated'].add(best)
                    render_tree()

            current_fig[0] = fig
            canvas_w = FigureCanvasTkAgg(fig, master=chart_frame)
            canvas_w.mpl_connect('button_press_event', on_click)
            canvas_w.draw()
            canvas_w.get_tk_widget().pack(fill='both', expand=True)

        def _on_gene(v):
            state['gene']    = v
            state['outgroup'] = TEXTS["branch_outgroup_none"]
            state['rotated']  = set()
            outgroup_combo.set(TEXTS["branch_outgroup_none"])
            render_tree()

        def _on_outgroup(v):
            state['outgroup'] = v
            state['rotated']  = set()
            render_tree()

        gene_combo.configure(command=_on_gene)
        outgroup_combo.configure(command=_on_outgroup)
        if valid_genes:
            render_tree()
    
    def _create_export_tab(self, parent):
        """Aba de exportação"""
        main_frame = ctk.CTkScrollableFrame(parent, fg_color=self.COLORS['bg_feed'],
                                           corner_radius=10)
        main_frame.pack(fill='both', expand=True, padx=15, pady=15)
        
        header_frame = ctk.CTkFrame(main_frame, fg_color='transparent')
        header_frame.pack(fill='x', pady=(0, 30))
        
        ctk.CTkLabel(header_frame, text=TEXTS["export_tab_title"],
                    font=("Roboto", 18, "bold"),
                    text_color=self.COLORS['text_primary']).pack()

        _export_callbacks = [
            self._export_excel,
            self._export_csv,
            self._export_charts,
            self._export_html,
        ]
        _export_colors = [
            self.COLORS['success'],
            self.COLORS['accent_blue'],
            self.COLORS['warning'],
            self.COLORS['accent_cyan'],
        ]
        export_options = [
            (title, desc, cmd, color)
            for (title, desc), cmd, color
            in zip(TEXTS["export_options"], _export_callbacks, _export_colors)
        ]

        for title, desc, command, color in export_options:
            card = ctk.CTkFrame(main_frame, fg_color=self.COLORS['bg_card'],
                                corner_radius=12, border_width=1,
                                border_color=self.COLORS['border'])
            card.pack(fill='x', pady=8)

            # Top strip
            ctk.CTkFrame(card, fg_color=color, height=3,
                         corner_radius=2).pack(fill='x')

            row = ctk.CTkFrame(card, fg_color='transparent')
            row.pack(fill='x', padx=20, pady=14)

            txt = ctk.CTkFrame(row, fg_color='transparent')
            txt.pack(side="left", fill='both', expand=True)
            ctk.CTkLabel(txt, text=title, font=("Roboto", 12, "bold"),
                         text_color=color, anchor='w').pack(anchor='w')
            ctk.CTkLabel(txt, text=desc, font=("Roboto", 9),
                         text_color=self.COLORS['text_muted'], anchor='w').pack(anchor='w', pady=(2, 0))

            ctk.CTkButton(row, text=TEXTS["export_btn"], width=130, height=34,
                          fg_color=self.COLORS['border'],
                          hover_color=color,
                          text_color=color,
                          border_width=1, border_color=color,
                          font=("Roboto", 10, "bold"),
                          corner_radius=8,
                          command=command).pack(side="right")
    
    # ═══════════════════════════════════════════════════════════
    # MÉTODOS AUXILIARES
    # ═══════════════════════════════════════════════════════════
    
    def _detect_positive_selection(self) -> dict:
        """
        Detecta seleção positiva baseado em:
        1. ω > 1.0 (evidência de seleção positiva)
        2. p-valor < 0.05 (significância estatística do LRT)
        
        Apenas M2a e M8 retornam dados de sítios sob seleção
        """
        positive_data = {}
        
        # Definir comparações de modelos (alternativo vs nulo)
        # Apenas M2a e M8 têm dados de sítios BEB/NEB
        comparisons = {
            'M2a': ('M2a_omega', 'lrt_M1a_vs_M2a', 2),  # (coluna_omega, coluna_lrt, df)
            'M8': ('M8_omega', 'lrt_M7_vs_M8', 2),
        }
        
        for idx, row in self.df.iterrows():
            gene_signals = {}
            
            # Verificar cada modelo
            for model_name, (omega_col, lrt_col, df) in comparisons.items():
                # Verificar se as colunas existem
                if omega_col not in self.df.columns or lrt_col not in self.df.columns:
                    continue
                
                omega_val = row[omega_col]
                lrt_val = row[lrt_col]
                
                # Verificar se omega > 1
                if pd.isna(omega_val) or omega_val <= 1.0:
                    continue
                
                # Verificar se LRT é válido
                if pd.isna(lrt_val):
                    continue
                
                # Calcular p-valor
                p_val = 1 - stats.chi2.cdf(lrt_val, df=df) if lrt_val > 0 else 1.0
                
                # Se significante, adicionar aos sinais
                if p_val < 0.05:
                    gene_signals[model_name] = {
                        'omega': omega_val,
                        'p_value': p_val,
                        'lrt': lrt_val
                    }
            
            # Se encontrou sinais, adicionar ao resultado
            if gene_signals:
                positive_data[row['Gene']] = gene_signals
        
        return positive_data
    
    def _count_models(self) -> str:
        """Conta modelos únicos"""
        model_cols = [col for col in self.df.columns if '_lnL' in col or '_omega' in col]
        unique_models = set()
        
        for col in model_cols:
            # Extrair nome do modelo removendo sufixos (lnL, omega, np, time, stops)
            # Ex: M8_omega -> M8, Branch-site_omega -> Branch-site, Branch_lnL -> Branch
            model_name = col
            for suffix in ['_lnL', '_omega', '_np', '_time', '_stops']:
                if suffix in model_name:
                    model_name = model_name.replace(suffix, '')
                    break
            
            if model_name:
                unique_models.add(model_name)
        
        return str(len(unique_models))
    
    def _calc_avg_omega(self) -> float:
        """Calcula ω médio de todos os modelos"""
        omega_cols = [col for col in self.df.columns if '_omega' in col]
        if not omega_cols:
            return 0.0
        
        omegas = []
        for col in omega_cols:
            vals = pd.to_numeric(self.df[col], errors='coerce').dropna()
            omegas.extend(vals.tolist())
        
        return np.mean(omegas) if omegas else 0.0
    
    def _get_available_lrt_columns(self):
        """Retorna (comparisons, descriptions):
           comparisons  = {label → col_name}
           descriptions = {col_name → null-hyp text shown below the combo}
        """
        KNOWN = {
            'lrt_M0_vs_M1a': (
                'M0 → M1a   (neutralidade)',
                'H₀  M0 — taxa ω única para todos os sítios  ·  '
                'H₁  M1a — ω₀ < 1 e ω₁ = 1 (neutralidade quase-neutra)   ·   df = 2   ·  '
                'Pré-teste; M1a vs M2a é o teste principal de seleção positiva',
            ),
            'lrt_M1a_vs_M2a': (
                'M1a → M2a   (sítios positivos)',
                'H₀  M1a — apenas purificação/neutralidade (ω ≤ 1)  ·  '
                'H₁  M2a — sítios com ω > 1   ·   df = 2   ·   '
                'Detecta seleção positiva em sítios ao longo de todos os ramos',
            ),
            'lrt_M7_vs_M8': (
                'M7 → M8   (Beta + ω > 1)',
                'H₀  M7 — distribuição Beta restrita a 0 < ω < 1  ·  '
                'H₁  M8 — Beta + classe com ω ≥ 1   ·   df = 2',
            ),
            'lrt_M0_vs_Branch': (
                'M0 → Branch   (ramos livres)',
                'H₀  M0 — uma única taxa ω para todos os ramos  ·  '
                'H₁  Branch — ω independente por grupo marcado   ·  '
                'df = n° de grupos foreground marcados (#1, #2 …)   ·  '
                'Ex: 1 marca → df=1 (χ²crit=3.84)  |  5 marcas → df=5 (χ²crit=11.07)   ·  '
                '[!]  Requer > 200 pb',
            ),
            'lrt_Branch-site_null_vs_Branch-site': (
                'Branch-site null → Branch-site   (seleção episódica)',
                'H₀  Branch-site null — ω ≤ 1 no foreground  ·  '
                'H₁  Branch-site — sítios com ω > 1 no foreground   ·   df = 1   ·  '
                '[!]  Requer > 200 pb',
            ),
            'lrt_M0_vs_Branch-site': (
                'M0 → Branch-site   (seleção episódica alt.)',
                'H₀  M0 — taxa única  ·  '
                'H₁  Branch-site — seleção episódica no foreground   ·   df = 2   ·  '
                '[!]  Requer > 200 pb',
            ),
        }

        comparisons  = {}   # label → col
        descriptions = {}   # col   → description text

        for col in self.df.columns:
            if not col.startswith('lrt_'):
                continue
            if col in KNOWN:
                label, desc = KNOWN[col]
            else:
                label = (col.replace('lrt_', '')
                            .replace('_vs_', ' → ')
                            .replace('_', ' '))
                label = ' '.join(
                    w.upper() if w.lower() in ('m0', 'm1a', 'm2a', 'm7', 'm8') else w
                    for w in label.split()
                )
                desc = ''
            comparisons[label]  = col
            descriptions[col]   = desc

        print(f"[INFO] LRT columns encontradas: {comparisons}")
        return comparisons, descriptions
    
    def _render_lrt_table(self, parent, lrt_col: str, comparison_name: str):
        """Renderiza tabela LRT com estatísticas e omegas recuperados
        Para Branch/BranchSite, exibe múltiplos omegas por tag"""
        # Parse model names reliably from lrt_col (e.g. lrt_M0_vs_Branch),
        # not from comparison_name (which may use → and extra text).
        col_parts = lrt_col.replace('lrt_', '').split('_vs_')
        if len(col_parts) != 2:
            ctk.CTkLabel(parent, text=TEXTS["viewer_lrt_parse_err"]).pack()
            return

        model1_raw = col_parts[0].lower()
        model2_raw = col_parts[1].lower()
        
        # Identificar qual é o modelo alternativo (com mais parâmetros)
        # M0 vs M1a -> M1a é alternativo
        # M1a vs M2a -> M2a é alternativo
        # M7 vs M8 -> M8 é alternativo
        # M0 vs Branch -> Branch é alternativo
        
        model_hierarchy = {
            'm0': 0, 'm1a': 1, 'm2a': 2, 'm7': 1, 'm8': 2, 'branch': 1
        }
        
        if model_hierarchy.get(model2_raw, 2) > model_hierarchy.get(model1_raw, 0):
            alternative_model = model2_raw
        else:
            alternative_model = model1_raw
        
        # Capitalizar corretamente
        # Para Branch-site usar nomes padronizados
        if alternative_model == 'branch':
            alt_display = 'Branch'
        elif alternative_model.lower() == 'branch-site':
            alt_display = 'Branch-site'
        elif alternative_model.lower() == 'branch-site_null':
            alt_display = 'Branch-site_null'
        else:
            alt_display = alternative_model.upper()
        
        is_branch_model = alternative_model == 'branch'
        is_branchsite_model = 'branch-site' in alternative_model.lower()

        # ── 200 bp notice (Branch / Branch-site) ──────────────────────
        if is_branch_model or is_branchsite_model:
            notice = ctk.CTkFrame(parent, fg_color='#1c1408', corner_radius=8,
                                  border_width=1,
                                  border_color=self.COLORS['warning'])
            notice.pack(fill='x', padx=8, pady=(4, 8))
            ctk.CTkLabel(
                notice,
                text=TEXTS["lrt_branch_warning"],
                font=("Roboto", 9),
                text_color=self.COLORS['warning'],
                wraplength=860,
                justify='left',
            ).pack(padx=14, pady=8, anchor='w')

        # ── Table header ──────────────────────────────────────────────
        header_frame = ctk.CTkFrame(parent, fg_color='#1a1a26', corner_radius=6)
        header_frame.pack(fill='x', padx=8, pady=(4, 2))

        if is_branchsite_model:
            headers    = ["Gene", "Class", "Proportion", "Background ω", "Foreground ω", "2Δℓ", "p-valor", "Sig."]
            col_widths = [200, 80, 120, 120, 120, 100, 120, 60]
        elif is_branch_model:
            headers    = ["Gene", "ω (Branch Tags)", "2Δℓ", "p-valor", "Sig."]
            col_widths = [250, 400, 100, 120, 60]
        else:
            headers    = ["Gene", f"ω ({alt_display})", "2Δℓ", "p-valor", "Sig."]
            col_widths = [260, 120, 100, 130, 60]

        for i, (h, width) in enumerate(zip(headers, col_widths)):
            ctk.CTkLabel(header_frame, text=h, font=("Roboto", 11, "bold"),
                         text_color=self.COLORS['accent_blue'], width=width).grid(
                             row=0, column=i, padx=5, pady=9, sticky="w")
        
        # Rows
        row_count = 0
        for idx, (_, row) in enumerate(self.df.iterrows()):
            lrt_val = row[lrt_col]
            
            if pd.isna(lrt_val):
                continue
            
            gene = row['Gene']
            
            # Para Branch/BranchSite, extrair múltiplos omegas por tag
            if is_branch_model:
                try:
                    from src.backend.sites_parser import SitesParser
                    from pathlib import Path
                    
                    results_file = self._find_results_file(gene, alt_display)
                    if results_file:
                        omegas_by_tag = SitesParser.extract_omega_by_tags(results_file)
                        if omegas_by_tag:
                            # Mapear tags para nomes mais informativos
                            # background -> "Background", #1 -> "#1", foreground -> "Foreground"
                            omega_items = []
                            
                            # Ordenar: background primeiro, depois #1, #2, etc., depois foreground
                            def sort_tags(item):
                                tag, val = item
                                if tag == 'background':
                                    return (0, tag)
                                elif tag == 'foreground':
                                    return (2, tag)
                                elif tag.startswith('#'):
                                    try:
                                        return (1, int(tag[1:]))
                                    except Exception:
                                        return (1.5, tag)
                                else:
                                    return (3, tag)
                            
                            for tag, val in sorted(omegas_by_tag.items(), key=sort_tags):
                                # Formatar com indicador de placeholder
                                if val == 999.0:
                                    # 999 é placeholder do CODEML (sem dados para essa tag)
                                    display_val = "N/A"
                                    display_tag = tag.replace('background', 'Background').replace('foreground', 'Foreground')
                                    if tag.startswith('#'):
                                        display_tag = tag
                                    omega_items.append(f"{display_tag}: {display_val}")
                                else:
                                    # Valor real
                                    if tag == 'background':
                                        display_tag = 'Background'
                                    elif tag == 'foreground':
                                        display_tag = 'Foreground'
                                    elif tag.startswith('#'):
                                        display_tag = tag
                                    else:
                                        display_tag = tag.replace('_', ' ').title()
                                    
                                    omega_items.append(f"{display_tag}: {val:.4f}")
                            
                            # Formatar com quebra de linha se houver muitos valores
                            if len(omega_items) <= 2:
                                omega_str = " | ".join(omega_items)
                            else:
                                omega_str = "\n".join(omega_items)
                        else:
                            # Fallback para omega global
                            omega = SitesParser.extract_omega_robust(results_file)
                            omega_str = f"{omega:.4f}" if omega else "N/A"
                    else:
                        omega_str = "N/A"
                except Exception as e:
                    omega_str = "N/A"
                    
                # Para cálculo de p-valor, usar omega global
                omega = None
                try:
                    results_file = self._find_results_file(gene, alt_display)
                    if results_file:
                        omega = SitesParser.extract_omega_robust(results_file)
                except Exception:
                    pass
            else:
                # Buscar omega do modelo alternativo - com fallback para arquivo
                omega_col = f"{alt_display}_omega"
                omega = row.get(omega_col, np.nan)
                
                # Se omega está faltando, tentar extrair do arquivo
                if pd.isna(omega):
                    try:
                        from src.backend.sites_parser import SitesParser
                        from pathlib import Path
                        
                        # Procurar arquivo de resultados com suporte a múltiplas variações
                        results_file = self._find_results_file(gene, alt_display)
                        if results_file:
                            omega = SitesParser.extract_omega_robust(results_file)
                    except Exception:
                        pass

                omega_str = f"{omega:.4f}" if pd.notna(omega) else "N/A"
            
            # Para Branch-site, preparar dados de site classes
            branchsite_class_data = None
            if is_branchsite_model:
                try:
                    branchsite_class_data = {}
                    for class_name in ['0', '1', '2a', '2b']:
                        prop_col = f"Branch-site_class{class_name}_prop"
                        bg_w_col = f"Branch-site_class{class_name}_bg_w"
                        fg_w_col = f"Branch-site_class{class_name}_fg_w"
                        
                        prop = row.get(prop_col, np.nan)
                        bg_w = row.get(bg_w_col, np.nan)
                        fg_w = row.get(fg_w_col, np.nan)
                        
                        if pd.notna(prop) and pd.notna(bg_w) and pd.notna(fg_w):
                            branchsite_class_data[class_name] = {
                                'prop': prop,
                                'bg_w': bg_w,
                                'fg_w': fg_w
                            }
                except Exception:
                    branchsite_class_data = None
            
            # Calcular p-valor com df correto por tipo de comparação:
            #   M2a / M8        : chi2(df=2)
            #   Branch          : chi2(df = n_grupos = abs(np_Branch - np_M0))
            #   Branch-site     : mistura 50:50 chi2(0)+chi2(1) → p = 0.5*chi2.sf(x, 1)
            #   demais (M1a etc): chi2(df=1)  [M0 vs M1a nunca chega aqui; já é tratado]
            if alternative_model in ['m2a', 'm8']:
                df_chi2 = 2
                p_val = 1 - stats.chi2.cdf(lrt_val, df=2) if lrt_val > 0 else 1.0
            elif is_branch_model:
                # df = (np_Branch − np_M0) − (ntime_Branch − ntime_M0)
                # Necessário pois M0 usa árvore não-enraizada e Branch usa enraizada.
                m0_np_val        = row.get('M0_np', np.nan)
                branch_np_val    = row.get('Branch_np', np.nan)
                m0_ntime_val     = row.get('M0_ntime', np.nan)
                branch_ntime_val = row.get('Branch_ntime', np.nan)
                if pd.notna(m0_np_val) and pd.notna(branch_np_val):
                    raw_df = abs(int(branch_np_val) - int(m0_np_val))
                    if pd.notna(m0_ntime_val) and pd.notna(branch_ntime_val):
                        df_chi2 = max(1, raw_df - (int(branch_ntime_val) - int(m0_ntime_val)))
                    else:
                        df_chi2 = max(1, raw_df)
                else:
                    df_chi2 = 1
                p_val = 1 - stats.chi2.cdf(lrt_val, df=df_chi2) if lrt_val > 0 else 1.0
            elif is_branchsite_model:
                # Distribuição nula: mistura 50:50 chi2(0)+chi2(1)
                # P(2Δl > x) = 0.5 * P(chi2(1) > x)  para x > 0
                # Valor crítico α=0.05: 2.706   α=0.01: 5.412
                p_val = (0.5 * stats.chi2.sf(lrt_val, df=1)) if lrt_val > 0 else 1.0
                df_chi2 = 1  # for display only
            else:
                p_val = 1 - stats.chi2.cdf(lrt_val, df=1) if lrt_val > 0 else 1.0
                df_chi2 = 1
            is_sig = p_val < 0.05
            
            p_val_str = self._fmt_pval(p_val)
            
            # ── Row color based on significance + omega ───────────────
            is_strong = is_sig and pd.notna(omega) and omega > 1.0 if not (is_branch_model or is_branchsite_model) else False
            row_idx = row_count  # for alternating

            if is_strong:
                bg_color    = '#0b2016'
                border_color = '#10b981'
            elif is_sig:
                bg_color    = '#11112a'
                border_color = self.COLORS['accent_blue']
            else:
                bg_color    = self.COLORS['bg_card'] if row_idx % 2 == 0 else self.COLORS['bg_sidebar']
                border_color = self.COLORS['border']
            
            # Para Branch-site, renderizar uma linha por site class
            if is_branchsite_model and branchsite_class_data:
                for class_idx, class_name in enumerate(['0', '1', '2a', '2b']):
                    if class_name not in branchsite_class_data:
                        continue
                    
                    class_info = branchsite_class_data[class_name]
                    
                    row_frame = ctk.CTkFrame(parent, 
                                            fg_color=bg_color,
                                            corner_radius=4, border_width=1,
                                            border_color=border_color)
                    row_frame.pack(fill='x', padx=8, pady=4)
                    
                    sig_text = "* Sim" if is_sig else "—"

                    vals = [
                        gene if class_idx == 0 else "",
                        f"Class {class_name}",
                        f"{class_info['prop']:.4f}",
                        f"{class_info['bg_w']:.4f}",
                        f"{class_info['fg_w']:.4f}",
                        f"{lrt_val:.4f}" if class_idx == 0 else "",
                        p_val_str if class_idx == 0 else "",
                        sig_text if class_idx == 0 else ""
                    ]

                    for i, (v, width) in enumerate(zip(vals, col_widths)):
                        if i == 7 and is_sig:
                            color = self.COLORS['success_light']
                        elif i == 7:
                            color = self.COLORS['text_muted']
                        elif i == 6 and is_sig:
                            color = self.COLORS['success_light']
                        elif i == 6:
                            color = self.COLORS['text_secondary']
                        else:
                            color = self.COLORS['text_secondary']
                        
                        label = ctk.CTkLabel(row_frame, text=v, font=("Roboto", 11),
                                   text_color=color, width=width)
                        label.grid(row=0, column=i, padx=5, pady=8, sticky="w")

                row_count += 1
            else:
                # Renderização padrão para outros modelos
                row_frame = ctk.CTkFrame(parent,
                                        fg_color=bg_color,
                                        corner_radius=4, border_width=1,
                                        border_color=border_color)
                row_frame.pack(fill='x', padx=8, pady=4)

                # Destaque especial para omega > 1 E significante (apenas para não-Branch)
                if is_strong:
                    sig_text = "** pos"
                elif is_sig:
                    sig_text = "* sig"
                else:
                    sig_text = "—"

                # Para Branch: mostrar df calculado na célula 2Δℓ para facilitar
                # leitura do threshold (ex: "12.345 (df=5)")
                lrt_display = (f"{lrt_val:.4f} (df={df_chi2})"
                               if is_branch_model else f"{lrt_val:.4f}")
                vals = [gene, omega_str, lrt_display, p_val_str, sig_text]

                for i, (v, width) in enumerate(zip(vals, col_widths)):
                    if i == 4 and is_strong:
                        color = self.COLORS['success_light']
                    elif i == 4 and is_sig:
                        color = self.COLORS['accent_blue']
                    elif i == 4:
                        color = self.COLORS['text_muted']
                    elif i == 3 and is_sig:
                        color = self.COLORS['success_light'] if is_strong else self.COLORS['accent_blue']
                    else:
                        color = self.COLORS['text_secondary']

                    if i == 1 and is_branch_model and '\n' in str(v):
                        label = ctk.CTkLabel(row_frame, text=v, font=("Roboto", 10),
                                             text_color=color, width=width, justify="left")
                    else:
                        label = ctk.CTkLabel(row_frame, text=v, font=("Roboto", 11),
                                             text_color=color, width=width)

                    label.grid(row=0, column=i, padx=5, pady=8, sticky="nw")
            
            row_count += 1
        
        if row_count == 0:
            ctk.CTkLabel(parent, text=TEXTS["lrt_no_data_for_comparison"],
                        font=("Roboto", 11),
                        text_color=self.COLORS['warning']).pack(pady=30)
        else:
            # ── Footer ───────────────────────────────────────────────
            footer = ctk.CTkFrame(parent, fg_color=self.COLORS['bg_sidebar'], corner_radius=6)
            footer.pack(fill='x', padx=8, pady=(10, 8))

            # Recalcular sig_count com a mesma distribuição usada nas linhas
            def _footer_pval(x):
                if not (pd.notna(x) and x > 0):
                    return 1.0
                if is_branchsite_model:
                    return 0.5 * stats.chi2.sf(x, df=1)
                return 1 - stats.chi2.cdf(x, df=df_chi2)

            sig_count = sum(1 for _, row in self.df.iterrows()
                            if pd.notna(row.get(lrt_col)) and
                            _footer_pval(row[lrt_col]) < 0.05)

            df_display_footer = "mixture(0,1)" if is_branchsite_model else str(df_chi2)
            footer_text = TEXTS["lrt_footer_template"].format(
                total=row_count, sig=sig_count, df=df_display_footer
            )
            ctk.CTkLabel(footer, text=footer_text, font=("Roboto", 10),
                         text_color=self.COLORS['text_tertiary']).pack(pady=8, padx=12)
    
    # ═══════════════════════════════════════════════════════════
    # EXPORTAÇÃO
    # ═══════════════════════════════════════════════════════════

    # ── df_chi2 per comparison ────────────────────────────────
    _DF_CHI2 = {
        'lrt_M0_vs_M1a':                           2,
        'lrt_M1a_vs_M2a':                          2,
        'lrt_M7_vs_M8':                            2,
        'lrt_M0_vs_Branch':                        1,   # approximation
        'lrt_Branch-site_null_vs_Branch-site':     1,
        'lrt_M0_vs_Branch-site':                   2,
    }

    def _build_export_df(self, lrt_col: str) -> pd.DataFrame:
        """Build a clean, filtered DataFrame for one LRT comparison.

        • Only genes where the LRT value is not NaN (model was run).
        • Shows Gene, relevant ω columns, 2Δℓ, p-valor, Sig.
        • No internal columns (_np, _time, _stops, _lnL).
        """
        # Filter: only rows that ran this model
        lrt_series = pd.to_numeric(self.df[lrt_col], errors='coerce')
        mask = lrt_series.notna()
        if not mask.any():
            return pd.DataFrame()
        df_f = self.df[mask].copy()
        lrt_vals = lrt_series[mask]

        # df for chi2
        df_chi2 = self._DF_CHI2.get(lrt_col, 1)

        # Para M0 vs Branch: df per-row via ntime correction
        # df = (np_Branch − np_M0) − (ntime_Branch − ntime_M0)
        branch_df_series = None
        if 'vs_Branch' in lrt_col and 'site' not in lrt_col.lower():
            if ('Branch_np' in df_f.columns and 'M0_np' in df_f.columns):
                raw_dfs = (
                    pd.to_numeric(df_f['Branch_np'], errors='coerce') -
                    pd.to_numeric(df_f['M0_np'],    errors='coerce')
                ).abs()
                if 'Branch_ntime' in df_f.columns and 'M0_ntime' in df_f.columns:
                    ntime_diff = (
                        pd.to_numeric(df_f['Branch_ntime'], errors='coerce') -
                        pd.to_numeric(df_f['M0_ntime'],    errors='coerce')
                    )
                    branch_df_series = (raw_dfs - ntime_diff).clip(lower=1)
                else:
                    branch_df_series = raw_dfs.clip(lower=1)

        if branch_df_series is not None:
            p_vals = pd.Series([
                float(1 - stats.chi2.cdf(x, df=int(d))) if pd.notna(x) and x > 0 else 1.0
                for x, d in zip(lrt_vals, branch_df_series)
            ], index=lrt_vals.index)
        else:
            p_vals = lrt_vals.apply(
                lambda x: float(1 - stats.chi2.cdf(x, df=df_chi2)) if pd.notna(x) and x > 0 else 1.0
            )

        # Build output columns
        out = pd.DataFrame({'Gene': df_f['Gene'].values})

        parts = lrt_col.replace('lrt_', '').split('_vs_')
        model_names = parts if len(parts) == 2 else []

        # ── ω columns ────────────────────────────────────────────────────
        # For M0 vs Branch: skip generic Branch_omega (replaced by per-tag below).
        # For all other models: add the generic omega column normally.
        _is_branch_lrt = (lrt_col == 'lrt_M0_vs_Branch')

        for mn in model_names:
            if _is_branch_lrt and mn == 'Branch':
                continue  # per-tag columns added below
            omega_col = f'{mn}_omega'
            if omega_col in df_f.columns:
                out[f'ω ({mn})'] = pd.to_numeric(df_f[omega_col].values, errors='coerce').round(4)

        # ── Per-tag ω columns for Branch model ───────────────────────────
        # PAML's "w (dN/dS) for branches:" lists groups as [bg, #1, #2, ...]
        # matching the user's branch labels in order of first appearance in tree.
        if _is_branch_lrt:
            from src.backend.sites_parser import SitesParser as _SP

            # 1) Try to get per-tag omegas from TSV columns (post-regeneration)
            _tag_cols_in_tsv = sorted(
                [c for c in df_f.columns
                 if re.match(r'^Branch_(background|#\d+)_omega$', c)],
                key=lambda c: (c != 'Branch_background_omega',
                               int(re.search(r'#(\d+)', c).group(1))
                               if re.search(r'#(\d+)', c) else 0)
            )
            if _tag_cols_in_tsv:
                for col in _tag_cols_in_tsv:
                    tag = col.replace('Branch_', '').replace('_omega', '')
                    label = 'ω (bg)' if tag == 'background' else f'ω ({tag})'
                    out[label] = pd.to_numeric(df_f[col].values, errors='coerce').round(4)
            else:
                # 2) Fallback: parse on-the-fly from result files
                _tag_data: dict = {}   # tag -> {gene -> omega}
                for gene in df_f['Gene'].values:
                    rf = self.output_folder / 'Branch' / f"{gene}_Branch_results.txt"
                    if not rf.exists():
                        continue
                    try:
                        for tag, omega in _SP.extract_omega_by_tags(rf).items():
                            _tag_data.setdefault(tag, {})[gene] = omega
                    except Exception:
                        pass

                # Sort: background first, then #1, #2, ... by number
                sorted_tags = sorted(
                    _tag_data.keys(),
                    key=lambda t: (t != 'background',
                                   int(re.search(r'(\d+)', t).group(1))
                                   if re.search(r'(\d+)', t) else 0)
                )
                for tag in sorted_tags:
                    label = 'ω (bg)' if tag == 'background' else f'ω ({tag})'
                    out[label] = [
                        round(float(_tag_data[tag][g]), 4) if g in _tag_data[tag] else float('nan')
                        for g in df_f['Gene'].values
                    ]

        # ── Other per-tag ω columns (non-Branch models, e.g. Branch-site) ─
        tag_re = re.compile(r'^([A-Za-z0-9\-]+)_(.+)_omega$')
        for col in df_f.columns:
            m_col = tag_re.match(col)
            if m_col and col not in [f'{mn}_omega' for mn in model_names]:
                model_part, tag = m_col.group(1), m_col.group(2)
                # Skip Branch per-tag columns — already handled above
                if model_part == 'Branch' and re.match(r'^(background|#\d+)$', tag):
                    continue
                if any(mn.lower() == model_part.lower() for mn in model_names):
                    out[f'ω ({model_part}/{tag})'] = (
                        pd.to_numeric(df_f[col].values, errors='coerce').round(4)
                    )

        out['2Δℓ']         = lrt_vals.values.round(4)
        out['p-value']     = p_vals.values.round(8)
        out['Sig. p<0.05'] = p_vals.apply(lambda p: 'yes' if p < 0.05 else 'no').values
        out['Sig. p<0.01'] = p_vals.apply(lambda p: 'yes' if p < 0.01 else 'no').values

        # ── BEB positive sites (M2a and M8 only) ─────────────────────────
        # Format: "32 R* (8.200 ± 2.238); 91 G** (8.444 ± 1.804)"
        # Separator is ";" to avoid conflicts with CSV field delimiters.
        _SITES_MODELS = {'M1a_vs_M2a': 'M2a', 'M7_vs_M8': 'M8'}
        _lrt_key = lrt_col.replace('lrt_', '')
        if _lrt_key in _SITES_MODELS:
            sites_model = _SITES_MODELS[_lrt_key]
            from src.backend.sites_parser import SitesParser as _SP2
            sites_col = []
            for gene in df_f['Gene'].values:
                rf = self.output_folder / sites_model / f"{gene}_{sites_model}_results.txt"
                if not rf.exists():
                    sites_col.append('')
                    continue
                try:
                    beb_df = _SP2.parse_sites_from_file(rf, method='BEB')
                    if beb_df.empty:
                        sites_col.append('')
                        continue
                    sig = beb_df[beb_df['pr_w_gt_1'] >= 0.95].sort_values('position')
                    if sig.empty:
                        sites_col.append('')
                        continue
                    parts_list = []
                    for _, sr in sig.iterrows():
                        star = sr['significance'] if sr['significance'] else (
                            '**' if sr['pr_w_gt_1'] >= 0.99 else '*')
                        parts_list.append(
                            f"{int(sr['position'])} {sr['amino_acid']}{star} "
                            f"({sr['post_mean']:.3f} ± {sr['post_se']:.3f})"
                        )
                    sites_col.append('; '.join(parts_list))  # ";" avoids CSV conflicts
                except Exception:
                    sites_col.append('')
            out['Positive Sites (BEB)'] = sites_col

        return out.reset_index(drop=True)

    def _export_excel(self):
        """Exporta para Excel — uma aba por modelo LRT, genes sem dados omitidos"""
        filepath = filedialog.asksaveasfilename(
            parent=self,
            defaultextension=".xlsx",
            filetypes=[("Excel", "*.xlsx")]
        )
        if not filepath:
            return
        try:
            import openpyxl
            from openpyxl.styles import PatternFill, Font, Alignment, Border, Side
            from openpyxl.utils import get_column_letter
        except ImportError:
            # fallback: plain pandas export (no formatting)
            lrt_cols = [c for c in self.df.columns if c.startswith('lrt_')]
            if not lrt_cols:
                messagebox.showwarning(TEXTS["msg_warning"], TEXTS["msg_no_lrt"], parent=self)
                return
            with pd.ExcelWriter(filepath, engine='openpyxl') as writer:
                for lrt_col in lrt_cols:
                    sheet_name = (lrt_col.replace('lrt_', '')
                                         .replace('_vs_', ' vs ')
                                         .replace('_', ' ')[:31])
                    df_out = self._build_export_df(lrt_col)
                    if not df_out.empty:
                        df_out.to_excel(writer, sheet_name=sheet_name, index=False)
            messagebox.showinfo(TEXTS["msg_success"], TEXTS["msg_exported_to"].format(path=filepath), parent=self)
            return

        lrt_cols = [c for c in self.df.columns if c.startswith('lrt_')]
        if not lrt_cols:
            messagebox.showwarning(TEXTS["msg_warning"], TEXTS["msg_no_lrt"], parent=self)
            return

        # ── known sheet labels ──────────────────────────────────────────
        SHEET_LABELS = {
            'lrt_M0_vs_M1a':                       'M0 vs M1a',
            'lrt_M1a_vs_M2a':                      'M1a vs M2a',
            'lrt_M7_vs_M8':                        'M7 vs M8',
            'lrt_M0_vs_Branch':                    'M0 vs Branch',
            'lrt_Branch-site_null_vs_Branch-site': 'Branch-site',
            'lrt_M0_vs_Branch-site':               'M0 vs Branch-site',
        }

        # Style helpers — clean light-background professional theme
        HEADER_FILL   = PatternFill('solid', fgColor='1F3864')  # deep navy
        SIG01_FILL    = PatternFill('solid', fgColor='D6F0E8')  # pale teal  (p<0.01)
        SIG05_FILL    = PatternFill('solid', fgColor='EBF5E0')  # pale green (p<0.05 only)
        ALT_FILL      = PatternFill('solid', fgColor='F5F7FB')  # very light blue-gray
        HEADER_FONT   = Font(bold=True, color='FFFFFF', size=11)
        SIG01_FONT    = Font(color='0E5E4A', size=10, bold=True)
        SIG05_FONT    = Font(color='2D6A1F', size=10)
        NORMAL_FONT   = Font(color='1A1A2E', size=10)
        CENTER        = Alignment(horizontal='center', vertical='center', wrap_text=True)
        thin          = Side(style='thin', color='CCCCCC')
        border        = Border(bottom=thin, left=thin, right=thin)

        try:
            with pd.ExcelWriter(filepath, engine='openpyxl') as writer:
                sheets_written = 0
                for lrt_col in lrt_cols:
                    df_out = self._build_export_df(lrt_col)
                    if df_out.empty:
                        continue
                    sheet_name = SHEET_LABELS.get(lrt_col,
                                    lrt_col.replace('lrt_', '').replace('_vs_', ' vs ')
                                           .replace('_', ' ')[:31])
                    df_out.to_excel(writer, sheet_name=sheet_name, index=False)
                    ws = writer.sheets[sheet_name]

                    # ── style header row ──────────────────────────────
                    for cell in ws[1]:
                        cell.fill      = HEADER_FILL
                        cell.font      = HEADER_FONT
                        cell.alignment = CENTER

                    # ── style data rows ───────────────────────────────
                    sig05_col_idx = sig01_col_idx = sites_col_idx = None
                    for i, col in enumerate(df_out.columns, 1):
                        if col == 'Sig. p<0.05':
                            sig05_col_idx = i
                        elif col == 'Sig. p<0.01':
                            sig01_col_idx = i
                        elif col == 'Positive Sites (BEB)':
                            sites_col_idx = i

                    for row_idx, row in enumerate(ws.iter_rows(min_row=2), 1):
                        is_sig01 = (sig01_col_idx and
                                    row[sig01_col_idx - 1].value == 'yes')
                        is_sig05 = (sig05_col_idx and
                                    row[sig05_col_idx - 1].value == 'yes')
                        if is_sig01:
                            fill = SIG01_FILL
                            font = SIG01_FONT
                        elif is_sig05:
                            fill = SIG05_FILL
                            font = SIG05_FONT
                        else:
                            fill = ALT_FILL if row_idx % 2 == 0 else None
                            font = NORMAL_FONT
                        for col_i, cell in enumerate(row, 1):
                            cell.font      = font
                            cell.border    = border
                            # sites column: left-align, wrap
                            if col_i == sites_col_idx:
                                cell.alignment = Alignment(
                                    horizontal='left', vertical='top', wrap_text=True)
                            else:
                                cell.alignment = CENTER
                            if fill:
                                cell.fill = fill

                    # ── auto-fit column widths ────────────────────────
                    for col_cells in ws.columns:
                        header_val = col_cells[0].value or ''
                        # Sites column: fixed wide + row height
                        if header_val == 'Positive Sites (BEB)':
                            ws.column_dimensions[
                                get_column_letter(col_cells[0].column)
                            ].width = 60
                        else:
                            max_len = max(
                                len(str(c.value)) if c.value is not None else 0
                                for c in col_cells
                            )
                            ws.column_dimensions[
                                get_column_letter(col_cells[0].column)
                            ].width = min(max_len + 4, 35)

                    # Set row heights: header taller, data rows auto
                    ws.row_dimensions[1].height = 22
                    for r in range(2, ws.max_row + 1):
                        ws.row_dimensions[r].height = 18

                    sheets_written += 1

                # ── Summary sheet ─────────────────────────────────────
                summary_rows = []
                for lrt_col in lrt_cols:
                    df_out = self._build_export_df(lrt_col)
                    if df_out.empty:
                        continue
                    n_sig05 = (df_out['Sig. p<0.05'] == 'yes').sum()
                    n_sig01 = (df_out['Sig. p<0.01'] == 'yes').sum()
                    sheet_name = SHEET_LABELS.get(lrt_col,
                                    lrt_col.replace('lrt_', '').replace('_vs_', ' vs ')
                                           .replace('_', ' '))
                    summary_rows.append({
                        'Comparison':         sheet_name,
                        'Genes analyzed':     len(df_out),
                        'Significant (p<0.05)': int(n_sig05),
                        'Significant (p<0.01)': int(n_sig01),
                        '% Sig. (p<0.05)':    f"{100*n_sig05/max(len(df_out),1):.1f}%",
                    })
                if summary_rows:
                    pd.DataFrame(summary_rows).to_excel(
                        writer, sheet_name='Summary', index=False)
                    ws_r = writer.sheets['Summary']
                    for cell in ws_r[1]:
                        cell.fill = HEADER_FILL
                        cell.font = HEADER_FONT
                        cell.alignment = CENTER
                    for col_cells in ws_r.columns:
                        max_len = max(
                            len(str(c.value)) if c.value is not None else 0
                            for c in col_cells
                        )
                        ws_r.column_dimensions[
                            get_column_letter(col_cells[0].column)
                        ].width = min(max_len + 4, 30)

            messagebox.showinfo(TEXTS["msg_success"],
                TEXTS["msg_excel_exported"].format(n=sheets_written, path=filepath),
                parent=self)
        except Exception as e:
            messagebox.showerror(TEXTS["msg_error"], TEXTS["msg_excel_err"].format(error=e), parent=self)

    def _export_csv(self):
        """Exporta para CSV — um arquivo por modelo LRT, genes sem dados omitidos"""
        lrt_cols = [c for c in self.df.columns if c.startswith('lrt_')]
        if not lrt_cols:
            messagebox.showwarning(TEXTS["msg_warning"], TEXTS["msg_no_lrt"], parent=self)
            return

        SHEET_LABELS = {
            'lrt_M0_vs_M1a':                       'M0_vs_M1a',
            'lrt_M1a_vs_M2a':                      'M1a_vs_M2a',
            'lrt_M7_vs_M8':                        'M7_vs_M8',
            'lrt_M0_vs_Branch':                    'M0_vs_Branch',
            'lrt_Branch-site_null_vs_Branch-site': 'Branch-site',
            'lrt_M0_vs_Branch-site':               'M0_vs_Branch-site',
        }

        # Ask for base path (files will be named <base>_<model>.csv)
        base_path = filedialog.asksaveasfilename(
            parent=self,
            title=TEXTS["dialog_save_csv"],
            defaultextension=".csv",
            filetypes=[("CSV", "*.csv")]
        )
        if not base_path:
            return

        from pathlib import Path as _P
        base = _P(base_path).with_suffix('')   # strip .csv if added
        files_written = []
        try:
            for lrt_col in lrt_cols:
                df_out = self._build_export_df(lrt_col)
                if df_out.empty:
                    continue
                label = SHEET_LABELS.get(lrt_col,
                    lrt_col.replace('lrt_', '').replace('_vs_', '_vs_').replace(' ', '_'))
                out_path = str(base) + f'_{label}.csv'
                df_out.to_csv(out_path, index=False, sep=',', encoding='utf-8-sig')
                files_written.append(out_path)

            if files_written:
                messagebox.showinfo(TEXTS["msg_success"],
                    TEXTS["msg_csv_exported"].format(n=len(files_written), files="\n".join(files_written)),
                    parent=self)
            else:
                messagebox.showwarning(TEXTS["msg_warning"], TEXTS["msg_no_csv_data"], parent=self)
        except Exception as e:
            messagebox.showerror(TEXTS["msg_error"], TEXTS["msg_csv_err"].format(error=e), parent=self)

    def _export_charts(self):
        """Exporta gráficos"""
        filepath = filedialog.asksaveasfilename(
            parent=self,
            defaultextension=".png",
            filetypes=[("PNG", "*.png"), ("PDF", "*.pdf")]
        )
        if not filepath:
            return
        
        try:
            fig, axes = plt.subplots(2, 2, figsize=(14, 10))
            fig.patch.set_facecolor('#0f0f0f')
            
            # Gráfico 1: Distribuição de ω
            omega_cols = [col for col in self.df.columns if '_omega' in col]
            if omega_cols:
                omega_data = []
                for col in omega_cols:
                    vals = pd.to_numeric(self.df[col], errors='coerce').dropna()
                    omega_data.extend(vals.tolist())
                
                if omega_data:
                    axes[0, 0].hist(omega_data, bins=30, color='#10b981', alpha=0.7, edgecolor='white')
                    axes[0, 0].set_title(TEXTS["chart_omega_dist"], color='white', fontsize=12)
                    axes[0, 0].set_xlabel('ω', color='white')
                    axes[0, 0].set_ylabel(TEXTS["chart_freq"], color='white')
                    axes[0, 0].set_facecolor('#1e1e1e')
                    axes[0, 0].tick_params(colors='white')
            
            # Gráfico 2: LRT values
            lrt_cols = [col for col in self.df.columns if col.startswith('lrt_')]
            if lrt_cols:
                lrt_data = pd.to_numeric(self.df[lrt_cols[0]], errors='coerce').dropna()
                if not lrt_data.empty:
                    axes[0, 1].hist(lrt_data, bins=20, color='#3b82f6', alpha=0.7, edgecolor='white')
                    axes[0, 1].set_title('2Δℓ Distribution', color='white', fontsize=12)
                    axes[0, 1].set_xlabel('2Δℓ', color='white')
                    axes[0, 1].set_ylabel(TEXTS["chart_freq"], color='white')
                    axes[0, 1].set_facecolor('#1e1e1e')
                    axes[0, 1].tick_params(colors='white')
            
            # Gráfico 3: Genes com seleção positiva
            positive_genes = self._detect_positive_selection()
            if positive_genes:
                gene_names = list(positive_genes.keys())[:10]
                gene_counts = [len(positive_genes[g]) for g in gene_names]
                axes[1, 0].barh(gene_names, gene_counts, color='#10b981', alpha=0.8)
                axes[1, 0].set_title('Top 10 Genes com Seleção Positiva', color='white', fontsize=12)
                axes[1, 0].set_xlabel('Nº de Sinais', color='white')
                axes[1, 0].set_facecolor('#1e1e1e')
                axes[1, 0].tick_params(colors='white')
            
            # Texto de resumo
            axes[1, 1].axis('off')
            summary_text = f"""
            RESUMO DA ANALISE

            Total de Genes: {len(self.df)}
            Genes com Selecao Positiva: {len(positive_genes)}
            omega Medio: {self._calc_avg_omega():.3f}
            Modelos Analisados: {self._count_models()}
            """
            axes[1, 1].text(0.1, 0.5, summary_text, color='white', fontsize=11,
                          verticalalignment='center', family='monospace',
                          bbox=dict(boxstyle='round', facecolor='#1e1e1e', alpha=0.8))
            
            plt.tight_layout()
            plt.savefig(filepath, dpi=300, facecolor='#0f0f0f')
            plt.close()
            
            messagebox.showinfo(TEXTS["msg_success"], TEXTS["msg_exported_to"].format(path=filepath), parent=self)
        except Exception as e:
            messagebox.showerror(TEXTS["msg_error"], TEXTS["msg_export_err"].format(error=e), parent=self)
    
    def _export_html(self):
        """Exporta relatório HTML — uma seção por modelo LRT, genes sem dados omitidos"""
        filepath = filedialog.asksaveasfilename(
            parent=self,
            defaultextension=".html",
            filetypes=[("HTML", "*.html")]
        )
        if not filepath:
            return

        SHEET_LABELS = {
            'lrt_M0_vs_M1a':                       ('M0 → M1a',  'Nearly neutral pre-test (M1a vs M0)'),
            'lrt_M1a_vs_M2a':                      ('M1a → M2a', 'Sítios positivos (M2a vs M1a)'),
            'lrt_M7_vs_M8':                        ('M7 → M8',   'Beta + ω > 1 (M8 vs M7)'),
            'lrt_M0_vs_Branch':                    ('M0 → Branch','Ramos livres (Branch vs M0)'),
            'lrt_Branch-site_null_vs_Branch-site': ('Branch-site','Seleção episódica no foreground'),
            'lrt_M0_vs_Branch-site':               ('M0 → Branch-site','Seleção episódica alt.'),
        }
        BRANCH_WARNING = (
            '<div class="warn">[!] Branch e Branch-site requerem sequencias &gt; 200 pb '
            'para estimativas confiaveis de ω.</div>'
        )

        try:
            positive_genes = self._detect_positive_selection()
            lrt_cols = [c for c in self.df.columns if c.startswith('lrt_')]
            now_str  = pd.Timestamp.now().strftime('%d/%m/%Y às %H:%M')

            # ── Build per-model section HTML ──────────────────────────
            sections_html = ''
            for lrt_col in lrt_cols:
                df_out = self._build_export_df(lrt_col)
                if df_out.empty:
                    continue
                label, subtitle = SHEET_LABELS.get(lrt_col, (lrt_col, ''))
                is_branch = 'branch' in lrt_col.lower()

                # table header
                th_cells = ''.join(f'<th>{c}</th>' for c in df_out.columns)
                # table rows
                tr_rows = ''
                for _, row in df_out.iterrows():
                    sig = row.get('Sig. p<0.05', '—') == 'sim'
                    row_cls = ' class="sig"' if sig else ''
                    cells = ''
                    for col_name, val in row.items():
                        cell_cls = ''
                        if 'ω' in col_name and isinstance(val, float) and val > 1.0:
                            cell_cls = ' class="pos"'
                        if pd.isna(val):
                            display = '—'
                        elif isinstance(val, float):
                            display = f'{val:.5g}'
                        else:
                            display = str(val)
                        cells += f'<td{cell_cls}>{display}</td>'
                    tr_rows += f'<tr{row_cls}>{cells}</tr>\n'

                n_sig = (df_out['Sig. p<0.05'] == 'sim').sum() if 'Sig. p<0.05' in df_out.columns else 0
                warn_html = BRANCH_WARNING if is_branch else ''

                sections_html += f"""
        <section>
          <h2>{label}</h2>
          <p class="subtitle">{subtitle}</p>
          {warn_html}
          <p class="meta">{len(df_out)} genes analisados &nbsp;·&nbsp; {n_sig} significantes (p &lt; 0.05)</p>
          <div style="overflow-x:auto">
          <table>
            <thead><tr>{th_cells}</tr></thead>
            <tbody>{tr_rows}</tbody>
          </table>
          </div>
        </section>
"""

            # ── Positive-selection cards ──────────────────────────────
            pos_html = ''
            if positive_genes:
                for gene, signals in positive_genes.items():
                    signals_inner = ''.join(
                        f'<div class="signal">{st}: ω = {sd["omega"]:.4f}, '
                        f'p = {self._fmt_pval(sd["p_value"])}</div>'
                        for st, sd in signals.items()
                    )
                    pos_html += (f'<div class="gene-card">'
                                 f'<div class="gene-name">{gene}</div>'
                                 f'{signals_inner}</div>\n')
            else:
                pos_html = ('<p style="color:#f59e0b;text-align:center;padding:20px">'
                            'Nenhum gene com selecao positiva detectado</p>')

            html_content = f"""<!DOCTYPE html>
<html lang="pt-BR">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>EasyPAML — Relatório de Análise</title>
<style>
*{{margin:0;padding:0;box-sizing:border-box}}
body{{font-family:'Segoe UI',sans-serif;background:#0c0c0f;color:#e0e0e0;padding:36px 20px}}
.container{{max-width:1280px;margin:0 auto;background:#16161a;border-radius:14px;padding:40px;
            box-shadow:0 8px 32px rgba(0,0,0,.5)}}
h1{{color:#6366f1;font-size:28px;margin-bottom:6px;text-align:center}}
h2{{color:#22d3ee;font-size:20px;margin:40px 0 8px;border-bottom:1px solid #222}}
.subtitle{{color:#9898a6;font-size:13px;margin-bottom:6px}}
.meta{{color:#5e5e6e;font-size:12px;margin-bottom:12px}}
.stats-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:16px;margin:24px 0}}
.stat-card{{background:#1e1e24;padding:18px;border-radius:10px;border:1px solid #222;text-align:center}}
.stat-label{{font-size:11px;color:#9898a6;text-transform:uppercase;margin-bottom:6px}}
.stat-value{{font-size:28px;font-weight:700;color:#22c55e}}
table{{width:100%;border-collapse:collapse;margin:8px 0;font-size:12px}}
th{{background:#1e3a5f;color:#fff;padding:10px 12px;text-align:left;font-weight:600}}
td{{padding:9px 12px;border-bottom:1px solid #222;color:#ccccd8}}
tr:hover td{{background:#1e1e24}}
tr.sig td{{background:#0b2016;color:#6ee7b7}}
td.pos{{color:#22d3ee;font-weight:700}}
.gene-card{{background:#2a2a2a;border:1px solid #10b981;border-radius:8px;padding:16px;margin:10px 0}}
.gene-name{{font-size:15px;font-weight:700;color:#10b981;margin-bottom:8px}}
.signal{{background:#1e1e1e;padding:6px 10px;border-radius:5px;margin:4px 0;font-size:13px;color:#a7f3d0}}
.warn{{background:#1c1408;border:1px solid #f59e0b;border-radius:6px;padding:8px 14px;
       margin:8px 0 12px;color:#fbbf24;font-size:12px}}
section{{margin-bottom:48px}}
.footer{{text-align:center;margin-top:48px;padding-top:24px;border-top:1px solid #222;
         color:#5e5e6e;font-size:12px}}
</style>
</head>
<body>
<div class="container">
  <h1>Relatorio de Analise — EasyPAML</h1>
  <p style="text-align:center;color:#5e5e6e;margin-top:6px">Gerado em {now_str}</p>

  <h2 style="margin-top:28px">Estatisticas Gerais</h2>
  <div class="stats-grid">
    <div class="stat-card"><div class="stat-label">Total de Genes</div>
      <div class="stat-value">{len(self.df)}</div></div>
    <div class="stat-card"><div class="stat-label">Modelos Rodados</div>
      <div class="stat-value">{self._count_models()}</div></div>
    <div class="stat-card"><div class="stat-label">Seleção Positiva</div>
      <div class="stat-value">{len(positive_genes)}</div></div>
    <div class="stat-card"><div class="stat-label">ω Médio</div>
      <div class="stat-value">{self._calc_avg_omega():.3f}</div></div>
  </div>

  <h2>Genes com Selecao Positiva</h2>
  {pos_html}

  <h2>Resultados por Modelo</h2>
  {sections_html}

  <div class="footer">
    <p>Relatorio gerado pelo <strong>EasyPAML</strong> — Pipeline CODEML / PAML</p>
  </div>
</div>
</body>
</html>
"""
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(html_content)
            messagebox.showinfo(TEXTS["msg_success"], TEXTS["msg_html_exported"].format(path=filepath), parent=self)
        except Exception as e:
            messagebox.showerror(TEXTS["msg_error"], TEXTS["msg_html_err"].format(error=e), parent=self)