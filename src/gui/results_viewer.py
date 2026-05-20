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
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy import stats
from tkinter import filedialog, messagebox
import re
import sys


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
        self.title("🧬 EasyPAML - Painel de Análise")
        self.geometry("1600x1000")
        self.attributes("-topmost", True)
        self.grab_set()
        
        self.configure(fg_color=self.COLORS['bg_dark'])
        self.output_folder = output_folder
        self.df = None
        self.tag_columns = {}
        
        if not self._load_data():
            self._show_error("Arquivo analysis_summary.tsv não encontrado!")
            return
        
        self._extract_tag_columns()
        self.setup_ui()
    
    def _load_data(self) -> bool:
        """Carrega dados do TSV com tratamento robusto"""
        tsv_file = self.output_folder / "analysis_summary.tsv"
        if not tsv_file.exists():
            return False
        
        try:
            self.df = pd.read_csv(tsv_file, sep='\t')
            self.df = self.df.replace(['NA', 'nan', '', 'None'], np.nan)
            
            numeric_cols = [col for col in self.df.columns if col != 'Gene']
            for col in numeric_cols:
                self.df[col] = pd.to_numeric(self.df[col], errors='coerce')
            
            # Tentar recuperar omegas faltantes dos arquivos de resultados
            self._recover_missing_omegas()
            
            print(f"✅ Dados carregados: {len(self.df)} genes")
            print(f"📊 Colunas: {list(self.df.columns)}")
            return True
        except Exception as e:
            print(f"❌ Erro ao carregar dados: {e}")
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
            
            print(f"\n🔍 Recuperando omegas faltantes para {model_name}...")
            
            for idx in missing_rows:
                gene_name = self.df.loc[idx, 'Gene']
                
                # Procurar arquivo de resultados com suporte a múltiplas variações
                results_file = self._find_results_file(gene_name, model_name)
                
                if results_file:
                    try:
                        omega = SitesParser.extract_omega_robust(results_file)
                        if omega is not None:
                            self.df.loc[idx, omega_col] = omega
                            print(f"  ✓ {gene_name} ({model_name}): ω = {omega:.4f}")
                        else:
                            print(f"  ✗ {gene_name} ({model_name}): não foi possível extrair")
                    except Exception as e:
                        print(f"  ✗ {gene_name} ({model_name}): erro - {e}")
                else:
                    print(f"  ℹ {gene_name} ({model_name}): arquivo não encontrado")
    
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
        
        ctk.CTkLabel(error_frame, text="❌", font=("Roboto", 48)).pack(pady=20)
        ctk.CTkLabel(error_frame, text=message,
                    font=("Roboto", 14, "bold"), 
                    text_color=self.COLORS['danger']).pack(pady=10)
        ctk.CTkLabel(error_frame, text="Execute uma análise para gerar resultados.",
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

        ctk.CTkLabel(left, text="🧬", font=("Roboto", 24)).pack(side="left", padx=(0, 12))
        title_block = ctk.CTkFrame(left, fg_color='transparent')
        title_block.pack(side="left", fill='y', pady=14)
        ctk.CTkLabel(title_block, text="EasyPAML  —  Resultados",
                     font=("Roboto", 16, "bold"),
                     text_color=self.COLORS['text_primary']).pack(anchor="w")
        ctk.CTkLabel(title_block, text="Análise de Seleção Positiva  ·  CODEML / PAML",
                     font=("Roboto", 9),
                     text_color=self.COLORS['text_muted']).pack(anchor="w")

        right = ctk.CTkFrame(header, fg_color='transparent')
        right.pack(side="right", padx=24, pady=16, fill='y')
        ctk.CTkLabel(right, text=f"{len(self.df)} genes carregados",
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
        
        tabs.add("LRT & P-values")
        tabs.add("Global ω > 1")
        tabs.add("Positive Sites")

        # Verificar se há dados de Branch-site para adicionar aba especial
        branchsite_cols = [col for col in self.df.columns if 'Branch-site_class' in col]
        if branchsite_cols:
            tabs.add("Branch-site Classes")

        tabs.add("Branch Analysis")
        tabs.add("Export")

        self._create_lrt_stats_tab(tabs.tab("LRT & P-values"))
        self._create_positive_selection_tab(tabs.tab("Global ω > 1"))
        self._create_sites_tab(tabs.tab("Positive Sites"))

        if branchsite_cols:
            self._create_branchsite_class_tab(tabs.tab("Branch-site Classes"))

        self._create_tree_tab(tabs.tab("Branch Analysis"))
        self._create_export_tab(tabs.tab("Export"))
    
    def _create_stats_panel(self, parent):
        """Painel com estatísticas gerais — cards premium"""
        positive_genes = self._detect_positive_selection()

        stats_data = [
            ("Total de Genes",    str(len(self.df)),            self.COLORS['accent_blue'],  "🧬"),
            ("Modelos Rodados",   self._count_models(),         self.COLORS['accent_cyan'],  "📊"),
            ("Seleção Positiva",  str(len(positive_genes)),     self.COLORS['success'],      "✅"),
            ("ω Médio",           f"{self._calc_avg_omega():.3f}", self.COLORS['warning'],   "📈"),
        ]

        for label, value, color, icon in stats_data:
            card = ctk.CTkFrame(parent, fg_color=self.COLORS['bg_card'],
                                corner_radius=10, border_width=1,
                                border_color=self.COLORS['border'])
            card.pack(side="left", fill="both", expand=True, padx=5)

            inner = ctk.CTkFrame(card, fg_color='transparent')
            inner.pack(fill='both', expand=True, padx=16, pady=14)

            ctk.CTkLabel(inner, text=f"{icon}  {label}",
                         font=("Roboto", 10),
                         text_color=self.COLORS['text_tertiary'],
                         anchor='w').pack(anchor='w')
            ctk.CTkLabel(inner, text=value,
                         font=("Roboto", 28, "bold"),
                         text_color=color,
                         anchor='w').pack(anchor='w', pady=(4, 0))
    
    def _create_lrt_stats_tab(self, parent):
        """Aba de Tabela LRT com p-valores"""
        ctrl_frame = ctk.CTkFrame(parent, fg_color=self.COLORS['bg_card'],
                                 corner_radius=8, height=60)
        ctrl_frame.pack(fill='x', padx=10, pady=10)
        ctrl_frame.pack_propagate(False)
        
        ctk.CTkLabel(ctrl_frame, text="Comparação:", 
                    font=("Roboto", 11, "bold")).pack(side='left', padx=15, pady=10)
        
        comparisons = self._get_available_lrt_columns()
        if not comparisons:
            ctk.CTkLabel(parent, text="⚠️ Sem comparações LRT disponíveis",
                        font=("Roboto", 12),
                        text_color=self.COLORS['warning']).pack(pady=50)
            return
        
        comp_combo = ctk.CTkComboBox(ctrl_frame, values=list(comparisons.keys()), width=300)
        comp_combo.pack(side='left', padx=(0, 15))
        comp_combo.set(list(comparisons.keys())[0])
        
        table_frame = ctk.CTkScrollableFrame(parent, fg_color=self.COLORS['bg_feed'],
                                            corner_radius=8)
        table_frame.pack(fill='both', expand=True, padx=10, pady=10)
        
        def update_lrt_table(*args):
            for widget in table_frame.winfo_children():
                widget.destroy()
            
            selected_comp = comp_combo.get()
            col_name = comparisons[selected_comp]
            self._render_lrt_table(table_frame, col_name, selected_comp)
        
        comp_combo.configure(command=update_lrt_table)
        update_lrt_table()
    
    def _create_positive_selection_tab(self, parent):
        """Genes com ω global > 1 (seleção positiva ao nível do gene inteiro)"""
        # ── Info banner ────────────────────────────────────────────────
        info = ctk.CTkFrame(parent, fg_color='#1a1a26', corner_radius=8)
        info.pack(fill='x', padx=10, pady=(10, 2))

        ctk.CTkLabel(info,
                     text="🌐  Seleção Global — ω > 1 no gene inteiro",
                     font=("Roboto", 11, "bold"),
                     text_color=self.COLORS['accent_blue']).pack(side="left", padx=14, pady=(10, 2))

        ctk.CTkLabel(info,
                     text="Critério: ω médio do modelo M2a ou M8 > 1.0  AND  LRT p < 0.05  ·  "
                          "Diferente de seleção em sítios específicos (aba Positive Sites)",
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
            ctk.CTkLabel(empty, text="Nenhum sinal de seleção positiva detectado",
                         font=("Roboto", 13, "bold"),
                         text_color=self.COLORS['text_tertiary']).pack()
            ctk.CTkLabel(empty, text="Critério: ω > 1.0  AND  p-valor < 0.05",
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

            ctk.CTkLabel(content, text=f"🧬  {gene_name}",
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
            ctk.CTkLabel(badge, text="★ Positivo",
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
        ctk.CTkLabel(line1, text="Modelo:", font=("Roboto", 10, "bold")).pack(side='left', padx=(0, 8))
        
        model_combo = ctk.CTkComboBox(line1, values=['M8', 'M2a', 'Branch-site'], width=150)
        model_combo.pack(side='left', padx=(0, 30))
        model_combo.set('M8')
        
        # 2º: Gene (será atualizado quando modelo mudar)
        ctk.CTkLabel(line1, text="Gene:", font=("Roboto", 10, "bold")).pack(side='left', padx=(0, 8))
        
        gene_combo = ctk.CTkComboBox(line1, values=[], width=200)
        gene_combo.pack(side='left', padx=(0, 30))
        
        # 3º: Análise
        ctk.CTkLabel(line1, text="Análise:", font=("Roboto", 10, "bold")).pack(side='left', padx=(0, 8))
        
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
        
        ctk.CTkLabel(line2, text="Filtrar Pr(w>1) ≥", font=("Roboto", 10, "bold")).pack(side='left', padx=(0, 8))
        
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
        
        btn_update = ctk.CTkButton(line2, text="🔄 Atualizar", width=120,
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
                        text=f"❌ Arquivo não encontrado: {display_name}",
                        font=("Roboto", 11),
                        text_color=self.COLORS['warning']).pack(pady=50)
            return
        
        # Parse com a classe SitesParser (se disponível)
        try:
            from backend.sites_parser import SitesParser
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
                        text=f"❌ Erro ao parsear arquivo:\n{str(e)}",
                        font=("Roboto", 11),
                        text_color=self.COLORS['danger']).pack(pady=50)
            return
        
        if df_filtered.empty:
            ctk.CTkLabel(parent, 
                        text=f"ℹ️ Nenhum sítio acima do limiar Pr(w>1) ≥ {p_threshold}",
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
                     text=f"🧬  {gene_name}",
                     font=("Roboto", 13, "bold"),
                     text_color=self.COLORS['text_primary']).pack(anchor="w")
        ctk.CTkLabel(banner_left,
                     text=f"Modelo: {model_name}  ·  Análise: {method}  ·  {omega_text}",
                     font=("Roboto", 9),
                     text_color=self.COLORS['text_tertiary']).pack(anchor="w")

        badge = ctk.CTkFrame(banner, fg_color=self.COLORS['accent_blue'],
                             corner_radius=8, width=80)
        badge.pack(side="right", padx=14, pady=10)
        badge.pack_propagate(False)
        ctk.CTkLabel(badge, text=f"{len(df_filtered)} sítios",
                     font=("Roboto", 11, "bold"),
                     text_color="white").pack(expand=True)

        # ── Table header ─────────────────────────────────────────────
        th = ctk.CTkFrame(parent, fg_color='#1a1a26', corner_radius=6)
        th.pack(fill='x', padx=8, pady=(0, 2))

        headers = [("Posição", 70), ("AA", 55), ("Pr(ω>1)", 100), ("Sig", 55),
                   ("ω (média)", 105), ("ω − SE", 105), ("ω + SE", 105)]

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
                sig_text    = "★★"
                sig_color   = '#10b981'
            elif is_95:
                bg_color    = '#0d1a10'
                text_color  = '#a7f3d0'
                border_color = '#059669'
                sig_text    = "★"
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
            f"★★ (p≥0.99): {sig99}",
            f"★ (p≥0.95): {sig95}",
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
                from backend.sites_parser import SitesParser
                omega_global = SitesParser.extract_omega_robust(str(filepath))
            except:
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
            print(f"❌ Erro no parser manual: {e}")
        
        return df_sites, omega_global
    
    def _create_branchsite_class_tab(self, parent):
        """Aba de visualizacao estruturada de classes de Branch-site"""
        ctrl_frame = ctk.CTkFrame(parent, fg_color=self.COLORS['bg_card'],
                                 corner_radius=8, height=60)
        ctrl_frame.pack(fill='x', padx=10, pady=10)
        ctrl_frame.pack_propagate(False)
        
        ctk.CTkLabel(ctrl_frame, text="Gene:", font=("Roboto", 11, "bold")).pack(side='left', padx=15, pady=10)
        
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
                ctk.CTkLabel(table_frame, text="Gene not found", 
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
        
        ctk.CTkLabel(header_frame, text=f"Branch-site Model Classes - {gene}",
                    font=("Roboto", 12, "bold"),
                    text_color=self.COLORS['accent_cyan']).pack(pady=8)
        
        # Table header
        table_header_frame = ctk.CTkFrame(parent, fg_color=self.COLORS['bg_card_hover'],
                                         corner_radius=6)
        table_header_frame.pack(fill='x', padx=8, pady=(0, 4))
        
        headers = [("Site Class", 120), ("Proportion", 150), ("Background w", 150), ("Foreground w", 150)]
        
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
            bg_w_str = "See file"
            
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
        
        footer_text = "Foreground w values: high values (> 1.0) indicate positive selection on the foreground branch for that site class"
        ctk.CTkLabel(footer_frame, text=footer_text,
                    font=("Roboto", 9),
                    text_color=self.COLORS['text_tertiary'],
                    wraplength=400).pack(pady=8, padx=8)
    
    def _create_tree_tab(self, parent):
        """Branch Analysis — ω por ramo/tag do Branch model com visualização colorida"""
        branch_data = self.tag_columns.get('Branch', {})
        has_branch_data = bool(branch_data.get('omega'))

        # ── Info banner ────────────────────────────────────────────────
        info = ctk.CTkFrame(parent, fg_color='#1a1a26', corner_radius=8)
        info.pack(fill='x', padx=10, pady=(10, 4))
        ctk.CTkLabel(info,
                     text="🌿  Branch Model — ω por Ramo/Tag",
                     font=("Roboto", 11, "bold"),
                     text_color=self.COLORS['accent_blue']).pack(side="left", padx=14, pady=(10, 2))
        ctk.CTkLabel(info,
                     text="Cada barra representa um ramo etiquetado na árvore  ·  "
                          "ω > 1 indica seleção positiva naquele ramo específico",
                     font=("Roboto", 9),
                     text_color=self.COLORS['text_tertiary']).pack(side="left", padx=(0, 14), pady=(10, 2))

        if not has_branch_data:
            empty = ctk.CTkFrame(parent, fg_color='transparent')
            empty.pack(expand=True)
            ctk.CTkLabel(empty, text="🌿", font=("Roboto", 40),
                         text_color=self.COLORS['text_muted']).pack(pady=(50, 8))
            ctk.CTkLabel(empty, text="Sem dados do Branch Model",
                         font=("Roboto", 14, "bold"),
                         text_color=self.COLORS['text_tertiary']).pack()
            ctk.CTkLabel(empty,
                         text="Execute o modelo Branch com uma árvore etiquetada para visualizar ω por ramo.",
                         font=("Roboto", 10),
                         text_color=self.COLORS['text_muted'],
                         wraplength=420).pack(pady=(6, 0))
            return

        # ── Gene selector ──────────────────────────────────────────────
        ctrl = ctk.CTkFrame(parent, fg_color=self.COLORS['bg_card'],
                            corner_radius=8, height=52)
        ctrl.pack(fill='x', padx=10, pady=(0, 4))
        ctrl.pack_propagate(False)

        ctk.CTkLabel(ctrl, text="Gene:", font=("Roboto", 11, "bold")).pack(side='left', padx=15)
        gene_list = self.df['Gene'].tolist()
        gene_combo = ctk.CTkComboBox(ctrl, values=gene_list, width=280)
        gene_combo.pack(side='left', padx=(0, 15))
        if gene_list:
            gene_combo.set(gene_list[0])

        # Summary info label (updated per gene)
        lrt_col = 'lrt_M0_vs_Branch'
        has_lrt = lrt_col in self.df.columns
        info_lbl = ctk.CTkLabel(ctrl, text="", font=("Roboto", 9),
                                text_color=self.COLORS['text_tertiary'])
        info_lbl.pack(side='left', padx=10)

        # ── Chart area ─────────────────────────────────────────────────
        chart_frame = ctk.CTkFrame(parent, fg_color=self.COLORS['bg_feed'], corner_radius=8)
        chart_frame.pack(fill='both', expand=True, padx=10, pady=(0, 10))

        def _omega_color(w: float) -> str:
            if w > 2.0:   return '#ef4444'   # bright red — strong positive
            if w > 1.5:   return '#f97316'   # orange — moderate positive
            if w > 1.0:   return '#fbbf24'   # amber — slight positive
            if w > 0.7:   return '#6366f1'   # indigo — nearly neutral
            return         '#22d3ee'          # cyan — purifying selection

        def render_branch(gene_name: str):
            for w in chart_frame.winfo_children():
                w.destroy()

            row_data = self.df[self.df['Gene'] == gene_name]
            if row_data.empty:
                ctk.CTkLabel(chart_frame, text="Gene não encontrado",
                             text_color=self.COLORS['text_tertiary']).pack(expand=True)
                return

            row = row_data.iloc[0]
            omega_cols = branch_data.get('omega', {})
            tags = sorted(omega_cols.keys(), key=lambda t: (t != 'background', t))
            omegas = []
            for tag in tags:
                col = omega_cols[tag]
                val = row[col] if col in self.df.columns else np.nan
                omegas.append(float(val) if pd.notna(val) else 0.0)

            if not omegas:
                ctk.CTkLabel(chart_frame, text="Sem dados de ω para este gene",
                             text_color=self.COLORS['text_tertiary']).pack(expand=True)
                return

            # Update LRT info label
            if has_lrt:
                lrt_val = row.get(lrt_col, np.nan)
                if pd.notna(lrt_val):
                    p = 1 - stats.chi2.cdf(lrt_val, df=1) if lrt_val > 0 else 1.0
                    sig = "  ★ p < 0.05" if p < 0.05 else ""
                    info_lbl.configure(
                        text=f"LRT (M0 vs Branch): 2Δℓ = {lrt_val:.3f}  ·  p = {self._fmt_pval(p)}{sig}",
                        text_color=self.COLORS['success'] if p < 0.05 else self.COLORS['text_tertiary']
                    )

            # Build labels
            tag_labels = []
            for t in tags:
                if t == 'background':   tag_labels.append('Background (bg)')
                elif t.startswith('#'): tag_labels.append(f'Ramo {t}')
                else:                   tag_labels.append(t.replace('_', ' ').capitalize())

            colors = [_omega_color(w) for w in omegas]
            n = len(tags)
            fig_h = max(3.2, n * 0.75)

            fig, ax = plt.subplots(figsize=(8.5, fig_h), facecolor='#111115')
            ax.set_facecolor('#16161a')

            y_pos = list(range(n))
            bars = ax.barh(y_pos, omegas, color=colors, alpha=0.88,
                           height=0.54, edgecolor='none')

            # ω = 1 reference line
            x_max = max(max(omegas) * 1.25, 2.5)
            ax.axvline(x=1.0, color='#9898a6', linestyle='--', linewidth=1.5, alpha=0.7)
            ax.text(1.02, n - 0.15, 'neutro', color='#9898a6', fontsize=8)

            # Value labels inside/outside bars
            for bar, omega in zip(bars, omegas):
                offset = x_max * 0.015
                ax.text(omega + offset, bar.get_y() + bar.get_height() / 2,
                        f'{omega:.3f}', va='center', ha='left',
                        color='#ededef', fontsize=9, fontweight='bold')

            ax.set_yticks(y_pos)
            ax.set_yticklabels(tag_labels, color='#ededef', fontsize=10)
            ax.set_xlabel('ω (dN/dS)', color='#9898a6', fontsize=10)
            ax.set_title(f'Perfil de ω por Ramo  —  {gene_name}',
                         color='#ededef', fontsize=12, fontweight='bold', pad=12)
            ax.tick_params(colors='#9898a6', length=3)
            ax.set_xlim(0, x_max)

            # Legend
            legend_els = [
                mpatches.Patch(facecolor='#22d3ee', label='ω < 0.7  Purificação'),
                mpatches.Patch(facecolor='#6366f1', label='0.7–1.0  Neutro'),
                mpatches.Patch(facecolor='#fbbf24', label='ω > 1.0  Seleção positiva'),
                mpatches.Patch(facecolor='#ef4444', label='ω > 2.0  Seleção forte'),
            ]
            ax.legend(handles=legend_els, loc='lower right', framealpha=0.35,
                      facecolor='#1e1e24', edgecolor='#32323e',
                      labelcolor='#9898a6', fontsize=8)

            for spine in ax.spines.values():
                spine.set_color('#222228')
            ax.grid(axis='x', color='#222228', linewidth=0.5, alpha=0.7)

            plt.tight_layout(pad=1.4)

            canvas_w = FigureCanvasTkAgg(fig, master=chart_frame)
            canvas_w.draw()
            canvas_w.get_tk_widget().pack(fill='both', expand=True)
            plt.close(fig)

        gene_combo.configure(command=lambda v: render_branch(v))
        if gene_list:
            render_branch(gene_list[0])
    
    def _create_export_tab(self, parent):
        """Aba de exportação"""
        main_frame = ctk.CTkScrollableFrame(parent, fg_color=self.COLORS['bg_feed'],
                                           corner_radius=10)
        main_frame.pack(fill='both', expand=True, padx=15, pady=15)
        
        header_frame = ctk.CTkFrame(main_frame, fg_color='transparent')
        header_frame.pack(fill='x', pady=(0, 30))
        
        ctk.CTkLabel(header_frame, text="💾 Exportar Resultados",
                    font=("Roboto", 18, "bold"),
                    text_color=self.COLORS['text_primary']).pack()
        
        export_options = [
            ("📊 Excel (.xlsx)", "Tabela com formatação profissional", 
             self._export_excel, self.COLORS['success']),
            ("📄 CSV", "Formato universal compatível",
             self._export_csv, self.COLORS['accent_blue']),
            ("📈 Gráficos (PNG)", "Exportar gráficos em alta resolução",
             self._export_charts, self.COLORS['warning']),
            ("📋 Relatório (HTML)", "Relatório completo interativo",
             self._export_html, self.COLORS['accent_cyan']),
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

            ctk.CTkButton(row, text="Exportar →", width=130, height=34,
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
    
    def _get_available_lrt_columns(self) -> dict:
        """Retorna comparações LRT disponíveis"""
        lrt_cols = {}
        for col in self.df.columns:
            if col.startswith('lrt_'):
                # Formatar label de forma legível
                label = col.replace('lrt_', '').replace('_vs_', ' vs ').replace('_', ' ')
                # Capitalizar nomes de modelos
                label = ' '.join(word.upper() if word.lower() in ['m0', 'm1a', 'm2a', 'm7', 'm8'] 
                                else word.capitalize() for word in label.split())
                lrt_cols[label] = col
        
        print(f"🔍 LRT columns encontradas: {lrt_cols}")
        return lrt_cols
    
    def _render_lrt_table(self, parent, lrt_col: str, comparison_name: str):
        """Renderiza tabela LRT com estatísticas e omegas recuperados
        Para Branch/BranchSite, exibe múltiplos omegas por tag"""
        # Parse do nome da comparação
        # Ex: "M1a vs M2a" ou "M0 Vs Branch"
        parts = comparison_name.lower().split(' vs ')
        if len(parts) != 2:
            ctk.CTkLabel(parent, text="❌ Erro ao parsear comparação").pack()
            return
        
        model1_raw, model2_raw = parts[0].strip(), parts[1].strip()
        
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
                                    except:
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
                except:
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
                    except:
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
                except:
                    branchsite_class_data = None
            
            # Calcular p-valor (df=2 para site models com classes, df=1 para outros)
            df_chi2 = 2 if alternative_model in ['m2a', 'm8', 'branch-site'] else 1
            p_val = 1 - stats.chi2.cdf(lrt_val, df=df_chi2) if lrt_val > 0 else 1.0
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
                    
                    sig_text = "★ Sim" if is_sig else "—"

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

                vals = [gene, omega_str, f"{lrt_val:.4f}", p_val_str, sig_text]

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
            ctk.CTkLabel(parent, text="⚠️ Nenhum dado LRT disponível para esta comparação",
                        font=("Roboto", 11),
                        text_color=self.COLORS['warning']).pack(pady=30)
        else:
            # ── Footer ───────────────────────────────────────────────
            footer = ctk.CTkFrame(parent, fg_color=self.COLORS['bg_sidebar'], corner_radius=6)
            footer.pack(fill='x', padx=8, pady=(10, 8))

            sig_count = sum(1 for _, row in self.df.iterrows()
                            if pd.notna(row.get(lrt_col)) and
                            (1 - stats.chi2.cdf(row[lrt_col], df=df_chi2) < 0.05))

            footer_text = (
                f"Total: {row_count} genes  ·  "
                f"Significantes p < 0.05: {sig_count}  ·  "
                f"df = {df_chi2}  ·  "
                f"★ = significante  ·  ★ ω>1 = seleção positiva confirmada"
            )
            ctk.CTkLabel(footer, text=footer_text, font=("Roboto", 10),
                         text_color=self.COLORS['text_tertiary']).pack(pady=8, padx=12)
    
    # ═══════════════════════════════════════════════════════════
    # EXPORTAÇÃO
    # ═══════════════════════════════════════════════════════════
    
    def _export_excel(self):
        """Exporta para Excel"""
        filepath = filedialog.asksaveasfilename(
            parent=self,
            defaultextension=".xlsx",
            filetypes=[("Excel", "*.xlsx")]
        )
        if filepath:
            try:
                self.df.to_excel(filepath, sheet_name='Resultados', index=False)
                messagebox.showinfo("Sucesso", f"Exportado para:\n{filepath}", parent=self)
            except Exception as e:
                messagebox.showerror("Erro", f"Erro ao exportar: {e}", parent=self)

    def _export_csv(self):
        """Exporta para CSV"""
        filepath = filedialog.asksaveasfilename(
            parent=self,
            defaultextension=".csv",
            filetypes=[("CSV", "*.csv")]
        )
        if filepath:
            try:
                self.df.to_csv(filepath, index=False, sep='\t')
                messagebox.showinfo("Sucesso", f"Exportado para:\n{filepath}", parent=self)
            except Exception as e:
                messagebox.showerror("Erro", f"Erro ao exportar: {e}", parent=self)

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
                    axes[0, 0].set_title('Distribuição de ω', color='white', fontsize=12)
                    axes[0, 0].set_xlabel('ω', color='white')
                    axes[0, 0].set_ylabel('Frequência', color='white')
                    axes[0, 0].set_facecolor('#1e1e1e')
                    axes[0, 0].tick_params(colors='white')
            
            # Gráfico 2: LRT values
            lrt_cols = [col for col in self.df.columns if col.startswith('lrt_')]
            if lrt_cols:
                lrt_data = pd.to_numeric(self.df[lrt_cols[0]], errors='coerce').dropna()
                if not lrt_data.empty:
                    axes[0, 1].hist(lrt_data, bins=20, color='#3b82f6', alpha=0.7, edgecolor='white')
                    axes[0, 1].set_title('Distribuição de 2Δℓ', color='white', fontsize=12)
                    axes[0, 1].set_xlabel('2Δℓ', color='white')
                    axes[0, 1].set_ylabel('Frequência', color='white')
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
            📊 RESUMO DA ANÁLISE
            
            Total de Genes: {len(self.df)}
            Genes com Seleção Positiva: {len(positive_genes)}
            ω Médio: {self._calc_avg_omega():.3f}
            Modelos Analisados: {self._count_models()}
            """
            axes[1, 1].text(0.1, 0.5, summary_text, color='white', fontsize=11,
                          verticalalignment='center', family='monospace',
                          bbox=dict(boxstyle='round', facecolor='#1e1e1e', alpha=0.8))
            
            plt.tight_layout()
            plt.savefig(filepath, dpi=300, facecolor='#0f0f0f')
            plt.close()
            
            messagebox.showinfo("Sucesso", f"Graficos exportados em:\n{filepath}", parent=self)
        except Exception as e:
            messagebox.showerror("Erro", f"Erro ao exportar: {e}", parent=self)
    
    def _export_html(self):
        """Exporta relatório HTML interativo"""
        filepath = filedialog.asksaveasfilename(
            parent=self,
            defaultextension=".html",
            filetypes=[("HTML", "*.html")]
        )
        if not filepath:
            return
        
        try:
            positive_genes = self._detect_positive_selection()
            
            html_content = f"""
<!DOCTYPE html>
<html lang="pt-BR">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EasyPAML - Relatório de Análise</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: linear-gradient(135deg, #0f0f0f 0%, #1a1a1a 100%);
            color: #e0e0e0;
            padding: 40px 20px;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: #1e1e1e;
            border-radius: 16px;
            padding: 40px;
            box-shadow: 0 8px 32px rgba(0,0,0,0.4);
        }}
        h1 {{
            color: #3b82f6;
            font-size: 32px;
            margin-bottom: 10px;
            text-align: center;
        }}
        h2 {{
            color: #06b6d4;
            font-size: 24px;
            margin-top: 40px;
            margin-bottom: 20px;
            border-bottom: 2px solid #2a2a2a;
            padding-bottom: 10px;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 30px 0;
        }}
        .stat-card {{
            background: #2a2a2a;
            padding: 20px;
            border-radius: 12px;
            border: 1px solid #3a3a3a;
            text-align: center;
        }}
        .stat-label {{
            font-size: 12px;
            color: #a8a8a8;
            text-transform: uppercase;
            margin-bottom: 8px;
        }}
        .stat-value {{
            font-size: 32px;
            font-weight: bold;
            color: #10b981;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background: #2a2a2a;
            border-radius: 8px;
            overflow: hidden;
        }}
        th {{
            background: #3b82f6;
            color: white;
            padding: 15px;
            text-align: left;
            font-weight: 600;
        }}
        td {{
            padding: 12px 15px;
            border-bottom: 1px solid #3a3a3a;
        }}
        tr:hover {{
            background: #333;
        }}
        .positive {{
            color: #10b981;
            font-weight: bold;
        }}
        .footer {{
            text-align: center;
            margin-top: 50px;
            padding-top: 30px;
            border-top: 1px solid #2a2a2a;
            color: #a8a8a8;
            font-size: 14px;
        }}
        .gene-card {{
            background: #2a2a2a;
            border: 2px solid #10b981;
            border-radius: 10px;
            padding: 20px;
            margin: 15px 0;
        }}
        .gene-name {{
            font-size: 18px;
            font-weight: bold;
            color: #10b981;
            margin-bottom: 10px;
        }}
        .signal {{
            background: #1e1e1e;
            padding: 8px 12px;
            border-radius: 6px;
            margin: 5px 0;
            font-size: 14px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 Relatório de Análise - EasyPAML</h1>
        <p style="text-align: center; color: #a8a8a8; margin-top: 10px;">
            Gerado em {pd.Timestamp.now().strftime('%d/%m/%Y às %H:%M')}
        </p>
        
        <h2>📊 Estatísticas Gerais</h2>
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-label">Total de Genes</div>
                <div class="stat-value">{len(self.df)}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Modelos Analisados</div>
                <div class="stat-value">{self._count_models()}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Seleção Positiva</div>
                <div class="stat-value">{len(positive_genes)}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">ω Médio</div>
                <div class="stat-value">{self._calc_avg_omega():.3f}</div>
            </div>
        </div>
        
        <h2>✅ Genes com Seleção Positiva</h2>
"""
            
            if positive_genes:
                for gene, signals in positive_genes.items():
                    html_content += f"""
        <div class="gene-card">
            <div class="gene-name">🧬 {gene}</div>
"""
                    for signal_type, data in signals.items():
                        html_content += f"""
            <div class="signal">
                {signal_type}: ω = {data['omega']:.4f}, p-valor = {self._fmt_pval(data['p_value'])}
            </div>
"""
                    html_content += "        </div>\n"
            else:
                html_content += """
        <p style="color: #f59e0b; text-align: center; padding: 30px;">
            ⚠️ Nenhum gene com seleção positiva detectado (ω > 1.0 AND p < 0.05)
        </p>
"""
            
            html_content += f"""
        <h2>📋 Tabela Completa de Resultados</h2>
        <div style="overflow-x: auto;">
        <table>
            <thead>
                <tr>
                    <th>Gene</th>
"""
            
            # Adicionar headers dinamicamente
            for col in self.df.columns:
                if col != 'Gene':
                    html_content += f"                    <th>{col}</th>\n"
            
            html_content += """
                </tr>
            </thead>
            <tbody>
"""
            
            # Adicionar linhas
            for _, row in self.df.iterrows():
                html_content += "                <tr>\n"
                html_content += f"                    <td><strong>{row['Gene']}</strong></td>\n"
                
                for col in self.df.columns:
                    if col != 'Gene':
                        val = row[col]
                        if pd.notna(val):
                            if isinstance(val, float):
                                formatted_val = f"{val:.4f}"
                            else:
                                formatted_val = str(val)
                        else:
                            formatted_val = "N/A"
                        
                        # Destacar ω > 1
                        if '_omega' in col and pd.notna(val) and val > 1.0:
                            html_content += f'                    <td class="positive">{formatted_val}</td>\n'
                        else:
                            html_content += f"                    <td>{formatted_val}</td>\n"
                
                html_content += "                </tr>\n"
            
            html_content += """
            </tbody>
        </table>
        </div>
        
        <div class="footer">
            <p>🧬 Relatório gerado pelo <strong>EasyPAML</strong></p>
            <p>Pipeline de análise evolutiva com CODEML (PAML)</p>
        </div>
    </div>
</body>
</html>
"""
            
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(html_content)
            
            messagebox.showinfo("Sucesso", f"Relatorio HTML exportado em:\n{filepath}", parent=self)
        except Exception as e:
            messagebox.showerror("Erro", f"Erro ao exportar HTML: {e}", parent=self)