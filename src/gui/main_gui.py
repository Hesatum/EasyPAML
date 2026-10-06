import customtkinter as ctk
from tkinter import filedialog
from tkinter import simpledialog, Canvas
from pathlib import Path
import threading
import traceback
import io
import sys
import os
import signal
import platform as _platform

# Ajuste de caminho para importação do backend
sys.path.insert(0, str(Path(__file__).parent.parent))

# ── Compatibilidade de plataforma ─────────────────────────────────────────────
_ON_LINUX  = _platform.system() == "Linux"
_ON_WIN    = _platform.system() == "Windows"
# Fonte sans-serif: Roboto no Windows/Mac, DejaVu Sans no Linux
_FONT_UI   = "Roboto" if not _ON_LINUX else "DejaVu Sans"
# Fonte monoespaçada: Cascadia Code no Windows, DejaVu Sans Mono em outros
_FONT_MONO = "Cascadia Code" if _ON_WIN else "DejaVu Sans Mono"

from backend.codeml_backend import CodemlBatchAnalysis
from backend import messages as backend_messages
from backend.ctl_params import (CODONFREQ_OPTIONS, DEFAULT_CODONFREQ, DEFAULT_CTL_PARAMS,
                                build_ctl_text, codonfreq_label, parse_codonfreq_label)
from backend.preflight import discover_per_gene_trees, group_by_gene, list_alignment_files, run_preflight
from .results_viewer import ResultsViewerWindow
from .gui_texts import TEXTS, set_language, get_language, tr
from .ui_helpers import (PALETTE, PreflightDialog, ask_yes_no, disable_mouse_wheel,
                         fit_to_screen, hover_tint, mix, open_folder, show_about, show_message)

try:
    from Bio import Phylo
except Exception:
    Phylo = None

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt

ctk.set_appearance_mode("System")
ctk.set_default_color_theme("blue")

class StdoutRedirect:
    """print() de bibliotecas/partes antigas vira mensagem 'debug' (escondida
    até o usuário marcar "Mostrar detalhes técnicos")."""
    def __init__(self, append_func):
        self.append = append_func
    def write(self, s):
        for line in str(s).splitlines():
            if line.strip():
                self.append(line, 'debug')
    def flush(self): pass

class ModelConfigWindow(ctk.CTkToplevel):
    """Edita os parâmetros do .ctl de um modelo (em memória). CodonFreq é uma
    lista suspensa com número e nome de cada opção; embaixo, a prévia de
    TODOS os parâmetros que irão para o .ctl."""

    COLORS = {
        'bg_dark':        '#0c0c0e',
        'bg_card':        '#16161a',
        'text_primary':   '#ededef',
        'text_secondary': '#a3a3b8',
        'success':        '#15803d',
        'success_hover':  '#166534',
    }

    _NUMERIC_FIELDS = (('NSsites', 'cfg_field_nssites'), ('model', 'cfg_field_model'),
                       ('fix_omega', 'cfg_field_fix_omega'), ('omega', 'cfg_field_omega'),
                       ('ncatG', 'cfg_field_ncatg'), ('kappa', 'cfg_field_kappa'))

    def __init__(self, parent, model_code: str, initial: dict):
        super().__init__(parent)
        self.title(TEXTS["model_config_header"].format(model_code=model_code))
        fit_to_screen(self, 560, 720, min_w=480, min_h=480)
        self.parent = parent
        self.model_code = model_code
        self.entries = {}
        self.configure(fg_color=self.COLORS['bg_dark'])
        self.transient(parent)
        self.after(50, self._grab)
        self.bind("<Escape>", lambda e: self.destroy())

        ctk.CTkLabel(self, text=TEXTS["model_config_header"].format(model_code=model_code),
                     font=(_FONT_UI, 16, "bold"),
                     text_color=self.COLORS['text_primary']).pack(anchor='w', padx=20, pady=(18, 8))

        form = ctk.CTkScrollableFrame(self, fg_color=self.COLORS['bg_card'], corner_radius=10,
                                      border_width=1, border_color='#2a2a2a')
        form.pack(fill="both", expand=True, padx=20, pady=(0, 12))

        current = dict(initial)
        current.update(parent.custom_model_params.get(model_code, {}))
        defaults = parent._global_ctl_options()

        def _label(key_text):
            name, hint = TEXTS[key_text]
            ctk.CTkLabel(form, text=name, font=(_FONT_UI, 13, "bold"),
                         text_color=self.COLORS['text_primary']).pack(anchor="w", padx=12, pady=(10, 0))
            ctk.CTkLabel(form, text=hint, font=(_FONT_UI, 12), wraplength=460, justify='left',
                         text_color=self.COLORS['text_secondary']).pack(anchor="w", padx=12, pady=(0, 4))

        for key, text_key in self._NUMERIC_FIELDS[:4]:
            _label(text_key)
            ent = ctk.CTkEntry(form, fg_color=self.COLORS['bg_dark'], border_color='#32323e',
                               text_color=self.COLORS['text_primary'])
            ent.pack(fill="x", padx=12, pady=(0, 4))
            ent.insert(0, str(current.get(key, '')))
            ent.bind('<KeyRelease>', lambda e: self._refresh_preview())
            self.entries[key] = ent

        _label('cfg_field_codonfreq')
        cf_values = [codonfreq_label(v) for v, _, _ in CODONFREQ_OPTIONS]
        self.codonfreq_menu = ctk.CTkOptionMenu(
            form, values=cf_values, command=lambda _v: self._refresh_preview(),
            fg_color='#262632', button_color='#3a3a4e', text_color=self.COLORS['text_primary'])
        self.codonfreq_menu.set(codonfreq_label(current.get('CodonFreq', defaults['CodonFreq'])))
        self.codonfreq_menu.pack(fill='x', padx=12, pady=(0, 4))

        for key, text_key in self._NUMERIC_FIELDS[4:]:
            _label(text_key)
            ent = ctk.CTkEntry(form, fg_color=self.COLORS['bg_dark'], border_color='#32323e',
                               text_color=self.COLORS['text_primary'])
            ent.pack(fill="x", padx=12, pady=(0, 4))
            ent.insert(0, str(current.get(key, defaults.get(key, ''))))
            ent.bind('<KeyRelease>', lambda e: self._refresh_preview())
            self.entries[key] = ent

        ctk.CTkLabel(form, text=TEXTS["cfg_preview"], font=(_FONT_UI, 12, "bold"),
                     text_color=self.COLORS['text_secondary']).pack(anchor='w', padx=12, pady=(12, 2))
        self.preview = ctk.CTkTextbox(form, height=260, font=(_FONT_MONO, 11),
                                      fg_color=self.COLORS['bg_dark'],
                                      text_color=self.COLORS['text_secondary'])
        self.preview.pack(fill='x', padx=12, pady=(0, 12))
        self._refresh_preview()

        btn_frame = ctk.CTkFrame(self, fg_color='transparent')
        btn_frame.pack(fill="x", padx=20, pady=(0, 18))
        ctk.CTkButton(btn_frame, text=TEXTS["model_config_btn_cancel"], fg_color='#3d3d4a',
                      hover_color='#4d4d5a', command=self.destroy, text_color='#ffffff',
                      font=(_FONT_UI, 13, "bold"), corner_radius=6).pack(side="left", fill="x",
                                                                          expand=True, padx=(0, 8))
        ctk.CTkButton(btn_frame, text=TEXTS["model_config_btn_save"], fg_color=self.COLORS['success'],
                      hover_color=self.COLORS['success_hover'], text_color='#ffffff',
                      command=self._on_save, font=(_FONT_UI, 13, "bold"),
                      corner_radius=6).pack(side="left", fill="x", expand=True)

    def _grab(self):
        try:
            self.grab_set()
            self.focus_force()
        except Exception:
            pass

    @staticmethod
    def _parse(val: str):
        val = val.strip()
        try:
            return float(val) if any(c in val for c in '.eE') else int(val)
        except ValueError:
            return val

    def _values(self) -> dict:
        out = {k: self._parse(ent.get()) for k, ent in self.entries.items() if ent.get().strip() != ''}
        out['CodonFreq'] = parse_codonfreq_label(self.codonfreq_menu.get())
        return out

    def _refresh_preview(self):
        params = dict(DEFAULT_CTL_PARAMS)
        params.update(self.parent._global_ctl_options())
        params.update(self._values())
        params['seqfile'] = f"GENE_{self.model_code}_seq.fasta"
        params['treefile'] = f"GENE_{self.model_code}_tree.nwk"
        params['outfile'] = f"GENE_{self.model_code}_results.txt"
        self.preview.configure(state='normal')
        self.preview.delete('1.0', 'end')
        self.preview.insert('end', build_ctl_text(params))
        self.preview.configure(state='disabled')

    def _on_save(self):
        self.parent.custom_model_params[self.model_code] = self._values()
        self.parent.model_ctl_labels[self.model_code].configure(
            text=TEXTS["model_status_configured"], text_color="#fbbf24")
        self.parent.append_log(TEXTS["cfg_saved"].format(model=self.model_code) + "\n", 'ok')
        self.destroy()


class TreeLabelWindow(ctk.CTkToplevel):
    """Janela para marcar ramos com cladograma retangular biologicamente correto - Premium Styling"""
    
    # Cores para consistent styling
    BG_DARK    = '#0c0c0e'
    BG_SIDEBAR = '#111115'
    BG_CARD    = '#16161a'
    TEXT_PRIMARY   = '#ededef'
    TEXT_SECONDARY = '#9898a6'
    ACCENT_BLUE  = '#6366f1'
    ACCENT_PINK  = '#f472b6'
    SUCCESS      = '#22c55e'
    SUCCESS_HOVER = '#16a34a'
    DANGER       = '#f87171'
    DANGER_HOVER = '#ef4444'
    
    def __init__(self, parent, tree_path: Path | None, mode: str = 'branchsite'):
        super().__init__(parent)
        self.title(tr("Marcar ramos", "Label branches") + f" - {mode}")
        self.geometry("1400x850")
        self.parent = parent
        self.mode = mode
        
        self.configure(fg_color=self.BG_DARK)
        self.attributes("-topmost", True)
        self.grab_set()
        
        self.tree_path = tree_path
        self.clade_tags = {}
        self.clade_positions = {}
        self.scatter_objects = []
        self.marked_clades = {}
        
        if Phylo is None:
            ctk.CTkLabel(self, text=TEXTS["tree_err_no_biopython"],
                        text_color=self.TEXT_SECONDARY).pack(padx=20, pady=20)
            return

        # Layout premium
        left_frame = ctk.CTkFrame(self, width=280, fg_color=self.BG_SIDEBAR,
                                 border_width=1, border_color=self.BG_CARD)
        left_frame.pack(side='left', fill='y', padx=0, pady=0)
        left_frame.pack_propagate(False)
        
        plot_frame = ctk.CTkFrame(self, fg_color=self.BG_DARK)
        plot_frame.pack(side='left', fill='both', expand=True, padx=10, pady=10)

        # Sidebar premium
        title_label = ctk.CTkLabel(left_frame, text=TEXTS["tree_labeler_sidebar_title"],
                                  font=(_FONT_UI, 14, "bold"),
                                  text_color=self.TEXT_PRIMARY)
        title_label.pack(pady=(15, 10), padx=15)

        if self.mode == 'branchsite':
            instructions = TEXTS["tree_labeler_instructions_branchsite"]
        else:
            instructions = TEXTS["tree_labeler_instructions_branch"]
        
        inst_label = ctk.CTkLabel(left_frame, text=instructions, wraplength=240, justify="left", 
                                 font=(_FONT_UI, 12), text_color=self.TEXT_SECONDARY)
        inst_label.pack(pady=10, padx=15)

        legend_header = ctk.CTkLabel(left_frame, text=TEXTS["tree_labeler_legend_title"], font=(_FONT_UI, 12, "bold"),
                                    text_color=self.ACCENT_BLUE)
        legend_header.pack(pady=(20, 5), padx=15)
        
        self.legend_frame = ctk.CTkFrame(left_frame, fg_color=self.BG_CARD, corner_radius=8)
        self.legend_frame.pack(fill='both', expand=True, padx=15, pady=5)

        btn_frame = ctk.CTkFrame(left_frame, fg_color='transparent')
        btn_frame.pack(side='bottom', fill='x', padx=12, pady=15)
        
        ctk.CTkButton(btn_frame, text=TEXTS["tree_labeler_btn_save"], fg_color=self.SUCCESS, hover_color=self.SUCCESS_HOVER,
                     command=self._on_save, height=40, font=(_FONT_UI, 12, "bold"),
                     text_color=self.TEXT_PRIMARY, corner_radius=6).pack(fill='x', pady=5)
        ctk.CTkButton(btn_frame, text=TEXTS["tree_labeler_btn_cancel"], fg_color='#3d3d3d', hover_color='#4d4d4d',
                     command=self.destroy, height=40, font=(_FONT_UI, 12, "bold"),
                     text_color=self.TEXT_PRIMARY, corner_radius=6).pack(fill='x', pady=5)

        # Gráfico
        if self.tree_path is None:
            ctk.CTkLabel(plot_frame, text=TEXTS["tree_err_no_tree"], font=("Arial", 14)).pack(pady=50)
            return

        try:
            self.tree = Phylo.read(str(self.tree_path), 'newick')
        except Exception as e:
            ctk.CTkLabel(plot_frame, text=TEXTS["tree_err_load"].format(error=e), font=("Arial", 12)).pack(pady=50)
            return

        # Matplotlib
        self.fig = Figure(figsize=(12, 10), dpi=110)
        self.fig.patch.set_facecolor('#0a0a0a')
        self.ax = self.fig.add_subplot(111)
        self.ax.set_facecolor('#0a0a0a')
        
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        for spine in self.ax.spines.values():
            spine.set_visible(False)

        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill='both', expand=True)
        self.canvas.mpl_connect('pick_event', self._on_pick)

        self._compute_rectangular_layout()
        self._draw_tree()
        self._refresh_legend()

    def _compute_rectangular_layout(self):
        """Layout em FORMATO DE CHAVES respeitando branch lengths"""
        def _collect_terminals_ordered(clade, acc):
            if clade.is_terminal():
                acc.append(clade)
            else:
                for c in clade.clades:
                    _collect_terminals_ordered(c, acc)
        
        terminals = []
        _collect_terminals_ordered(self.tree.root, terminals)
        
        Y_SPACING = 3.0
        terminal_y_map = {term: float(idx) * Y_SPACING for idx, term in enumerate(terminals)}
        
        depths = {}
        
        def calc_depth_with_lengths(clade, accumulated_depth=0.0):
            depths[clade] = accumulated_depth

            for child in clade.clades:
                # Cladograma: passo uniforme independente dos branch lengths reais.
                # Branch lengths reais distorcem muito a visualização quando as
                # sequências têm taxas de evolução heterogêneas — para o propósito
                # de etiquetagem, topologia importa; comprimento, não.
                calc_depth_with_lengths(child, accumulated_depth + 1.0)
        
        calc_depth_with_lengths(self.tree.root, 0.0)

        # Profundidade máxima entre todos os terminais — todos serão alinhados neste X.
        # Isso cria o visual clássico de cladograma retangular onde todos os tips
        # ficam na mesma coluna da direita, independente de sua profundidade real.
        max_depth = max(
            (depths[t] for t in terminals if t in depths),
            default=1.0
        )
        # Guardado como atributo para que _draw_tree possa usar o offset de label correto
        self._cladogram_depth = max_depth

        for clade in self.tree.find_clades(order='postorder'):
            if clade.is_terminal():
                # Alinhar todos os tips na coluna máxima
                x = max_depth
                y = terminal_y_map[clade]
            else:
                # Nós internos: profundidade real (posição topológica na árvore)
                x = depths.get(clade, 0.0)
                child_ys = [self.clade_positions[child][1] for child in clade.clades
                            if child in self.clade_positions]
                y = sum(child_ys) / len(child_ys) if child_ys else 0.0

            self.clade_positions[clade] = (x, y)

    def _draw_tree(self):
        """Desenha cladograma colorindo APENAS do nó marcado até os tips
        
        Lógica: A cor flui como uma tubulação, do nó marcado até seus tips.
        Se um nó descendente também foi marcado, aquela cor SOBREPÕE a anterior.
        """
        self.ax.clear()
        self.ax.set_facecolor('#0a0a0a')
        
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        for spine in self.ax.spines.values():
            spine.set_visible(False)
        
        self.scatter_objects.clear()

        def get_tag_for_branch(clade):
            """
            Retorna a tag que deve colorir este clado.
            
            Lógica (como tubulação, flui para BAIXO):
            1. Se o próprio clado foi marcado → usa sua tag
            2. Senão, procura o nó marcado mais PRÓXIMO na cadeia ancestral direta
            3. Senão, retorna None (cinza)
            
            "Mais próximo" = primeiro nó marcado quando sobe na árvore
            """
            # Primeiro: este clado foi marcado?
            if clade in self.marked_clades:
                return self.marked_clades[clade]
            
            # Segundo: percorrer a cadeia de ancestrais e achar o primeiro marcado
            # Biopython não tem parent direto, então vamos comparar com todos os marcados
            try:
                # Para cada nó marcado, verificar se é ancestral
                closest_tag = None
                min_distance = float('inf')
                
                for marked_clade, tag in self.marked_clades.items():
                    # Verificar se marked_clade é ancestral de clade
                    try:
                        all_descendants = list(marked_clade.find_clades())
                        if clade in all_descendants:
                            # Calcular distância: número de nós entre marked_clade e clade
                            distance = len(all_descendants) - len(list(clade.find_clades()))
                            if distance < min_distance:
                                min_distance = distance
                                closest_tag = tag
                    except Exception:
                        pass

                return closest_tag
            except Exception:
                pass
            
            return None

        def is_descendant_of_marked(clade):
            """Verifica se este clado é descendente de um nó marcado"""
            tag = get_tag_for_branch(clade)
            return (tag is not None, tag)

        # === DESENHAR LINHAS ===
        for clade in self.tree.find_clades():
            if clade.is_terminal():
                continue
            
            if clade not in self.clade_positions:
                continue
            
            x_parent, y_parent = self.clade_positions[clade]
            
            children_positions = []
            for child in clade.clades:
                if child in self.clade_positions:
                    children_positions.append((child, self.clade_positions[child]))
            
            if not children_positions:
                continue
            
            child_ys = [pos[1] for _, pos in children_positions]
            y_min = min(child_ys)
            y_max = max(child_ys)
            
            # LINHA VERTICAL: usar a tag do nó marcado mais próximo na subárvore
            parent_tag = get_tag_for_branch(clade)
            
            if parent_tag:
                vertical_color = self._get_tag_color(parent_tag)
                vertical_width = 2.8
                vertical_alpha = 1.0
            else:
                vertical_color = '#555555'
                vertical_width = 1.2
                vertical_alpha = 0.7
            
            if abs(y_max - y_min) > 0.01:
                self.ax.plot([x_parent, x_parent], [y_min, y_max],
                           color=vertical_color, linewidth=vertical_width, 
                           zorder=1, alpha=vertical_alpha, solid_capstyle='round')
            
            for child, (x_child, y_child) in children_positions:
                # LINHA HORIZONTAL: cada filho herda a tag do seu próprio ramo
                branch_tag = get_tag_for_branch(child)
                
                if branch_tag:
                    branch_color = self._get_tag_color(branch_tag)
                    branch_width = 2.8
                    branch_alpha = 1.0
                else:
                    branch_color = '#555555'
                    branch_width = 1.2
                    branch_alpha = 0.7

                self.ax.plot([x_parent, x_child], [y_child, y_child],
                           color=branch_color, linewidth=branch_width, 
                           zorder=2, alpha=branch_alpha, solid_capstyle='round')

        # === MARCADORES ===
        for clade in self.tree.find_clades():
            if clade not in self.clade_positions:
                continue
            
            x, y = self.clade_positions[clade]
            
            is_marked, tag = is_descendant_of_marked(clade)
            color = self._get_tag_color(tag) if is_marked else '#777777'
            
            if clade.is_terminal():
                size = 90 if is_marked else 45
                edge_width = 2.2 if is_marked else 1.0
            else:
                size = 55 if is_marked else 28
                edge_width = 1.8 if is_marked else 0.8
            
            edge_color = '#ffffff' if is_marked else '#999999'
            
            scatter = self.ax.scatter(
                [x], [y],
                s=size,
                c=[color],
                edgecolors=edge_color,
                linewidths=edge_width,
                picker=12,
                zorder=100,
                alpha=0.95
            )
            
            scatter.clade_obj = clade
            self.scatter_objects.append((scatter, clade))

        # === LABELS ===
        for term in self.tree.get_terminals():
            if term not in self.clade_positions:
                continue
            
            x, y = self.clade_positions[term]
            name = term.name or "terminal"
            
            is_marked, tag = is_descendant_of_marked(term)
            if is_marked:
                text_color = self._get_tag_color(tag)
                weight = 'bold'
                fontsize = 11
            else:
                text_color = '#e5e5e5'
                weight = 'normal'
                fontsize = 10
            
            # Offset proporcional à profundidade do cladograma para não sobrepor o círculo
            label_offset = getattr(self, '_cladogram_depth', 10) * 0.03
            self.ax.text(x + label_offset, y, name,
                        va='center', ha='left',
                        fontsize=fontsize,
                        color=text_color,
                        weight=weight,
                        zorder=20,
                        family='monospace')

        # === LIMITES ===
        if self.clade_positions:
            all_x = [x for x, y in self.clade_positions.values()]
            all_y = [y for x, y in self.clade_positions.values()]
            
            min_x, max_x = min(all_x), max(all_x)
            min_y, max_y = min(all_y), max(all_y)
            
            max_label_len = 0
            try:
                max_label_len = max(len(t.name or '') for t in self.tree.get_terminals())
            except Exception:
                max_label_len = 20
            
            x_range = max_x - min_x if max_x > min_x else 0.01
            right_margin = max(x_range * 0.35, 0.003) + (max_label_len * 0.0012)
            
            self.ax.set_xlim(min_x - x_range * 0.03, max_x + right_margin)
            self.ax.set_ylim(min_y - 2.5, max_y + 2.5)

        self.canvas.draw_idle()

    def _on_pick(self, event):
        """Callback quando nódulo é clicado"""
        artist = event.artist
        clade = getattr(artist, 'clade_obj', None)
        if clade is None:
            return

        if self.mode == 'branchsite':
            current_tag = self.marked_clades.get(clade)
            
            if current_tag == '#1':
                self._remove_tag_recursively(clade)
                self.parent.append_log(tr("Marca #1 removida de ", "Tag #1 removed from ") + f"{self._get_clade_name(clade)}\n")
            else:
                self._remove_tag_recursively(clade)
                self._apply_tag_recursively(clade, '#1')
                self.marked_clades[clade] = '#1'
                self.parent.append_log(tr("Marca #1 aplicada a ", "Tag #1 applied to ") + f"{self._get_clade_name(clade)}\n")
        
        else:
            current_tag = self.marked_clades.get(clade)
            clade_name = self._get_clade_name(clade)
            
            if current_tag:
                response = simpledialog.askstring(
                    TEXTS["tag_dialog_edit_title"],
                    TEXTS["tag_dialog_edit_prompt"].format(tag=current_tag),
                    parent=self
                )
                
                if response:
                    response = response.strip().lower()
                    if response == 'remover':
                        self._remove_tag_recursively(clade)
                        self.parent.append_log(tr(f"Marca {current_tag} removida de {clade_name}", f"Tag {current_tag} removed from {clade_name}") + "\n")
                    elif response.isdigit():
                        self._remove_tag_recursively(clade)
                        new_tag = f"#{response}"
                        self._apply_tag_recursively(clade, new_tag)
                        self.marked_clades[clade] = new_tag
                        self.parent.append_log(tr(f"Marca alterada para {new_tag} em {clade_name}", f"Tag changed to {new_tag} on {clade_name}") + "\n")
            else:
                response = simpledialog.askstring(
                    TEXTS["tag_dialog_new_title"],
                    TEXTS["tag_dialog_new_prompt"],
                    parent=self
                )
                
                if response and response.strip().isdigit():
                    new_tag = f"#{response.strip()}"
                    self._remove_tag_recursively(clade)
                    self._apply_tag_recursively(clade, new_tag)
                    self.marked_clades[clade] = new_tag
                    self.parent.append_log(tr(f"Marca {new_tag} aplicada a {clade_name}", f"Tag {new_tag} applied to {clade_name}") + "\n")

        self._draw_tree()
        self._refresh_legend()

    def _apply_tag_recursively(self, clade, tag: str):
        """Aplica tag aos terminais descendentes"""
        clades_to_remove = []
        for marked_clade in list(self.marked_clades.keys()):
            if marked_clade in clade.find_clades():
                clades_to_remove.append(marked_clade)
        
        for old_marked in clades_to_remove:
            self.marked_clades.pop(old_marked, None)
        
        for terminal in clade.get_terminals():
            self.clade_tags[terminal] = tag

    def _remove_tag_recursively(self, clade):
        """Remove tag dos terminais descendentes"""
        self.marked_clades.pop(clade, None)
        
        for terminal in clade.get_terminals():
            self.clade_tags.pop(terminal, None)
        
        for desc in list(self.marked_clades.keys()):
            if desc in clade.find_clades() and desc != clade:
                self.marked_clades.pop(desc, None)

    def _get_clade_name(self, clade) -> str:
        """Nome legível do clado"""
        if clade.is_terminal():
            return clade.name or "terminal"
        else:
            terminals = clade.get_terminals()
            if len(terminals) <= 3:
                names = [t.name or "?" for t in terminals[:3]]
                return f"clado({', '.join(names)})"
            else:
                return f"clado({len(terminals)} terminais)"

    def _get_tag_color(self, tag: str) -> str:
        """Cor baseada no número da tag"""
        if not tag or not tag.startswith('#'):
            return '#888888'
        
        try:
            num = int(tag.replace('#', ''))
            cmap = plt.get_cmap('tab20')
            rgba = cmap(num % 20)
            return matplotlib.colors.to_hex(rgba)
        except Exception:
            return '#ff0000'

    def _delete_tag(self, tag: str):
        """Remove todas as ocorrências de uma tag da árvore e redesenha."""
        for clade in [c for c, t in list(self.clade_tags.items()) if t == tag]:
            self.clade_tags.pop(clade, None)
        for clade in [c for c, t in list(self.marked_clades.items()) if t == tag]:
            self.marked_clades.pop(clade, None)
        self._draw_tree()
        self._refresh_legend()
        self.parent.append_log(tr(f"Marca {tag} removida.", f"Tag {tag} removed.") + "\n")

    def _refresh_legend(self):
        """Atualiza legenda com tags ativas e botão de exclusão por tag."""
        for widget in self.legend_frame.winfo_children():
            widget.destroy()

        active_tags = sorted(set(self.clade_tags.values()))

        if not active_tags:
            ctk.CTkLabel(self.legend_frame, text=TEXTS["tree_labeler_no_tags"],
                         text_color="#888888", font=(_FONT_UI, 13, "italic")
                         ).pack(anchor='w', padx=15, pady=10)
        else:
            for tag in active_tags:
                color = self._get_tag_color(tag)

                row_frame = ctk.CTkFrame(self.legend_frame, fg_color="transparent")
                row_frame.pack(fill='x', padx=8, pady=3)

                color_box = Canvas(row_frame, width=22, height=16, highlightthickness=0)
                try:
                    color_box.configure(bg=self.legend_frame.cget('fg_color')[1])
                except Exception:
                    color_box.configure(bg='#16161a')
                color_box.create_rectangle(2, 2, 20, 14, fill=color, outline='white', width=1)
                color_box.pack(side='left', padx=(4, 4))

                tag_label = ctk.CTkLabel(row_frame, text=tag,
                                         font=(_FONT_UI, 12, "bold"),
                                         text_color=color)
                tag_label.pack(side='left', padx=(0, 4))

                count = sum(1 for t in self.clade_tags.values() if t == tag)
                ctk.CTkLabel(row_frame, text=f"({count})",
                             font=(_FONT_UI, 12),
                             text_color="#888888").pack(side='left')

                del_btn = ctk.CTkButton(
                    row_frame, text="X", width=22, height=22,
                    fg_color='transparent',
                    hover_color=self.DANGER,
                    text_color='#888888',
                    border_width=0,
                    corner_radius=4,
                    font=(_FONT_UI, 13, "bold"),
                    command=lambda t=tag: self._delete_tag(t)
                )
                del_btn.pack(side='right', padx=(0, 4))

    def _on_save(self):
        """Salva árvore etiquetada em formato Newick"""
        import re

        for c in self.tree.find_clades():
            if getattr(c, 'name', None):
                cleaned = re.sub(r"\s*#\d+\b", "", str(c.name)).strip()
                c.name = cleaned if cleaned != '' else None

        applied_terminals = 0

        for clade, tag in list(self.clade_tags.items()):
            if not tag:
                continue

            try:
                current_name = clade.name or ''
                base = re.sub(r"\s*#\d+\b", "", str(current_name)).strip()

                if clade.is_terminal():
                    new_name = f"{base}{tag}" if base else f"{tag}"
                    if clade.name != new_name:
                        clade.name = new_name
                        applied_terminals += 1
            except Exception:
                continue

        try:
            sio = io.StringIO()
            Phylo.write(self.tree, sio, 'newick')
            newick_str = sio.getvalue().strip()
            newick_str = re.sub(r"\s+", " ", newick_str)

            if self.mode == 'branchsite':
                self.parent.tree_branchsite_labeled = newick_str
                self.parent.append_log(tr(f"[OK] Árvore branch-site salva ({applied_terminals} terminais).", f"[OK] Branch-site tree saved ({applied_terminals} tips).") + "\n")
            else:
                self.parent.tree_branch_labeled = newick_str
                self.parent.append_log(tr(f"[OK] Árvore do modelo Branch salva ({applied_terminals} terminais).", f"[OK] Branch tree saved ({applied_terminals} tips).") + "\n")

                branchsite_version = re.sub(r"\s*#(?!1)\d+\b", "", newick_str)
                branchsite_version = re.sub(r"\s+", " ", branchsite_version).strip()
                self.parent.tree_branchsite_labeled = branchsite_version

        except Exception as e:
            self.parent.append_log(tr("Erro ao gerar o Newick: ", "Error generating Newick: ") + f"{e}\n", "error")

        self.destroy()


class App(ctk.CTk):
    # ═══════════════════════════════════════════════════════════════════════
    # PALETA DE CORES PREMIUM - Estilo YouTube/Instagram Dark
    # ═══════════════════════════════════════════════════════════════════════
    COLORS = {
        # Backgrounds — near-black, layered dark
        'bg_darkest':     '#07070a',
        'bg_dark':        '#0d0d11',
        'bg_sidebar':     '#111118',
        'bg_card':        '#17171f',
        'bg_card_hover':  '#20202c',
        'bg_feed':        '#1c1c26',   # elevated surface — lighter than bg_card
        'bg_input':       '#1e1e2a',

        # Text hierarchy — higher contrast than before
        'text_primary':   '#eeeef2',
        'text_secondary': '#a0a0b4',
        'text_tertiary':  '#8e8ea4',   # contraste >= 4,5:1 nos cartões
        'text_muted':     '#8a8aa0',

        # Primary accent — indigo
        'accent_blue':        '#6366f1',
        'accent_blue_hover':  '#4f46e5',
        'accent_blue_light':  '#818cf8',

        # Secondary accents
        'accent_cyan':        '#22d3ee',
        'accent_cyan_hover':  '#06b6d4',
        'accent_purple':      '#a78bfa',
        'accent_purple_hover':'#7c3aed',
        'accent_pink':        '#f472b6',
        'accent_pink_hover':  '#db2777',

        # Status
        'success':        '#22c55e',
        'success_hover':  '#16a34a',
        'success_light':  '#86efac',
        'warning':        '#f59e0b',
        'warning_hover':  '#d97706',
        'danger':         '#f87171',
        'danger_hover':   '#ef4444',
        'info':           '#22d3ee',
        'info_hover':     '#06b6d4',

        # Borders — slightly more visible
        'border':         '#262632',
        'border_hover':   '#3a3a4e',
    }
    
    def __init__(self):
        super().__init__()
        from backend.version import __version__
        self.title(f"EasyPAML {__version__}")
        fit_to_screen(self, 1400, 850)
        
        # Aplicar tema escuro profissional
        ctk.set_appearance_mode("dark")
        ctk.set_default_color_theme("blue")
        
        # Configure janela principal
        self.configure(fg_color=self.COLORS['bg_dark'])

        # Inicializar backend
        self.codeml_backend = CodemlBatchAnalysis()

        self.custom_model_params = {}
        self.last_run_summary = None
        self.show_details_var = ctk.BooleanVar(value=False)
        self._detail_lines = []
        self.input_folder = None
        self.per_gene_trees = {}
        self.tree_file = None
        self.output_folder = None
        self.analysis_thread = None
        self.analysis_instance = None
        self.pause_event = threading.Event()
        self.stop_event = threading.Event()
        
        # Opção para incluir modelos neutros automaticamente
        self.include_neutral_models = ctk.BooleanVar(value=True)
        self.include_m8a_var = ctk.BooleanVar(value=True)

        # Auto-detect CPU cores; user can adjust via slider
        _max_cores = CodemlBatchAnalysis.available_cores()
        self.cores_var = ctk.IntVar(value=_max_cores)
        self.ignore_stop_codons_var = ctk.BooleanVar(value=False)
        self.auto_prune_tree_var = ctk.BooleanVar(value=True)

        self.tree_branch_labeled = None
        self.tree_branchsite_labeled = None

        # ═══ SIDEBAR ═══
        self.sidebar = ctk.CTkFrame(self, width=290, corner_radius=0,
                                    fg_color=self.COLORS['bg_sidebar'],
                                    border_width=1, border_color=self.COLORS['bg_card_hover'])
        self.sidebar.pack(side="left", fill="y", padx=0, pady=0)
        self.sidebar.pack_propagate(False)

        # ── Logo ─────────────────────────────────────────────────────
        logo_frame = ctk.CTkFrame(self.sidebar, fg_color='transparent')
        logo_frame.pack(fill="x", padx=16, pady=(20, 4))

        badge_row = ctk.CTkFrame(logo_frame, fg_color='transparent')
        badge_row.pack(anchor='w')

        badge = ctk.CTkFrame(badge_row, fg_color=self.COLORS['accent_blue'],
                             width=38, height=38, corner_radius=10)
        badge.pack(side='left', padx=(0, 11))
        badge.pack_propagate(False)
        ctk.CTkLabel(badge, text="EP", font=(_FONT_UI, 13, "bold"),
                     text_color='#ffffff').pack(expand=True)

        title_col = ctk.CTkFrame(badge_row, fg_color='transparent')
        title_col.pack(side='left', anchor='center')
        ctk.CTkLabel(title_col, text=TEXTS["app_sidebar_title"],
                     font=(_FONT_UI, 17, "bold"),
                     text_color=self.COLORS['text_primary']).pack(anchor='w')
        ctk.CTkLabel(title_col, text=TEXTS["app_sidebar_subtitle"],
                     font=(_FONT_UI, 13),
                     text_color=self.COLORS['text_tertiary']).pack(anchor='w')

        ctk.CTkFrame(logo_frame, fg_color=self.COLORS['border'],
                     height=1, corner_radius=0).pack(fill='x', pady=(14, 0))

        # ── Scrollable section container ─────────────────────────────
        _sb = ctk.CTkScrollableFrame(self.sidebar, fg_color='transparent',
                                     scrollbar_button_color='#52526a',
                                     scrollbar_button_hover_color='#6366f1')
        _sb.pack(fill='both', expand=True, padx=0, pady=(4, 0))

        # ── Rodapé de idioma (fixo, fora do scroll) ──────────────────────
        self._build_lang_footer()

        def _sec(parent, title, icon=""):
            """Labeled card section with inner content frame."""
            card = ctk.CTkFrame(parent, fg_color=self.COLORS['bg_card'],
                                corner_radius=10, border_width=1,
                                border_color=self.COLORS['border'])
            card.pack(fill='x', padx=10, pady=(0, 8))
            hdr = ctk.CTkFrame(card, fg_color='transparent')
            hdr.pack(fill='x', padx=14, pady=(10, 0))
            # Small accent dot before section title
            ctk.CTkFrame(hdr, fg_color=self.COLORS['accent_blue'],
                         width=3, height=13, corner_radius=2).pack(side='left', padx=(0, 7), anchor='center')
            ctk.CTkLabel(hdr,
                         text=title,
                         font=(_FONT_UI, 12, "bold"),
                         text_color=self.COLORS['text_secondary']).pack(side='left', anchor='w')
            ctk.CTkFrame(card, fg_color=self.COLORS['border'],
                         height=1, corner_radius=0).pack(fill='x', padx=12, pady=(6, 0))
            inner = ctk.CTkFrame(card, fg_color='transparent')
            inner.pack(fill='x', padx=10, pady=(10, 12))
            return inner

        def _obtn(parent, text, cmd, color, **kw):
            """Outline-style button — uniform border, accent on hover."""
            return ctk.CTkButton(
                parent, text=text, command=cmd,
                fg_color=self.COLORS['bg_card_hover'],
                hover_color=color,
                text_color=self.COLORS['text_primary'],
                border_width=1, border_color=self.COLORS['border_hover'],
                corner_radius=8, **kw)

        # ── Arquivos ─────────────────────────────────────────────────
        fi = _sec(_sb, TEXTS["section_files"])

        self.btn_input = _obtn(fi, TEXTS["btn_input_folder"],
                               self.select_input_folder,
                               self.COLORS['accent_blue'],
                               font=(_FONT_UI, 13, "bold"), height=36)
        self.btn_input.pack(fill='x', pady=(0, 2))
        self.label_input = ctk.CTkLabel(fi, text=TEXTS["label_not_selected"],
                                        font=(_FONT_UI, 13),
                                        wraplength=230,
                                        text_color=self.COLORS['text_tertiary'])
        self.label_input.pack(anchor='w', padx=4, pady=(0, 8))

        self.btn_tree = _obtn(fi, TEXTS["btn_tree_file"],
                              self.select_tree_file,
                              self.COLORS['accent_blue'],
                              font=(_FONT_UI, 13, "bold"), height=36)
        self.btn_tree.pack(fill='x', pady=(0, 2))
        self.label_tree = ctk.CTkLabel(fi, text=TEXTS["label_not_selected"],
                                       font=(_FONT_UI, 13),
                                       wraplength=230,
                                       text_color=self.COLORS['text_tertiary'])
        self.label_tree.pack(anchor='w', padx=4, pady=(0, 8))

        self.btn_output = _obtn(fi, TEXTS["btn_output_folder"],
                                self.select_output_folder,
                                self.COLORS['accent_blue'],
                                font=(_FONT_UI, 13, "bold"), height=36)
        self.btn_output.pack(fill='x', pady=(0, 2))
        self.label_output = ctk.CTkLabel(fi, text=TEXTS["label_not_selected"],
                                         font=(_FONT_UI, 13),
                                         wraplength=230,
                                         text_color=self.COLORS['text_tertiary'])
        self.label_output.pack(anchor='w', padx=4)

        # ── Resultados ───────────────────────────────────────────────
        ri = _sec(_sb, TEXTS["section_results"])

        self.btn_results = _obtn(ri, TEXTS["btn_view_results"],
                                  self._open_results_viewer,
                                  self.COLORS['accent_blue'],
                                  font=(_FONT_UI, 13, "bold"), height=36)
        self.btn_results.pack(fill='x', pady=(0, 6))
        self.btn_results.configure(state="disabled")

        self.btn_update_results = _obtn(ri, TEXTS["btn_update_results"],
                                         self._update_results_files,
                                         self.COLORS['success'],
                                         font=(_FONT_UI, 13, "bold"), height=36)
        self.btn_update_results.pack(fill='x', pady=(0, 2))
        self.label_update_results = ctk.CTkLabel(ri,
                                                  text=TEXTS["label_update_results_hint"],
                                                  font=(_FONT_UI, 13, "italic"),
                                                  wraplength=230,
                                                  text_color=self.COLORS['text_muted'])
        self.label_update_results.pack(anchor='w', padx=4)

        # ── Configurações ────────────────────────────────────────────
        ci = _sec(_sb, TEXTS["section_config"])

        ctk.CTkLabel(ci, text=TEXTS["label_omega_initial"],
                     font=(_FONT_UI, 13, "bold"), anchor='w',
                     text_color=self.COLORS['text_secondary']).pack(anchor='w')
        self.omega_label = ctk.CTkLabel(ci, text="")  # kept for compat, unused
        self.entry_omega = ctk.CTkEntry(ci, placeholder_text="0.5",
                                        fg_color=self.COLORS['bg_card_hover'],
                                        border_color=self.COLORS['border_hover'],
                                        border_width=1,
                                        corner_radius=8,
                                        text_color=self.COLORS['text_primary'],
                                        height=32)
        self.entry_omega.insert(0, "0.5")
        self.entry_omega.pack(fill='x', pady=(4, 10))

        # Frequências de códons (global; a janela "editar" de cada modelo pode mudar)
        row_cf = ctk.CTkFrame(ci, fg_color='transparent')
        row_cf.pack(fill='x')
        ctk.CTkLabel(row_cf, text=TEXTS["label_codonfreq"],
                     font=(_FONT_UI, 13, "bold"), anchor='w',
                     text_color=self.COLORS['text_secondary']).pack(side='left')
        ctk.CTkButton(row_cf, text="?", width=22, height=22, corner_radius=11,
                      font=(_FONT_UI, 12, "bold"), fg_color=self.COLORS['border'],
                      hover_color=self.COLORS['border_hover'],
                      text_color=self.COLORS['text_primary'],
                      command=lambda: self._show_help(
                          TEXTS["label_codonfreq"], TEXTS["label_codonfreq_hint"])
                      ).pack(side='right')
        self.codonfreq_var = ctk.StringVar(value=codonfreq_label(DEFAULT_CODONFREQ))
        ctk.CTkOptionMenu(ci, variable=self.codonfreq_var,
                          values=[codonfreq_label(v) for v, _, _ in CODONFREQ_OPTIONS],
                          fg_color=self.COLORS['bg_card_hover'], button_color=self.COLORS['border_hover'],
                          text_color=self.COLORS['text_primary'], height=30
                          ).pack(fill='x', pady=(4, 10))

        ctk.CTkLabel(ci, text=TEXTS["label_ncatg"],
                     font=(_FONT_UI, 13, "bold"), anchor='w',
                     text_color=self.COLORS['text_secondary']).pack(anchor='w')
        self.entry_ncatg = ctk.CTkEntry(ci, fg_color=self.COLORS['bg_card_hover'],
                                        border_color=self.COLORS['border_hover'], border_width=1,
                                        corner_radius=8, text_color=self.COLORS['text_primary'],
                                        height=32)
        self.entry_ncatg.insert(0, str(DEFAULT_CTL_PARAMS['ncatG']))
        self.entry_ncatg.pack(fill='x', pady=(4, 10))

        row_to = ctk.CTkFrame(ci, fg_color='transparent')
        row_to.pack(fill='x')
        ctk.CTkLabel(row_to, text=TEXTS["label_timeout"],
                     font=(_FONT_UI, 13, "bold"), anchor='w',
                     text_color=self.COLORS['text_secondary']).pack(side='left')
        ctk.CTkButton(row_to, text="?", width=22, height=22, corner_radius=11,
                      font=(_FONT_UI, 12, "bold"), fg_color=self.COLORS['border'],
                      hover_color=self.COLORS['border_hover'],
                      text_color=self.COLORS['text_primary'],
                      command=lambda: self._show_help(
                          TEXTS["label_timeout"], TEXTS["label_timeout_hint"])
                      ).pack(side='right')
        self.entry_timeout = ctk.CTkEntry(ci, placeholder_text=TEXTS["label_timeout_auto"],
                                          fg_color=self.COLORS['bg_card_hover'],
                                          border_color=self.COLORS['border_hover'], border_width=1,
                                          corner_radius=8, text_color=self.COLORS['text_primary'],
                                          height=32)
        self.entry_timeout.pack(fill='x', pady=(4, 10))

        # Remover gaps toggle
        self.cleandata_var = ctk.BooleanVar(value=True)
        row_g = ctk.CTkFrame(ci, fg_color='transparent')
        row_g.pack(fill='x', pady=(0, 8))
        ctk.CTkLabel(row_g, text=TEXTS["label_remove_gaps"],
                     font=(_FONT_UI, 12),
                     text_color=self.COLORS['text_secondary']).pack(side='left')
        self.cb_cleandata = ctk.CTkSwitch(
            row_g, text="",
            variable=self.cleandata_var,
            onvalue=True, offvalue=False,
            switch_width=36, switch_height=18,
            progress_color=self.COLORS['success'],
            button_color='#f0fdf4',
            button_hover_color='#dcfce7',
            fg_color=self.COLORS['border'])
        self.cb_cleandata.pack(side='right')
        ctk.CTkButton(row_g, text="?", width=18, height=18, corner_radius=9,
                      font=(_FONT_UI, 12, "bold"), fg_color=self.COLORS['border'],
                      hover_color=self.COLORS['border_hover'],
                      text_color=self.COLORS['text_primary'],
                      command=lambda: self._show_help(
                          TEXTS["label_remove_gaps"], TEXTS["label_remove_gaps_hint"])
                      ).pack(side='right', padx=(0, 4))

        # CPU slider
        ctk.CTkLabel(ci, text=TEXTS["label_cpus"],
                     font=(_FONT_UI, 13, "bold"), anchor='w',
                     text_color=self.COLORS['text_secondary']).pack(anchor='w', pady=(0, 4))
        cores_row = ctk.CTkFrame(ci, fg_color='transparent')
        cores_row.pack(fill='x', pady=(0, 2))
        _max = CodemlBatchAnalysis.available_cores()
        self.cores_slider = ctk.CTkSlider(
            cores_row, from_=1, to=max(2, _max),
            number_of_steps=max(1, _max - 1),
            variable=self.cores_var,
            command=self._update_cores_label,
            button_color=self.COLORS['accent_blue'],
            progress_color=self.COLORS['accent_blue'])
        self.cores_slider.pack(fill='x', side='left', expand=True)
        disable_mouse_wheel(self.cores_slider)   # rolar o painel não muda o valor
        self.cores_disp = ctk.CTkLabel(
            cores_row, text=f"{_max}×",
            font=(_FONT_UI, 12, "bold"),
            text_color=self.COLORS['accent_blue_light'],
            width=40)
        self.cores_disp.pack(side='right', padx=(6, 0))
        ctk.CTkLabel(ci, text=TEXTS["label_cpus_detected"].format(n=_max),
                     font=(_FONT_UI, 13),
                     text_color=self.COLORS['text_muted']).pack(anchor='w')

        # Ignorar Stop Codons toggle
        row_stops = ctk.CTkFrame(ci, fg_color='transparent')
        row_stops.pack(fill='x', pady=(10, 0))
        ctk.CTkLabel(row_stops, text=TEXTS["label_ignore_stops"],
                     font=(_FONT_UI, 12),
                     text_color=self.COLORS['text_secondary']).pack(side='left')
        self.cb_ignore_stops = ctk.CTkSwitch(
            row_stops, text="",
            variable=self.ignore_stop_codons_var,
            onvalue=True, offvalue=False,
            switch_width=36, switch_height=18,
            progress_color=self.COLORS['accent_cyan'],
            button_color='#f0fdff',
            button_hover_color='#cffafe',
            fg_color=self.COLORS['border'])
        self.cb_ignore_stops.pack(side='right')
        ctk.CTkButton(row_stops, text="?", width=18, height=18, corner_radius=9,
                      font=(_FONT_UI, 12, "bold"), fg_color=self.COLORS['border'],
                      hover_color=self.COLORS['border_hover'],
                      text_color=self.COLORS['text_primary'],
                      command=lambda: self._show_help(
                          TEXTS["label_ignore_stops"], TEXTS["label_ignore_stops_hint"])
                      ).pack(side='right', padx=(0, 4))

        # Incluir M8a como nulo automático do M8
        row_m8a = ctk.CTkFrame(ci, fg_color='transparent')
        row_m8a.pack(fill='x', pady=(10, 0))
        ctk.CTkLabel(row_m8a, text=TEXTS["label_include_m8a"],
                     font=(_FONT_UI, 12),
                     text_color=self.COLORS['text_secondary']).pack(side='left')
        ctk.CTkSwitch(row_m8a, text="", variable=self.include_m8a_var,
                      onvalue=True, offvalue=False, switch_width=36, switch_height=18,
                      progress_color=self.COLORS['accent_cyan'], button_color='#f0fdff',
                      button_hover_color='#cffafe', fg_color=self.COLORS['border']).pack(side='right')
        ctk.CTkButton(row_m8a, text="?", width=18, height=18, corner_radius=9,
                      font=(_FONT_UI, 12, "bold"), fg_color=self.COLORS['border'],
                      hover_color=self.COLORS['border_hover'],
                      text_color=self.COLORS['text_primary'],
                      command=lambda: self._show_help(
                          TEXTS["label_include_m8a"], TEXTS["label_include_m8a_hint"])
                      ).pack(side='right', padx=(0, 4))

        # Poda automática de árvore toggle
        row_prune = ctk.CTkFrame(ci, fg_color='transparent')
        row_prune.pack(fill='x', pady=(10, 0))
        ctk.CTkLabel(row_prune, text=TEXTS["label_auto_prune"],
                     font=(_FONT_UI, 12),
                     text_color=self.COLORS['text_secondary']).pack(side='left')
        self.cb_auto_prune = ctk.CTkSwitch(
            row_prune, text="",
            variable=self.auto_prune_tree_var,
            onvalue=True, offvalue=False,
            switch_width=36, switch_height=18,
            progress_color=self.COLORS['accent_cyan'],
            button_color='#f0fdff',
            button_hover_color='#cffafe',
            fg_color=self.COLORS['border'])
        self.cb_auto_prune.pack(side='right')
        ctk.CTkButton(row_prune, text="?", width=18, height=18, corner_radius=9,
                      font=(_FONT_UI, 12, "bold"), fg_color=self.COLORS['border'],
                      hover_color=self.COLORS['border_hover'],
                      text_color=self.COLORS['text_primary'],
                      command=lambda: self._show_help(
                          TEXTS["label_auto_prune"], TEXTS["label_auto_prune_hint"])
                      ).pack(side='right', padx=(0, 4))

        self.main_frame = ctk.CTkFrame(self, fg_color=self.COLORS['bg_dark'])
        self.main_frame.pack(side="right", fill="both", expand=True, padx=20, pady=20)

        self.files_hint = ctk.CTkLabel(self.main_frame, text=TEXTS["hint_select_files"],
                                       font=(_FONT_UI, 13, "bold"), anchor='w', justify='left',
                                       text_color=self.COLORS['warning'], wraplength=900)
        self.files_hint.pack(fill='x', pady=(0, 6))
        self.tabs = ctk.CTkTabview(self.main_frame, fg_color=self.COLORS['bg_card'],
                                   segmented_button_fg_color=self.COLORS['bg_card'],
                                   segmented_button_selected_color=PALETTE['accent_fill'],
                                   segmented_button_selected_hover_color=self.COLORS['accent_blue_hover'],
                                   text_color='#ffffff',
                                   corner_radius=10)
        self.tabs.pack(fill="both", expand=True, padx=0, pady=(0, 15))
        self.tabs.add(TEXTS["tab_site_models"])
        self.tabs.add(TEXTS["tab_branch_model"])
        self.tabs.add(TEXTS["tab_branchsite"])

        self.model_vars = {}
        self.model_ctl_labels = {}
        self.model_checkboxes = {}
        self.model_gear_buttons = {}

        self._setup_model_list()

        self.ctrl_frame = ctk.CTkFrame(self.main_frame, fg_color=self.COLORS['bg_card'],
                                       corner_radius=10, border_width=1,
                                       border_color=self.COLORS['border'])
        self.ctrl_frame.pack(fill="x", padx=0, pady=(0, 15))

        # ── Status + neutral-models row ───────────────────────────────
        status_bar = ctk.CTkFrame(self.ctrl_frame, fg_color=self.COLORS['bg_sidebar'], corner_radius=8)
        status_bar.pack(fill="x", padx=12, pady=(10, 6))

        self.status_indicator = ctk.CTkLabel(
            status_bar, text=TEXTS["status_ready"],
            font=(_FONT_UI, 13, "bold"),
            text_color=self.COLORS['text_tertiary']
        )
        self.status_indicator.pack(side="left", padx=(14, 20), pady=8)

        self.stop_label = ctk.CTkLabel(
            status_bar, text=TEXTS["status_stops_template"].format(n=0),
            font=(_FONT_UI, 13, "bold"),
            text_color=self.COLORS['danger']
        )
        self.stop_label.pack(side="left")

        neutral_row = ctk.CTkFrame(status_bar, fg_color='transparent')
        neutral_row.pack(side="right", padx=(0, 6))
        ctk.CTkLabel(neutral_row, text=TEXTS["label_neutral_models"],
                     font=(_FONT_UI, 13),
                     text_color=self.COLORS['text_secondary']).pack(side='left', padx=(0, 8))
        neutral_sw = ctk.CTkSwitch(
            neutral_row, text="",
            variable=self.include_neutral_models,
            onvalue=True, offvalue=False,
            switch_width=34, switch_height=17,
            progress_color=self.COLORS['accent_blue'],
            button_color='#f0f0ff',
            button_hover_color='#e0e0ff',
            fg_color=self.COLORS['border'])
        neutral_sw.pack(side='right')

        help_btn = ctk.CTkButton(
            status_bar, text="?", width=26, height=26,
            font=(_FONT_UI, 12, "bold"),
            fg_color=self.COLORS['border'],
            hover_color=self.COLORS['accent_blue'],
            text_color=self.COLORS['accent_blue_light'],
            border_width=1, border_color=self.COLORS['accent_blue'],
            command=self._show_neutral_models_info,
            corner_radius=6
        )
        help_btn.pack(side="right", padx=(0, 12))

        # Barra de progresso "gene X de N"
        prog_row = ctk.CTkFrame(self.ctrl_frame, fg_color='transparent')
        prog_row.pack(fill='x', padx=14, pady=(0, 2))
        self.progress_bar = ctk.CTkProgressBar(prog_row, height=10,
                                               progress_color=self.COLORS['success'])
        self.progress_bar.pack(side='left', fill='x', expand=True, padx=(0, 10))
        self.progress_bar.set(0)
        self.progress_label = ctk.CTkLabel(prog_row, text=TEXTS["progress_idle"],
                                           font=(_FONT_UI, 12),
                                           text_color=self.COLORS['text_secondary'])
        self.progress_label.pack(side='right')

        # Botões de controle
        btn_frame = ctk.CTkFrame(self.ctrl_frame, fg_color='transparent')
        btn_frame.pack(fill="x", padx=0, pady=(0, 10))

        def _action_btn(parent, text, cmd, color):
            # hover = tom escuro da cor, para o texto (na cor) continuar legível
            return ctk.CTkButton(
                parent, text=text, command=cmd,
                fg_color=self.COLORS['bg_card'],
                hover_color=hover_tint(color, self.COLORS['bg_card']),
                text_color=color,
                text_color_disabled=self.COLORS['text_tertiary'],
                border_width=1, border_color=color,
                font=(_FONT_UI, 13, "bold"),
                height=44, corner_radius=8)

        self.btn_run = _action_btn(btn_frame, TEXTS["btn_run"],
                                   self.start_analysis, self.COLORS['success'])
        self.btn_run.pack(side="left", fill="both", expand=True, padx=(12, 4), pady=10)

        self.btn_pause = _action_btn(btn_frame, TEXTS["btn_pause"],
                                     self._toggle_pause, self.COLORS['warning'])
        self.btn_pause.pack(side="left", fill="both", expand=True, padx=4, pady=10)

        self.btn_stop = _action_btn(btn_frame, TEXTS["btn_stop"],
                                    self._stop_analysis, self.COLORS['danger'])
        self.btn_stop.pack(side="left", fill="both", expand=True, padx=4, pady=10)

        self.btn_open_output = _action_btn(btn_frame, TEXTS["btn_open_output"],
                                           self._open_output_folder, self.COLORS['info'])
        self.btn_open_output.pack(side="left", fill="both", expand=True, padx=(4, 12), pady=10)

        # ═══ LOG FRAME COM HEADER ═══
        log_container = ctk.CTkFrame(self.main_frame, fg_color='transparent')
        log_container.pack(fill="both", expand=True, padx=0, pady=0)
        
        log_header = ctk.CTkFrame(log_container, fg_color='transparent')
        log_header.pack(fill="x", padx=0, pady=(10, 0))
        
        ctk.CTkLabel(log_header, text=TEXTS["log_header_title"], font=(_FONT_UI, 13, "bold"),
                    text_color=self.COLORS['text_secondary']).pack(side="left", padx=0)

        ctk.CTkLabel(log_header, text="·",
                    font=(_FONT_UI, 13), text_color=self.COLORS['text_muted']).pack(side="left", padx=8)

        ctk.CTkLabel(log_header, text=TEXTS["log_header_subtitle"],
                    font=(_FONT_UI, 13), text_color=self.COLORS['text_muted']).pack(side="left", padx=0)
        ctk.CTkCheckBox(log_header, text=TEXTS["chk_show_details"], variable=self.show_details_var,
                        command=self._toggle_details, font=(_FONT_UI, 12),
                        text_color=self.COLORS['text_secondary'], checkbox_width=18,
                        checkbox_height=18).pack(side="right")

        self.log = ctk.CTkTextbox(log_container, font=(_FONT_MONO, 12), wrap='word',
                                  fg_color=self.COLORS['bg_dark'],
                                  text_color=self.COLORS['text_secondary'],
                                  border_color=self.COLORS['border'],
                                  border_width=1,
                                  corner_radius=8)
        self.log.pack(fill="both", expand=True, padx=0, pady=(8, 0))

        self.log.tag_config("success", foreground=self.COLORS['success_light'])
        self.log.tag_config("error", foreground=self.COLORS['danger'])
        self.log.tag_config("warning", foreground=self.COLORS['warning'])
        self.log.tag_config("info", foreground=self.COLORS['accent_cyan'])
        self.log.tag_config("header", foreground=self.COLORS['accent_blue_light'])
        self.log.tag_config("debug", foreground=self.COLORS['text_tertiary'], elide=True)

        # Mensagem inicial — editar em gui_texts.py › "log_welcome"
        self.log.insert("end", TEXTS["log_welcome"])
        backend_messages.set_language(get_language())

        self._update_models_state()
        self._poll_stop_count()
        self.bind_all("<Control-q>", lambda e: self.destroy())

    def _global_ctl_options(self) -> dict:
        """Opções globais do .ctl escolhidas em Configurações."""
        try:
            ncatg = int(self.entry_ncatg.get().strip() or DEFAULT_CTL_PARAMS['ncatG'])
        except (ValueError, AttributeError):
            ncatg = DEFAULT_CTL_PARAMS['ncatG']
        try:
            cf = parse_codonfreq_label(self.codonfreq_var.get())
        except (ValueError, AttributeError):
            cf = DEFAULT_CODONFREQ
        return {'CodonFreq': cf, 'ncatG': ncatg, 'kappa': DEFAULT_CTL_PARAMS['kappa']}

    def _open_output_folder(self):
        if not self.output_folder:
            self.append_log(TEXTS["log_no_output_folder"], 'warn')
            return
        self.output_folder.mkdir(parents=True, exist_ok=True)
        open_folder(self.output_folder)

    def _toggle_details(self):
        """Mostra/esconde as mensagens de depuração no log."""
        show = self.show_details_var.get()
        try:
            self.log.tag_config("debug", elide=not show)
        except Exception:
            pass

    def _show_help(self, title: str, body: str) -> None:
        win = ctk.CTkToplevel(self)
        win.title(title)
        win.resizable(False, False)
        win.transient(self)
        win.after(50, lambda: (win.grab_set(), win.focus_set()))
        win.bind("<Escape>", lambda e: win.destroy())
        ctk.CTkLabel(win, text=title, font=(_FONT_UI, 14, "bold"),
                     text_color=self.COLORS['text_primary']).pack(padx=18, pady=(16, 4), anchor='w')
        ctk.CTkLabel(win, text=body, font=(_FONT_UI, 13), wraplength=420,
                     justify='left',
                     text_color=self.COLORS['text_secondary']).pack(padx=18, pady=(0, 12))
        ctk.CTkButton(win, text="OK", width=80, command=win.destroy, text_color='#ffffff',
                      fg_color=PALETTE['accent_fill']).pack(pady=(0, 14))

    def _update_cores_label(self, value=None):
        n = int(self.cores_var.get())
        self.cores_disp.configure(text=f"{n}×")

    def _setup_model_list(self):
        MODEL_META = {
            'M0': {
                'color': '#3b82f6',
                'desc': 'Um único ω para todo o gene. Não detecta variação entre sítios. '
                        'Usado como baseline e como modelo nulo para o Branch Model.'
            },
            'M1a': {
                'color': '#06b6d4',
                'desc': 'Permite purificação (0 < ω < 1) e neutralidade (ω = 1) — sem seleção positiva. '
                        'Modelo nulo (null) para M2a no teste LRT.'
            },
            'M2a': {
                'color': '#10b981',
                'desc': 'Acrescenta classe com ω > 1 ao M1a. LRT M2a vs M1a indica seleção positiva por sítio. '
                        'BEB/NEB identificam os sítios sob seleção.'
            },
            'M7': {
                'color': '#8b5cf6',
                'desc': 'ω segue distribuição Beta contínua — todos os sítios têm ω ≤ 1. '
                        'Modelo nulo mais flexível para comparar com M8.'
            },
            'M8': {
                'color': '#ec4899',
                'desc': 'Beta(p,q) + classe discreta com ω > 1. LRT M8 vs M7 é o teste mais '
                        'robusto para seleção positiva por sítio. Requer comparação com M7.'
            },
            'M8a': {
                'color': '#f472b6',
                'desc': '',
            },
            'Branch': {
                'color': '#f59e0b',
                'desc': 'Estima ω independente por ramo etiquetado. LRT com M0 testa se há '
                        'pressão seletiva diferente nas linhagens marcadas.'
            },
            'Branch-site': {
                'color': '#ef4444',
                'desc': 'Detecta seleção positiva em sítios específicos do ramo foreground (#1). '
                        'Combina variação por sítio e por linhagem — o teste mais poderoso.'
            },
            'Branch-site_null': {
                'color': '#6d6d6d',
                'desc': 'Versão restrita do Branch-site com ω₂=1 fixado. '
                        'Adicionado automaticamente como null para o LRT do Branch-site.'
            },
        }

        models = {
            TEXTS["tab_site_models"]:  ['M0', 'M1a', 'M2a', 'M7', 'M8', 'M8a'],
            TEXTS["tab_branch_model"]: ['Branch'],
            TEXTS["tab_branchsite"]:   ['Branch-site', 'Branch-site_null'],
        }

        for tab_name, codes in models.items():
            outer = ctk.CTkFrame(
                self.tabs.tab(tab_name),
                fg_color='transparent'
            )
            outer.pack(fill="x", padx=8, pady=8)
            outer.grid_columnconfigure(0, weight=1)
            outer.grid_columnconfigure(1, weight=1)

            for idx, code in enumerate(codes):
                meta   = MODEL_META.get(code, {'color': self.COLORS['accent_blue'], 'desc': ''})
                accent = meta['color']

                row = idx // 2
                col = idx % 2

                # ── Compact card (2-column grid) ──────────────────────────────
                card = ctk.CTkFrame(outer, fg_color=self.COLORS['bg_feed'],
                                    corner_radius=8, border_width=1,
                                    border_color=self.COLORS['border'],
                                    height=64)
                card.pack_propagate(False)
                card.grid(row=row, column=col, sticky='ew', pady=3, padx=3)

                # Left accent bar
                accent_bar = ctk.CTkFrame(card, fg_color=accent, width=4, corner_radius=2)
                accent_bar.pack(side="left", fill="y", padx=(5, 6), pady=5)
                accent_bar.pack_propagate(False)

                # Right icon buttons (stacked)
                btn_col = ctk.CTkFrame(card, fg_color='transparent')
                btn_col.pack(side="right", padx=(0, 4), pady=4)

                info_btn = ctk.CTkButton(
                    btn_col, text="?", width=22, height=22,
                    fg_color=self.COLORS['bg_card_hover'],
                    hover_color=hover_tint(accent, self.COLORS['bg_card_hover']),
                    text_color=accent,
                    corner_radius=4,
                    border_width=1, border_color=self.COLORS['border_hover'],
                    font=(_FONT_UI, 13, "bold"),
                    command=lambda c=code: self._show_model_info(c)
                )
                info_btn.pack(pady=(0, 2))

                gear = ctk.CTkButton(
                    btn_col, text=TEXTS["cfg_btn"], width=52, height=22,
                    fg_color=self.COLORS['bg_card'],
                    hover_color=hover_tint(self.COLORS['accent_purple'], self.COLORS['bg_card']),
                    text_color=self.COLORS['text_primary'],
                    text_color_disabled=self.COLORS['text_tertiary'],
                    corner_radius=4,
                    border_width=1, border_color=self.COLORS['border_hover'],
                    font=(_FONT_UI, 11),
                    command=lambda c=code: self._open_config_window(c)
                )
                gear.pack()

                # Content area (toggle + status)
                content = ctk.CTkFrame(card, fg_color='transparent')
                content.pack(side="left", fill="both", expand=True, pady=4)

                display_name = self.codeml_backend.MODEL_CONFIGS[code].get('display_name', code)
                var = ctk.BooleanVar(value=False)
                cb = ctk.CTkSwitch(
                    content,
                    text=f"  {display_name}",
                    variable=var,
                    onvalue=True, offvalue=False,
                    command=self._update_models_state,
                    switch_width=32, switch_height=16,
                    progress_color=accent,
                    button_color='#f0f0ff',
                    button_hover_color='white',
                    fg_color=self.COLORS['border'],
                    text_color=self.COLORS['text_primary'],
                    font=(_FONT_UI, 12, "bold")
                )
                cb.pack(anchor='w')

                lbl = ctk.CTkLabel(content, text=TEXTS["model_desc"].get(code, TEXTS["model_status_default"]),
                                   font=(_FONT_UI, 11), anchor='w', justify='left',
                                   wraplength=330,
                                   text_color=self.COLORS['text_tertiary'])
                lbl.pack(anchor='w', padx=(2, 0))

                self.model_vars[code]         = var
                self.model_ctl_labels[code]   = lbl
                self.model_checkboxes[code]   = cb
                self.model_gear_buttons[code] = gear
        
        # ═══ BRANCH: Botão de etiquetagem ═══
        branch_tab = self.tabs.tab(TEXTS["tab_branch_model"])

        self.btn_label_branch = ctk.CTkButton(
            branch_tab,
            text=TEXTS["btn_label_branch"],
            fg_color=self.COLORS['bg_card_hover'],
            hover_color=self.COLORS['accent_blue'],
            command=lambda: self._open_tree_labeler(mode='branch'),
            height=40,
            font=(_FONT_UI, 13, "bold"),
            text_color=self.COLORS['accent_blue_light'],
            border_width=1, border_color=self.COLORS['border_hover'],
            corner_radius=8
        )
        self.btn_label_branch.pack(fill='x', padx=12, pady=(16, 12))

        # ═══ BRANCHSITE: Botão de etiquetagem ═══
        branchsite_tab = self.tabs.tab(TEXTS["tab_branchsite"])

        self.btn_label_branchsite = ctk.CTkButton(
            branchsite_tab,
            text=TEXTS["btn_label_branchsite"],
            fg_color=self.COLORS['bg_card_hover'],
            hover_color=self.COLORS['accent_blue'],
            command=lambda: self._open_tree_labeler(mode='branchsite'),
            height=40,
            font=(_FONT_UI, 13, "bold"),
            text_color=self.COLORS['accent_blue_light'],
            border_width=1, border_color=self.COLORS['border_hover'],
            corner_radius=8
        )
        self.btn_label_branchsite.pack(fill='x', padx=12, pady=(16, 12))

    def _open_config_window(self, code):
        default_config = CodemlBatchAnalysis.MODEL_CONFIGS.get(code, {})
        ModelConfigWindow(self, code, default_config)

    def _open_tree_labeler(self, mode: str = 'branchsite'):
        if self.tree_file is None:
            self.append_log(TEXTS["log_no_tree_selected"])
            return

        try:
            TreeLabelWindow(self, self.tree_file, mode=mode)
        except Exception as e:
            self.append_log(tr("Erro ao abrir a marcação de ramos: ", "Error opening branch labelling: ") + f"{e}", "error")
            self.append_log(traceback.format_exc(), "debug")

    def select_input_folder(self):
        start = str(self.input_folder) if self.input_folder else str(Path.home())
        path = filedialog.askdirectory(initialdir=start, title=TEXTS["btn_input_folder"])
        if path:
            self.input_folder = Path(path)
            files = list_alignment_files(self.input_folder)
            chosen, _ = group_by_gene(files)
            self.per_gene_trees = discover_per_gene_trees(self.input_folder, genes=set(chosen))
            if files:
                names = ", ".join(f.name for f in files[:4]) + (" …" if len(files) > 4 else "")
                text = TEXTS["label_found_alignments"].format(n=len(files), names=names)
                if self.per_gene_trees:
                    text += "\n" + TEXTS["label_per_gene_trees"].format(n=len(self.per_gene_trees),
                                                                       total=len(chosen))
                color = self.COLORS['text_secondary']
            else:
                text = TEXTS["label_no_alignments"]
                color = self.COLORS['warning']
            self.label_input.configure(text=f"{self.input_folder.name}\n{text}", text_color=color)
            self._update_models_state()

    def select_tree_file(self):
        start = str(self.tree_file.parent) if self.tree_file else (
            str(self.input_folder) if self.input_folder else str(Path.home()))
        path = filedialog.askopenfilename(
            initialdir=start, title=TEXTS["btn_tree_file"],
            filetypes=[('Newick', '*.nwk *.tree *.tre *.newick *.txt'), ('*', '*.*')])
        if path:
            self.tree_file = Path(path)
            self.label_tree.configure(text=str(self.tree_file.name),
                                      text_color=self.COLORS['text_secondary'])
            self._update_models_state()

    def select_output_folder(self):
        """Aceita uma pasta que ainda não existe (digitada no diálogo) e a cria."""
        start = str(self.output_folder.parent if self.output_folder else
                    (self.input_folder.parent if self.input_folder else Path.home()))
        path = filedialog.askdirectory(initialdir=start, mustexist=False,
                                       title=TEXTS["dialog_choose_output"])
        if not path:
            return
        folder = Path(path)
        created = not folder.exists()
        try:
            folder.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            show_message(self, TEXTS["msg_error"],
                         TEXTS["msg_output_folder_error"].format(error=exc), 'error')
            return
        self.output_folder = folder
        label = (TEXTS["label_output_created"].format(name=folder.name) if created else folder.name)
        self.label_output.configure(text=label, text_color=self.COLORS['text_secondary'])
        self._update_models_state()

    def _show_neutral_models_info(self):
        """Pares de modelos nulo/alternativo do LRT (texto em gui_texts › lrt_pairs)."""
        win = ctk.CTkToplevel(self)
        win.title(TEXTS["neutral_window_title"])
        fit_to_screen(win, 760, 680, min_w=520, min_h=400)
        win.transient(self)
        win.after(50, lambda: (win.grab_set(), win.focus_set()))
        win.bind("<Escape>", lambda e: win.destroy())
        win.configure(fg_color=self.COLORS['bg_dark'])

        hdr = ctk.CTkFrame(win, fg_color='transparent')
        hdr.pack(fill='x', padx=18, pady=(16, 4))
        ctk.CTkLabel(hdr, text=TEXTS["neutral_header"], font=(_FONT_UI, 15, "bold"),
                     text_color=self.COLORS['accent_blue_light']).pack(anchor='w')
        ctk.CTkLabel(hdr, text=TEXTS["neutral_intro"], font=(_FONT_UI, 13), wraplength=700,
                     justify='left', text_color=self.COLORS['text_secondary']).pack(anchor='w', pady=(4, 0))

        scroll = ctk.CTkScrollableFrame(win, fg_color='transparent')
        scroll.pack(fill='both', expand=True, padx=14, pady=6)
        for null, alt, color, title, test, detail, note in TEXTS["lrt_pairs"]:
            card = ctk.CTkFrame(scroll, fg_color=self.COLORS['bg_card'], corner_radius=10,
                                border_width=1, border_color=self.COLORS['border'])
            card.pack(fill='x', padx=6, pady=5)
            ctk.CTkFrame(card, fg_color=color, width=5, corner_radius=2).pack(
                side='left', fill='y', padx=(6, 10), pady=8)
            content = ctk.CTkFrame(card, fg_color='transparent')
            content.pack(side='left', fill='both', expand=True, pady=10, padx=(0, 10))
            row = ctk.CTkFrame(content, fg_color='transparent')
            row.pack(fill='x', anchor='w')
            ctk.CTkLabel(row, text=title, font=(_FONT_UI, 13, "bold"), text_color=color).pack(side='left')
            ctk.CTkLabel(row, text="   " + TEXTS["neutral_null_alt"].format(null=null, alt=alt),
                         font=(_FONT_UI, 12), text_color=self.COLORS['text_secondary']).pack(side='left')
            ctk.CTkLabel(content, text=test, font=(_FONT_UI, 13, "bold"),
                         text_color=self.COLORS['text_primary'], anchor='w').pack(anchor='w', pady=(4, 2))
            ctk.CTkLabel(content, text=detail, font=(_FONT_UI, 12), wraplength=600, justify='left',
                         text_color=self.COLORS['text_secondary'], anchor='w').pack(anchor='w', pady=(0, 3))
            ctk.CTkLabel(content, text=TEXTS["neutral_note"].format(note=note), font=(_FONT_UI, 12, "italic"),
                         wraplength=600, justify='left', text_color=self.COLORS['text_tertiary'],
                         anchor='w').pack(anchor='w')
        ctk.CTkLabel(win, text=TEXTS["neutral_footer"], font=(_FONT_UI, 12),
                     text_color=self.COLORS['success'], wraplength=700).pack(padx=18, pady=12)

    def _show_model_info(self, model_code: str):
        """Mostra informações detalhadas sobre um modelo específico"""
        # Obter informações do modelo
        model_info = self.codeml_backend.MODEL_INFO.get(model_code)
        if not model_info:
            return
        
        # Criar janela de informações
        info_window = ctk.CTkToplevel(self)
        info_window.title(TEXTS["model_config_header"].format(model_code=model_code))
        fit_to_screen(info_window, 750, 600, min_w=480, min_h=360)
        info_window.transient(self)
        info_window.after(50, lambda: (info_window.grab_set(), info_window.focus_set()))
        info_window.bind("<Escape>", lambda e: info_window.destroy())
        
        # Header com nome do modelo
        display_name = self.codeml_backend.MODEL_CONFIGS[model_code].get('display_name', model_code)
        header = ctk.CTkLabel(info_window, 
                             text=f"{model_code} - {model_info.get('full_name', '')}",
                             font=(_FONT_UI, 13, "bold"),
                             text_color=self.COLORS['accent_blue_light'])
        header.pack(padx=15, pady=15)
        
        # Scrollable frame para conteúdo
        scroll_frame = ctk.CTkScrollableFrame(info_window, fg_color=self.COLORS['bg_feed'],
                                             corner_radius=8)
        scroll_frame.pack(fill="both", expand=True, padx=12, pady=(0, 12))
        
        # Tipo de teste
        test_type_label = ctk.CTkLabel(scroll_frame, text=TEXTS["model_info_test_type"],
                                      font=(_FONT_UI, 13, "bold"),
                                      text_color=self.COLORS['accent_purple'])
        test_type_label.pack(anchor="w", padx=8, pady=(8, 2))
        
        test_type_value = ctk.CTkLabel(scroll_frame, text=model_info.get('test_type', ''),
                                      font=(_FONT_UI, 12),
                                      text_color=self.COLORS['text_primary'],
                                      wraplength=700, justify="left")
        test_type_value.pack(anchor="w", padx=25, pady=(0, 8))
        
        # Parâmetros
        params_label = ctk.CTkLabel(scroll_frame, text=TEXTS["model_info_params"],
                                   font=(_FONT_UI, 13, "bold"),
                                   text_color=self.COLORS['accent_purple'])
        params_label.pack(anchor="w", padx=8, pady=(8, 2))
        
        params_value = ctk.CTkLabel(scroll_frame, text=model_info.get('parameters', ''),
                                   font=(_FONT_UI, 12),
                                   text_color=self.COLORS['text_primary'],
                                   wraplength=700, justify="left")
        params_value.pack(anchor="w", padx=25, pady=(0, 8))
        
        # Propósito
        purpose_label = ctk.CTkLabel(scroll_frame, text=TEXTS["model_info_purpose"],
                                    font=(_FONT_UI, 13, "bold"),
                                    text_color=self.COLORS['accent_cyan'])
        purpose_label.pack(anchor="w", padx=8, pady=(8, 2))
        
        purpose_value = ctk.CTkLabel(scroll_frame, text=model_info.get('purpose', ''),
                                    font=(_FONT_UI, 12),
                                    text_color=self.COLORS['text_primary'],
                                    wraplength=700, justify="left")
        purpose_value.pack(anchor="w", padx=25, pady=(0, 8))
        
        # Interpretação
        interp_label = ctk.CTkLabel(scroll_frame, text=TEXTS["model_info_interpretation"],
                                   font=(_FONT_UI, 13, "bold"),
                                   text_color=self.COLORS['success'])
        interp_label.pack(anchor="w", padx=8, pady=(8, 2))

        interp_value = ctk.CTkLabel(scroll_frame, text=model_info.get('interpretation', ''),
                                   font=(_FONT_UI, 12),
                                   text_color=self.COLORS['text_primary'],
                                   wraplength=700, justify="left")
        interp_value.pack(anchor="w", padx=25, pady=(0, 8))

        # Caso de uso
        use_case_label = ctk.CTkLabel(scroll_frame, text=TEXTS["model_info_use_case"],
                                     font=(_FONT_UI, 13, "bold"),
                                     text_color=self.COLORS['warning'])
        use_case_label.pack(anchor="w", padx=8, pady=(8, 2))

        use_case_value = ctk.CTkLabel(scroll_frame, text=model_info.get('use_case', ''),
                                     font=(_FONT_UI, 12),
                                     text_color=self.COLORS['text_primary'],
                                     wraplength=700, justify="left")
        use_case_value.pack(anchor="w", padx=25, pady=(0, 8))

        # Referências
        refs_label = ctk.CTkLabel(scroll_frame, text=TEXTS["model_info_references"],
                                 font=(_FONT_UI, 13, "bold"),
                                 text_color='#eab308')   # amarelo — não tem par no COLORS
        refs_label.pack(anchor="w", padx=8, pady=(8, 2))
        
        refs_value = ctk.CTkLabel(scroll_frame, text=model_info.get('references', ''),
                                 font=(_FONT_UI, 12, "italic"),
                                 text_color=self.COLORS['text_secondary'],
                                 wraplength=700, justify="left")
        refs_value.pack(anchor="w", padx=25, pady=(0, 8))

    # ── Idioma ───────────────────────────────────────────────────────────────

    @staticmethod
    def _lang_pref_path() -> Path:
        """Caminho do arquivo de preferência de idioma."""
        return Path.home() / '.easypml_lang'

    @staticmethod
    def load_language_pref() -> None:
        """Idioma inicial: o salvo pelo usuário; sem preferência salva, o do
        sistema (português -> PT; qualquer outro -> EN, o padrão)."""
        lang = None
        try:
            p = App._lang_pref_path()
            if p.exists():
                lang = p.read_text(encoding='utf-8').strip()
        except Exception:
            lang = None
        if lang not in ('pt', 'en'):
            lang = backend_messages.system_language()
        set_language(lang)
        backend_messages.set_language(lang)

    def _save_language_pref(self, lang: str) -> None:
        try:
            self._lang_pref_path().write_text(lang, encoding='utf-8')
        except Exception:
            pass

    def _build_lang_footer(self) -> None:
        """Rodapé fixo da sidebar com botões PT / EN."""
        footer = ctk.CTkFrame(self.sidebar,
                              fg_color=self.COLORS['bg_card'],
                              corner_radius=0,
                              border_width=1,
                              border_color=self.COLORS['border'])
        footer.pack(fill='x', side='bottom', padx=0, pady=0)

        inner = ctk.CTkFrame(footer, fg_color='transparent')
        inner.pack(fill='x', padx=10, pady=6)


        current = get_language()

        def _make_btn(code: str, label: str):
            active = (code == current)
            btn = ctk.CTkButton(
                inner, text=label, width=46, height=24,
                font=(_FONT_UI, 12, "bold" if active else "normal"),
                corner_radius=6,
                fg_color=self.COLORS['accent_blue'] if active else self.COLORS['border'],
                hover_color=self.COLORS['accent_blue_hover'],
                text_color='#ffffff',
                command=lambda c=code: self._switch_language(c)
            )
            btn.pack(side='left', padx=2)

        _make_btn('pt', 'PT')
        _make_btn('en', 'EN')
        ctk.CTkButton(inner, text=TEXTS["btn_about"], width=70, height=24,
                      font=(_FONT_UI, 12), corner_radius=6,
                      fg_color=self.COLORS['border'], hover_color=self.COLORS['border_hover'],
                      text_color='#ffffff', command=lambda: show_about(self)).pack(side='right', padx=2)

    def _switch_language(self, lang: str) -> None:
        """Salva preferência e reinicia o app para aplicar o idioma."""
        import subprocess
        current = get_language()
        if lang == current:
            return
        self._save_language_pref(lang)
        set_language(lang)

        restart = ask_yes_no(self, TEXTS['lang_restart_title'],
                             f"{TEXTS['lang_switch_message']}\n\n{TEXTS['lang_switch_confirm']}")
        if restart:
            try:
                subprocess.Popen([sys.executable] + sys.argv)
                self.destroy()
            except Exception as exc:
                show_message(self, "EasyPAML", TEXTS["lang_switch_err"].format(error=exc), 'error')

    def _open_results_viewer(self):
        if not self.output_folder:
            self.append_log(TEXTS["log_no_output_folder"])
            return
        try:
            ResultsViewerWindow(self, self.output_folder)
        except Exception as e:
            self.append_log(tr("Erro ao abrir o painel de resultados: ", "Error opening the results panel: ") + f"{e}", "error")
            self.append_log(traceback.format_exc())

    def _regenerate_summary_files(self):
        """Abre diálogo para selecionar pasta e regenera os 3 arquivos de síntese"""
        results_folder = filedialog.askdirectory(
            title=TEXTS["dialog_select_results_folder"],
            initialdir=str(Path.home() / "Desktop")
        )
        
        if not results_folder:
            return
        
        results_folder = Path(results_folder)
        
        self.append_log(TEXTS["log_updating_results"], "header")
        self.append_log(tr("Pasta de resultados: ", "Results folder: ") + f"{results_folder}", "info")
        
        # Executar em thread separada para não travar GUI
        def _update_thread():
            try:
                self.append_log(TEXTS["log_detecting_models"])

                # Descobrir quais modelos estão presentes
                models = set()
                for item in results_folder.iterdir():
                    if item.is_dir() and item.name not in ['reports']:
                        models.add(item.name)

                models = sorted(models)
                self.append_log(tr("Modelos encontrados: ", "Models found: ") + ", ".join(models), "info")

                # Determinar comparações disponíveis
                self.append_log(TEXTS["log_lrt_comparisons"])
                comparisons = []

                if 'M0' in models and 'M1a' in models:
                    comparisons.append("M0 vs M1a")
                    self.append_log(TEXTS["log_lrt_M0_M1a"])

                if 'M1a' in models and 'M2a' in models:
                    comparisons.append("M1a vs M2a")
                    self.append_log(TEXTS["log_lrt_M1a_M2a"])

                if 'M7' in models and 'M8' in models:
                    comparisons.append("M7 vs M8")
                    self.append_log(TEXTS["log_lrt_M7_M8"])

                if 'M0' in models and 'Branch' in models:
                    comparisons.append("M0 vs Branch")
                    self.append_log(TEXTS["log_lrt_M0_Branch"])

                if 'Branch-site_null' in models and 'Branch-site' in models:
                    comparisons.append("Branch-site_null vs Branch-site")
                    self.append_log(TEXTS["log_lrt_BranchSite"])

                self.append_log(TEXTS["log_lrt_total"].format(n=len(comparisons)))

                # Regenerar arquivos
                self.append_log(TEXTS["log_regenerating"])
                generated_files = CodemlBatchAnalysis.regenerate_summary_files(results_folder)
                
                if generated_files:
                    self.append_log(TEXTS["log_update_done"])
                    for file_type, file_path in generated_files.items():
                        filepath = Path(file_path)
                        size = filepath.stat().st_size if filepath.exists() else 0
                        self.append_log(f"  [OK] {file_type:25s} | {size:,} bytes", "ok")
                    # Atualiza pasta de saída para a pasta selecionada e habilita o botão
                    self.output_folder = results_folder
                    self.after(0, self._update_models_state)
                else:
                    self.append_log(tr("Nenhum arquivo foi gerado.", "No file was generated."), "error")
            
            except Exception as e:
                self.append_log(str(e), "error")
                self.append_log(traceback.format_exc(), "debug")
        
        update_thread = threading.Thread(target=_update_thread, daemon=True)
        update_thread.start()

    def _update_results_files(self):
        """Alias para _regenerate_summary_files com novo nome"""
        self._regenerate_summary_files()

    def _all_genes_have_trees(self) -> bool:
        if not self.input_folder or not self.per_gene_trees:
            return False
        chosen, _ = group_by_gene(list_alignment_files(self.input_folder))
        return bool(chosen) and all(g in self.per_gene_trees for g in chosen)

    def _update_models_state(self):
        # árvore: o arquivo escolhido OU uma árvore por gene (GENE.nwk) para todos
        has_tree = bool(self.tree_file) or self._all_genes_have_trees()
        enabled = all([self.input_folder, has_tree, self.output_folder])
        state = "normal" if enabled else "disabled"
        try:
            if enabled:
                self.files_hint.pack_forget()
            elif not self.files_hint.winfo_ismapped():
                self.files_hint.pack(fill='x', pady=(0, 6), before=self.tabs)
        except Exception:
            pass
        
        for btn in self.model_gear_buttons.values(): 
            btn.configure(state=state)
        for cb in self.model_checkboxes.values(): 
            cb.configure(state=state)
        
        try:
            branch_selected = self.model_vars.get('Branch', ctk.BooleanVar()).get()
            branchsite_selected = self.model_vars.get('Branch-site', ctk.BooleanVar()).get()
            
            self.btn_label_branch.configure(
                state="normal" if (enabled and branch_selected) else "disabled"
            )
            self.btn_label_branchsite.configure(
                state="normal" if (enabled and branchsite_selected) else "disabled"
            )
        except Exception:
            pass

        try:
            if self.output_folder and (self.output_folder / 'analysis_summary.tsv').exists():
                self.btn_results.configure(state='normal')
            else:
                self.btn_results.configure(state='disabled')
        except Exception:
            try:
                self.btn_results.configure(state='disabled')
            except Exception:
                pass

    _LEVEL_TAGS = {'ok': 'success', 'error': 'error', 'warn': 'warning', 'info': None,
                   'header': 'header', 'debug': 'debug'}

    def append_log(self, text: str, level: str = None):
        """Uma mensagem por linha. level: ok | error | warn | info | header | debug.
        Sem level, a cor é deduzida do texto (mensagens antigas). Mensagens
        'debug' só aparecem com "Mostrar detalhes técnicos"."""
        def _append():
            if level is not None:
                tag = self._LEVEL_TAGS.get(level)
            elif any(x in text for x in ("[OK]", "SUCESSO", "CONCLUÍD")):
                tag = "success"
            elif any(x in text for x in ("[Erro]", "ERRO", "FALHOU", "FAILED", "Error")):
                tag = "error"
            elif any(x in text for x in ("[!]", "AVISO", "WARNING")):
                tag = "warning"
            else:
                tag = None
            line = text if text.endswith("\n") else text + "\n"
            tags = (tag,) if tag else ()
            if level == 'debug':
                tags = ('debug',)
            self.log.insert("end", line, tags)
            self.log.see("end")
        self.after(0, _append)

    def _backend_log(self, level: str, text: str):
        """log_callback do backend (chamado de threads de trabalho)."""
        self.append_log(text, level)

    def _backend_progress(self, done: int, total: int, gene: str = ''):
        def _upd():
            frac = (done / total) if total else 0
            self.progress_bar.set(frac)
            self.progress_label.configure(text=TEXTS["progress_template"].format(
                done=done, total=total))
        self.after(0, _upd)

    def _set_pause_button(self, paused: bool):
        if paused:
            self.btn_pause.configure(text=TEXTS["btn_resume"], fg_color=PALETTE['success_fill'],
                                     hover_color=mix(PALETTE['success_fill'], '#000000', 0.2),
                                     text_color='#ffffff', border_color=PALETTE['success_fill'])
        else:
            self.btn_pause.configure(text=TEXTS["btn_pause"], fg_color=self.COLORS['bg_card'],
                                     hover_color=hover_tint(self.COLORS['warning'], self.COLORS['bg_card']),
                                     text_color=self.COLORS['warning'], border_color=self.COLORS['warning'])

    def _toggle_pause(self):
        if not self.pause_event:
            return
        if not self.analysis_thread or not self.analysis_thread.is_alive():
            return
        if self.pause_event.is_set():
            self.pause_event.clear()
            self._suspend_active_codeml()
            self._set_pause_button(True)
            self.status_indicator.configure(text=TEXTS["status_paused"],
                                            text_color=self.COLORS['warning'])
        else:
            self._resume_active_codeml()
            self.pause_event.set()
            self._set_pause_button(False)
            self.status_indicator.configure(text=TEXTS["status_running"],
                                            text_color=self.COLORS['success'])

    def _suspend_active_codeml(self):
        """Suspende os processos codeml ativos (psutil)."""
        try:
            import psutil
            procs = list(getattr(self.analysis_instance, '_active_processes', []) if self.analysis_instance else [])
            for proc in procs:
                try:
                    psutil.Process(proc.pid).suspend()
                except Exception:
                    pass
        except ImportError:
            pass  # psutil indisponível — pausa só entre genes

    def _resume_active_codeml(self):
        try:
            import psutil
            procs = list(getattr(self.analysis_instance, '_active_processes', []) if self.analysis_instance else [])
            for proc in procs:
                try:
                    psutil.Process(proc.pid).resume()
                except Exception:
                    pass
        except ImportError:
            pass

    def _stop_analysis(self):
        if not self.analysis_thread or not self.analysis_thread.is_alive():
            return
        # 1. backend para de iniciar genes/modelos; 2. encerra e recolhe os codeml
        self.stop_event.set()
        if self.pause_event:
            self.pause_event.set()
        self._resume_active_codeml()
        self.status_indicator.configure(text=TEXTS["status_stopped"], text_color=self.COLORS['danger'])
        if self.analysis_instance:
            try:
                self.analysis_instance.stop_all_processes()
            except Exception as e:
                self.append_log(f"{e}", 'debug')

    def _selected_models(self):
        selected = [k for k, v in self.model_vars.items() if v.get()]
        if self.include_neutral_models.get():
            selected = CodemlBatchAnalysis.auto_complete_null_models(
                selected, include_neutral=True, include_m8a=bool(self.include_m8a_var.get()))
        return selected

    def start_analysis(self):
        """Iniciar: valida os dados ANTES de rodar e mostra um diálogo com os
        problemas encontrados ("corrigir e voltar" / "continuar mesmo assim")."""
        if self.analysis_thread and self.analysis_thread.is_alive():
            return
        original = [k for k, v in self.model_vars.items() if v.get()]
        if not original:
            self.append_log(tr("Selecione pelo menos um modelo.", "Select at least one model."), 'warn')
            return
        needs_branchsite = any('Branch-site' in m for m in original)
        if needs_branchsite and not self.tree_branchsite_labeled:
            self.append_log(("Branch-site precisa de um ramo marcado: use "
                             f"'{TEXTS['btn_label_branchsite']}'.") if get_language() == 'pt' else
                            (f"Branch-site needs a labelled branch: use '{TEXTS['btn_label_branchsite']}'."),
                            'error')
            return
        if 'Branch' in original and not self.tree_branch_labeled:
            self.append_log(("O modelo Branch precisa de ramos marcados: use "
                             f"'{TEXTS['btn_label_branch']}'.") if get_language() == 'pt' else
                            (f"The Branch model needs labelled branches: use '{TEXTS['btn_label_branch']}'."),
                            'error')
            return
        selected = self._selected_models()
        added = sorted(set(selected) - set(original))
        if added:
            self.append_log(("Modelos nulos adicionados automaticamente: " if get_language() == 'pt'
                             else "Null models added automatically: ") + ", ".join(added), 'info')

        self.btn_run.configure(state="disabled")
        self.status_indicator.configure(text=TEXTS["preflight_running"], text_color=self.COLORS['info'])
        ignore = bool(self.ignore_stop_codons_var.get())
        prune = bool(self.auto_prune_tree_var.get())

        def _check():
            try:
                report = run_preflight(self.input_folder, self.tree_file, auto_prune=prune,
                                       ignore_stop_codons=ignore,
                                       per_gene_trees=self.per_gene_trees)
                err = None
            except Exception as exc:
                report, err = None, exc
            self.after(0, lambda: self._after_preflight(selected, report, err))
        threading.Thread(target=_check, daemon=True).start()

    def _after_preflight(self, selected, report, error):
        self.status_indicator.configure(text=TEXTS["status_ready"], text_color=self.COLORS['text_tertiary'])
        ignore_stops = bool(self.ignore_stop_codons_var.get())
        if error is not None:
            self.append_log(f"{error}", 'error')
        elif report is not None and report.has_problems:
            for line in report.format_text(get_language(), include_info=False).splitlines():
                self.append_log(line, 'warn')
            choice = PreflightDialog(self, report).show()
            if choice != 'continue':
                self.btn_run.configure(state="normal")
                return
            # "Continuar mesmo assim": o codeml trata stop codons como dado ausente
            if any(i.kind == 'stop_codon' for i in report.issues):
                ignore_stops = True
        self._launch(selected, ignore_stops)

    def _launch(self, selected, ignore_stops: bool):
        self.stop_event.clear()
        self.pause_event.set()
        self._set_pause_button(False)
        self.progress_bar.set(0)
        self.progress_label.configure(text=TEXTS["progress_template"].format(done=0, total='…'))
        self.status_indicator.configure(text=TEXTS["status_running"], text_color=self.COLORS['success'])
        self.analysis_thread = threading.Thread(target=self._run_thread, args=(selected, ignore_stops),
                                                daemon=True)
        self.analysis_thread.start()

    def _timeout_seconds(self) -> int:
        """Minutos digitados em Configurações -> segundos; vazio ou inválido = 0 (automático)."""
        try:
            minutes = float((self.entry_timeout.get() or '').replace(',', '.'))
        except ValueError:
            return 0
        return int(minutes * 60) if minutes > 0 else 0

    def _run_thread(self, selected, ignore_stops=False):
        old_stdout = sys.stdout
        sys.stdout = StdoutRedirect(self.append_log)
        summary = None
        try:
            analysis = CodemlBatchAnalysis()
            self.analysis_instance = analysis
            try:
                omega = float(self.entry_omega.get() or 0.5)
            except ValueError:
                omega = 0.5
            opts = self._global_ctl_options()
            analysis.config = {
                'input_folder': self.input_folder,
                'tree_file': self.tree_file,
                'per_gene_trees': dict(self.per_gene_trees),
                'output_folder': self.output_folder,
                'models': selected,
                'custom_model_params': self.custom_model_params,
                'labeled_tree_content': self.tree_branch_labeled,
                'labeled_tree_branchsite': self.tree_branchsite_labeled,
                'omega': omega,
                'CodonFreq': opts['CodonFreq'],
                'ncatG': opts['ncatG'],
                'cleandata': int(self.cleandata_var.get()),
                'timeout': self._timeout_seconds(),
                'idle_timeout': 300,
                'run_lrt': True,
                'n_workers': int(self.cores_var.get()),
                'ignore_stop_codons': ignore_stops,
                'auto_prune_tree': self.auto_prune_tree_var.get(),
                'pause_event': self.pause_event,
                'stop_event': self.stop_event,
                'log_callback': self._backend_log,
                'progress_callback': self._backend_progress,
                'interface': 'gui',
            }
            self.append_log(TEXTS["log_analysis_start"], 'header')
            summary = analysis.run_batch_analysis()
        except Exception as e:
            self.append_log(f"{e}", 'error')
            self.append_log(traceback.format_exc(), 'debug')
        finally:
            sys.stdout = old_stdout
            self.analysis_instance = None
            self.last_run_summary = summary
            self.after(0, lambda: self._on_run_finished(summary))

    def _on_run_finished(self, summary):
        self.btn_run.configure(state="normal")
        self.status_indicator.configure(text=TEXTS["status_ready"], text_color=self.COLORS['text_tertiary'])
        self._set_pause_button(False)
        self._update_models_state()
        if not summary:
            return
        total, done = summary.get('total', 0), summary.get('ok', 0) + summary.get('failed', 0)
        self.progress_label.configure(text=TEXTS["progress_done"].format(done=done, total=total))
        if summary.get('failed'):
            items = "\n".join(f"• {g}: {r}" for g, r in sorted(summary['failures'].items())[:30])
            if len(summary['failures']) > 30:
                items += f"\n… (+{len(summary['failures']) - 30})"
            show_message(self, TEXTS["failures_title"], TEXTS["failures_text"].format(
                failed=summary['failed'], total=total, items=items,
                path=Path(summary['output_folder']) / 'genes_status.tsv'), 'error')
        # Abre os resultados ao terminar (como o README promete), se algum gene rodou
        if summary.get('ok') and not summary.get('stopped'):
            self._open_results_viewer()

    def _poll_stop_count(self):
        if self.analysis_instance:
            cnt = getattr(self.analysis_instance, 'current_stop_count', 0)
            self.stop_label.configure(text=TEXTS["status_stops_template"].format(n=cnt))
        self.after(1000, self._poll_stop_count)


if __name__ == "__main__":
    app = App()
    app.mainloop()
