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
from .results_viewer import ResultsViewerWindow
from .gui_texts import TEXTS

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
    def __init__(self, append_func):
        self.append = append_func
    def write(self, s):
        if s and s.strip() != "": self.append(s)
    def flush(self): pass

class ModelConfigWindow(ctk.CTkToplevel):
    """Janela para editar os parâmetros do .ctl em memória - Estilo Premium"""
    
    COLORS = {
        'bg_dark':        '#0c0c0e',
        'bg_card':        '#16161a',
        'text_primary':   '#ededef',
        'text_secondary': '#9898a6',
        'accent_blue':    '#6366f1',
        'accent_blue_hover': '#4f46e5',
        'success':        '#22c55e',
        'success_hover':  '#16a34a',
        'danger':         '#f87171',
    }
    
    def __init__(self, parent, model_code: str, initial: dict):
        super().__init__(parent)
        self.title(f"Configurar {model_code}")
        self.geometry("500x500")
        self.parent = parent
        self.model_code = model_code
        self.entries = {}
        
        self.configure(fg_color=self.COLORS['bg_dark'])
        self.attributes("-topmost", True)
        self.grab_set()

        # ═══ HEADER ═══
        header = ctk.CTkFrame(self, fg_color='transparent')
        header.pack(fill="x", padx=20, pady=(20, 10))
        
        ctk.CTkLabel(header, text=TEXTS["model_config_header"].format(model_code=model_code),
                    font=(_FONT_UI, 16, "bold"),
                    text_color=self.COLORS['text_primary']).pack(side="left")
        
        # ═══ FORM FRAME ═══
        form_frame = ctk.CTkScrollableFrame(self, fg_color=self.COLORS['bg_card'],
                                           corner_radius=10, border_width=1,
                                           border_color='#2a2a2a')
        form_frame.pack(fill="both", expand=True, padx=20, pady=(0, 15))

        fields = [
            ('NSsites', 'NSsites', '0=M0  1=M1a  2=M2a  7=M7  8=M8'),
            ('model',   'model (CODEML)', '0=NSsites (sítios)  2=Branch (ramos)'),
            ('fix_omega', 'fix_omega', '0=estima ω livremente  1=fixa ω'),
            ('omega',   'omega inicial', 'Valor inicial de dN/dS (ex: 0.5)'),
            ('CodonFreq', 'CodonFreq', '7=F3×4 (recomendado)'),
            ('description', 'Descrição', 'Texto descritivo'),
        ]

        for name, label, placeholder in fields:
            # Label com ícone
            lbl_frame = ctk.CTkFrame(form_frame, fg_color='transparent')
            lbl_frame.pack(fill="x", padx=12, pady=(12, 4))
            
            ctk.CTkLabel(lbl_frame, text=f"{label}:", font=(_FONT_UI, 11, "bold"),
                        text_color=self.COLORS['text_primary']).pack(anchor="w")

            ctk.CTkLabel(lbl_frame, text=placeholder, font=(_FONT_UI, 9),
                        text_color='#999999').pack(anchor="w", pady=(0, 4))
            
            # Entry com estilo premium
            ent = ctk.CTkEntry(form_frame, placeholder_text=placeholder,
                              fg_color=self.COLORS['bg_dark'],
                              border_color='#32323e',
                              border_width=1,
                              corner_radius=8,
                              text_color=self.COLORS['text_primary'],
                              placeholder_text_color='#666666')
            ent.pack(fill="x", padx=12, pady=(0, 8))
            
            val = self._initial_val(name, initial)
            ent.insert(0, str(val))
            self.entries[name] = ent

        # ═══ BUTTONS ═══
        btn_frame = ctk.CTkFrame(self, fg_color='transparent')
        btn_frame.pack(fill="x", padx=20, pady=(0, 20))
        
        btn_cancel = ctk.CTkButton(btn_frame, text=TEXTS["model_config_btn_cancel"],
                                   fg_color='#3d3d3d',
                                   hover_color='#4d4d4d', command=self.destroy,
                                   font=(_FONT_UI, 11, "bold"), corner_radius=6)
        btn_cancel.pack(side="left", fill="x", expand=True, padx=(0, 8))

        btn_save = ctk.CTkButton(btn_frame, text=TEXTS["model_config_btn_save"],
                                 fg_color=self.COLORS['success'],
                                 hover_color=self.COLORS['success_hover'],
                                 command=self._on_save, font=(_FONT_UI, 11, "bold"), corner_radius=6)
        btn_save.pack(side="left", fill="x", expand=True)

    def _initial_val(self, name, initial):
        if name in self.parent.custom_model_params.get(self.model_code, {}):
            return self.parent.custom_model_params[self.model_code][name]
        return initial.get(name, "")

    def _on_save(self):
        out = {}
        for k, ent in self.entries.items():
            val = ent.get().strip()
            try:
                if '.' in val:
                    out[k] = float(val)
                else:
                    out[k] = int(val)
            except Exception:
                out[k] = val  # valor mantido como string se não for numérico
        
        self.parent.custom_model_params[self.model_code] = out
        self.parent.model_ctl_labels[self.model_code].configure(
            text=TEXTS["model_status_configured"], text_color="#fbbf24")
        self.parent.append_log(f"[OK] Parâmetros de {self.model_code} salvos com sucesso!\n")
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
        self.title(f"Etiquetar Arvore - {mode.upper()}")
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
            ctk.CTkLabel(self, text="[!] Biopython não instalado. Execute: pip install biopython",
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
                                 font=(_FONT_UI, 10), text_color=self.TEXT_SECONDARY)
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
            ctk.CTkLabel(plot_frame, text="[!] Nenhuma arvore selecionada.", font=("Arial", 14)).pack(pady=50)
            return

        try:
            self.tree = Phylo.read(str(self.tree_path), 'newick')
        except Exception as e:
            ctk.CTkLabel(plot_frame, text=f"[Erro] Erro ao carregar arvore:\n{e}", font=("Arial", 12)).pack(pady=50)
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
                child_length = child.branch_length if child.branch_length else 1.0
                calc_depth_with_lengths(child, accumulated_depth + child_length)
        
        calc_depth_with_lengths(self.tree.root, 0.0)
        
        for clade in self.tree.find_clades(order='postorder'):
            x = depths.get(clade, 0.0)
            
            if clade.is_terminal():
                y = terminal_y_map[clade]
            else:
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
            
            self.ax.text(x + 0.002, y, name, 
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
                self.parent.append_log(f"[-] Tag #1 removida de {self._get_clade_name(clade)}\n")
            else:
                self._remove_tag_recursively(clade)
                self._apply_tag_recursively(clade, '#1')
                self.marked_clades[clade] = '#1'
                self.parent.append_log(f"[+] Tag #1 aplicada a {self._get_clade_name(clade)}\n")
        
        else:
            current_tag = self.marked_clades.get(clade)
            clade_name = self._get_clade_name(clade)
            
            if current_tag:
                response = simpledialog.askstring(
                    "Editar Tag",
                    f"Ramo atual: {current_tag}\n\nDigite novo número ou 'remover':",
                    parent=self
                )
                
                if response:
                    response = response.strip().lower()
                    if response == 'remover':
                        self._remove_tag_recursively(clade)
                        self.parent.append_log(f"[-] Tag {current_tag} removida de {clade_name}\n")
                    elif response.isdigit():
                        self._remove_tag_recursively(clade)
                        new_tag = f"#{response}"
                        self._apply_tag_recursively(clade, new_tag)
                        self.marked_clades[clade] = new_tag
                        self.parent.append_log(f"[edit] Tag alterada para {new_tag} em {clade_name}\n")
            else:
                response = simpledialog.askstring(
                    "Número da Tag",
                    f"Digite o número da tag:\n(ex: 1 para #1, 2 para #2)",
                    parent=self
                )
                
                if response and response.strip().isdigit():
                    new_tag = f"#{response.strip()}"
                    self._remove_tag_recursively(clade)
                    self._apply_tag_recursively(clade, new_tag)
                    self.marked_clades[clade] = new_tag
                    self.parent.append_log(f"[tag] Tag {new_tag} aplicada a {clade_name}\n")

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
        self.parent.append_log(f"Tag {tag} removida.\n")

    def _refresh_legend(self):
        """Atualiza legenda com tags ativas e botão de exclusão por tag."""
        for widget in self.legend_frame.winfo_children():
            widget.destroy()

        active_tags = sorted(set(self.clade_tags.values()))

        if not active_tags:
            ctk.CTkLabel(self.legend_frame, text=TEXTS["tree_labeler_no_tags"],
                         text_color="#888888", font=(_FONT_UI, 11, "italic")
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
                             font=(_FONT_UI, 10),
                             text_color="#888888").pack(side='left')

                del_btn = ctk.CTkButton(
                    row_frame, text="X", width=22, height=22,
                    fg_color='transparent',
                    hover_color=self.DANGER,
                    text_color='#888888',
                    border_width=0,
                    corner_radius=4,
                    font=(_FONT_UI, 9, "bold"),
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
                self.parent.append_log(f"[OK] Arvore Branch-Site salva ({applied_terminals} terminais).\n")
            else:
                self.parent.tree_branch_labeled = newick_str
                self.parent.append_log(f"[OK] Arvore Branch salva ({applied_terminals} terminais).\n")

                branchsite_version = re.sub(r"\s*#(?!1)\d+\b", "", newick_str)
                branchsite_version = re.sub(r"\s+", " ", branchsite_version).strip()
                self.parent.tree_branchsite_labeled = branchsite_version

        except Exception as e:
            self.parent.append_log(f"[Erro] Erro ao gerar Newick: {e}\n")

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
        'text_tertiary':  '#686880',
        'text_muted':     '#424252',

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
        self.title("EasyPAML")
        self.geometry("1400x850")
        
        # Aplicar tema escuro profissional
        ctk.set_appearance_mode("dark")
        ctk.set_default_color_theme("blue")
        
        # Configure janela principal
        self.configure(fg_color=self.COLORS['bg_dark'])

        # Inicializar backend
        self.codeml_backend = CodemlBatchAnalysis()

        self.custom_model_params = {} 
        self.input_folder = None
        self.tree_file = None
        self.output_folder = None
        self.analysis_thread = None
        self.analysis_instance = None
        self.pause_event = threading.Event()
        self.manual_event = threading.Event()
        self.manual_all_event = threading.Event()
        self.stop_event = threading.Event()
        
        # Opção para incluir modelos neutros automaticamente
        self.include_neutral_models = ctk.BooleanVar(value=True)

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
                     font=(_FONT_UI, 9),
                     text_color=self.COLORS['text_tertiary']).pack(anchor='w')

        ctk.CTkFrame(logo_frame, fg_color=self.COLORS['border'],
                     height=1, corner_radius=0).pack(fill='x', pady=(14, 0))

        # ── Scrollable section container ─────────────────────────────
        _sb = ctk.CTkScrollableFrame(self.sidebar, fg_color='transparent',
                                     scrollbar_button_color=self.COLORS['border'],
                                     scrollbar_button_hover_color=self.COLORS['border_hover'])
        _sb.pack(fill='both', expand=True, padx=0, pady=(4, 4))

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
                         font=(_FONT_UI, 10, "bold"),
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
                               font=(_FONT_UI, 11, "bold"), height=36)
        self.btn_input.pack(fill='x', pady=(0, 2))
        self.label_input = ctk.CTkLabel(fi, text=TEXTS["label_not_selected"],
                                        font=(_FONT_UI, 9),
                                        wraplength=230,
                                        text_color=self.COLORS['text_tertiary'])
        self.label_input.pack(anchor='w', padx=4, pady=(0, 8))

        self.btn_tree = _obtn(fi, TEXTS["btn_tree_file"],
                              self.select_tree_file,
                              self.COLORS['accent_blue'],
                              font=(_FONT_UI, 11, "bold"), height=36)
        self.btn_tree.pack(fill='x', pady=(0, 2))
        self.label_tree = ctk.CTkLabel(fi, text=TEXTS["label_not_selected"],
                                       font=(_FONT_UI, 9),
                                       wraplength=230,
                                       text_color=self.COLORS['text_tertiary'])
        self.label_tree.pack(anchor='w', padx=4, pady=(0, 8))

        self.btn_output = _obtn(fi, TEXTS["btn_output_folder"],
                                self.select_output_folder,
                                self.COLORS['accent_blue'],
                                font=(_FONT_UI, 11, "bold"), height=36)
        self.btn_output.pack(fill='x', pady=(0, 2))
        self.label_output = ctk.CTkLabel(fi, text=TEXTS["label_not_selected"],
                                         font=(_FONT_UI, 9),
                                         wraplength=230,
                                         text_color=self.COLORS['text_tertiary'])
        self.label_output.pack(anchor='w', padx=4)

        # ── Resultados ───────────────────────────────────────────────
        ri = _sec(_sb, TEXTS["section_results"])

        self.btn_results = _obtn(ri, TEXTS["btn_view_results"],
                                  self._open_results_viewer,
                                  self.COLORS['accent_blue'],
                                  font=(_FONT_UI, 11, "bold"), height=36)
        self.btn_results.pack(fill='x', pady=(0, 6))
        self.btn_results.configure(state="disabled")

        self.btn_update_results = _obtn(ri, TEXTS["btn_update_results"],
                                         self._update_results_files,
                                         self.COLORS['success'],
                                         font=(_FONT_UI, 11, "bold"), height=36)
        self.btn_update_results.pack(fill='x', pady=(0, 2))
        self.label_update_results = ctk.CTkLabel(ri,
                                                  text=TEXTS["label_update_results_hint"],
                                                  font=(_FONT_UI, 9, "italic"),
                                                  wraplength=230,
                                                  text_color=self.COLORS['text_muted'])
        self.label_update_results.pack(anchor='w', padx=4)

        # ── Configurações ────────────────────────────────────────────
        ci = _sec(_sb, TEXTS["section_config"])

        ctk.CTkLabel(ci, text=TEXTS["label_omega_initial"],
                     font=(_FONT_UI, 9, "bold"), anchor='w',
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

        # Remover gaps toggle
        self.cleandata_var = ctk.BooleanVar(value=True)
        row_g = ctk.CTkFrame(ci, fg_color='transparent')
        row_g.pack(fill='x', pady=(0, 8))
        ctk.CTkLabel(row_g, text=TEXTS["label_remove_gaps"],
                     font=(_FONT_UI, 10),
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
                      font=(_FONT_UI, 9, "bold"), fg_color=self.COLORS['border'],
                      hover_color=self.COLORS['border_hover'],
                      text_color=self.COLORS['text_muted'],
                      command=lambda: self._show_help(
                          TEXTS["label_remove_gaps"], TEXTS["label_remove_gaps_hint"])
                      ).pack(side='right', padx=(0, 4))

        # CPU slider
        ctk.CTkLabel(ci, text=TEXTS["label_cpus"],
                     font=(_FONT_UI, 9, "bold"), anchor='w',
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
        self.cores_disp = ctk.CTkLabel(
            cores_row, text=f"{_max}×",
            font=(_FONT_UI, 10, "bold"),
            text_color=self.COLORS['accent_blue'],
            width=32)
        self.cores_disp.pack(side='right', padx=(6, 0))
        ctk.CTkLabel(ci, text=f"(detectado: {_max} núcleos)",
                     font=(_FONT_UI, 9),
                     text_color=self.COLORS['text_muted']).pack(anchor='w')

        # ── Modo Heurístico (fix_kappa do M0) ────────────────────────────────
        # fix_kappa=1: fixa κ (ts/tv) no valor estimado pelo M0, em vez de
        # re-estimá-lo livremente em cada modelo.  Economiza ~20-30 % de
        # iterações por modelo complexo (M1a, M2a, M7, M8, Branch-site).
        # NOTA: é uma aproximação — κ real de M1a/M2a pode diferir levemente do M0.
        # O LRT permanece válido desde que AMBOS os modelos do par usem o mesmo κ.
        self.heuristic_mode_var = ctk.BooleanVar(value=False)
        row_h = ctk.CTkFrame(ci, fg_color='transparent')
        row_h.pack(fill='x', pady=(10, 0))
        ctk.CTkLabel(row_h, text=TEXTS["label_heuristic_mode"],
                     font=(_FONT_UI, 10),
                     text_color=self.COLORS['text_secondary']).pack(side='left')
        self.cb_heuristic = ctk.CTkSwitch(
            row_h, text="",
            variable=self.heuristic_mode_var,
            onvalue=True, offvalue=False,
            switch_width=36, switch_height=18,
            progress_color='#f59e0b',   # âmbar — indica modo de aproximação
            button_color='#fffbeb',
            button_hover_color='#fef3c7',
            fg_color=self.COLORS['border'])
        self.cb_heuristic.pack(side='right')
        ctk.CTkButton(row_h, text="?", width=18, height=18, corner_radius=9,
                      font=(_FONT_UI, 9, "bold"), fg_color=self.COLORS['border'],
                      hover_color=self.COLORS['border_hover'],
                      text_color=self.COLORS['text_muted'],
                      command=lambda: self._show_help(
                          TEXTS["label_heuristic_mode"], TEXTS["label_heuristic_hint"])
                      ).pack(side='right', padx=(0, 4))

        # Ignorar Stop Codons toggle
        row_stops = ctk.CTkFrame(ci, fg_color='transparent')
        row_stops.pack(fill='x', pady=(10, 0))
        ctk.CTkLabel(row_stops, text=TEXTS["label_ignore_stops"],
                     font=(_FONT_UI, 10),
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
                      font=(_FONT_UI, 9, "bold"), fg_color=self.COLORS['border'],
                      hover_color=self.COLORS['border_hover'],
                      text_color=self.COLORS['text_muted'],
                      command=lambda: self._show_help(
                          TEXTS["label_ignore_stops"], TEXTS["label_ignore_stops_hint"])
                      ).pack(side='right', padx=(0, 4))

        # Poda automática de árvore toggle
        row_prune = ctk.CTkFrame(ci, fg_color='transparent')
        row_prune.pack(fill='x', pady=(10, 0))
        ctk.CTkLabel(row_prune, text=TEXTS["label_auto_prune"],
                     font=(_FONT_UI, 10),
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
                      font=(_FONT_UI, 9, "bold"), fg_color=self.COLORS['border'],
                      hover_color=self.COLORS['border_hover'],
                      text_color=self.COLORS['text_muted'],
                      command=lambda: self._show_help(
                          TEXTS["label_auto_prune"], TEXTS["label_auto_prune_hint"])
                      ).pack(side='right', padx=(0, 4))

        self.main_frame = ctk.CTkFrame(self, fg_color=self.COLORS['bg_dark'])
        self.main_frame.pack(side="right", fill="both", expand=True, padx=20, pady=20)

        self.tabs = ctk.CTkTabview(self.main_frame, fg_color=self.COLORS['bg_card'],
                                   segmented_button_fg_color=self.COLORS['bg_card'],
                                   segmented_button_selected_color=self.COLORS['accent_blue'],
                                   text_color=self.COLORS['text_primary'],
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
            font=(_FONT_UI, 11, "bold"),
            text_color=self.COLORS['text_tertiary']
        )
        self.status_indicator.pack(side="left", padx=(14, 20), pady=8)

        self.stop_label = ctk.CTkLabel(
            status_bar, text=TEXTS["status_stops_template"].format(n=0),
            font=(_FONT_UI, 11, "bold"),
            text_color=self.COLORS['danger']
        )
        self.stop_label.pack(side="left")

        neutral_row = ctk.CTkFrame(status_bar, fg_color='transparent')
        neutral_row.pack(side="right", padx=(0, 6))
        ctk.CTkLabel(neutral_row, text=TEXTS["label_neutral_models"],
                     font=(_FONT_UI, 9),
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
            text_color=self.COLORS['accent_blue'],
            border_width=1, border_color=self.COLORS['accent_blue'],
            command=self._show_neutral_models_info,
            corner_radius=6
        )
        help_btn.pack(side="right", padx=(0, 12))

        # Botões de controle
        btn_frame = ctk.CTkFrame(self.ctrl_frame, fg_color='transparent')
        btn_frame.pack(fill="x", padx=0, pady=(0, 10))

        def _action_btn(parent, text, cmd, color):
            return ctk.CTkButton(
                parent, text=text, command=cmd,
                fg_color=self.COLORS['bg_card'],
                hover_color=color,
                text_color=color,
                border_width=1, border_color=color,
                font=(_FONT_UI, 12, "bold"),
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

        self.btn_resolve_stops = _action_btn(btn_frame, TEXTS["btn_resolve_stops"],
                                              self._continue_all_for_gene, self.COLORS['info'])
        self.btn_resolve_stops.pack(side="left", fill="both", expand=True, padx=(4, 12), pady=10)

        # ═══ LOG FRAME COM HEADER ═══
        log_container = ctk.CTkFrame(self.main_frame, fg_color='transparent')
        log_container.pack(fill="both", expand=True, padx=0, pady=0)
        
        log_header = ctk.CTkFrame(log_container, fg_color='transparent')
        log_header.pack(fill="x", padx=0, pady=(10, 0))
        
        ctk.CTkLabel(log_header, text=TEXTS["log_header_title"], font=(_FONT_UI, 11, "bold"),
                    text_color=self.COLORS['text_secondary']).pack(side="left", padx=0)

        ctk.CTkLabel(log_header, text="·",
                    font=(_FONT_UI, 13), text_color=self.COLORS['text_muted']).pack(side="left", padx=8)

        ctk.CTkLabel(log_header, text=TEXTS["log_header_subtitle"],
                    font=(_FONT_UI, 9), text_color=self.COLORS['text_muted']).pack(side="left", padx=0)

        self.log = ctk.CTkTextbox(log_container, font=(_FONT_MONO, 11),
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
        self.log.tag_config("header", foreground=self.COLORS['accent_blue'])

        # Welcome message — editar em gui_texts.py › "log_welcome"
        self.log.insert("end", TEXTS["log_welcome"])

        self._update_models_state()
        self._poll_stop_count()

    def _show_help(self, title: str, body: str) -> None:
        win = ctk.CTkToplevel(self)
        win.title(title)
        win.resizable(False, False)
        win.grab_set()
        win.focus_set()
        ctk.CTkLabel(win, text=title, font=(_FONT_UI, 11, "bold"),
                     text_color=self.COLORS['text_primary']).pack(padx=18, pady=(16, 4), anchor='w')
        ctk.CTkLabel(win, text=body, font=(_FONT_UI, 10), wraplength=290,
                     justify='left',
                     text_color=self.COLORS['text_secondary']).pack(padx=18, pady=(0, 12))
        ctk.CTkButton(win, text="OK", width=80, command=win.destroy,
                      fg_color=self.COLORS['accent_blue']).pack(pady=(0, 14))

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
            TEXTS["tab_site_models"]:  ['M0', 'M1a', 'M2a', 'M7', 'M8'],
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
                                    height=58)
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
                    btn_col, text="?", width=20, height=20,
                    fg_color=self.COLORS['bg_card_hover'],
                    hover_color=accent,
                    text_color=accent,
                    corner_radius=4,
                    border_width=1, border_color=self.COLORS['border_hover'],
                    font=(_FONT_UI, 9, "bold"),
                    command=lambda c=code: self._show_model_info(c)
                )
                info_btn.pack(pady=(0, 2))

                gear = ctk.CTkButton(
                    btn_col, text="cfg", width=20, height=20,
                    fg_color=self.COLORS['bg_card'],
                    hover_color=self.COLORS['accent_purple'],
                    text_color=self.COLORS['text_secondary'],
                    corner_radius=4,
                    border_width=1, border_color=self.COLORS['border_hover'],
                    font=(_FONT_UI, 8),
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

                lbl = ctk.CTkLabel(content, text=TEXTS["model_status_default"],
                                   font=(_FONT_UI, 9),
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
            font=(_FONT_UI, 11, "bold"),
            text_color=self.COLORS['accent_blue'],
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
            font=(_FONT_UI, 11, "bold"),
            text_color=self.COLORS['accent_blue'],
            border_width=1, border_color=self.COLORS['border_hover'],
            corner_radius=8
        )
        self.btn_label_branchsite.pack(fill='x', padx=12, pady=(16, 12))

    def _open_config_window(self, code):
        default_config = CodemlBatchAnalysis.MODEL_CONFIGS.get(code, {})
        ModelConfigWindow(self, code, default_config)

    def _open_tree_labeler(self, mode: str = 'branchsite'):
        if self.tree_file is None:
            self.append_log("[!] Selecione uma arvore (.nwk) primeiro.\n")
            return

        try:
            TreeLabelWindow(self, self.tree_file, mode=mode)
        except Exception as e:
            self.append_log(f"[Erro] Erro ao abrir TreeLabelWindow: {e}\n")
            self.append_log(f"{traceback.format_exc()}\n")

    def select_input_folder(self):
        path = filedialog.askdirectory()
        if path:
            self.input_folder = Path(path)
            self.label_input.configure(text=str(self.input_folder.name))
            self._update_models_state()

    def select_tree_file(self):
        path = filedialog.askopenfilename(filetypes=[('Tree files', '*.tree *.tre *.nwk *.txt'), ('All files', '*.*')])
        if path:
            self.tree_file = Path(path)
            self.label_tree.configure(text=str(self.tree_file.name))
            self._update_models_state()

    def select_output_folder(self):
        path = filedialog.askdirectory()
        if path:
            self.output_folder = Path(path)
            self.label_output.configure(text=str(self.output_folder.name))
            self._update_models_state()

    def _show_neutral_models_info(self):
        """Mostra informações sobre todos os pares de modelos nulos/alternativos para o LRT"""
        info_window = ctk.CTkToplevel(self)
        info_window.title("Modelos Nulos – Comparações LRT")
        info_window.geometry("740x640")
        info_window.attributes("-topmost", True)
        info_window.grab_set()
        info_window.configure(fg_color=self.COLORS['bg_dark'])

        # ── Header ──────────────────────────────────────────────────────────
        hdr = ctk.CTkFrame(info_window, fg_color='transparent')
        hdr.pack(fill='x', padx=18, pady=(16, 4))
        ctk.CTkLabel(hdr, text="Modelos Nulos e Comparacoes LRT",
                     font=(_FONT_UI, 14, "bold"),
                     text_color=self.COLORS['accent_blue']).pack(anchor='w')
        ctk.CTkLabel(hdr,
                     text="Quando ativado, o EasyPAML adiciona automaticamente o modelo nulo de cada par LRT "
                          "— sem precisar marcá-lo manualmente. Cada comparação usa o Teste da Razão de "
                          "Verossimilhança (LRT): 2ΔlnL comparado ao χ² com os graus de liberdade corretos.",
                     font=(_FONT_UI, 10),
                     wraplength=690, justify='left',
                     text_color=self.COLORS['text_secondary']).pack(anchor='w', pady=(4, 0))

        ctk.CTkFrame(info_window, fg_color=self.COLORS['border'], height=1
                     ).pack(fill='x', padx=18, pady=(10, 0))

        # ── Scrollable content ───────────────────────────────────────────────
        scroll = ctk.CTkScrollableFrame(info_window, fg_color='transparent')
        scroll.pack(fill='both', expand=True, padx=14, pady=6)

        LRT_PAIRS = [
            {
                'null': 'M1a', 'alt': 'M2a', 'color': '#10b981',
                'title': 'M2a  vs  M1a',
                'test': 'Testa seleção positiva por sítio — Site Models',
                'detail': (
                    '2 graus de liberdade (χ²). M1a permite apenas purificação (ω < 1) e '
                    'neutralidade (ω = 1). M2a acrescenta uma classe com ω > 1 (seleção '
                    'positiva). Se 2ΔlnL > 5.99 (α=0.05), há evidência de seleção positiva; '
                    'os sítios são identificados por BEB/NEB.'
                ),
                'note': 'M1a não precisa de fix_omega=1 — o CODEML restringe ω₁=1 internamente via NSsites=1.',
            },
            {
                'null': 'M7', 'alt': 'M8', 'color': '#ec4899',
                'title': 'M8  vs  M7',
                'test': 'Testa seleção positiva com distribuição Beta — teste mais robusto',
                'detail': (
                    '2 graus de liberdade (χ²). M7 restringe toda a distribuição de ω ao '
                    'intervalo (0, 1) via Beta(p, q). M8 acrescenta uma classe discreta com '
                    'ω > 1. Este par é o mais recomendado para identificar seleção positiva '
                    'por sítio, pois a distribuição Beta é mais biológica que classes discretas.'
                ),
                'note': 'M7 não precisa de fix_omega=1.',
            },
            {
                'null': 'M0', 'alt': 'Branch', 'color': '#f59e0b',
                'title': 'Branch  vs  M0',
                'test': 'Testa variação de ω entre linhagens — Branch Model',
                'detail': (
                    'Os graus de liberdade dependem do número de ramos etiquetados. '
                    'M0 usa um único ω global para todos os ramos e sítios. O Branch Model '
                    'estima ω independente para cada ramo ou grupo marcado. Rejeitar M0 indica '
                    'que a pressão seletiva varia entre as linhagens analisadas.'
                ),
                'note': 'M0 não precisa de fix_omega=1 — ω é estimado livremente como baseline.',
            },
            {
                'null': 'Branch-site_null', 'alt': 'Branch-site', 'color': '#ef4444',
                'title': 'Branch-site  vs  Branch-site_null',
                'test': 'Testa seleção episódica em sítios do ramo foreground (#1)',
                'detail': (
                    'Distribuição MISTA 50:50 (χ²₀ + χ²₁) — não o χ² convencional! '
                    'Valor crítico: 2.706 (α=0.05) e 5.412 (α=0.01). '
                    'Branch-site_null fixa ω₂=1 no foreground (Yang et al. 2005, Zhang et al. 2005). '
                    'Este é o teste mais poderoso para detectar seleção positiva episódica em '
                    'linhagens específicas, combinando variação por sítio e por ramo.'
                ),
                'note': 'Branch-site_null requer fix_omega=1 e omega=1.0 no .ctl — configurado automaticamente.',
            },
        ]

        for info in LRT_PAIRS:
            card = ctk.CTkFrame(scroll, fg_color=self.COLORS['bg_card'],
                                corner_radius=10, border_width=1,
                                border_color=self.COLORS['border'])
            card.pack(fill='x', padx=6, pady=5)

            # Color bar
            ctk.CTkFrame(card, fg_color=info['color'], width=5,
                         corner_radius=2).pack(side='left', fill='y', padx=(6, 10), pady=8)

            content = ctk.CTkFrame(card, fg_color='transparent')
            content.pack(side='left', fill='both', expand=True, pady=10, padx=(0, 10))

            # Title row
            title_row = ctk.CTkFrame(content, fg_color='transparent')
            title_row.pack(fill='x', anchor='w')
            ctk.CTkLabel(title_row, text=info['title'],
                         font=(_FONT_UI, 12, "bold"),
                         text_color=info['color']).pack(side='left')
            ctk.CTkLabel(title_row,
                         text=f"   null: {info['null']}  →  alternativo: {info['alt']}",
                         font=(_FONT_UI, 10),
                         text_color=self.COLORS['text_secondary']).pack(side='left')

            # Test label
            ctk.CTkLabel(content, text=f"  {info['test']}",
                         font=(_FONT_UI, 10, "bold"),
                         text_color=self.COLORS['text_primary'],
                         anchor='w').pack(anchor='w', pady=(4, 2))

            # Detail
            ctk.CTkLabel(content, text=info['detail'],
                         font=(_FONT_UI, 10),
                         text_color=self.COLORS['text_secondary'],
                         wraplength=590, justify='left',
                         anchor='w').pack(anchor='w', pady=(0, 3))

            # Implementation note
            ctk.CTkLabel(content, text=f"  Nota: {info['note']}",
                         font=(_FONT_UI, 9, "italic"),
                         text_color=self.COLORS['text_tertiary'],
                         wraplength=590, justify='left',
                         anchor='w').pack(anchor='w')

        # ── Footer ───────────────────────────────────────────────────────────
        ctk.CTkFrame(info_window, fg_color=self.COLORS['border'], height=1
                     ).pack(fill='x', padx=18)
        ctk.CTkLabel(info_window,
                     text="[OK]  Com a opcao ATIVADA, o modelo nulo de cada par selecionado e "
                          "adicionado automaticamente — voce nao precisa marca-lo.",
                     font=(_FONT_UI, 10),
                     text_color=self.COLORS['success'],
                     wraplength=690).pack(padx=18, pady=12)

    def _show_model_info(self, model_code: str):
        """Mostra informações detalhadas sobre um modelo específico"""
        # Obter informações do modelo
        model_info = self.codeml_backend.MODEL_INFO.get(model_code)
        if not model_info:
            return
        
        # Criar janela de informações
        info_window = ctk.CTkToplevel(self)
        info_window.title(f"Modelo: {model_code}")
        info_window.geometry("750x600")
        info_window.attributes("-topmost", True)
        info_window.grab_set()
        
        # Header com nome do modelo
        display_name = self.codeml_backend.MODEL_CONFIGS[model_code].get('display_name', model_code)
        header = ctk.CTkLabel(info_window, 
                             text=f"{model_code} - {model_info.get('full_name', '')}",
                             font=(_FONT_UI, 13, "bold"),
                             text_color=self.COLORS['accent_blue'])
        header.pack(padx=15, pady=15)
        
        # Scrollable frame para conteúdo
        scroll_frame = ctk.CTkScrollableFrame(info_window, fg_color=self.COLORS['bg_feed'],
                                             corner_radius=8)
        scroll_frame.pack(fill="both", expand=True, padx=12, pady=(0, 12))
        
        # Tipo de teste
        test_type_label = ctk.CTkLabel(scroll_frame, text=TEXTS["model_info_test_type"],
                                      font=(_FONT_UI, 11, "bold"),
                                      text_color=self.COLORS['accent_purple'])
        test_type_label.pack(anchor="w", padx=8, pady=(8, 2))
        
        test_type_value = ctk.CTkLabel(scroll_frame, text=model_info.get('test_type', ''),
                                      font=(_FONT_UI, 10),
                                      text_color=self.COLORS['text_primary'],
                                      wraplength=700, justify="left")
        test_type_value.pack(anchor="w", padx=25, pady=(0, 8))
        
        # Parâmetros
        params_label = ctk.CTkLabel(scroll_frame, text=TEXTS["model_info_params"],
                                   font=(_FONT_UI, 11, "bold"),
                                   text_color=self.COLORS['accent_purple'])
        params_label.pack(anchor="w", padx=8, pady=(8, 2))
        
        params_value = ctk.CTkLabel(scroll_frame, text=model_info.get('parameters', ''),
                                   font=(_FONT_UI, 10),
                                   text_color=self.COLORS['text_primary'],
                                   wraplength=700, justify="left")
        params_value.pack(anchor="w", padx=25, pady=(0, 8))
        
        # Propósito
        purpose_label = ctk.CTkLabel(scroll_frame, text=TEXTS["model_info_purpose"],
                                    font=(_FONT_UI, 11, "bold"),
                                    text_color=self.COLORS['accent_cyan'])
        purpose_label.pack(anchor="w", padx=8, pady=(8, 2))
        
        purpose_value = ctk.CTkLabel(scroll_frame, text=model_info.get('purpose', ''),
                                    font=(_FONT_UI, 10),
                                    text_color=self.COLORS['text_primary'],
                                    wraplength=700, justify="left")
        purpose_value.pack(anchor="w", padx=25, pady=(0, 8))
        
        # Interpretação
        interp_label = ctk.CTkLabel(scroll_frame, text=TEXTS["model_info_interpretation"],
                                   font=(_FONT_UI, 11, "bold"),
                                   text_color=self.COLORS['success'])
        interp_label.pack(anchor="w", padx=8, pady=(8, 2))

        interp_value = ctk.CTkLabel(scroll_frame, text=model_info.get('interpretation', ''),
                                   font=(_FONT_UI, 10),
                                   text_color=self.COLORS['text_primary'],
                                   wraplength=700, justify="left")
        interp_value.pack(anchor="w", padx=25, pady=(0, 8))

        # Caso de uso
        use_case_label = ctk.CTkLabel(scroll_frame, text=TEXTS["model_info_use_case"],
                                     font=(_FONT_UI, 11, "bold"),
                                     text_color=self.COLORS['warning'])
        use_case_label.pack(anchor="w", padx=8, pady=(8, 2))

        use_case_value = ctk.CTkLabel(scroll_frame, text=model_info.get('use_case', ''),
                                     font=(_FONT_UI, 10),
                                     text_color=self.COLORS['text_primary'],
                                     wraplength=700, justify="left")
        use_case_value.pack(anchor="w", padx=25, pady=(0, 8))

        # Referências
        refs_label = ctk.CTkLabel(scroll_frame, text=TEXTS["model_info_references"],
                                 font=(_FONT_UI, 11, "bold"),
                                 text_color='#eab308')   # amarelo — não tem par no COLORS
        refs_label.pack(anchor="w", padx=8, pady=(8, 2))
        
        refs_value = ctk.CTkLabel(scroll_frame, text=model_info.get('references', ''),
                                 font=(_FONT_UI, 10, "italic"),
                                 text_color=self.COLORS['text_secondary'],
                                 wraplength=700, justify="left")
        refs_value.pack(anchor="w", padx=25, pady=(0, 8))

    def _open_results_viewer(self):
        if not self.output_folder:
            self.append_log("[!] Selecione uma pasta de saida primeiro.\n")
            return
        try:
            ResultsViewerWindow(self, self.output_folder)
        except Exception as e:
            self.append_log(f"[Erro] Erro ao abrir visualizador: {e}\n")
            self.append_log(traceback.format_exc())

    def _regenerate_summary_files(self):
        """Abre diálogo para selecionar pasta e regenera os 3 arquivos de síntese"""
        results_folder = filedialog.askdirectory(
            title="Selecione a pasta com resultados para atualizar síntese",
            initialdir=str(Path.home() / "Desktop")
        )
        
        if not results_folder:
            return
        
        results_folder = Path(results_folder)
        
        self.append_log("\n" + "="*80 + "\n")
        self.append_log(">> ATUALIZANDO RESULTADOS\n")
        self.append_log("="*80 + "\n")
        self.append_log(f"Pasta de resultados: {results_folder}\n\n")
        
        # Executar em thread separada para não travar GUI
        def _update_thread():
            try:
                self.append_log("Detectando modelos presentes...\n")
                
                # Descobrir quais modelos estão presentes
                models = set()
                for item in results_folder.iterdir():
                    if item.is_dir() and item.name not in ['reports']:
                        models.add(item.name)
                
                models = sorted(models)
                self.append_log(f"[OK] Modelos encontrados: {', '.join(models)}\n\n")
                
                # Determinar comparações disponíveis
                self.append_log("Determinando comparações para LRT:\n")
                comparisons = []
                
                if 'M0' in models and 'M1a' in models:
                    comparisons.append("M0 vs M1a")
                    self.append_log("  • M0 (null) vs M1a (alt) - Variação de ω entre sítios\n")
                
                if 'M1a' in models and 'M2a' in models:
                    comparisons.append("M1a vs M2a")
                    self.append_log("  • M1a (null) vs M2a (alt) - Seleção positiva\n")
                
                if 'M7' in models and 'M8' in models:
                    comparisons.append("M7 vs M8")
                    self.append_log("  • M7 (null) vs M8 (alt) - Seleção positiva (Beta)\n")
                
                if 'M0' in models and 'Branch' in models:
                    comparisons.append("M0 vs Branch")
                    self.append_log("  • M0 (null) vs Branch (alt) - Seleção por ramo\n")
                
                if 'Branch-site_null' in models and 'Branch-site' in models:
                    comparisons.append("Branch-site_null vs Branch-site")
                    self.append_log("  • Branch-site_null (null) vs Branch-site (alt) - Seleção branch-site\n")
                
                self.append_log(f"\nTotal de {len(comparisons)} comparações encontradas.\n\n")
                
                # Regenerar arquivos
                self.append_log("Regenerando arquivos de síntese...\n")
                generated_files = CodemlBatchAnalysis.regenerate_summary_files(results_folder)
                
                if generated_files:
                    self.append_log("\n[OK] ATUALIZACAO CONCLUIDA COM SUCESSO!\n")
                    self.append_log("="*80 + "\n")
                    for file_type, file_path in generated_files.items():
                        filepath = Path(file_path)
                        size = filepath.stat().st_size if filepath.exists() else 0
                        self.append_log(f"  [OK] {file_type:25s} | {size:,} bytes\n")
                    self.append_log("="*80 + "\n")
                    # Atualiza pasta de saída para a pasta selecionada e habilita o botão
                    self.output_folder = results_folder
                    self.after(0, self._update_models_state)
                else:
                    self.append_log("\n[Erro] Nenhum arquivo foi gerado.\n")
            
            except Exception as e:
                self.append_log(f"\n[Erro] ERRO: {str(e)}\n")
                self.append_log(traceback.format_exc())
        
        update_thread = threading.Thread(target=_update_thread, daemon=True)
        update_thread.start()

    def _update_results_files(self):
        """Alias para _regenerate_summary_files com novo nome"""
        self._regenerate_summary_files()

    def _update_models_state(self):
        enabled = all([self.input_folder, self.tree_file, self.output_folder])
        state = "normal" if enabled else "disabled"
        
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

    def append_log(self, text: str):
        """Adiciona mensagem ao log com tag de cor apropriada e feedback visual premium"""
        def _append():
            tag = None
            
            # Detectar tipo de mensagem por palavra-chave
            if any(x in text for x in ["[OK]", "SUCESSO", "COMPLETADO", "FINALIZADO", "OK", "SALVO"]):
                tag = "success"
            elif any(x in text for x in ["[Erro]", "ERRO", "TIMEOUT", "FALHA", "PROBLEMA"]):
                tag = "error"
            elif any(x in text for x in ["[!]", "[AVISO]", "AVISO", "CUIDADO", "ATENÇÃO"]):
                tag = "warning"
            elif any(x in text for x in ["[+]", "[-]", "[tag]", "[edit]", "INFO", "INICIANDO", "PROCESSANDO"]):
                tag = "info"
            elif any(x in text for x in ["═", "───", "╔", "╚", "║"]):
                tag = "header"
            
            self.log.insert("end", text, tag)
            self.log.see("end")  # Auto-scroll para o final
        
        self.after(0, _append)

    def _poll_stop_count(self):
        if self.analysis_instance:
            cnt = getattr(self.analysis_instance, 'current_stop_count', 0)
            self.stop_label.configure(text=TEXTS["status_stops_template"].format(n=cnt))
        self.after(1000, self._poll_stop_count)

    def _continue_all_for_gene(self):
        if self.manual_all_event:
            self.manual_all_event.set()
            if self.manual_event: 
                self.manual_event.set()
            self.append_log("▶ Resolvendo todos os stops deste Exon automaticamente...\n")

    def _toggle_pause(self):
        if not self.pause_event:
            return
        if not self.analysis_thread or not self.analysis_thread.is_alive():
            return

        if self.pause_event.is_set():
            # Pausar: bloquear próximos genes/modelos + suspender processos em curso
            self.pause_event.clear()
            self._suspend_active_codeml()
            self.btn_pause.configure(text=TEXTS["btn_resume"], fg_color="#10b981")
            self.status_indicator.configure(text=TEXTS["status_paused"],
                                            text_color=self.COLORS['warning'])
            self.append_log("|| Analise pausada (CODEML suspenso).\n")
        else:
            # Retomar: desbloquear threads + retomar processos suspensos
            self._resume_active_codeml()
            self.pause_event.set()
            self.btn_pause.configure(text=TEXTS["btn_pause"],
                                     fg_color=self.COLORS['warning'])
            self.status_indicator.configure(text=TEXTS["status_running"],
                                            text_color=self.COLORS['success'])
            self.append_log(">> Analise retomada.\n")

    def _suspend_active_codeml(self):
        """Suspende todos os processos CODEML ativos usando psutil."""
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
        """Retoma todos os processos CODEML suspensos usando psutil."""
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
            self.append_log("[!] Nenhuma analise em execucao.\n")
            return

        self.append_log("[STOP] PARANDO ANALISE...\n")
        # 1. Sinalizar para o backend parar de iniciar novos genes/modelos
        self.stop_event.set()
        self.status_indicator.configure(text=TEXTS["status_stopped"], text_color=self.COLORS['danger'])

        # 2. Matar TODOS os processos CODEML ativos imediatamente
        if self.analysis_instance:
            try:
                procs = list(getattr(self.analysis_instance, '_active_processes', []))
                if procs:
                    self.append_log(f">> Terminando {len(procs)} processo(s) CODEML...\n")
                for proc in procs:
                    try:
                        if proc.poll() is None:
                            proc.terminate()
                            proc.wait(timeout=2)
                    except Exception:
                        try:
                            proc.kill()
                        except Exception:
                            pass
            except Exception as e:
                self.append_log(f"[!] Erro ao terminar processos: {e}\n")

        # 3. Desbloquear qualquer evento de pausa/manual para não travar threads
        if self.pause_event:
            self.pause_event.set()
        if self.manual_event:
            self.manual_event.set()
        if self.manual_all_event:
            self.manual_all_event.set()

        self.btn_run.configure(state="normal", fg_color=self.COLORS['success'],
                               text_color=self.COLORS['text_primary'])
        self.append_log("[STOP] Analise interrompida pelo usuario.\n")

    def start_analysis(self):
        selected = [k for k, v in self.model_vars.items() if v.get()]
        if not selected:
            self.append_log("[!] Selecione pelo menos um modelo.\n")
            return
        
        # ═══ AUTO-COMPLETAR MODELOS NULOS ═══
        original_selected = selected.copy()
        include_neutral = self.include_neutral_models.get()
        selected = CodemlBatchAnalysis.auto_complete_null_models(selected, include_neutral=include_neutral)
        
        # Informar ao usuário quais modelos foram auto-adicionados
        if len(selected) > len(original_selected):
            added = set(selected) - set(original_selected)
            self.append_log(f"[OK] Modelos auto-adicionados: {', '.join(sorted(added))}\n")
            self.append_log(f"   (Necessários para comparação LRT automática)\n\n")
        
        self.stop_event.clear()
        self.btn_run.configure(state="disabled")
        self.pause_event.set()
        self.manual_all_event.clear()
        self.manual_event.clear()
        
        # Feedback visual: mudar cor do botão e status
        self.btn_run.configure(fg_color=self.COLORS['bg_card'],
                               text_color=self.COLORS['text_muted'])
        self.status_indicator.configure(text=TEXTS["status_running"], text_color=self.COLORS['success'])

        self.analysis_thread = threading.Thread(target=self._run_thread, args=(selected,), daemon=True)
        self.analysis_thread.start()

    def _run_thread(self, selected):
        sys.stdout = StdoutRedirect(self.append_log)
        try:
            analysis = CodemlBatchAnalysis()
            self.analysis_instance = analysis
            
            needs_branchsite = any('Branch-site' in m for m in selected)
            needs_branch = 'Branch' in selected
            
            # ═══ VALIDAÇÕES ═══
            if needs_branchsite and not self.tree_branchsite_labeled:
                self.append_log("\n[Erro] ERRO: Modelo Branch-Site sem arvore etiquetada!\n")
                self.append_log("   -> Use 'Marcar Branch-site' na aba Branch-Site.\n\n")
                self.btn_run.configure(state="normal", 
                                      fg_color=self.COLORS['success'],
                                      text_color=self.COLORS['text_primary'])
                return
            
            if needs_branch and not self.tree_branch_labeled:
                self.append_log("\n[Erro] ERRO: Modelo Branch sem arvore etiquetada!\n")
                self.append_log("   -> Use 'Marcar Ramos (Multiplas Tags)' na aba Branch Model.\n\n")
                self.btn_run.configure(state="normal",
                                      fg_color=self.COLORS['success'],
                                      text_color=self.COLORS['text_primary'])
                return

            # ═══ PREPARAR CONFIGURAÇÃO ═══
            # Auto-completar modelos nulos se checkbox estiver ativado
            models_to_run = selected
            if self.include_neutral_models.get():
                models_to_run = CodemlBatchAnalysis.auto_complete_null_models(
                    selected, 
                    include_neutral=True
                )
                if models_to_run != selected:
                    self.append_log(f"\n[OK] Auto-adicionados modelos nulos: {', '.join(set(models_to_run) - set(selected))}\n")
            
            analysis.config = {
                'input_folder': self.input_folder,
                'tree_file': self.tree_file,
                'output_folder': self.output_folder,
                'models': models_to_run,
                'custom_model_params': self.custom_model_params,
                'labeled_tree_content': self.tree_branch_labeled,
                'labeled_tree_branchsite': self.tree_branchsite_labeled,
                'omega': float(self.entry_omega.get() or 0.5),
                'cleandata': int(self.cleandata_var.get()),
                'timeout': 1600,
                'run_lrt': True,
                'n_workers': int(self.cores_var.get()),
                'heuristic_mode': self.heuristic_mode_var.get(),
                'auto_continue_stop_codons': self.ignore_stop_codons_var.get(),
                'auto_prune_tree': self.auto_prune_tree_var.get(),
                'pause_event': self.pause_event,
                'manual_continue_event': self.manual_event,
                'manual_continue_all_event': self.manual_all_event,
                'stop_event': self.stop_event
            }
            
            self.append_log("╔" + "═"*58 + "╗\n")
            self.append_log("║     INICIANDO ANALISE DE SELECAO POSITIVA     " + " "*10 + "║\n")
            self.append_log("╚" + "═"*58 + "╝\n\n")
            self.append_log("▶ INICIANDO ANALISE BATCH\n")
            self.append_log("="*60 + "\n")
            
            analysis.run_batch_analysis()
            
        except Exception as e:
            self.append_log(f"[Erro] Erro critico: {e}\n")
            self.append_log(f"{traceback.format_exc()}\n")
        finally:
            self.analysis_instance = None

            if self.stop_event.is_set():
                self.append_log("\n" + "="*60 + "\n")
                self.append_log("■ ANALISE INTERROMPIDA\n")
                self.append_log("="*60 + "\n")
            else:
                self.append_log("\n" + "="*60 + "\n")
                self.append_log("[OK] ANALISE CONCLUIDA\n")
                self.append_log("="*60 + "\n")

            def _reset_ui():
                self.btn_run.configure(
                    state="normal",
                    fg_color=self.COLORS['success'],
                    text_color=self.COLORS['text_primary']
                )
                self.status_indicator.configure(
                    text=TEXTS["status_ready"],
                    text_color=self.COLORS['text_tertiary']
                )
                self.btn_pause.configure(
                    text=TEXTS["btn_pause"],
                    fg_color=self.COLORS['warning']
                )
                self._update_models_state()

            self.after(0, _reset_ui)


if __name__ == "__main__":
    app = App()
    app.mainloop()
