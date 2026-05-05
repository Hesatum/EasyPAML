import customtkinter as ctk
from tkinter import filedialog
from tkinter import simpledialog, Canvas
from pathlib import Path
import threading
import io
import sys
import os
import signal

# Importar backend - ajuste de caminho relativo
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.codeml_backend import CodemlBatchAnalysis
from .results_viewer import ResultsViewerWindow

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
        self.title(f"⚙️ Configurar {model_code}")
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
        
        ctk.CTkLabel(header, text=f"Modelo: {model_code}", font=("Roboto", 16, "bold"),
                    text_color=self.COLORS['text_primary']).pack(side="left")
        
        # ═══ FORM FRAME ═══
        form_frame = ctk.CTkScrollableFrame(self, fg_color=self.COLORS['bg_card'],
                                           corner_radius=10, border_width=1,
                                           border_color='#2a2a2a')
        form_frame.pack(fill="both", expand=True, padx=20, pady=(0, 15))

        fields = [
            ('description', 'Descrição', 'Texto descritivo'),
            ('model', 'Model', 'Tipo de modelo (1-6)'),
            ('NSsites', 'NSsites', 'Modelos de sítios'),
            ('CodonFreq', 'Codon Freq', 'Frequência de códons'),
            ('fix_omega', 'Fix Omega', '0=livre, 1=fixo'),
            ('omega', 'Omega Inicial', 'Razão dN/dS')
        ]

        for name, label, placeholder in fields:
            # Label com ícone
            lbl_frame = ctk.CTkFrame(form_frame, fg_color='transparent')
            lbl_frame.pack(fill="x", padx=12, pady=(12, 4))
            
            ctk.CTkLabel(lbl_frame, text=f"📝 {label}:", font=("Roboto", 11, "bold"),
                        text_color=self.COLORS['text_primary']).pack(anchor="w")
            
            ctk.CTkLabel(lbl_frame, text=placeholder, font=("Roboto", 9),
                        text_color='#999999').pack(anchor="w", pady=(0, 4))
            
            # Entry com estilo premium
            ent = ctk.CTkEntry(form_frame, placeholder_text=placeholder,
                              fg_color=self.COLORS['bg_dark'],
                              border_color=self.COLORS['accent_blue'],
                              border_width=2,
                              corner_radius=6,
                              text_color=self.COLORS['text_primary'],
                              placeholder_text_color='#666666')
            ent.pack(fill="x", padx=12, pady=(0, 8))
            
            val = self._initial_val(name, initial)
            ent.insert(0, str(val))
            self.entries[name] = ent

        # ═══ BUTTONS ═══
        btn_frame = ctk.CTkFrame(self, fg_color='transparent')
        btn_frame.pack(fill="x", padx=20, pady=(0, 20))
        
        btn_cancel = ctk.CTkButton(btn_frame, text="✕ Cancelar", fg_color='#3d3d3d',
                                   hover_color='#4d4d4d', command=self.destroy,
                                   font=("Roboto", 11, "bold"), corner_radius=6)
        btn_cancel.pack(side="left", fill="x", expand=True, padx=(0, 8))
        
        btn_save = ctk.CTkButton(btn_frame, text="✓ Salvar", fg_color=self.COLORS['success'],
                                hover_color=self.COLORS['success_hover'],
                                command=self._on_save, font=("Roboto", 11, "bold"), corner_radius=6)
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
                if '.' in val: out[k] = float(val)
                else: out[k] = int(val)
            except: out[k] = val
        
        self.parent.custom_model_params[self.model_code] = out
        self.parent.model_ctl_labels[self.model_code].configure(text="Configurado ✨", text_color="#fbbf24")
        self.parent.append_log(f"✅ Parâmetros de {self.model_code} salvos com sucesso!\n")
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
        self.title(f"🏷️ Etiquetar Árvore - {mode.upper()}")
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
            ctk.CTkLabel(self, text="⚠️ Biopython não instalado. Execute: pip install biopython",
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
        title_label = ctk.CTkLabel(left_frame, text="📋 Instruções", 
                                  font=("Roboto", 14, "bold"),
                                  text_color=self.TEXT_PRIMARY)
        title_label.pack(pady=(15, 10), padx=15)
        
        if self.mode == 'branchsite':
            instructions = "• Clique nos CÍRCULOS\n  para marcar/desmarcar\n  foreground (#1).\n\n• Apenas tag #1\n  é permitida.\n\n• Cor: Vermelho"
        else:
            instructions = "• Clique nos CÍRCULOS\n  para atribuir tags.\n\n• Digite número\n  da tag (1, 2, 3...).\n\n• Cada número\n  recebe cor única."
        
        inst_label = ctk.CTkLabel(left_frame, text=instructions, wraplength=240, justify="left", 
                                 font=("Roboto", 10), text_color=self.TEXT_SECONDARY)
        inst_label.pack(pady=10, padx=15)

        legend_header = ctk.CTkLabel(left_frame, text="🎨 Tags Ativas", font=("Roboto", 12, "bold"),
                                    text_color=self.ACCENT_BLUE)
        legend_header.pack(pady=(20, 5), padx=15)
        
        self.legend_frame = ctk.CTkFrame(left_frame, fg_color=self.BG_CARD, corner_radius=8)
        self.legend_frame.pack(fill='both', expand=True, padx=15, pady=5)

        btn_frame = ctk.CTkFrame(left_frame, fg_color='transparent')
        btn_frame.pack(side='bottom', fill='x', padx=12, pady=15)
        
        ctk.CTkButton(btn_frame, text='✓ SALVAR', fg_color=self.SUCCESS, hover_color=self.SUCCESS_HOVER, 
                     command=self._on_save, height=40, font=("Roboto", 12, "bold"),
                     text_color=self.TEXT_PRIMARY, corner_radius=6).pack(fill='x', pady=5)
        ctk.CTkButton(btn_frame, text='✕ CANCELAR', fg_color='#3d3d3d', hover_color='#4d4d4d',
                     command=self.destroy, height=40, font=("Roboto", 12, "bold"),
                     text_color=self.TEXT_PRIMARY, corner_radius=6).pack(fill='x', pady=5)

        # Gráfico
        if self.tree_path is None:
            ctk.CTkLabel(plot_frame, text="⚠️ Nenhuma árvore selecionada.", font=("Arial", 14)).pack(pady=50)
            return

        try:
            self.tree = Phylo.read(str(self.tree_path), 'newick')
        except Exception as e:
            ctk.CTkLabel(plot_frame, text=f"❌ Erro ao carregar árvore:\n{e}", font=("Arial", 12)).pack(pady=50)
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
                            
                            # Quanto menor a distância, mais próximo está
                            if distance < min_distance:
                                min_distance = distance
                                closest_tag = tag
                    except:
                        pass
                
                return closest_tag
            except:
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
            except:
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
                self.parent.append_log(f"🔴 Tag #1 removida de {self._get_clade_name(clade)}\n")
            else:
                self._remove_tag_recursively(clade)
                self._apply_tag_recursively(clade, '#1')
                self.marked_clades[clade] = '#1'
                self.parent.append_log(f"🟢 Tag #1 aplicada a {self._get_clade_name(clade)}\n")
        
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
                        self.parent.append_log(f"🗑️ Tag {current_tag} removida de {clade_name}\n")
                    elif response.isdigit():
                        self._remove_tag_recursively(clade)
                        new_tag = f"#{response}"
                        self._apply_tag_recursively(clade, new_tag)
                        self.marked_clades[clade] = new_tag
                        self.parent.append_log(f"✏️ Tag alterada para {new_tag} em {clade_name}\n")
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
                    self.parent.append_log(f"🏷️ Tag {new_tag} aplicada a {clade_name}\n")

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
        except:
            return '#ff0000'

    def _refresh_legend(self):
        """Atualiza legenda com tags ativas"""
        for widget in self.legend_frame.winfo_children():
            widget.destroy()

        active_tags = sorted(set(self.clade_tags.values()))
        
        if not active_tags:
            no_tags = ctk.CTkLabel(self.legend_frame, text="Nenhuma tag ativa", 
                                  text_color="#888888", font=("Arial", 11, "italic"))
            no_tags.pack(anchor='w', padx=15, pady=10)
        else:
            for tag in active_tags:
                color = self._get_tag_color(tag)
                
                row_frame = ctk.CTkFrame(self.legend_frame, fg_color="transparent")
                row_frame.pack(fill='x', padx=10, pady=4)
                
                color_box = Canvas(row_frame, width=24, height=18, highlightthickness=0)
                try:
                    color_box.configure(bg=self.legend_frame.cget('fg_color')[1])
                except:
                    color_box.configure(bg='#2b2b2b')
                color_box.create_rectangle(2, 2, 22, 16, fill=color, outline='white', width=2)
                color_box.pack(side='left', padx=5)
                
                tag_label = ctk.CTkLabel(row_frame, text=tag, font=("Arial", 12, "bold"))
                tag_label.pack(side='left', padx=8)
                
                count = sum(1 for t in self.clade_tags.values() if t == tag)
                count_label = ctk.CTkLabel(row_frame, text=f"({count})", 
                                          font=("Arial", 10), text_color="#888888")
                count_label.pack(side='left')

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
                self.parent.append_log(f"✅ Árvore Branch-Site salva ({applied_terminals} terminais).\n")
            else:
                self.parent.tree_branch_labeled = newick_str
                self.parent.append_log(f"✅ Árvore Branch salva ({applied_terminals} terminais).\n")

                branchsite_version = re.sub(r"\s*#(?!1)\d+\b", "", newick_str)
                branchsite_version = re.sub(r"\s+", " ", branchsite_version).strip()
                self.parent.tree_branchsite_labeled = branchsite_version

        except Exception as e:
            self.parent.append_log(f"❌ Erro ao gerar Newick: {e}\n")

        self.destroy()


class App(ctk.CTk):
    # ═══════════════════════════════════════════════════════════════════════
    # PALETA DE CORES PREMIUM - Estilo YouTube/Instagram Dark
    # ═══════════════════════════════════════════════════════════════════════
    COLORS = {
        # Backgrounds — near-black, Linear/Discord inspired
        'bg_darkest':     '#080808',
        'bg_dark':        '#0c0c0e',
        'bg_sidebar':     '#111115',
        'bg_card':        '#16161a',
        'bg_card_hover':  '#1e1e24',
        'bg_feed':        '#111115',
        'bg_input':       '#1e1e24',

        # Text hierarchy
        'text_primary':   '#ededef',
        'text_secondary': '#9898a6',
        'text_tertiary':  '#5e5e6e',
        'text_muted':     '#3a3a48',

        # Primary accent — indigo (Linear-inspired)
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

        # Borders
        'border':         '#222228',
        'border_hover':   '#32323e',
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
        self.attributes("-alpha", 0.93)

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
        self.wgs_mode_var = ctk.BooleanVar(value=False)

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
        logo_frame.pack(fill="x", pady=(18, 12))
        ctk.CTkLabel(logo_frame, text="🧬", font=("Roboto", 34, "bold")).pack()
        ctk.CTkLabel(logo_frame, text="EasyPAML",
                     font=("Roboto", 18, "bold"),
                     text_color=self.COLORS['text_primary']).pack(pady=(4, 0))
        ctk.CTkLabel(logo_frame, text="Seleção Positiva",
                     font=("Roboto", 9),
                     text_color=self.COLORS['accent_blue']).pack()

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
            hdr.pack(fill='x', padx=12, pady=(8, 0))
            ctk.CTkLabel(hdr,
                         text=f"{icon}  {title}" if icon else title,
                         font=("Roboto", 8, "bold"),
                         text_color=self.COLORS['text_tertiary']).pack(anchor='w')
            ctk.CTkFrame(card, fg_color=self.COLORS['border'],
                         height=1, corner_radius=0).pack(fill='x', padx=12, pady=(5, 0))
            inner = ctk.CTkFrame(card, fg_color='transparent')
            inner.pack(fill='x', padx=10, pady=(8, 10))
            return inner

        def _obtn(parent, text, cmd, color, **kw):
            """Outline-style button."""
            return ctk.CTkButton(
                parent, text=text, command=cmd,
                fg_color=self.COLORS['bg_card_hover'],
                hover_color=color,
                text_color=color,
                border_width=1, border_color=color,
                corner_radius=8, **kw)

        # ── Arquivos ─────────────────────────────────────────────────
        fi = _sec(_sb, "ARQUIVOS", "📁")

        self.btn_input = _obtn(fi, "📁 Pasta .fas",
                               self.select_input_folder,
                               self.COLORS['accent_blue'],
                               font=("Roboto", 10, "bold"), height=32)
        self.btn_input.pack(fill='x', pady=(0, 2))
        self.label_input = ctk.CTkLabel(fi, text="Não selecionado",
                                        font=("Roboto", 8),
                                        wraplength=230,
                                        text_color=self.COLORS['text_muted'])
        self.label_input.pack(anchor='w', padx=4, pady=(0, 6))

        self.btn_tree = _obtn(fi, "🌳 Árvore (.nwk)",
                              self.select_tree_file,
                              self.COLORS['accent_cyan'],
                              font=("Roboto", 10, "bold"), height=32)
        self.btn_tree.pack(fill='x', pady=(0, 2))
        self.label_tree = ctk.CTkLabel(fi, text="Não selecionado",
                                       font=("Roboto", 8),
                                       wraplength=230,
                                       text_color=self.COLORS['text_muted'])
        self.label_tree.pack(anchor='w', padx=4, pady=(0, 6))

        self.btn_output = _obtn(fi, "💾 Pasta Saída",
                                self.select_output_folder,
                                self.COLORS['accent_purple'],
                                font=("Roboto", 10, "bold"), height=32)
        self.btn_output.pack(fill='x', pady=(0, 2))
        self.label_output = ctk.CTkLabel(fi, text="Não selecionado",
                                         font=("Roboto", 8),
                                         wraplength=230,
                                         text_color=self.COLORS['text_muted'])
        self.label_output.pack(anchor='w', padx=4)

        # ── Resultados ───────────────────────────────────────────────
        ri = _sec(_sb, "RESULTADOS", "📊")

        self.btn_results = _obtn(ri, "📊 Ver Resultados",
                                  self._open_results_viewer,
                                  self.COLORS['info'],
                                  font=("Roboto", 10, "bold"), height=32)
        self.btn_results.pack(fill='x', pady=(0, 6))
        self.btn_results.configure(state="disabled")

        self.btn_update_results = _obtn(ri, "⚡ Atualizar Resultados",
                                         self._update_results_files,
                                         '#10b981',
                                         font=("Roboto", 10, "bold"), height=32)
        self.btn_update_results.pack(fill='x', pady=(0, 2))
        self.label_update_results = ctk.CTkLabel(ri,
                                                  text="Atualizar arquivos de análise",
                                                  font=("Roboto", 8, "italic"),
                                                  wraplength=230,
                                                  text_color=self.COLORS['text_muted'])
        self.label_update_results.pack(anchor='w', padx=4)

        # ── Configurações ────────────────────────────────────────────
        ci = _sec(_sb, "CONFIGURAÇÕES", "⚙")

        ctk.CTkLabel(ci, text="dN/dS Inicial (ω):",
                     font=("Roboto", 9, "bold"), anchor='w',
                     text_color=self.COLORS['text_secondary']).pack(anchor='w')
        self.omega_label = ctk.CTkLabel(ci, text="")  # kept for compat, unused
        self.entry_omega = ctk.CTkEntry(ci, placeholder_text="0.5",
                                        fg_color=self.COLORS['bg_card_hover'],
                                        border_color=self.COLORS['border_hover'],
                                        text_color=self.COLORS['text_primary'],
                                        height=32)
        self.entry_omega.insert(0, "0.5")
        self.entry_omega.pack(fill='x', pady=(4, 10))

        # Remover gaps toggle
        self.cleandata_var = ctk.BooleanVar(value=True)
        row_g = ctk.CTkFrame(ci, fg_color='transparent')
        row_g.pack(fill='x', pady=(0, 8))
        ctk.CTkLabel(row_g, text="Remover gaps",
                     font=("Roboto", 10),
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

        # CPU slider
        ctk.CTkLabel(ci, text="⚡ CPUs (paralelismo):",
                     font=("Roboto", 9, "bold"), anchor='w',
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
            font=("Roboto", 10, "bold"),
            text_color=self.COLORS['accent_blue'],
            width=32)
        self.cores_disp.pack(side='right', padx=(6, 0))
        ctk.CTkLabel(ci, text=f"(detectado: {_max} núcleos)",
                     font=("Roboto", 8),
                     text_color=self.COLORS['text_muted']).pack(anchor='w')

        # Modo WGS toggle
        row_wgs = ctk.CTkFrame(ci, fg_color='transparent')
        row_wgs.pack(fill='x', pady=(10, 0))
        ctk.CTkLabel(row_wgs, text="Modo WGS (ndata)",
                     font=("Roboto", 10),
                     text_color=self.COLORS['text_secondary']).pack(side='left')
        self.cb_wgs = ctk.CTkSwitch(
            row_wgs, text="",
            variable=self.wgs_mode_var,
            onvalue=True, offvalue=False,
            switch_width=36, switch_height=18,
            progress_color=self.COLORS['accent_cyan'],
            button_color='#f0fdff',
            button_hover_color='#cffafe',
            fg_color=self.COLORS['border'])
        self.cb_wgs.pack(side='right')

        self.main_frame = ctk.CTkFrame(self, fg_color=self.COLORS['bg_dark'])
        self.main_frame.pack(side="right", fill="both", expand=True, padx=20, pady=20)

        self.tabs = ctk.CTkTabview(self.main_frame, fg_color=self.COLORS['bg_card'],
                                   segmented_button_fg_color=self.COLORS['bg_card'],
                                   segmented_button_selected_color=self.COLORS['accent_blue'],
                                   text_color=self.COLORS['text_primary'],
                                   corner_radius=10)
        self.tabs.pack(fill="both", expand=True, padx=0, pady=(0, 15))
        self.tabs.add("Site Models")
        self.tabs.add("Branch Model")
        self.tabs.add("Branch-Site")

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
            status_bar, text="● Pronto",
            font=("Roboto", 10, "bold"),
            text_color=self.COLORS['text_muted']
        )
        self.status_indicator.pack(side="left", padx=(14, 20), pady=8)

        self.stop_label = ctk.CTkLabel(
            status_bar, text="◆ Stops: 0",
            font=("Roboto", 10, "bold"),
            text_color=self.COLORS['danger']
        )
        self.stop_label.pack(side="left")

        neutral_row = ctk.CTkFrame(status_bar, fg_color='transparent')
        neutral_row.pack(side="right", padx=(0, 6))
        ctk.CTkLabel(neutral_row, text="Modelos nulos automáticos",
                     font=("Roboto", 9),
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
            font=("Roboto", 12, "bold"),
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
                font=("Roboto", 11, "bold"),
                height=42, corner_radius=8)

        self.btn_run = _action_btn(btn_frame, "▶  INICIAR",
                                   self.start_analysis, self.COLORS['success'])
        self.btn_run.pack(side="left", fill="both", expand=True, padx=(12, 4), pady=10)

        self.btn_pause = _action_btn(btn_frame, "⏸  PAUSAR",
                                     self._toggle_pause, self.COLORS['warning'])
        self.btn_pause.pack(side="left", fill="both", expand=True, padx=4, pady=10)

        self.btn_stop = _action_btn(btn_frame, "⏹  PARAR",
                                    self._stop_analysis, self.COLORS['danger'])
        self.btn_stop.pack(side="left", fill="both", expand=True, padx=4, pady=10)

        self.btn_resolve_stops = _action_btn(btn_frame, "⏭  STOPS",
                                              self._continue_all_for_gene, self.COLORS['info'])
        self.btn_resolve_stops.pack(side="left", fill="both", expand=True, padx=(4, 12), pady=10)

        # ═══ LOG FRAME COM HEADER ═══
        log_container = ctk.CTkFrame(self.main_frame, fg_color='transparent')
        log_container.pack(fill="both", expand=True, padx=0, pady=0)
        
        log_header = ctk.CTkFrame(log_container, fg_color='transparent')
        log_header.pack(fill="x", padx=0, pady=(10, 0))
        
        ctk.CTkLabel(log_header, text="📋 LOG DE EXECUÇÃO", font=("Roboto", 11, "bold"),
                    text_color=self.COLORS['accent_blue']).pack(side="left", padx=0)
        
        ctk.CTkLabel(log_header, text="Feedback em tempo real da análise",
                    font=("Roboto", 9), text_color=self.COLORS['text_tertiary']).pack(side="left", padx=10)

        self.log = ctk.CTkTextbox(log_container, font=("Cascadia Code", 9),
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

        # Welcome message
        self.log.insert("end",
            "  ╔══════════════════════════════════════════════╗\n"
            "  ║   EasyPAML  ·  PAML / CODEML interface      ║\n"
            "  ╚══════════════════════════════════════════════╝\n\n"
            "  1. Selecione a pasta de alinhamentos (.fas)\n"
            "  2. Escolha a árvore filogenética (.nwk)\n"
            "  3. Defina a pasta de saída\n"
            "  4. Marque os modelos e clique  ▶ INICIAR\n\n"
            "  ─────────────────────────────────────────────────\n\n"
        )

        self._update_models_state()
        self._poll_stop_count()

    def _update_cores_label(self, value=None):
        n = int(self.cores_var.get())
        self.cores_disp.configure(text=f"{n}×")

    def _setup_model_list(self):
        MODEL_META = {
            'M0':               {'color': '#3b82f6', 'desc': 'ω único para todos os sítios  ·  modelo baseline'},
            'M1a':              {'color': '#06b6d4', 'desc': 'Dois sítios: purificação (0<ω<1) e neutro (ω=1)'},
            'M2a':              {'color': '#10b981', 'desc': 'Três sítios: adiciona classe com ω > 1  ·  BEB/NEB'},
            'M7':               {'color': '#8b5cf6', 'desc': 'Distribuição Beta de ω contida em (0, 1)'},
            'M8':               {'color': '#ec4899', 'desc': 'Beta + sítios sob seleção positiva (ω > 1)  ·  BEB/NEB'},
            'Branch':           {'color': '#f59e0b', 'desc': 'ω livre por ramo  ·  requer árvore etiquetada'},
            'Branch-site':      {'color': '#ef4444', 'desc': 'Seleção positiva episódica no ramo alvo'},
            'Branch-site_null': {'color': '#6d6d6d', 'desc': 'Modelo nulo para o teste LRT Branch-site'},
        }

        models = {
            "Site Models": ['M0', 'M1a', 'M2a', 'M7', 'M8'],
            "Branch Model": ['Branch'],
            "Branch-Site": ['Branch-site', 'Branch-site_null']
        }

        for tab_name, codes in models.items():
            scroll_frame = ctk.CTkScrollableFrame(
                self.tabs.tab(tab_name),
                fg_color=self.COLORS['bg_card'],
                label_text="",
                label_font=("Roboto", 11, "bold"),
                label_text_color=self.COLORS['accent_blue']
            )
            scroll_frame.pack(fill="both", expand=True, padx=8, pady=8)

            for code in codes:
                meta   = MODEL_META.get(code, {'color': self.COLORS['accent_blue'], 'desc': ''})
                accent = meta['color']
                desc   = meta['desc']

                # ── Outer card ────────────────────────────────────────────────
                card = ctk.CTkFrame(scroll_frame, fg_color=self.COLORS['bg_feed'],
                                    corner_radius=10, border_width=1,
                                    border_color=self.COLORS['border'])
                card.pack(fill="x", pady=3, padx=4)

                # Left accent bar (slim, 4px)
                accent_bar = ctk.CTkFrame(card, fg_color=accent, width=4, corner_radius=2)
                accent_bar.pack(side="left", fill="y", padx=(5, 8), pady=6)
                accent_bar.pack_propagate(False)

                # Text block
                text_block = ctk.CTkFrame(card, fg_color='transparent')
                text_block.pack(side="left", fill="both", expand=True, pady=7)

                display_name = self.codeml_backend.MODEL_CONFIGS[code].get('display_name', code)
                var = ctk.BooleanVar(value=False)
                sw_row = ctk.CTkFrame(text_block, fg_color='transparent')
                sw_row.pack(fill='x', anchor='w')
                cb = ctk.CTkSwitch(
                    sw_row,
                    text=f"  {display_name}",
                    variable=var,
                    onvalue=True, offvalue=False,
                    command=self._update_models_state,
                    switch_width=36, switch_height=18,
                    progress_color=accent,
                    button_color='#f0f0ff',
                    button_hover_color='white',
                    fg_color=self.COLORS['border'],
                    text_color=self.COLORS['text_primary'],
                    font=("Roboto", 11, "bold")
                )
                cb.pack(side='left', anchor='w', pady=(2, 0))

                if desc:
                    ctk.CTkLabel(
                        text_block, text=f"   {desc}",
                        font=("Roboto", 9),
                        text_color=self.COLORS['text_tertiary'],
                        anchor='w',
                        justify='left',
                        wraplength=260
                    ).pack(anchor="w", pady=(1, 4), fill='x')

                # Right controls
                controls = ctk.CTkFrame(card, fg_color='transparent')
                controls.pack(side="right", padx=8, pady=7)

                lbl = ctk.CTkLabel(controls, text="padrão", font=("Roboto", 8),
                                   text_color=self.COLORS['text_tertiary'])
                lbl.pack(side="left", padx=(0, 8))

                info_btn = ctk.CTkButton(
                    controls, text="?", width=26, height=26,
                    fg_color=self.COLORS['bg_card_hover'],
                    hover_color=accent,
                    text_color=accent,
                    corner_radius=6,
                    border_width=1, border_color=self.COLORS['border_hover'],
                    font=("Roboto", 11, "bold"),
                    command=lambda c=code: self._show_model_info(c)
                )
                info_btn.pack(side="left", padx=(0, 4))

                gear = ctk.CTkButton(
                    controls, text="⚙", width=26, height=26,
                    fg_color=self.COLORS['bg_card'],
                    hover_color=self.COLORS['accent_purple'],
                    text_color=self.COLORS['text_secondary'],
                    corner_radius=6,
                    border_width=1, border_color=self.COLORS['border_hover'],
                    font=("Roboto", 11),
                    command=lambda c=code: self._open_config_window(c)
                )
                gear.pack(side="left")

                self.model_vars[code]         = var
                self.model_ctl_labels[code]   = lbl
                self.model_checkboxes[code]   = cb
                self.model_gear_buttons[code] = gear
        
        # ═══ BRANCH: Botão de etiquetagem ═══
        branch_tab = self.tabs.tab("Branch Model")
        
        self.btn_label_branch = ctk.CTkButton(
            branch_tab,
            text="🏷️ Marcar Ramos (Múltiplas Tags)",
            fg_color=self.COLORS['bg_card'],
            hover_color=self.COLORS['accent_blue'],
            command=lambda: self._open_tree_labeler(mode='branch'),
            height=44,
            font=("Roboto", 12, "bold"),
            text_color=self.COLORS['accent_blue'],
            border_width=1, border_color=self.COLORS['accent_blue'],
            corner_radius=8
        )
        self.btn_label_branch.pack(fill='x', padx=12, pady=(20, 12))

        # ═══ BRANCHSITE: Botão de etiquetagem ═══
        branchsite_tab = self.tabs.tab("Branch-Site")
        
        self.btn_label_branchsite = ctk.CTkButton(
            branchsite_tab,
            text="🏷️ Marcar Branch-site",
            fg_color=self.COLORS['bg_card'],
            hover_color=self.COLORS['accent_pink'],
            command=lambda: self._open_tree_labeler(mode='branchsite'),
            height=44,
            font=("Roboto", 12, "bold"),
            text_color=self.COLORS['accent_pink'],
            border_width=1, border_color=self.COLORS['accent_pink'],
            corner_radius=8
        )
        self.btn_label_branchsite.pack(fill='x', padx=12, pady=(20, 12))

    def _open_config_window(self, code):
        default_config = CodemlBatchAnalysis.MODEL_CONFIGS.get(code, {})
        ModelConfigWindow(self, code, default_config)

    def _open_tree_labeler(self, mode: str = 'branchsite'):
        if self.tree_file is None:
            self.append_log("⚠️ Selecione uma árvore (.nwk) primeiro.\n")
            return
        
        try:
            TreeLabelWindow(self, self.tree_file, mode=mode)
        except Exception as e:
            import traceback
            self.append_log(f"❌ Erro ao abrir TreeLabelWindow: {e}\n")
            self.append_log(f"{traceback.format_exc()}\n")

    def select_input_folder(self):
        path = filedialog.askdirectory()
        if path:
            self.input_folder = Path(path)
            self.label_input.configure(text=str(self.input_folder.name))
            self._update_models_state()

    def select_tree_file(self):
        path = filedialog.askopenfilename(filetypes=[('Tree files', '*.tree *.txt *.nwk')])
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
        """Mostra informações sobre modelos nulos e suas configurações"""
        # Criar janela de informações
        info_window = ctk.CTkToplevel(self)
        info_window.title("Configuração de Modelos Nulos")
        info_window.geometry("700x500")
        info_window.attributes("-topmost", True)
        info_window.grab_set()
        
        # Header
        header = ctk.CTkLabel(info_window, text="📊 Modelos Nulos - Configuração Automática",
                             font=("Roboto", 12, "bold"),
                             text_color=self.COLORS['accent_blue'])
        header.pack(padx=15, pady=15)
        
        # Scrollable frame
        scroll_frame = ctk.CTkScrollableFrame(info_window, fg_color=self.COLORS['bg_feed'],
                                             corner_radius=8)
        scroll_frame.pack(fill="both", expand=True, padx=12, pady=(0, 12))
        
        # Informações sobre cada modelo nulo
        for null_model, config in self.codeml_backend.NEUTRAL_MODELS.items():
            alt_model = config.get('corresponding_alternative', 'N/A')
            reason = config.get('reason', '')
            
            # Card para cada modelo
            card = ctk.CTkFrame(scroll_frame, fg_color=self.COLORS['bg_card'],
                               corner_radius=8, border_width=1,
                               border_color=self.COLORS['bg_card_hover'])
            card.pack(fill="x", padx=8, pady=6)
            
            # Título do modelo
            title = ctk.CTkLabel(card, text=f"🔗 {null_model}", font=("Roboto", 11, "bold"),
                                text_color=self.COLORS['accent_blue'])
            title.pack(anchor="w", padx=12, pady=(8, 4))
            
            # Alternativa correspondente
            alt_text = ctk.CTkLabel(card, 
                                   text=f"▶ Usado como nulo para: {alt_model}",
                                   font=("Roboto", 10),
                                   text_color=self.COLORS['text_secondary'])
            alt_text.pack(anchor="w", padx=20, pady=2)
            
            # Configuração
            config_text = ctk.CTkLabel(card,
                                      text=f"⚙️ Configuração: fix_omega=1, omega=1.0 (FIXADO)",
                                      font=("Roboto", 10),
                                      text_color=self.COLORS['success'])
            config_text.pack(anchor="w", padx=20, pady=2)
            
            # Motivo/Citação
            reason_text = ctk.CTkLabel(card, text=f"📚 {reason}",
                                      font=("Roboto", 9),
                                      text_color=self.COLORS['text_tertiary'],
                                      wraplength=600, justify="left")
            reason_text.pack(anchor="w", padx=20, pady=(2, 8))
        
        # Footer com instrução
        footer = ctk.CTkLabel(info_window, 
                             text="✅ Quando a opção está ATIVADA, estes modelos nulos são adicionados automaticamente",
                             font=("Roboto", 9),
                             text_color=self.COLORS['text_tertiary'],
                             wraplength=650, justify="center")
        footer.pack(padx=12, pady=12)

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
                             text=f"📊 {model_code} - {model_info.get('full_name', '')}",
                             font=("Roboto", 13, "bold"),
                             text_color=self.COLORS['accent_blue'])
        header.pack(padx=15, pady=15)
        
        # Scrollable frame para conteúdo
        scroll_frame = ctk.CTkScrollableFrame(info_window, fg_color=self.COLORS['bg_feed'],
                                             corner_radius=8)
        scroll_frame.pack(fill="both", expand=True, padx=12, pady=(0, 12))
        
        # Tipo de teste
        test_type_label = ctk.CTkLabel(scroll_frame, text="📋 Tipo de Teste:",
                                      font=("Roboto", 11, "bold"),
                                      text_color=self.COLORS['accent_purple'])
        test_type_label.pack(anchor="w", padx=8, pady=(8, 2))
        
        test_type_value = ctk.CTkLabel(scroll_frame, text=model_info.get('test_type', ''),
                                      font=("Roboto", 10),
                                      text_color=self.COLORS['text_primary'],
                                      wraplength=700, justify="left")
        test_type_value.pack(anchor="w", padx=25, pady=(0, 8))
        
        # Parâmetros
        params_label = ctk.CTkLabel(scroll_frame, text="⚙️ Parâmetros:",
                                   font=("Roboto", 11, "bold"),
                                   text_color=self.COLORS['accent_purple'])
        params_label.pack(anchor="w", padx=8, pady=(8, 2))
        
        params_value = ctk.CTkLabel(scroll_frame, text=model_info.get('parameters', ''),
                                   font=("Roboto", 10),
                                   text_color=self.COLORS['text_primary'],
                                   wraplength=700, justify="left")
        params_value.pack(anchor="w", padx=25, pady=(0, 8))
        
        # Propósito
        purpose_label = ctk.CTkLabel(scroll_frame, text="🎯 Propósito:",
                                    font=("Roboto", 11, "bold"),
                                    text_color=self.COLORS['accent_cyan'])
        purpose_label.pack(anchor="w", padx=8, pady=(8, 2))
        
        purpose_value = ctk.CTkLabel(scroll_frame, text=model_info.get('purpose', ''),
                                    font=("Roboto", 10),
                                    text_color=self.COLORS['text_primary'],
                                    wraplength=700, justify="left")
        purpose_value.pack(anchor="w", padx=25, pady=(0, 8))
        
        # Interpretação
        interp_label = ctk.CTkLabel(scroll_frame, text="💡 Interpretação:",
                                   font=("Roboto", 11, "bold"),
                                   text_color=self.COLORS['accent_green'])
        interp_label.pack(anchor="w", padx=8, pady=(8, 2))
        
        interp_value = ctk.CTkLabel(scroll_frame, text=model_info.get('interpretation', ''),
                                   font=("Roboto", 10),
                                   text_color=self.COLORS['text_primary'],
                                   wraplength=700, justify="left")
        interp_value.pack(anchor="w", padx=25, pady=(0, 8))
        
        # Caso de uso
        use_case_label = ctk.CTkLabel(scroll_frame, text="✔️ Quando usar:",
                                     font=("Roboto", 11, "bold"),
                                     text_color=self.COLORS['accent_orange'])
        use_case_label.pack(anchor="w", padx=8, pady=(8, 2))
        
        use_case_value = ctk.CTkLabel(scroll_frame, text=model_info.get('use_case', ''),
                                     font=("Roboto", 10),
                                     text_color=self.COLORS['text_primary'],
                                     wraplength=700, justify="left")
        use_case_value.pack(anchor="w", padx=25, pady=(0, 8))
        
        # Referências
        refs_label = ctk.CTkLabel(scroll_frame, text="📚 Referências:",
                                 font=("Roboto", 11, "bold"),
                                 text_color=self.COLORS['accent_yellow'])
        refs_label.pack(anchor="w", padx=8, pady=(8, 2))
        
        refs_value = ctk.CTkLabel(scroll_frame, text=model_info.get('references', ''),
                                 font=("Roboto", 10, "italic"),
                                 text_color=self.COLORS['text_secondary'],
                                 wraplength=700, justify="left")
        refs_value.pack(anchor="w", padx=25, pady=(0, 8))

    def _open_results_viewer(self):
        if not self.output_folder:
            self.append_log("⚠️ Selecione uma pasta de saída primeiro.\n")
            return
        try:
            ResultsViewerWindow(self, self.output_folder)
        except Exception as e:
            import traceback
            self.append_log(f"❌ Erro ao abrir visualizador: {e}\n")
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
        self.append_log("⚡ ATUALIZANDO RESULTADOS\n")
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
                self.append_log(f"✓ Modelos encontrados: {', '.join(models)}\n\n")
                
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
                    self.append_log("\n✅ ATUALIZAÇÃO CONCLUÍDA COM SUCESSO!\n")
                    self.append_log("="*80 + "\n")
                    for file_type, file_path in generated_files.items():
                        filepath = Path(file_path)
                        size = filepath.stat().st_size if filepath.exists() else 0
                        self.append_log(f"  ✓ {file_type:25s} | {size:,} bytes\n")
                    self.append_log("="*80 + "\n")
                else:
                    self.append_log("\n❌ Nenhum arquivo foi gerado.\n")
            
            except Exception as e:
                self.append_log(f"\n❌ ERRO: {str(e)}\n")
                import traceback
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
            
            # Detectar tipo de mensagem por emoji/palavra-chave com mais precisão
            if any(x in text for x in ["✓", "✅", "SUCESSO", "COMPLETADO", "FINALIZADO", "OK", "SALVO"]):
                tag = "success"
            elif any(x in text for x in ["✗", "❌", "ERRO", "TIMEOUT", "FALHA", "PROBLEMA"]):
                tag = "error"
            elif any(x in text for x in ["⚠️", "AVISO", "CUIDADO", "ATENÇÃO"]):
                tag = "warning"
            elif any(x in text for x in ["ℹ️", "INFO", "INICIANDO", "PROCESSANDO", "📝", "🏷️", "🔴", "🟢"]):
                tag = "info"
            elif any(x in text for x in ["═", "───", "╔", "╚", "║"]):
                tag = "header"
            
            self.log.insert("end", text, tag)
            self.log.see("end")  # Auto-scroll para o final
        
        self.after(0, _append)

    def _poll_stop_count(self):
        if self.analysis_instance:
            cnt = getattr(self.analysis_instance, 'current_stop_count', 0)
            self.stop_label.configure(text=f"⚠️ Stops: {cnt}")
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
        if self.pause_event.is_set():
            self.pause_event.clear()
            self.btn_pause.configure(text="▶ RETOMAR", fg_color="#10b981")
            self.status_indicator.configure(text="● Pausada", text_color=self.COLORS['warning'])
            self.append_log("⏸ Análise pausada.\n")
        else:
            self.pause_event.set()
            self.btn_pause.configure(text="⏸ PAUSAR", fg_color=self.COLORS['warning'])
            self.status_indicator.configure(text="● Executando", text_color=self.COLORS['success'])
            self.append_log("▶ Análise retomada.\n")

    def _stop_analysis(self):
        if not self.analysis_thread or not self.analysis_thread.is_alive():
            self.append_log("⚠️ Nenhuma análise em execução.\n")
            return
        
        self.append_log("🛑 PARANDO ANÁLISE...\n")
        self.stop_event.set()
        self.status_indicator.configure(text="● Parada", text_color=self.COLORS['danger'])
        
        if self.analysis_instance:
            try:
                if hasattr(self.analysis_instance, 'current_process'):
                    process = self.analysis_instance.current_process
                    if process and process.poll() is None:
                        self.append_log("⚡ Terminando processo CODEML...\n")
                        try:
                            process.terminate()
                            process.wait(timeout=3)
                        except:
                            try:
                                process.kill()
                            except:
                                pass
            except Exception as e:
                self.append_log(f"⚠️ Erro ao terminar processo: {e}\n")
        
        if self.pause_event:
            self.pause_event.set()
        if self.manual_event:
            self.manual_event.set()
        if self.manual_all_event:
            self.manual_all_event.set()
        
        self.btn_run.configure(state="normal", fg_color=self.COLORS['success'], 
                              text_color=self.COLORS['text_primary'])
        self.append_log("🛑 Análise interrompida pelo usuário.\n")

    def start_analysis(self):
        selected = [k for k, v in self.model_vars.items() if v.get()]
        if not selected: 
            self.append_log("⚠️ Selecione pelo menos um modelo.\n")
            return
        
        # ═══ AUTO-COMPLETAR MODELOS NULOS ═══
        original_selected = selected.copy()
        include_neutral = self.include_neutral_models.get()
        selected = CodemlBatchAnalysis.auto_complete_null_models(selected, include_neutral=include_neutral)
        
        # Informar ao usuário quais modelos foram auto-adicionados
        if len(selected) > len(original_selected):
            added = set(selected) - set(original_selected)
            self.append_log(f"✅ Modelos auto-adicionados: {', '.join(sorted(added))}\n")
            self.append_log(f"   (Necessários para comparação LRT automática)\n\n")
        
        self.stop_event.clear()
        self.btn_run.configure(state="disabled")
        self.pause_event.set()
        self.manual_all_event.clear()
        self.manual_event.clear()
        
        # Feedback visual: mudar cor do botão e status
        self.btn_run.configure(fg_color=self.COLORS['bg_card'],
                               text_color=self.COLORS['text_muted'])
        self.status_indicator.configure(text="● Executando", text_color=self.COLORS['success'])

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
                self.append_log("\n❌ ERRO: Modelo Branch-Site sem árvore etiquetada!\n")
                self.append_log("   👉 Use 'Marcar Branch-site' na aba Branch-Site.\n\n")
                self.btn_run.configure(state="normal", 
                                      fg_color=self.COLORS['success'],
                                      text_color=self.COLORS['text_primary'])
                return
            
            if needs_branch and not self.tree_branch_labeled:
                self.append_log("\n❌ ERRO: Modelo Branch sem árvore etiquetada!\n")
                self.append_log("   👉 Use 'Marcar Ramos (Múltiplas Tags)' na aba Branch Model.\n\n")
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
                    self.append_log(f"\n✓ Auto-adicionados modelos nulos: {', '.join(set(models_to_run) - set(selected))}\n")
            
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
                'wgs_mode': self.wgs_mode_var.get(),
                'pause_event': self.pause_event,
                'manual_continue_event': self.manual_event,
                'manual_continue_all_event': self.manual_all_event,
                'stop_event': self.stop_event
            }
            
            self.append_log("╔" + "═"*58 + "╗\n")
            self.append_log("║  🧬 INICIANDO ANÁLISE DE SELEÇÃO POSITIVA  🧬" + " "*10 + "║\n")
            self.append_log("╚" + "═"*58 + "╝\n\n")
            self.append_log("🚀 INICIANDO ANÁLISE BATCH\n")
            self.append_log("="*60 + "\n")
            
            analysis.run_batch_analysis()
            
        except Exception as e:
            import traceback
            self.append_log(f"❌ Erro crítico: {e}\n")
            self.append_log(f"{traceback.format_exc()}\n")
        finally:
            self.btn_run.configure(state="normal")
            self.analysis_instance = None
            
            if self.stop_event.is_set():
                self.append_log("\n" + "="*60 + "\n")
                self.append_log("🛑 ANÁLISE INTERROMPIDA\n")
                self.append_log("="*60 + "\n")
            else:
                self.append_log("\n" + "="*60 + "\n")
                self.append_log("✅ ANÁLISE CONCLUÍDA\n")
                self.append_log("="*60 + "\n")


if __name__ == "__main__":
    app = App()
    app.mainloop()
