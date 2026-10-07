import customtkinter as ctk
from tkinter import Canvas
from pathlib import Path
import threading
import time
import traceback
import io
import sys
import os
import signal
import platform as _platform

sys.path.insert(0, str(Path(__file__).parent.parent))

# ── Platform ─────────────────────────────────────────────
_ON_LINUX  = _platform.system() == "Linux"
_ON_WIN    = _platform.system() == "Windows"
# sans-serif font: Roboto on Windows/macOS, DejaVu Sans on Linux
_FONT_UI   = "Roboto" if not _ON_LINUX else "DejaVu Sans"
# monospace font: Cascadia Code on Windows, DejaVu Sans Mono elsewhere
_FONT_MONO = "Cascadia Code" if _ON_WIN else "DejaVu Sans Mono"

from backend.codeml_backend import CodemlBatchAnalysis
from backend import messages as backend_messages
from backend.ctl_params import (CODONFREQ_OPTIONS, DEFAULT_CODONFREQ, DEFAULT_CTL_PARAMS,
                                build_ctl_text, codonfreq_label, parse_codonfreq_label)
from backend.lrt_stats import pairs_for as lrt_pairs_for
from backend.preflight import discover_per_gene_trees, group_by_gene, list_alignment_files, run_preflight
from .results_viewer import ResultsViewerWindow
from .gui_texts import TEXTS, set_language, get_language, tr
from .ui_helpers import (CURRENT_THEME, FONT_SIZE, PALETTE, RADIUS, SPACE, THEME_CHOICES, PreflightDialog,
                         ask_directory, ask_open_file, ask_string, ask_yes_no, os_error_text, save_theme_pref, system_theme,
                         disable_mouse_wheel, fit_to_screen, hover_tint, mix, open_folder,
                         show_about, show_message)

try:
    from Bio import Phylo
except Exception:
    Phylo = None

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt

ctk.set_default_color_theme("blue")

class StdoutRedirect:
    """print() output from libraries goes to the log as 'debug' messages, hidden
    until "Show technical details" is ticked."""
    def __init__(self, append_func):
        self.append = append_func
    def write(self, s):
        for line in str(s).splitlines():
            if line.strip():
                self.append(line, 'debug')
    def flush(self): pass

class ModelConfigWindow(ctk.CTkToplevel):
    """Edit one model's .ctl parameters, with a preview of the whole .ctl."""

    COLORS = {
        'bg_dark':        PALETTE['bg_window'],
        'bg_card':        PALETTE['bg_surface'],
        'bg_input':       PALETTE['bg_inset'],
        'bg_elevated':    PALETTE['bg_elevated'],
        'border':         PALETTE['control_border'],
        'text_primary':   PALETTE['text_primary'],
        'text_secondary': PALETTE['text_secondary'],
    }

    _NUMERIC_FIELDS = (('NSsites', 'cfg_field_nssites'), ('model', 'cfg_field_model'),
                       ('fix_omega', 'cfg_field_fix_omega'), ('omega', 'cfg_field_omega'),
                       ('ncatG', 'cfg_field_ncatg'), ('kappa', 'cfg_field_kappa'))

    def __init__(self, parent, model_code: str, initial: dict):
        super().__init__(parent)
        self.title(TEXTS["model_config_header"].format(model_code=model_code))
        fit_to_screen(self, 1100, 720, min_w=900, min_h=480)
        self.parent = parent
        self.model_code = model_code
        self.entries = {}
        self.configure(fg_color=self.COLORS['bg_dark'])
        self.transient(parent)
        self.after(50, self._grab)
        self.bind("<Escape>", lambda e: self.destroy())

        ctk.CTkLabel(self, text=TEXTS["model_config_header"].format(model_code=model_code),
                     font=(_FONT_UI, FONT_SIZE['lg'], "bold"),
                     text_color=self.COLORS['text_primary']).pack(anchor='w', padx=SPACE['xl'],
                                                                  pady=(SPACE['xl'], SPACE['md']))

        body = ctk.CTkFrame(self, fg_color='transparent')
        body.pack(fill="both", expand=True, padx=SPACE['xl'], pady=(0, SPACE['lg']))
        form = ctk.CTkScrollableFrame(body, fg_color=self.COLORS['bg_card'], width=440,
                                      corner_radius=RADIUS['panel'], border_width=0,
                                      scrollbar_button_color=self.COLORS['border'])
        form.pack(side='left', fill="y", padx=(0, SPACE['lg']))
        side = ctk.CTkFrame(body, fg_color='transparent')
        side.pack(side='left', fill='both', expand=True)

        current = dict(initial)
        current.update(parent.custom_model_params.get(model_code, {}))
        defaults = parent._global_ctl_options()

        def _label(key_text):
            name, hint = TEXTS[key_text]
            ctk.CTkLabel(form, text=name, font=(_FONT_UI, FONT_SIZE['md'], "bold"),
                         text_color=self.COLORS['text_primary']).pack(anchor="w", padx=SPACE['md'],
                                                                      pady=(SPACE['md'], 0))
            ctk.CTkLabel(form, text=hint, font=(_FONT_UI, FONT_SIZE['sm']), wraplength=400,
                         justify='left', text_color=self.COLORS['text_secondary']).pack(
                             anchor="w", padx=SPACE['md'], pady=(0, SPACE['xs']))

        for key, text_key in self._NUMERIC_FIELDS[:4]:
            _label(text_key)
            ent = ctk.CTkEntry(form, fg_color=self.COLORS['bg_input'], border_color=self.COLORS['border'],
                               border_width=1, corner_radius=RADIUS['field'], height=32,
                               text_color=self.COLORS['text_primary'])
            ent.pack(fill="x", padx=SPACE['md'], pady=(0, SPACE['xs']))
            ent.insert(0, str(current.get(key, '')))
            ent.bind('<KeyRelease>', lambda e: self._refresh_preview())
            self.entries[key] = ent

        _label('cfg_field_codonfreq')
        cf_values = [codonfreq_label(v) for v, _, _ in CODONFREQ_OPTIONS]
        self.codonfreq_menu = ctk.CTkOptionMenu(
            form, values=cf_values, command=lambda _v: self._refresh_preview(),
            fg_color=self.COLORS['bg_elevated'], button_color=self.COLORS['border'],
            button_hover_color=PALETTE['control_border_hover'], height=32,
            corner_radius=RADIUS['field'], text_color=self.COLORS['text_primary'])
        self.codonfreq_menu.set(codonfreq_label(current.get('CodonFreq', defaults['CodonFreq'])))
        self.codonfreq_menu.pack(fill='x', padx=SPACE['md'], pady=(0, SPACE['xs']))

        for key, text_key in self._NUMERIC_FIELDS[4:]:
            _label(text_key)
            ent = ctk.CTkEntry(form, fg_color=self.COLORS['bg_input'], border_color=self.COLORS['border'],
                               border_width=1, corner_radius=RADIUS['field'], height=32,
                               text_color=self.COLORS['text_primary'])
            ent.pack(fill="x", padx=SPACE['md'], pady=(0, SPACE['xs']))
            ent.insert(0, str(current.get(key, defaults.get(key, ''))))
            ent.bind('<KeyRelease>', lambda e: self._refresh_preview())
            self.entries[key] = ent

        ctk.CTkLabel(side, text=TEXTS["cfg_preview"], font=(_FONT_UI, FONT_SIZE['md'], "bold"),
                     text_color=self.COLORS['text_primary']).pack(anchor='w', pady=(0, SPACE['xs']))
        self.preview = ctk.CTkTextbox(side, font=(_FONT_MONO, FONT_SIZE['sm']),
                                      fg_color=self.COLORS['bg_input'], corner_radius=RADIUS['card'],
                                      text_color=self.COLORS['text_primary'], wrap='none')
        self.preview.pack(fill='both', expand=True)
        self._refresh_preview()

        btn_frame = ctk.CTkFrame(self, fg_color='transparent')
        btn_frame.pack(fill="x", padx=SPACE['xl'], pady=(0, SPACE['xl']))
        ctk.CTkButton(btn_frame, text=TEXTS["model_config_btn_cancel"], fg_color=PALETTE['neutral_fill'],
                      hover_color=mix(PALETTE['neutral_fill'], '#ffffff', 0.08), command=self.destroy,
                      text_color='#ffffff', height=36,
                      font=(_FONT_UI, FONT_SIZE['md'], "bold"), corner_radius=RADIUS['card']).pack(
                          side="left", fill="x", expand=True, padx=(0, SPACE['sm']))
        ctk.CTkButton(btn_frame, text=TEXTS["model_config_btn_save"], fg_color=PALETTE['accent_fill'],
                      hover_color=mix(PALETTE['accent_fill'], '#000000', 0.2), text_color='#ffffff',
                      command=self._on_save, font=(_FONT_UI, FONT_SIZE['md'], "bold"), height=36,
                      corner_radius=RADIUS['card']).pack(side="left", fill="x", expand=True)

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
            text=TEXTS["model_status_configured"], text_color=PALETTE['warning_text'])
        self.parent.append_log(TEXTS["cfg_saved"].format(model=self.model_code) + "\n", 'ok')
        self.destroy()


class TreeLabelWindow(ctk.CTkToplevel):
    """Window for labelling branches on a rectangular cladogram."""
    
    BG_DARK    = PALETTE['bg_window']
    BG_SIDEBAR = PALETTE['bg_panel']
    BG_CARD    = PALETTE['bg_surface']
    TEXT_PRIMARY   = PALETTE['text_primary']
    TEXT_SECONDARY = PALETTE['text_secondary']
    TEXT_TERTIARY  = PALETTE['text_tertiary']
    ACCENT_BLUE  = PALETTE['accent_text']
    ACCENT_PINK  = '#f472b6'
    SUCCESS      = PALETTE['success_fill']
    SUCCESS_HOVER = mix(PALETTE['success_fill'], '#000000', 0.2)
    DANGER       = PALETTE['danger_text']
    DANGER_HOVER = '#ef4444'
    
    def __init__(self, parent, tree_path: Path | None, mode: str = 'branchsite'):
        super().__init__(parent)
        self.title(tr("Marcar ramos", "Label branches") + f" - {mode}")
        fit_to_screen(self, 1400, 850, min_w=900, min_h=560)
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
                        text_color=self.TEXT_SECONDARY).pack(padx=SPACE['xl'], pady=SPACE['xl'])
            return

        left_frame = ctk.CTkFrame(self, width=280, fg_color=self.BG_SIDEBAR,
                                 corner_radius=0, border_width=0)
        left_frame.pack(side='left', fill='y', padx=0, pady=0)
        left_frame.pack_propagate(False)
        
        plot_frame = ctk.CTkFrame(self, fg_color=self.BG_DARK)
        plot_frame.pack(side='left', fill='both', expand=True, padx=SPACE['md'], pady=SPACE['md'])

        title_label = ctk.CTkLabel(left_frame, text=TEXTS["tree_labeler_sidebar_title"],
                                  font=(_FONT_UI, FONT_SIZE['lg'], "bold"),
                                  text_color=self.TEXT_PRIMARY)
        title_label.pack(pady=(SPACE['xl'], SPACE['sm']), padx=SPACE['lg'], anchor='w')

        if self.mode == 'branchsite':
            instructions = TEXTS["tree_labeler_instructions_branchsite"]
        else:
            instructions = TEXTS["tree_labeler_instructions_branch"]
        
        inst_label = ctk.CTkLabel(left_frame, text=instructions, wraplength=244, justify="left", 
                                 font=(_FONT_UI, FONT_SIZE['sm']), text_color=self.TEXT_SECONDARY)
        inst_label.pack(pady=(0, SPACE['sm']), padx=SPACE['lg'], anchor='w')

        legend_header = ctk.CTkLabel(left_frame, text=TEXTS["tree_labeler_legend_title"],
                                     font=(_FONT_UI, FONT_SIZE['xs'], "bold"),
                                     text_color=self.TEXT_SECONDARY)
        legend_header.pack(pady=(SPACE['lg'], SPACE['xs']), padx=SPACE['lg'], anchor='w')
        
        self.legend_frame = ctk.CTkFrame(left_frame, fg_color=self.BG_CARD, corner_radius=RADIUS['card'])
        self.legend_frame.pack(fill='both', expand=True, padx=SPACE['lg'], pady=SPACE['xs'])

        btn_frame = ctk.CTkFrame(left_frame, fg_color='transparent')
        btn_frame.pack(side='bottom', fill='x', padx=SPACE['lg'], pady=SPACE['lg'])
        
        ctk.CTkButton(btn_frame, text=TEXTS["tree_labeler_btn_save"], fg_color=self.SUCCESS, hover_color=self.SUCCESS_HOVER,
                     command=self._on_save, height=40, font=(_FONT_UI, FONT_SIZE['md'], "bold"),
                     text_color='#ffffff', corner_radius=RADIUS['card']).pack(fill='x', pady=(0, SPACE['sm']))
        ctk.CTkButton(btn_frame, text=TEXTS["tree_labeler_btn_cancel"], fg_color=PALETTE['neutral_fill'],
                     hover_color=mix(PALETTE['neutral_fill'], '#ffffff', 0.08),
                     command=self.destroy, height=40, font=(_FONT_UI, FONT_SIZE['md'], "bold"),
                     text_color='#ffffff', corner_radius=RADIUS['card']).pack(fill='x')

        if self.tree_path is None:
            ctk.CTkLabel(plot_frame, text=TEXTS["tree_err_no_tree"], font=(_FONT_UI, FONT_SIZE['lg']),
                         text_color=self.TEXT_SECONDARY).pack(pady=48)
            return

        try:
            self.tree = Phylo.read(str(self.tree_path), 'newick')
        except Exception as e:
            ctk.CTkLabel(plot_frame, text=TEXTS["tree_err_load"].format(error=e), font=(_FONT_UI, FONT_SIZE['sm']),
                         text_color=PALETTE['danger_text']).pack(pady=48)
            return

        self.fig = Figure(figsize=(12, 10), dpi=110)
        self.fig.patch.set_facecolor(self.BG_DARK)
        self.ax = self.fig.add_subplot(111)
        self.fig.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)
        self.ax.set_facecolor(self.BG_DARK)
        
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
                # for labelling only the topology matters; branch lengths are ignored
                calc_depth_with_lengths(child, accumulated_depth + 1.0)
        
        calc_depth_with_lengths(self.tree.root, 0.0)

        # every tip in the rightmost column
        max_depth = max(
            (depths[t] for t in terminals if t in depths),
            default=1.0
        )
        self._cladogram_depth = max_depth

        for clade in self.tree.find_clades(order='postorder'):
            if clade.is_terminal():
                x = max_depth
                y = terminal_y_map[clade]
            else:
                x = depths.get(clade, 0.0)
                child_ys = [self.clade_positions[child][1] for child in clade.clades
                            if child in self.clade_positions]
                y = sum(child_ys) / len(child_ys) if child_ys else 0.0

            self.clade_positions[clade] = (x, y)

    def _draw_tree(self):
        """Draw the cladogram, colouring from each labelled node down to its tips; a
        labelled descendant overrides its ancestor."""
        self.ax.clear()
        self.ax.set_facecolor(self.BG_DARK)
        
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        for spine in self.ax.spines.values():
            spine.set_visible(False)
        
        self.scatter_objects.clear()

        def get_tag_for_branch(clade):
            """Label that colours this clade: its own, else the nearest labelled ancestor,
            else None."""
            if clade in self.marked_clades:
                return self.marked_clades[clade]
            
            try:
                closest_tag = None
                min_distance = float('inf')
                
                for marked_clade, tag in self.marked_clades.items():
                    try:
                        all_descendants = list(marked_clade.find_clades())
                        if clade in all_descendants:
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
            """True if this clade descends from a labelled node."""
            tag = get_tag_for_branch(clade)
            return (tag is not None, tag)

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
            
            parent_tag = get_tag_for_branch(clade)
            
            if parent_tag:
                vertical_color = self._get_tag_color(parent_tag)
                vertical_width = 2.8
                vertical_alpha = 1.0
            else:
                vertical_color = PALETTE['canvas_line']
                vertical_width = 1.2
                vertical_alpha = 0.7
            
            if abs(y_max - y_min) > 0.01:
                self.ax.plot([x_parent, x_parent], [y_min, y_max],
                           color=vertical_color, linewidth=vertical_width, 
                           zorder=1, alpha=vertical_alpha, solid_capstyle='round')
            
            for child, (x_child, y_child) in children_positions:
                branch_tag = get_tag_for_branch(child)
                
                if branch_tag:
                    branch_color = self._get_tag_color(branch_tag)
                    branch_width = 2.8
                    branch_alpha = 1.0
                else:
                    branch_color = PALETTE['canvas_line']
                    branch_width = 1.2
                    branch_alpha = 0.7

                self.ax.plot([x_parent, x_child], [y_child, y_child],
                           color=branch_color, linewidth=branch_width, 
                           zorder=2, alpha=branch_alpha, solid_capstyle='round')

        for clade in self.tree.find_clades():
            if clade not in self.clade_positions:
                continue
            
            x, y = self.clade_positions[clade]
            
            is_marked, tag = is_descendant_of_marked(clade)
            color = self._get_tag_color(tag) if is_marked else PALETTE['canvas_line']
            
            if clade.is_terminal():
                size = 90 if is_marked else 45
                edge_width = 2.2 if is_marked else 1.0
            else:
                size = 55 if is_marked else 28
                edge_width = 1.8 if is_marked else 0.8
            
            edge_color = PALETTE['canvas_text'] if is_marked else PALETTE['canvas_line']
            
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
                text_color = PALETTE['canvas_text']
                weight = 'normal'
                fontsize = 10
            
            label_offset = getattr(self, '_cladogram_depth', 10) * 0.03
            self.ax.text(x + label_offset, y, name,
                        va='center', ha='left',
                        fontsize=fontsize,
                        color=text_color,
                        weight=weight,
                        zorder=20,
                        family='monospace')

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
            # right margin sized to the longest name so labels are not cut off
            right_margin = max(x_range * (0.25 + 0.012 * max_label_len), 0.003)
            
            self.ax.set_xlim(min_x - x_range * 0.03, max_x + right_margin)
            self.ax.set_ylim(min_y - 2.5, max_y + 2.5)

        self.canvas.draw_idle()

    def _on_pick(self, event):
        """Label the clicked node."""
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
                response = ask_string(self, TEXTS["tag_dialog_edit_title"],
                                      TEXTS["tag_dialog_edit_prompt"].format(tag=current_tag))
                
                if response:
                    response = response.strip().lower()
                    if response in ('remove', 'remover'):
                        self._remove_tag_recursively(clade)
                        self.parent.append_log(tr(f"Marca {current_tag} removida de {clade_name}", f"Tag {current_tag} removed from {clade_name}") + "\n")
                    elif response.isdigit():
                        self._remove_tag_recursively(clade)
                        new_tag = f"#{response}"
                        self._apply_tag_recursively(clade, new_tag)
                        self.marked_clades[clade] = new_tag
                        self.parent.append_log(tr(f"Marca alterada para {new_tag} em {clade_name}", f"Tag changed to {new_tag} on {clade_name}") + "\n")
            else:
                response = ask_string(self, TEXTS["tag_dialog_new_title"], TEXTS["tag_dialog_new_prompt"])
                
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
        """Readable name of a clade."""
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
        """Colour of a label number."""
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
        """Remove a label everywhere and redraw."""
        for clade in [c for c, t in list(self.clade_tags.items()) if t == tag]:
            self.clade_tags.pop(clade, None)
        for clade in [c for c, t in list(self.marked_clades.items()) if t == tag]:
            self.marked_clades.pop(clade, None)
        self._draw_tree()
        self._refresh_legend()
        self.parent.append_log(tr(f"Marca {tag} removida.", f"Tag {tag} removed.") + "\n")

    def _refresh_legend(self):
        """Refresh the legend of active labels."""
        for widget in self.legend_frame.winfo_children():
            widget.destroy()

        active_tags = sorted(set(self.clade_tags.values()))

        if not active_tags:
            ctk.CTkLabel(self.legend_frame, text=TEXTS["tree_labeler_no_tags"],
                         text_color=self.TEXT_TERTIARY, font=(_FONT_UI, FONT_SIZE['sm'])
                         ).pack(anchor='w', padx=SPACE['md'], pady=SPACE['md'])
        else:
            for tag in active_tags:
                color = self._get_tag_color(tag)

                row_frame = ctk.CTkFrame(self.legend_frame, fg_color="transparent")
                row_frame.pack(fill='x', padx=SPACE['sm'], pady=SPACE['xs'])

                color_box = Canvas(row_frame, width=22, height=16, highlightthickness=0)
                color_box.configure(bg=self.BG_CARD)
                color_box.create_rectangle(2, 2, 20, 14, fill=color, outline=PALETTE['canvas_line'], width=1)
                color_box.pack(side='left', padx=SPACE['xs'])

                tag_label = ctk.CTkLabel(row_frame, text=tag,
                                         font=(_FONT_UI, FONT_SIZE['sm'], "bold"),
                                         text_color=color)
                tag_label.pack(side='left', padx=(0, SPACE['xs']))

                count = sum(1 for t in self.clade_tags.values() if t == tag)
                ctk.CTkLabel(row_frame, text=f"({count})",
                             font=(_FONT_UI, FONT_SIZE['sm']),
                             text_color=self.TEXT_TERTIARY).pack(side='left')

                del_btn = ctk.CTkButton(
                    row_frame, text="X", width=24, height=24,
                    fg_color='transparent',
                    hover_color=hover_tint(self.DANGER, self.BG_CARD),
                    text_color=self.TEXT_TERTIARY,
                    border_width=0,
                    corner_radius=RADIUS['field'],
                    font=(_FONT_UI, FONT_SIZE['sm'], "bold"),
                    command=lambda t=tag: self._delete_tag(t)
                )
                del_btn.pack(side='right', padx=(0, SPACE['xs']))

    def _on_save(self):
        """Save the labelled tree as Newick."""
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
    COLORS = {
        'bg_darkest':     PALETTE['bg_window'],
        'bg_dark':        PALETTE['bg_window'],
        'bg_sidebar':     PALETTE['bg_panel'],
        'bg_card':        PALETTE['bg_surface'],
        'bg_feed':        PALETTE['bg_elevated'],
        'bg_card_hover':  PALETTE['bg_elevated'],
        'bg_hover':       PALETTE['bg_elevated_hover'],
        'bg_input':       PALETTE['bg_inset'],

        'text_primary':   PALETTE['text_primary'],
        'text_secondary': PALETTE['text_secondary'],
        'text_tertiary':  PALETTE['text_tertiary'],
        'text_muted':     PALETTE['text_muted'],

        'accent_blue':        PALETTE['accent_blue'],
        'accent_blue_hover':  mix(PALETTE['accent_fill'], '#000000', 0.2),
        'accent_blue_light':  PALETTE['accent_text'],
        'accent_fill':        PALETTE['accent_fill'],

        'accent_cyan':        PALETTE['accent_cyan'],
        'accent_cyan_hover':  PALETTE['info_fill'],
        'accent_purple':      PALETTE['accent_purple'],
        'accent_purple_hover':'#7c3aed',
        'accent_pink':        '#f472b6',
        'accent_pink_hover':  '#db2777',

        'success':        PALETTE['success_text'],
        'success_hover':  '#16a34a',
        'success_light':  PALETTE['success_light'],
        'warning':        PALETTE['warning_text'],
        'warning_hover':  '#d97706',
        'danger':         PALETTE['danger_text'],
        'danger_hover':   '#ef4444',
        'info':           PALETTE['info_text'],
        'info_hover':     '#06b6d4',

        'border':         PALETTE['control_border'],
        'border_hover':   PALETTE['control_border_hover'],
        'divider':        PALETTE['divider'],
    }

    SIDEBAR_WIDTH = 320

    def __init__(self):
        super().__init__()
        self.withdraw()          # shown complete at the end, not piece by piece
        from backend.version import __version__
        self.title(f"EasyPAML {__version__}")
        fit_to_screen(self, 1400, 850)
        
        
        self.configure(fg_color=self.COLORS['bg_dark'])

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
        
        self.include_neutral_models = ctk.BooleanVar(value=True)
        self._excluded_nulls = set()     # automatic nulls the user switched off

        _max_cores = CodemlBatchAnalysis.available_cores()
        # default: half the cores, so the computer stays usable
        self.cores_var = ctk.IntVar(value=max(1, _max_cores // 2))
        self.ignore_stop_codons_var = ctk.BooleanVar(value=False)
        self.auto_prune_tree_var = ctk.BooleanVar(value=True)

        self.tree_branch_labeled = None
        self.tree_branchsite_labeled = None

        self._visual_sig = None
        self._gene_count_cache = {}
        self._slots = {}
        self._tiles = {}
        self._tile_parts = {}
        self._auto_nulls = {}
        self._off_nulls = {}

        C = self.COLORS
        sp = SPACE
        fs = FONT_SIZE
        side_inner = self.SIDEBAR_WIDTH - 28 - sp['lg'] - sp['sm']

        self.sidebar = ctk.CTkFrame(self, width=self.SIDEBAR_WIDTH, corner_radius=0,
                                    fg_color=C['bg_sidebar'], border_width=0)
        self.sidebar.pack(side="left", fill="y", padx=0, pady=0)
        self.sidebar.pack_propagate(False)

        # ── Logo, name and version ──────────────────────
        logo_frame = ctk.CTkFrame(self.sidebar, fg_color='transparent')
        logo_frame.pack(fill="x", padx=sp['lg'], pady=(sp['lg'], 0))

        badge = ctk.CTkFrame(logo_frame, fg_color=C['accent_fill'],
                             width=32, height=32, corner_radius=RADIUS['card'])
        badge.pack(side='left', padx=(0, sp['md']))
        badge.pack_propagate(False)
        ctk.CTkLabel(badge, text="EP", font=(_FONT_UI, fs['sm'], "bold"),
                     text_color='#ffffff').pack(expand=True)

        title_col = ctk.CTkFrame(logo_frame, fg_color='transparent')
        title_col.pack(side='left', anchor='center')
        ctk.CTkLabel(title_col, text=TEXTS["app_sidebar_title"],
                     font=(_FONT_UI, fs['lg'], "bold"), height=20,
                     text_color=C['text_primary']).pack(anchor='w')
        ctk.CTkLabel(logo_frame, text=f"v{__version__}", font=(_FONT_UI, fs['xs']),
                     fg_color=C['bg_card_hover'], corner_radius=RADIUS['field'], height=20,
                     text_color=C['text_secondary']).pack(side='right', ipadx=sp['xs'])

        # ── Scrollable steps ──────────────────────────
        _sb = ctk.CTkScrollableFrame(self.sidebar, fg_color='transparent', corner_radius=0,
                                     scrollbar_button_color=C['border'],
                                     scrollbar_button_hover_color=C['border_hover'])
        _sb.pack(fill='both', expand=True, padx=0, pady=(sp['md'], 0))
        self._sidebar_scroll = _sb

        # ── Footer (outside the scroll area) ──────────────────────
        self._build_lang_footer()

        def _step(parent, number, title, first=False):
            """Step header: numbered circle (✓ when done) and title."""
            sec = ctk.CTkFrame(parent, fg_color='transparent')
            sec.pack(fill='x', padx=(sp['lg'], sp['sm']), pady=(sp['sm'] if first else sp['xl'], 0))
            hdr = ctk.CTkFrame(sec, fg_color='transparent')
            hdr.pack(fill='x')
            circle = ctk.CTkLabel(hdr, text=str(number) if number else "", width=20, height=20,
                                  corner_radius=10, fg_color=C['bg_card_hover'],
                                  font=(_FONT_UI, fs['xs'], "bold"), text_color=C['text_secondary'])
            if number:
                circle.pack(side='left', padx=(0, sp['sm']))
            ctk.CTkLabel(hdr, text=title.upper(), font=(_FONT_UI, fs['xs'], "bold"), height=20,
                         text_color=C['text_secondary'], anchor='w').pack(side='left')
            side = ctk.CTkLabel(hdr, text="", font=(_FONT_UI, fs['xs']), height=20,
                                text_color=C['text_tertiary'], anchor='e')
            side.pack(side='right')
            inner = ctk.CTkFrame(sec, fg_color='transparent')
            inner.pack(fill='x', pady=(sp['sm'], 0))
            return inner, circle, side

        def _obtn(parent, text, cmd, **kw):
            """Outlined secondary button."""
            return ctk.CTkButton(
                parent, text=text, command=cmd,
                fg_color=C['bg_card'],
                hover_color=C['bg_card_hover'],
                text_color=C['text_primary'],
                text_color_disabled=C['text_tertiary'],
                border_width=1, border_color=C['border'],
                font=(_FONT_UI, fs['md'], "bold"), height=36,
                corner_radius=RADIUS['card'], **kw)

        def _help(parent, title_key, hint_key, fill=None, command=None):
            """The '?' help button."""
            fill = fill or C['bg_card_hover']
            return ctk.CTkButton(
                parent, text="?", width=20, height=20, corner_radius=10,
                font=(_FONT_UI, fs['xs'], "bold"), fg_color=fill,
                hover_color=mix(fill, '#ffffff', 0.08),
                text_color=C['text_secondary'],
                command=command or (lambda: self._show_help(TEXTS[title_key], TEXTS[hint_key])))

        def _switch(parent, variable, color, **kw):
            return ctk.CTkSwitch(
                parent, text="", variable=variable, onvalue=True, offvalue=False,
                switch_width=36, switch_height=20, width=36,
                progress_color=color, button_color=PALETTE['switch_knob'],
                button_hover_color=PALETTE['switch_knob'], fg_color=PALETTE['switch_track'], **kw)

        def _slot(parent, title_key, command):
            """File card; the whole card opens the same chooser as its button."""
            card = ctk.CTkFrame(parent, fg_color=C['bg_card'], corner_radius=RADIUS['card'],
                                border_width=1, border_color=C['border'])
            card.pack(fill='x', pady=(0, sp['sm']))
            col = ctk.CTkFrame(card, fg_color='transparent')
            col.pack(side='left', fill='x', expand=True, pady=sp['md'], padx=sp['md'])
            title = ctk.CTkLabel(col, text=TEXTS[title_key], font=(_FONT_UI, fs['sm'], "bold"),
                                 anchor='w', height=16, text_color=C['text_primary'])
            title.pack(fill='x')
            row = ctk.CTkFrame(col, fg_color='transparent')
            row.pack(fill='x', pady=(sp['xs'], 0))
            btn = ctk.CTkButton(row, text=TEXTS["slot_choose"], command=command, width=64, height=24,
                                corner_radius=RADIUS['field'], font=(_FONT_UI, fs['xs'], "bold"),
                                fg_color='transparent', hover_color=C['bg_hover'],
                                text_color=C['text_secondary'], border_width=1,
                                border_color=C['border'])
            btn.pack(side='right', anchor='s')
            detail = ctk.CTkLabel(card, text=TEXTS["label_not_selected"])
            shown = ctk.CTkLabel(row, text=TEXTS["label_not_selected"], font=(_FONT_UI, fs['xs']),
                                 anchor='w', justify='left', text_color=C['text_tertiary'],
                                 wraplength=150)
            shown.pack(side='left', fill='x', expand=True)
            row.bind('<Configure>', lambda e, l=shown: l.configure(
                wraplength=max(100, e.width - 64 - sp['sm'])), add='+')
            self._make_clickable(card, (card, col, title, row, shown), command,
                                 lambda: self._slots and self._apply_visual_state(force=True))
            return card, title, btn, detail, shown

        # ── 1 · Data ────────────────────────────────────────────────
        fi, self._step1_circle, _ = _step(_sb, 1, TEXTS["step_data"], first=True)
        for key, title_key, cmd in (
                ('input', "btn_input_folder", self.select_input_folder),
                ('tree', "btn_tree_file", self.select_tree_file),
                ('output', "btn_output_folder", self.select_output_folder)):
            card, title, btn, detail, shown = _slot(fi, title_key, cmd)
            self._slots[key] = (card, title, btn, detail, shown)
            setattr(self, f"label_{key}", detail)
            setattr(self, f"btn_{key}", btn)

        # ── 2 · Models ─────────
        mi, self._step2_circle, _ = _step(_sb, 2, TEXTS["step_models"])
        self._step2_text = ctk.CTkLabel(mi, text=TEXTS["step_models_none"], font=(_FONT_UI, fs['sm']),
                                        anchor='w', text_color=C['text_tertiary'])
        self._step2_text.pack(fill='x', padx=(sp['xs'], 0))

        # ── 3 · Advanced settings ────────────────────────────────────────
        ci, _c3, side3 = _step(_sb, 3, TEXTS["step_settings"])
        side3.destroy()
        sec3, hdr3 = ci.master, ci.master.winfo_children()[0]
        ci.pack_forget()
        self._settings_body = ci
        srow = ctk.CTkFrame(sec3, fg_color='transparent')
        srow.pack(fill='x', pady=(sp['sm'], 0))
        self._settings_btn = ctk.CTkButton(srow, text=TEXTS["settings_show"], width=64, height=20,
                                           corner_radius=RADIUS['field'], fg_color='transparent',
                                           hover_color=C['bg_card_hover'], text_color=C['accent_blue_light'],
                                           border_width=1, border_color=C['border'],
                                           font=(_FONT_UI, fs['xs'], "bold"),
                                           command=self._toggle_settings)
        self._settings_summary = ctk.CTkLabel(srow, text="", font=(_FONT_UI, fs['sm']), anchor='w',
                                              justify='left', wraplength=side_inner - sp['xs'],
                                              text_color=C['text_tertiary'])
        self._settings_summary.pack(fill='x', padx=(sp['xs'], 0))
        self._settings_btn.pack(anchor='w', pady=(sp['sm'], 0))
        self._make_clickable(hdr3, [hdr3, srow, self._settings_summary] + [
            w for w in hdr3.winfo_children() if not isinstance(w, ctk.CTkButton)], self._toggle_settings)
        CTRL_W = 72
        ci.grid_columnconfigure(0, weight=1)
        ci.grid_columnconfigure(1, minsize=20 + sp['sm'])
        ci.grid_columnconfigure(2, minsize=CTRL_W)
        label_wrap = side_inner - CTRL_W - 20 - 2 * sp['sm']
        _row = [0]

        def _setting_label(text_key, gap):
            lbl = ctk.CTkLabel(ci, text=TEXTS[text_key], font=(_FONT_UI, fs['sm']),
                               anchor='w', justify='left', wraplength=label_wrap,
                               text_color=C['text_secondary'])
            lbl.grid(row=_row[0], column=0, sticky='w', pady=(gap, 0))
            return lbl

        def _place_help(widget, gap):
            widget.grid(row=_row[0], column=1, sticky='e', padx=(sp['sm'], sp['sm']), pady=(gap, 0))

        def _place_ctrl(widget, gap):
            widget.grid(row=_row[0], column=2, sticky='e', pady=(gap, 0))

        def _place_full(widget, gap=sp['sm']):
            _row[0] += 1
            widget.grid(row=_row[0], column=0, columnspan=3, sticky='ew', pady=(gap, 0))

        _setting_label("label_omega_initial", 0)
        self.entry_omega = ctk.CTkEntry(ci, placeholder_text="0.5", width=CTRL_W,
                                        fg_color=C['bg_input'],
                                        border_color=C['border'],
                                        border_width=1,
                                        corner_radius=RADIUS['field'],
                                        text_color=C['text_primary'],
                                        justify='right', height=32)
        self.entry_omega.insert(0, "0.5")
        _place_ctrl(self.entry_omega, 0)
        _place_help(_help(ci, "label_omega_initial", "label_omega_initial_hint"), 0)

        _row[0] += 1
        _setting_label("label_codonfreq", sp['md'])
        _place_help(_help(ci, "label_codonfreq", "label_codonfreq_hint"), sp['md'])
        self.codonfreq_var = ctk.StringVar(value=codonfreq_label(DEFAULT_CODONFREQ))
        _place_full(ctk.CTkOptionMenu(ci, variable=self.codonfreq_var,
                                      values=[codonfreq_label(v) for v, _, _ in CODONFREQ_OPTIONS],
                                      fg_color=C['bg_card_hover'], button_color=C['border'],
                                      button_hover_color=C['border_hover'],
                                      dropdown_fg_color=C['bg_card_hover'],
                                      dropdown_hover_color=C['bg_hover'],
                                      dropdown_text_color=C['text_primary'],
                                      text_color=C['text_primary'], height=32,
                                      corner_radius=RADIUS['field'],
                                      font=(_FONT_UI, fs['sm'])))

        _row[0] += 1
        _setting_label("label_ncatg", sp['md'])
        self.entry_ncatg = ctk.CTkEntry(ci, fg_color=C['bg_input'], width=CTRL_W,
                                        border_color=C['border'], border_width=1,
                                        corner_radius=RADIUS['field'], text_color=C['text_primary'],
                                        justify='right', height=32)
        self.entry_ncatg.insert(0, str(DEFAULT_CTL_PARAMS['ncatG']))
        _place_ctrl(self.entry_ncatg, sp['md'])
        _place_help(_help(ci, "label_ncatg", "label_ncatg_hint"), sp['md'])

        _row[0] += 1
        _setting_label("label_timeout", sp['md'])
        self.entry_timeout = ctk.CTkEntry(ci, placeholder_text=TEXTS["label_timeout_auto"],
                                          width=CTRL_W, fg_color=C['bg_input'],
                                          border_color=C['border'], border_width=1,
                                          corner_radius=RADIUS['field'],
                                          text_color=C['text_primary'],
                                          placeholder_text_color=C['text_muted'],
                                          justify='right', height=32)
        _place_ctrl(self.entry_timeout, sp['md'])
        _place_help(_help(ci, "label_timeout", "label_timeout_hint"), sp['md'])

        self.cleandata_var = ctk.BooleanVar(value=True)
        _row[0] += 1
        _setting_label("label_remove_gaps", sp['md'])
        self.cb_cleandata = _switch(ci, self.cleandata_var, C['success'])
        _place_ctrl(self.cb_cleandata, sp['md'])
        _place_help(_help(ci, "label_remove_gaps", "label_remove_gaps_hint"), sp['md'])

        _row[0] += 1
        _setting_label("label_cpus", sp['md'])
        _max = CodemlBatchAnalysis.available_cores()
        self.cores_disp = ctk.CTkLabel(
            ci, text=f"{int(self.cores_var.get())}×",
            font=(_FONT_UI, fs['sm'], "bold"),
            text_color=C['accent_blue_light'], anchor='e',
            width=CTRL_W)
        _place_ctrl(self.cores_disp, sp['md'])
        _place_help(_help(ci, "label_cpus", "label_cpus_hint"), sp['md'])
        self.cores_slider = ctk.CTkSlider(
            ci, from_=1, to=max(2, _max),
            number_of_steps=max(1, _max - 1),
            variable=self.cores_var,
            command=self._update_cores_label,
            height=16,
            fg_color=C['border'],
            button_color=C['accent_blue'],
            button_hover_color=C['accent_blue_light'],
            progress_color=C['accent_blue'])
        _place_full(self.cores_slider, sp['xs'])
        disable_mouse_wheel(self.cores_slider)
        _place_full(ctk.CTkLabel(ci, text=TEXTS["label_cpus_detected"].format(n=_max),
                                 font=(_FONT_UI, fs['xs']), anchor='w', height=16,
                                 text_color=C['text_muted']), sp['xs'])

        _row[0] += 1
        _setting_label("label_ignore_stops", sp['md'])
        self.cb_ignore_stops = _switch(ci, self.ignore_stop_codons_var, C['accent_cyan'])
        _place_ctrl(self.cb_ignore_stops, sp['md'])
        _place_help(_help(ci, "label_ignore_stops", "label_ignore_stops_hint"), sp['md'])

        _row[0] += 1
        _setting_label("label_auto_prune", sp['md'])
        self.cb_auto_prune = _switch(ci, self.auto_prune_tree_var, C['accent_cyan'])
        _place_ctrl(self.cb_auto_prune, sp['md'])
        _place_help(_help(ci, "label_auto_prune", "label_auto_prune_hint"), sp['md'])

        ctk.CTkFrame(_sb, fg_color='transparent', height=sp['xl']).pack(fill='x')

        # ── Main area ──
        self.main_frame = ctk.CTkFrame(self, fg_color=C['bg_dark'], corner_radius=0)
        self.main_frame.pack(side="right", fill="both", expand=True, padx=sp['xl'], pady=sp['lg'])

        # ── Header ────────────────
        top = ctk.CTkFrame(self.main_frame, fg_color='transparent')
        top.pack(fill='x', pady=(0, sp['md']))
        tcol = ctk.CTkFrame(top, fg_color='transparent')
        tcol.pack(side='left', fill='x', expand=True)
        ctk.CTkLabel(tcol, text=TEXTS["app_main_title"], font=(_FONT_UI, fs['xxl'], "bold"),
                     height=24, anchor='w', text_color=C['text_primary']).pack(fill='x')
        ctk.CTkLabel(tcol, text=TEXTS["app_main_subtitle"], font=(_FONT_UI, fs['sm']),
                     height=16, anchor='w', text_color=C['text_tertiary']).pack(fill='x')

        self.stop_label = ctk.CTkLabel(
            top, text=TEXTS["status_stops_template"].format(n="–"),
            font=(_FONT_UI, fs['sm'], "bold"), height=28, corner_radius=12,
            fg_color=C['bg_card'], text_color=C['danger'])
        self.stop_label.pack(side='right', ipadx=sp['sm'])
        self.status_indicator = ctk.CTkLabel(
            top, text=TEXTS["status_ready"],
            font=(_FONT_UI, fs['sm'], "bold"), height=28, corner_radius=12,
            fg_color=mix(C['bg_dark'], C['text_tertiary'], 0.16),
            text_color=C['text_tertiary']
        )
        self.status_indicator.pack(side='right', padx=(0, sp['sm']), ipadx=sp['sm'])

        self.files_hint = ctk.CTkLabel(self.main_frame, text="⚠  " + TEXTS["hint_select_files"],
                                       font=(_FONT_UI, fs['sm']), anchor='w', justify='left',
                                       height=32, corner_radius=RADIUS['card'],
                                       fg_color=mix(C['bg_dark'], C['warning'], 0.12),
                                       text_color=C['warning'], wraplength=900)
        self.files_hint.pack(fill='x', pady=(0, sp['md']), ipadx=sp['md'])

        # ── Model cards ──────────────────────────
        self.tabs = ctk.CTkTabview(self.main_frame, fg_color=C['bg_dark'],
                                   segmented_button_fg_color=C['bg_card'],
                                   segmented_button_selected_color=C['border_hover'],
                                   segmented_button_selected_hover_color=C['border_hover'],
                                   segmented_button_unselected_color=C['bg_card'],
                                   segmented_button_unselected_hover_color=C['bg_card_hover'],
                                   text_color=C['text_primary'], anchor='nw',
                                   corner_radius=RADIUS['card'], border_width=0, height=224)
        self.tabs.pack(fill="x", padx=0, pady=(0, sp['md']))
        self.tabs.add(TEXTS["tab_site_models"])
        self.tabs.add(TEXTS["tab_branch_model"])
        self.tabs.add(TEXTS["tab_branchsite"])

        models_hdr = ctk.CTkFrame(self.tabs, fg_color=C['bg_dark'], corner_radius=0)
        models_hdr.place(relx=1.0, x=0, y=0, anchor='ne')
        self._models_count = ctk.CTkLabel(models_hdr, text="", font=(_FONT_UI, fs['xs'], "bold"),
                                          height=28, text_color=C['text_tertiary'])
        self._models_count.pack(side='left', padx=(0, sp['lg']))
        ctk.CTkLabel(models_hdr, text=TEXTS["label_neutral_models"],
                     font=(_FONT_UI, fs['sm']), height=28,
                     text_color=C['text_secondary']).pack(side='left')
        help_btn = _help(models_hdr, None, None, command=self._show_neutral_models_info)
        help_btn.pack(side="left", padx=sp['sm'])
        neutral_sw = _switch(models_hdr, self.include_neutral_models, C['accent_blue'])
        neutral_sw.pack(side='left')

        self.model_vars = {}
        self.model_ctl_labels = {}
        self.model_checkboxes = {}
        self.model_gear_buttons = {}

        self._setup_model_list()

        # ── Run panel ─────
        self.ctrl_frame = ctk.CTkFrame(self.main_frame, fg_color=C['bg_card'],
                                       corner_radius=RADIUS['panel'], border_width=0)
        self.ctrl_frame.pack(fill="x", padx=0, pady=(0, sp['md']))

        run_row = ctk.CTkFrame(self.ctrl_frame, fg_color='transparent')
        run_row.pack(fill='x', padx=sp['lg'], pady=(sp['md'], sp['sm']))

        self.btn_run = ctk.CTkButton(
            run_row, text=TEXTS["btn_run"], command=self.start_analysis,
            fg_color=C['accent_fill'], hover_color=C['accent_blue_hover'],
            text_color='#ffffff',
            text_color_disabled=mix('#ffffff', C['accent_fill'], 0.45),
            border_width=0, font=(_FONT_UI, fs['lg'], "bold"),
            width=176, height=48, corner_radius=RADIUS['card'])
        self.btn_run.pack(side="left")

        def _action_btn(parent, text, cmd, color):
            return ctk.CTkButton(
                parent, text=text, command=cmd,
                fg_color='transparent',
                hover_color=hover_tint(color, C['bg_card']),
                text_color=color,
                text_color_disabled=C['text_tertiary'],
                border_width=1, border_color=C['border'],
                font=(_FONT_UI, fs['sm'], "bold"),
                height=32, corner_radius=RADIUS['card'])

        sum_col = ctk.CTkFrame(run_row, fg_color='transparent')
        sum_col.pack(side='left', fill='x', expand=True, padx=(sp['lg'], 0))
        self._summary_label = ctk.CTkLabel(sum_col, text="", font=(_FONT_UI, fs['md'], "bold"),
                                           anchor='w', justify='left', wraplength=640,
                                           text_color=C['text_primary'])
        self._summary_label.pack(fill='x')
        self._summary_models = ctk.CTkLabel(sum_col, text="", font=(_FONT_UI, fs['sm']),
                                            anchor='w', justify='left', wraplength=640,
                                            text_color=C['text_secondary'])
        self._summary_models.pack(fill='x', pady=(sp['xs'], 0))
        sum_col.bind('<Configure>', lambda e: [lbl.configure(wraplength=max(240, e.width - sp['sm']))
                                               for lbl in (self._summary_label, self._summary_models)],
                     add='+')

        prog_row = ctk.CTkFrame(self.ctrl_frame, fg_color='transparent')
        prog_row.pack(fill='x', padx=sp['lg'], pady=(0, sp['sm']))
        self.btn_open_output = _action_btn(prog_row, TEXTS["btn_open_output"],
                                           self._open_output_folder, C['text_primary'])
        self.btn_open_output.configure(hover_color=C['bg_card_hover'])
        self.btn_open_output.pack(side="right")
        self.btn_view_results_main = _action_btn(prog_row, TEXTS["btn_view_results"],
                                                 self._open_results_viewer, C['accent_blue_light'])
        self.btn_view_results_main.configure(hover_color=C['bg_card_hover'])
        self.btn_view_results_main.pack(side="right", padx=(0, sp['sm']))
        self.btn_stop = _action_btn(prog_row, TEXTS["btn_stop"],
                                    self._stop_analysis, C['danger'])
        self.btn_stop.pack(side="right", padx=(0, sp['sm']))
        self.btn_pause = _action_btn(prog_row, TEXTS["btn_pause"],
                                     self._toggle_pause, C['warning'])
        self.btn_pause.pack(side="right", padx=(0, sp['sm']))
        self.progress_label = ctk.CTkLabel(self.ctrl_frame, text=TEXTS["progress_idle"],
                                           font=(_FONT_UI, fs['sm']), anchor='w', height=18,
                                           text_color=C['text_secondary'])
        self.progress_label.pack(fill='x', padx=sp['lg'], pady=(0, sp['xs']))
        self.progress_bar = ctk.CTkProgressBar(self.ctrl_frame, height=8, corner_radius=4,
                                               fg_color=C['bg_card_hover'],
                                               progress_color=C['success'])
        self.progress_bar.pack(fill='x', padx=sp['lg'], pady=(0, sp['md']))
        self.progress_bar.set(0)

        log_container = ctk.CTkFrame(self.main_frame, fg_color=C['bg_input'],
                                     corner_radius=RADIUS['panel'], border_width=1,
                                     border_color=C['divider'])
        log_container.pack(fill="both", expand=True, padx=0, pady=0)
        self._log_container = log_container

        log_header = ctk.CTkFrame(log_container, fg_color='transparent', height=36)
        log_header.pack(fill="x", padx=sp['md'], pady=(sp['sm'], 0))

        ctk.CTkLabel(log_header, text=TEXTS["log_header_title"], font=(_FONT_UI, fs['xs'], "bold"),
                    text_color=C['text_secondary']).pack(side="left", padx=(sp['xs'], 0))

        ctk.CTkLabel(log_header, text="·",
                    font=(_FONT_UI, fs['sm']), text_color=C['text_muted']).pack(side="left", padx=sp['sm'])

        ctk.CTkLabel(log_header, text=TEXTS["log_header_subtitle"],
                    font=(_FONT_UI, fs['xs']), text_color=C['text_muted']).pack(side="left", padx=0)
        ctk.CTkCheckBox(log_header, text=TEXTS["chk_show_details"], variable=self.show_details_var,
                        command=self._toggle_details, font=(_FONT_UI, fs['xs']),
                        text_color=C['text_secondary'], checkbox_width=16,
                        checkbox_height=16, corner_radius=4, border_width=1,
                        border_color=C['border_hover'], fg_color=C['accent_fill'],
                        hover_color=C['bg_card_hover']).pack(side="right")
        self._log_toggle = ctk.CTkButton(log_header, text=TEXTS["log_collapse"], width=72, height=24,
                                         corner_radius=RADIUS['field'], fg_color='transparent',
                                         hover_color=C['bg_card_hover'], text_color=C['text_secondary'],
                                         border_width=1, border_color=C['border'],
                                         font=(_FONT_UI, fs['xs'], "bold"),
                                         command=self._toggle_log_panel)
        self._log_toggle.pack(side='right', padx=(0, sp['md']))

        self.log = ctk.CTkTextbox(log_container, font=(_FONT_MONO, fs['sm']), wrap='word',
                                  fg_color=C['bg_input'],
                                  text_color=C['text_secondary'],
                                  border_width=0,
                                  corner_radius=RADIUS['card'],
                                  spacing1=2, spacing3=2,
                                  scrollbar_button_color=C['border'],
                                  scrollbar_button_hover_color=C['border_hover'])
        self.log.pack(fill="both", expand=True, padx=sp['xs'], pady=(0, sp['xs']))

        self.log.tag_config("success", foreground=C['success_light'])
        self.log.tag_config("error", foreground=C['danger'])
        self.log.tag_config("warning", foreground=C['warning'])
        self.log.tag_config("info", foreground=C['accent_cyan'])
        self.log.tag_config("header", foreground=C['accent_blue_light'])
        self.log.tag_config("debug", foreground=C['text_tertiary'], elide=True)

        self.log.insert("end", TEXTS["log_welcome"])
        backend_messages.set_language(get_language())

        self._update_models_state()
        self._restore_state()
        self._poll_stop_count()
        self._refresh_visual_state()
        self.bind_all("<Control-q>", lambda e: self.destroy())
        self.update_idletasks()
        self.deiconify()

    # ── Appearance derived from the state ─────────

    def _make_clickable(self, frame, widgets, command, after=None):
        """Make a whole card clickable, with a hover tint."""
        def _enter(_e=None):
            frame._hovering = True
            base = getattr(frame, '_base_fg', None)
            if base and not getattr(frame, '_disabled', False):
                frame.configure(fg_color=mix(base, '#ffffff', 0.04))

        def _leave(_e=None):
            def _check():
                try:
                    x, y = frame.winfo_pointerxy()
                    fx, fy = frame.winfo_rootx(), frame.winfo_rooty()
                    inside = fx <= x < fx + frame.winfo_width() and fy <= y < fy + frame.winfo_height()
                except Exception:
                    inside = False
                if not inside:
                    frame._hovering = False
                    base = getattr(frame, '_base_fg', None)
                    if base:
                        frame.configure(fg_color=base)
            frame.after(16, _check)

        def _click(_e=None):
            command()
            if after:
                frame.after(50, after)

        for w in widgets:
            try:
                w.bind("<Button-1>", _click, add='+')
                w.bind("<Enter>", _enter, add='+')
                w.bind("<Leave>", _leave, add='+')
                w.configure(cursor='hand2')
            except Exception:
                pass

    def _toggle_settings(self):
        """Show or hide the advanced settings."""
        body = self._settings_body
        if body.winfo_ismapped():
            body.pack_forget()
            self._settings_btn.configure(text=TEXTS["settings_show"])
        else:
            body.pack(fill='x', pady=(SPACE['md'], 0))
            self._settings_btn.configure(text=TEXTS["settings_hide"])
            self.after(50, self._scroll_sidebar_to, self._settings_btn)

    def _scroll_sidebar_to(self, widget):
        """Scroll the sidebar so that widget is near the top."""
        try:
            canvas = self._sidebar_scroll._parent_canvas
            self.update_idletasks()
            inner = self._sidebar_scroll
            total = max(1, inner.winfo_height())
            y = widget.winfo_rooty() - inner.winfo_rooty() - 3 * SPACE['xl']
            canvas.yview_moveto(max(0.0, min(1.0, y / total)))
        except Exception:
            pass

    def _on_model_switch(self, code):
        var = self.model_vars[code]
        if code in self._auto_nulls and not var.get():
            self._excluded_nulls.add(code)       # automatic null switched off
        elif code in self._off_nulls and var.get():
            var.set(False)                       # back to automatic
            self._excluded_nulls.discard(code)
        elif var.get():
            self._excluded_nulls.discard(code)
        self._models_changed()

    def _tile_click(self, code, switch):
        if code in self._auto_nulls:
            self._excluded_nulls.add(code)
            self._models_changed()
        elif code in self._off_nulls:
            self._excluded_nulls.discard(code)
            self._models_changed()
        else:
            switch.toggle()

    def _models_changed(self):
        self._update_models_state()
        self._apply_visual_state(force=True)

    def _paint_tile(self, code):
        """Card colours: on, added automatically as a null, or off."""
        p = self._tile_parts.get(code)
        if not p:
            return
        C = self.COLORS
        card, cb, lbl, auto_lbl, accent, bg = (p['card'], p['cb'], p['lbl'], p['auto'],
                                               p['accent'], p['bg'])
        on = bool(p['var'].get())
        alt = None if on else self._auto_nulls.get(code)
        left_out = None if on or alt else self._off_nulls.get(code)
        try:
            if on:
                fill, border = mix(bg, accent, 0.16), mix(bg, accent, 0.7)
            elif alt:
                fill, border = mix(bg, accent, 0.06), mix(bg, accent, 0.4)
            else:
                fill, border = bg, C['border']
            card._base_fg = fill
            card.configure(fg_color=fill, border_color=border)
            if alt or left_out:
                if alt:
                    cb.configure(progress_color=mix(bg, accent, 0.55),
                                 button_color=mix(PALETTE['switch_knob'], bg, 0.25))
                    cb.select(from_variable_callback=True)
                    auto_lbl.configure(text=TEXTS["tile_auto"].format(alt=alt),
                                       text_color=C['text_secondary'])
                else:
                    cb.configure(progress_color=accent, fg_color=PALETTE['switch_track'],
                                 button_color=PALETTE['switch_knob'])
                    cb.deselect(from_variable_callback=True)
                    auto_lbl.configure(text=TEXTS["tile_auto_off"].format(alt=left_out, null=code),
                                       text_color=C['warning'])
                if lbl.winfo_manager():
                    lbl.pack_forget()
                if not auto_lbl.winfo_manager():
                    auto_lbl.pack(fill='x', padx=(SPACE['md'] + 4, SPACE['md']),
                                  pady=(SPACE['xs'], SPACE['sm']))
            else:
                cb.configure(progress_color=accent, fg_color=PALETTE['switch_track'],
                             button_color=PALETTE['switch_knob'])
                if not on:
                    cb.deselect(from_variable_callback=True)
                if auto_lbl.winfo_manager():
                    auto_lbl.pack_forget()
                if not lbl.winfo_manager():
                    lbl.pack(fill='x', padx=(SPACE['md'] + 4, SPACE['md']),
                             pady=(SPACE['xs'], SPACE['sm']))
                lbl.configure(text_color=C['text_secondary'] if on else C['text_tertiary'])
        except Exception:
            pass

    def _toggle_log_panel(self):
        """Collapse or expand the log."""
        if self.log.winfo_ismapped():
            self.log.pack_forget()
            self._log_container.pack_configure(expand=False)
            self._log_toggle.configure(text=TEXTS["log_expand"])
        else:
            self.log.pack(fill="both", expand=True, padx=SPACE['xs'], pady=(0, SPACE['xs']))
            self._log_container.pack_configure(expand=True)
            self._log_toggle.configure(text=TEXTS["log_collapse"])

    def _gene_count(self):
        folder = self.input_folder
        if not folder:
            return None
        key = str(folder)
        if key not in self._gene_count_cache:
            try:
                chosen, _ = group_by_gene(list_alignment_files(folder))
                self._gene_count_cache[key] = len(chosen)
            except Exception:
                self._gene_count_cache[key] = 0
        return self._gene_count_cache[key]

    def _refresh_visual_state(self):
        """Poll the state and refresh the appearance when it changes."""
        try:
            sig = (str(self.input_folder), str(self.tree_file), len(self.per_gene_trees),
                   str(self.output_folder),
                   tuple(k for k, v in self.model_vars.items() if v.get()),
                   bool(self.include_neutral_models.get()), tuple(sorted(self._excluded_nulls)),
                   int(self.cores_var.get()), str(self.status_indicator.cget('text_color')),
                   str(self.stop_label.cget('text')),
                   self.codonfreq_var.get(), self.entry_omega.get(), self.entry_ncatg.get(),
                   tuple(str(s[3].cget('text')) + str(s[3].cget('text_color')) for s in self._slots.values()))
            if sig != self._visual_sig:
                self._visual_sig = sig
                self._apply_visual_state()
        except Exception:
            pass
        self.after(400, self._refresh_visual_state)

    def _apply_visual_state(self, force=False):
        C = self.COLORS
        done = {'input': bool(self.input_folder),
                'tree': bool(self.tree_file) or bool(self.per_gene_trees),
                'output': bool(self.output_folder)}
        next_slot = next((k for k in ('input', 'tree', 'output') if not done[k]), None)
        for key, (card, title, btn, detail, shown) in self._slots.items():
            lines = []
            for line in str(detail.cget('text')).splitlines():
                lines.append(line.split(': ')[0] if ': ' in line else line)
            if key == 'tree' and not self.tree_file and self.per_gene_trees:
                lines = [TEXTS["slot_tree_per_gene"].format(n=len(self.per_gene_trees))]
            col = detail.cget('text_color')
            shown.configure(text="\n".join(lines[:3]),
                            text_color=col if done[key] or col == C['warning'] else C['text_tertiary'])
            if done[key]:
                card._base_fg = C['bg_card_hover']
                card.configure(fg_color=C['bg_card_hover'], border_color=C['bg_card_hover'])
                btn.configure(text=TEXTS["slot_change"], fg_color='transparent',
                              text_color=C['text_secondary'], border_color=C['border'],
                              hover_color=C['bg_hover'])
            else:
                card._base_fg = C['bg_card']
                card.configure(fg_color=C['bg_card'],
                               border_color=mix(C['bg_card'], C['accent_fill'], 0.7)
                               if key == next_slot else C['border'])
                if key == next_slot:
                    btn.configure(text=TEXTS["slot_choose"], fg_color=C['accent_fill'],
                                  text_color='#ffffff', border_color=C['accent_fill'],
                                  hover_color=C['accent_blue_hover'])
                else:
                    btn.configure(text=TEXTS["slot_choose"], fg_color='transparent',
                                  text_color=C['text_secondary'], border_color=C['border'],
                                  hover_color=C['bg_hover'])

        def _circle(lbl, ok, number):
            if ok:
                lbl.configure(text="✓", fg_color=PALETTE['success_fill'], text_color='#ffffff')
            else:
                lbl.configure(text=str(number), fg_color=C['bg_card_hover'], text_color=C['text_secondary'])

        data_ok = all(done.values())
        _circle(self._step1_circle, data_ok, 1)
        chosen = [k for k, v in self.model_vars.items() if v.get()]
        _circle(self._step2_circle, bool(chosen), 2)
        count = TEXTS["step_models_count"].format(n=len(chosen)) if chosen else ""
        self._models_count.configure(text=count)
        self._step2_text.configure(
            text=(count + ": " + ", ".join(chosen)) if chosen else TEXTS["step_models_none"],
            text_color=C['text_secondary'] if chosen else C['text_tertiary'])

        models = self._selected_models() if chosen else []
        order = list(self.model_vars)
        models = sorted(models, key=lambda m: order.index(m) if m in order else len(order))
        def _alt_of(m):
            for alt in chosen:
                nulls = CodemlBatchAnalysis.NULL_MODEL_PAIRS.get(alt, [])
                nulls = [nulls] if isinstance(nulls, str) else list(nulls)
                if m in nulls:
                    return alt
            return ", ".join(chosen)

        auto = {m: _alt_of(m) for m in models if m not in chosen}
        off = {}
        if chosen and self.include_neutral_models.get():
            every = CodemlBatchAnalysis.auto_complete_null_models(chosen, include_neutral=True)
            off = {m: _alt_of(m) for m in every if m not in chosen and m in self._excluded_nulls}
        self._auto_nulls, self._off_nulls = auto, off
        for code in self._tile_parts:
            self._paint_tile(code)
        if auto:
            self._step2_text.configure(text=self._step2_text.cget('text') + f" (+{len(auto)} auto)")

        n_genes = self._gene_count()
        parts = [TEXTS["run_summary_no_data"] if n_genes is None
                 else TEXTS["run_summary_genes"].format(n=n_genes)]
        if models:
            parts.append(TEXTS["run_summary_models"].format(n=len(models)) + ": " + ", ".join(models))
        parts.append(TEXTS["run_summary_cpus"].format(n=int(self.cores_var.get())))
        self._summary_label.configure(text="  ·  ".join(parts))
        pairs = lrt_pairs_for(models) if models else []
        if not models:
            tests, tcolor = TEXTS["run_summary_no_models"], C['text_tertiary']
        elif not pairs:
            tests, tcolor = TEXTS["run_tests_none"], C['warning']
        else:
            tests = TEXTS["run_tests_label"].format(
                tests="  ·  ".join(f"{alt} vs {null}" for null, alt in pairs))
            tcolor = C['text_secondary']
        self._summary_models.configure(text=tests, text_color=tcolor)

        try:
            cf = str(self.codonfreq_var.get()).split('= ')[-1]
            omega = self.entry_omega.get().strip() or '0.5'
            if get_language() == 'pt':
                omega = omega.replace('.', ',')
            self._settings_summary.configure(text=" · ".join((
                cf, f"ncatG {self.entry_ncatg.get().strip()}", f"ω {omega}",
                f"κ {DEFAULT_CTL_PARAMS['kappa']}",
                TEXTS["run_summary_cpus"].format(n=int(self.cores_var.get())))))
        except Exception:
            pass

        try:
            import re as _re
            n_stops = int((_re.findall(r'\d+', str(self.stop_label.cget('text'))) or ['0'])[-1])
            self.stop_label.configure(text_color=C['danger'] if n_stops else C['text_tertiary'])
        except Exception:
            pass

        try:
            col = self.status_indicator.cget('text_color')
            col = col if isinstance(col, str) and col.startswith('#') else C['text_tertiary']
            self.status_indicator.configure(fg_color=mix(C['bg_dark'], col, 0.16))
        except Exception:
            pass

    def _global_ctl_options(self) -> dict:
        """.ctl options chosen under Advanced settings."""
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
        """Show or hide debug messages in the log."""
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
        win.configure(fg_color=self.COLORS['bg_dark'])
        ctk.CTkLabel(win, text=title, font=(_FONT_UI, FONT_SIZE['lg'], "bold"),
                     text_color=self.COLORS['text_primary']).pack(
                         padx=SPACE['xl'], pady=(SPACE['xl'], SPACE['sm']), anchor='w')
        ctk.CTkLabel(win, text=body, font=(_FONT_UI, FONT_SIZE['md']), wraplength=600,
                     justify='left',
                     text_color=self.COLORS['text_secondary']).pack(
                         padx=SPACE['xl'], pady=(0, SPACE['lg']), anchor='w')
        ctk.CTkButton(win, text="OK", width=80, height=32, command=win.destroy, text_color='#ffffff',
                      fg_color=PALETTE['accent_fill'], hover_color=self.COLORS['accent_blue_hover'],
                      corner_radius=RADIUS['card'], font=(_FONT_UI, FONT_SIZE['md'], "bold")).pack(
                          padx=SPACE['xl'], pady=(0, SPACE['xl']), anchor='e')

    def _update_cores_label(self, value=None):
        n = int(self.cores_var.get())
        self.cores_disp.configure(text=f"{n}×")

    def _setup_model_list(self):
        MODEL_META = {
            'M0': {'color': '#3b82f6'},
            'M1a': {'color': '#06b6d4'},
            'M2a': {'color': '#10b981'},
            'M7': {'color': '#8b5cf6'},
            'M8': {'color': '#ec4899'},
            'M8a': {'color': '#f472b6'},
            'Branch': {'color': '#f59e0b'},
            'Branch-site': {'color': '#ef4444'},
            'Branch-site_null': {'color': '#6d6d6d'},
        }

        models = {
            TEXTS["tab_site_models"]:  ['M0', 'M1a', 'M2a', 'M7', 'M8', 'M8a'],
            TEXTS["tab_branch_model"]: ['Branch'],
            TEXTS["tab_branchsite"]:   ['Branch-site', 'Branch-site_null'],
        }

        C = self.COLORS
        sp = SPACE
        tile_bg = C['bg_card']
        # card order follows the LRT pairs (null before alternative)
        layout = {
            TEXTS["tab_site_models"]: [['M0', 'M1a', 'M2a'], ['M7', 'M8', 'M8a']],
            TEXTS["tab_branch_model"]: [['Branch']],
            TEXTS["tab_branchsite"]: [['Branch-site_null', 'Branch-site']],
        }

        def _tile(parent, code):
            meta = MODEL_META.get(code, {'color': C['accent_blue']})
            accent = meta['color']
            card = ctk.CTkFrame(parent, fg_color=tile_bg, corner_radius=RADIUS['card'],
                                border_width=1, border_color=C['border'])
            card._base_fg = tile_bg

            top = ctk.CTkFrame(card, fg_color='transparent')
            top.pack(fill='x', padx=(sp['md'], sp['sm']), pady=(sp['sm'], 0))

            display_name = self.codeml_backend.MODEL_CONFIGS[code].get('display_name', code)
            var = ctk.BooleanVar(value=False)
            cb = ctk.CTkSwitch(
                top,
                text=f"  {display_name}",
                variable=var,
                onvalue=True, offvalue=False,
                command=lambda c=code: self._on_model_switch(c),
                switch_width=32, switch_height=18,
                progress_color=accent,
                button_color=PALETTE['switch_knob'],
                button_hover_color=PALETTE['switch_knob'],
                fg_color=PALETTE['switch_track'],
                text_color=C['text_primary'],
                text_color_disabled=C['text_tertiary'],
                font=(_FONT_UI, FONT_SIZE['md'], "bold")
            )
            cb.pack(side='left')

            gear = ctk.CTkButton(
                top, text=TEXTS["cfg_btn"], width=48, height=20,
                fg_color='transparent',
                hover_color=C['bg_hover'],
                text_color=C['text_secondary'],
                text_color_disabled=C['text_tertiary'],
                corner_radius=RADIUS['field'],
                border_width=0,
                font=(_FONT_UI, FONT_SIZE['xs']),
                command=lambda c=code: self._open_config_window(c)
            )
            gear.pack(side='right')
            info_btn = ctk.CTkButton(
                top, text="?", width=20, height=20, corner_radius=10,
                fg_color='transparent',
                hover_color=C['bg_hover'],
                text_color=C['text_secondary'],
                border_width=0,
                font=(_FONT_UI, FONT_SIZE['xs'], "bold"),
                command=lambda c=code: self._show_model_info(c)
            )
            info_btn.pack(side='right', padx=(0, sp['xs']))

            lbl = ctk.CTkLabel(card, text=TEXTS["model_desc"].get(code, TEXTS["model_status_default"]),
                               font=(_FONT_UI, FONT_SIZE['xs']), anchor='nw', justify='left',
                               wraplength=240, height=32,
                               text_color=C['text_tertiary'])
            lbl.pack(fill='x', padx=(sp['md'] + 4, sp['md']), pady=(sp['xs'], sp['sm']))

            def _wrap(e, l=lbl):
                w = max(160, e.width - 2 * sp['md'] - 4)
                if abs(int(l.cget('wraplength')) - w) > 4:
                    l.configure(wraplength=w)
            card.bind('<Configure>', _wrap, add='+')

            auto_lbl = ctk.CTkLabel(card, text="", font=(_FONT_UI, FONT_SIZE['xs']), anchor='nw',
                                    justify='left', wraplength=240, height=32,
                                    text_color=C['text_secondary'])
            card.bind('<Configure>', lambda e, l=auto_lbl: l.configure(
                wraplength=max(160, e.width - 2 * sp['md'] - 4)), add='+')

            self._make_clickable(card, (card, top, lbl, auto_lbl),
                                 lambda c=code, sw=cb: self._tile_click(c, sw))

            self._tile_parts[code] = {'card': card, 'cb': cb, 'lbl': lbl, 'auto': auto_lbl,
                                      'accent': accent, 'var': var, 'bg': tile_bg}
            var.trace_add('write', lambda *_a, c=code: self._paint_tile(c))

            self.model_vars[code]         = var
            self.model_ctl_labels[code]   = lbl
            self.model_checkboxes[code]   = cb
            self.model_gear_buttons[code] = gear
            self._tiles[code] = card
            self._paint_tile(code)
            return card

        outers = {}
        for tab_name in layout:
            outer = ctk.CTkFrame(self.tabs.tab(tab_name), fg_color='transparent')
            outer.pack(fill="x", padx=0, pady=(sp['xs'], 0))
            for c in range(3):
                outer.grid_columnconfigure(c, weight=1, uniform='tiles')
            outers[tab_name] = outer
        for tab_name, codes in models.items():
            for code in codes:
                _tile(outers[tab_name], code)
        for tab_name, rows in layout.items():
            for r, codes in enumerate(rows):
                for c, code in enumerate(codes):
                    self._tiles[code].grid(row=r, column=c, sticky='nsew',
                                           padx=sp['xs'], pady=(0, sp['sm']))

        def _label_btn(tab, text, mode):
            btn = ctk.CTkButton(
                tab, text=text,
                fg_color='transparent',
                hover_color=C['bg_card_hover'],
                command=lambda: self._open_tree_labeler(mode=mode),
                height=40,
                font=(_FONT_UI, FONT_SIZE['md'], "bold"),
                text_color=C['accent_blue_light'],
                text_color_disabled=C['text_tertiary'],
                border_width=1, border_color=C['border'],
                corner_radius=RADIUS['card']
            )
            btn.pack(fill='x', padx=sp['xs'], pady=(sp['md'], 0))
            return btn

        self.btn_label_branch = _label_btn(self.tabs.tab(TEXTS["tab_branch_model"]),
                                           TEXTS["btn_label_branch"], 'branch')

        self.btn_label_branchsite = _label_btn(self.tabs.tab(TEXTS["tab_branchsite"]),
                                               TEXTS["btn_label_branchsite"], 'branchsite')

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
        path = ask_directory(self, TEXTS["btn_input_folder"], start)
        if path:
            self._use_input_folder(path)

    def _use_input_folder(self, path):
        self.input_folder = Path(path)
        self.stop_label.configure(text=TEXTS["status_stops_template"].format(n="–"))
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
        path = ask_open_file(self, TEXTS["btn_tree_file"], start,
                             filetypes=[('Newick', '*.nwk *.tree *.tre *.newick *.nh *.txt')])
        if path:
            self._use_tree_file(path)

    def _use_tree_file(self, path):
        self.tree_file = Path(path)
        self.label_tree.configure(text=str(self.tree_file.name),
                                  text_color=self.COLORS['text_secondary'])
        self._update_models_state()

    def select_output_folder(self):
        """Choose the output folder, creating it if needed."""
        start = str(self.output_folder.parent if self.output_folder else
                    (self.input_folder.parent if self.input_folder else Path.home()))
        path = ask_directory(self, TEXTS["dialog_choose_output"], start,
                             allow_new=True, must_exist=False)
        if path:
            self._use_output_folder(path)

    def _use_output_folder(self, path):
        folder = Path(path)
        created = not folder.exists()
        try:
            folder.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            show_message(self, TEXTS["msg_error"],
                         TEXTS["msg_output_folder_error"].format(error=os_error_text(exc)), 'error')
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
        hdr.pack(fill='x', padx=SPACE['xl'], pady=(SPACE['xl'], SPACE['sm']))
        ctk.CTkLabel(hdr, text=TEXTS["neutral_header"], font=(_FONT_UI, FONT_SIZE['lg'], "bold"),
                     text_color=self.COLORS['text_primary']).pack(anchor='w')
        ctk.CTkLabel(hdr, text=TEXTS["neutral_intro"], font=(_FONT_UI, FONT_SIZE['md']), wraplength=700,
                     justify='left', text_color=self.COLORS['text_secondary']).pack(
                         anchor='w', pady=(SPACE['xs'], 0))

        scroll = ctk.CTkScrollableFrame(win, fg_color='transparent',
                                        scrollbar_button_color=self.COLORS['border'])
        scroll.pack(fill='both', expand=True, padx=SPACE['lg'], pady=SPACE['sm'])
        for null, alt, color, title, test, detail, note in TEXTS["lrt_pairs"]:
            card = ctk.CTkFrame(scroll, fg_color=self.COLORS['bg_card'], corner_radius=RADIUS['card'],
                                border_width=0)
            card.pack(fill='x', padx=SPACE['xs'], pady=SPACE['xs'])
            ctk.CTkFrame(card, fg_color=color, width=4, height=8, corner_radius=2).pack(
                side='left', fill='y', padx=(SPACE['sm'], SPACE['md']), pady=SPACE['md'])
            content = ctk.CTkFrame(card, fg_color='transparent')
            content.pack(side='left', fill='both', expand=True, pady=SPACE['md'], padx=(0, SPACE['md']))
            row = ctk.CTkFrame(content, fg_color='transparent')
            row.pack(fill='x', anchor='w')
            ctk.CTkLabel(row, text=title, font=(_FONT_UI, FONT_SIZE['md'], "bold"),
                         text_color=self.COLORS['text_primary']).pack(side='left')
            ctk.CTkLabel(row, text=TEXTS["neutral_null_alt"].format(null=null, alt=alt),
                         font=(_FONT_UI, FONT_SIZE['sm']), text_color=self.COLORS['text_secondary']).pack(
                             side='left', padx=(SPACE['md'], 0))
            ctk.CTkLabel(content, text=test, font=(_FONT_UI, FONT_SIZE['md'], "bold"),
                         text_color=self.COLORS['accent_blue_light'], anchor='w').pack(
                             anchor='w', pady=(SPACE['xs'], SPACE['xs']))
            ctk.CTkLabel(content, text=detail, font=(_FONT_UI, FONT_SIZE['sm']), wraplength=600, justify='left',
                         text_color=self.COLORS['text_secondary'], anchor='w').pack(anchor='w', pady=(0, SPACE['xs']))
            ctk.CTkLabel(content, text=TEXTS["neutral_note"].format(note=note), font=(_FONT_UI, FONT_SIZE['sm'], "italic"),
                         wraplength=600, justify='left', text_color=self.COLORS['text_tertiary'],
                         anchor='w').pack(anchor='w')
        ctk.CTkLabel(win, text=TEXTS["neutral_footer"], font=(_FONT_UI, FONT_SIZE['sm']),
                     text_color=self.COLORS['text_secondary'], wraplength=700, justify='left').pack(
                         padx=SPACE['xl'], pady=(SPACE['sm'], SPACE['lg']), anchor='w')

    def _show_model_info(self, model_code: str):
        """Help window for one model."""
        info_src = (self.codeml_backend.MODEL_INFO_PT if get_language() == 'pt'
                    else self.codeml_backend.MODEL_INFO)
        model_info = info_src.get(model_code) or self.codeml_backend.MODEL_INFO.get(model_code)
        if not model_info:
            return
        
        info_window = ctk.CTkToplevel(self)
        info_window.title(TEXTS["model_config_header"].format(model_code=model_code))
        fit_to_screen(info_window, 750, 600, min_w=480, min_h=360)
        info_window.transient(self)
        info_window.after(50, lambda: (info_window.grab_set(), info_window.focus_set()))
        info_window.bind("<Escape>", lambda e: info_window.destroy())
        
        info_window.configure(fg_color=self.COLORS['bg_dark'])

        display_name = self.codeml_backend.MODEL_CONFIGS[model_code].get('display_name', model_code)
        header = ctk.CTkLabel(info_window, 
                             text=f"{model_code} - {model_info.get('full_name', '')}",
                             font=(_FONT_UI, FONT_SIZE['lg'], "bold"), anchor='w', justify='left',
                             wraplength=680,
                             text_color=self.COLORS['text_primary'])
        header.pack(fill='x', padx=SPACE['xl'], pady=(SPACE['xl'], SPACE['md']))
        
        scroll_frame = ctk.CTkScrollableFrame(info_window, fg_color=self.COLORS['bg_card'],
                                              corner_radius=RADIUS['panel'],
                                              scrollbar_button_color=self.COLORS['border'])
        scroll_frame.pack(fill="both", expand=True, padx=SPACE['xl'], pady=(0, SPACE['xl']))

        sections = (
            (TEXTS["model_info_test_type"], model_info.get('test_type', ''), self.COLORS['text_primary'], ''),
            (TEXTS["model_info_params"], model_info.get('parameters', ''), self.COLORS['text_primary'], ''),
            (TEXTS["model_info_purpose"], model_info.get('purpose', ''), self.COLORS['text_primary'], ''),
            (TEXTS["model_info_interpretation"], model_info.get('interpretation', ''),
             self.COLORS['text_primary'], ''),
            (TEXTS["model_info_use_case"], model_info.get('use_case', ''), self.COLORS['text_primary'], ''),
            (TEXTS["model_info_references"], model_info.get('references', ''),
             self.COLORS['text_secondary'], 'italic'),
        )
        for i, (title, value, color, style) in enumerate(sections):
            ctk.CTkLabel(scroll_frame, text=title, font=(_FONT_UI, FONT_SIZE['sm'], "bold"),
                         text_color=self.COLORS['accent_blue_light']).pack(
                             anchor="w", padx=SPACE['md'], pady=(SPACE['lg'] if i else SPACE['md'], SPACE['xs']))
            font = (_FONT_UI, FONT_SIZE['sm'], style) if style else (_FONT_UI, FONT_SIZE['sm'])
            ctk.CTkLabel(scroll_frame, text=value, font=font, text_color=color,
                         wraplength=640, justify="left").pack(anchor="w", padx=SPACE['md'], pady=0)

    # ── Language and theme ───────────────────────────────────────────────────────────────

    @staticmethod
    def _lang_pref_path() -> Path:
        """File that stores the chosen language."""
        return Path.home() / '.easypml_lang'

    @staticmethod
    def load_language_pref() -> None:
        """Saved language, or the system language (Portuguese -> PT, anything
        else -> EN)."""
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
        """Sidebar footer: language, About and theme."""
        C = self.COLORS
        footer = ctk.CTkFrame(self.sidebar, fg_color='transparent', corner_radius=0, border_width=0)
        footer.pack(fill='x', side='bottom', padx=0, pady=0)
        ctk.CTkFrame(footer, fg_color=C['divider'], height=1, corner_radius=0).pack(fill='x')

        inner = ctk.CTkFrame(footer, fg_color='transparent')
        inner.pack(fill='x', padx=SPACE['lg'], pady=SPACE['md'])

        current = get_language()

        def _make_btn(code: str, label: str):
            active = (code == current)
            btn = ctk.CTkButton(
                inner, text=label, width=44, height=28,
                font=(_FONT_UI, FONT_SIZE['sm'], "bold"),
                corner_radius=RADIUS['field'],
                fg_color=C['accent_fill'] if active else 'transparent',
                hover_color=C['accent_blue_hover'] if active else C['bg_card_hover'],
                text_color='#ffffff' if active else C['text_secondary'],
                command=lambda c=code: self._switch_language(c)
            )
            btn.pack(side='left', padx=(0, SPACE['xs']))

        _make_btn('pt', 'PT')
        _make_btn('en', 'EN')
        ctk.CTkButton(inner, text=TEXTS["btn_about"], width=72, height=28,
                      font=(_FONT_UI, FONT_SIZE['sm']), corner_radius=RADIUS['field'],
                      fg_color='transparent', hover_color=C['bg_card_hover'],
                      text_color=C['text_secondary'], command=lambda: show_about(self)).pack(side='right')

        row = ctk.CTkFrame(footer, fg_color='transparent')
        row.pack(fill='x', padx=SPACE['lg'], pady=(0, SPACE['md']))
        ctk.CTkLabel(row, text=TEXTS["theme_label"], font=(_FONT_UI, FONT_SIZE['sm']),
                     text_color=C['text_secondary']).pack(side='left', padx=(0, SPACE['sm']))
        names = {'system': TEXTS["theme_system"], 'light': TEXTS["theme_light"],
                 'dark': TEXTS["theme_dark"]}
        self._theme_names = names
        seg = ctk.CTkSegmentedButton(
            row, values=[names[k] for k in THEME_CHOICES], height=28,
            font=(_FONT_UI, FONT_SIZE['sm']), corner_radius=RADIUS['field'],
            fg_color=C['bg_card_hover'], selected_color=C['accent_fill'],
            selected_hover_color=C['accent_blue_hover'], unselected_color=C['bg_card_hover'],
            unselected_hover_color=C['bg_hover'], text_color=C['text_primary'],
            command=lambda label: self._switch_theme(
                next(k for k, v in names.items() if v == label)))
        seg.set(names[CURRENT_THEME['choice']])
        seg.pack(side='left', fill='x', expand=True)
        self._theme_seg = seg

    _RESTORE_ENV = 'EASYPAML_RESTORE'

    def _reopen(self) -> None:
        """Start a new window with the same folders and models, then close this one."""
        import json
        import subprocess
        state = {'input': str(self.input_folder or ''), 'tree': str(self.tree_file or ''),
                 'output': str(self.output_folder or ''),
                 'models': [k for k, v in self.model_vars.items() if v.get()],
                 'excluded': sorted(self._excluded_nulls)}
        env = dict(os.environ, **{self._RESTORE_ENV: json.dumps(state)})
        try:
            subprocess.Popen([sys.executable] + sys.argv, env=env)
            self.destroy()
        except Exception as exc:
            show_message(self, "EasyPAML", TEXTS["lang_switch_err"].format(error=exc), 'error')

    def _restore_state(self) -> None:
        import json
        try:
            state = json.loads(os.environ.pop(self._RESTORE_ENV, '') or '{}')
        except ValueError:
            return
        if state.get('input') and Path(state['input']).is_dir():
            self._use_input_folder(state['input'])
        if state.get('tree') and Path(state['tree']).is_file():
            self._use_tree_file(state['tree'])
        if state.get('output') and Path(state['output']).is_dir():
            self._use_output_folder(state['output'])
        for name in state.get('models', []):
            if name in self.model_vars:
                self.model_vars[name].set(True)
        self._excluded_nulls = set(state.get('excluded', []))
        self._update_models_state()

    def _analysis_running(self) -> bool:
        if self.analysis_thread and self.analysis_thread.is_alive():
            show_message(self, "EasyPAML", TEXTS["msg_wait_for_run"], 'warning')
            return True
        return False

    def _switch_theme(self, choice: str) -> None:
        """Save the chosen theme and reopen the program to apply it."""
        if choice == CURRENT_THEME['choice']:
            return
        save_theme_pref(choice)
        new_mode = system_theme() if choice == 'system' else choice
        if new_mode == CURRENT_THEME['mode']:
            CURRENT_THEME['choice'] = choice
            return
        if not self._analysis_running() and ask_yes_no(self, TEXTS['theme_label'], TEXTS['theme_restart']):
            self._reopen()

    def _switch_language(self, lang: str) -> None:
        """Save the chosen language and reopen the program to apply it."""
        current = get_language()
        if lang == current:
            return
        self._save_language_pref(lang)
        set_language(lang)

        if self._analysis_running():
            return
        if ask_yes_no(self, TEXTS['lang_restart_title'],
                      f"{TEXTS['lang_switch_message']}\n\n{TEXTS['lang_switch_confirm']}"):
            self._reopen()

    def _open_results_viewer(self):
        """Open the results panel; ask for a results folder when the current one has
        no results."""
        folder = self.output_folder
        if self.analysis_thread and self.analysis_thread.is_alive():
            show_message(self, TEXTS["btn_view_results"], TEXTS["msg_still_running"])
            return
        if not (folder and (Path(folder) / 'analysis_summary.tsv').exists()):
            chosen = ask_directory(self, TEXTS["dialog_select_results_folder"],
                                   folder or self.input_folder or Path.home())
            if not chosen:
                return
            folder = Path(chosen)
            if not (folder / 'analysis_summary.tsv').exists():
                show_message(self, TEXTS["btn_view_results"],
                             TEXTS["msg_not_results_folder"].format(path=folder), 'warning')
                return
        try:
            ResultsViewerWindow(self, folder)
        except Exception as e:
            self.append_log(tr("Erro ao abrir o painel de resultados: ", "Error opening the results panel: ") + f"{e}", "error")
            self.append_log(traceback.format_exc())

    def _all_genes_have_trees(self) -> bool:
        if not self.input_folder or not self.per_gene_trees:
            return False
        chosen, _ = group_by_gene(list_alignment_files(self.input_folder))
        return bool(chosen) and all(g in self.per_gene_trees for g in chosen)

    def _update_models_state(self):
        has_tree = bool(self.tree_file) or self._all_genes_have_trees()
        enabled = all([self.input_folder, has_tree, self.output_folder])
        state = "normal" if enabled else "disabled"
        try:
            if enabled:
                self.files_hint.pack_forget()
            elif not self.files_hint.winfo_ismapped():
                self.files_hint.pack(fill='x', pady=(0, SPACE['md']), ipadx=SPACE['md'], before=self.tabs)
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


    _LEVEL_TAGS = {'ok': 'success', 'error': 'error', 'warn': 'warning', 'info': None,
                   'header': 'header', 'debug': 'debug'}

    def append_log(self, text: str, level: str = None):
        """Add one line to the log. level: ok, error, warn, info, header or debug;
        without it the colour is guessed from the text. debug lines show only with
        "Show technical details"."""
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
            a = self.analysis_instance
            if a is not None and getattr(a, 'runs_total', 0) and not self.stop_event.is_set():
                self._show_live_progress(a)
        self.after(0, _upd)

    def _set_pause_button(self, paused: bool):
        C = self.COLORS
        if paused:
            self.btn_pause.configure(text=TEXTS["btn_resume"], fg_color=PALETTE['success_fill'],
                                     hover_color=mix(PALETTE['success_fill'], '#000000', 0.2),
                                     text_color='#ffffff', border_color=PALETTE['success_fill'])
        else:
            self.btn_pause.configure(text=TEXTS["btn_pause"], fg_color='transparent',
                                     hover_color=hover_tint(C['warning'], C['bg_card']),
                                     text_color=C['warning'], border_color=C['border'])

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
            pass  # without psutil, pausing happens between genes

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
        if not ask_yes_no(self, TEXTS["stop_confirm_title"], TEXTS["stop_confirm_text"],
                          yes=TEXTS["stop_confirm_yes"], no=TEXTS["stop_confirm_no"]):
            return
        # stop starting new runs, then end the running codeml processes
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
            chosen = set(selected)
            selected = [m for m in CodemlBatchAnalysis.auto_complete_null_models(selected, include_neutral=True)
                        if m in chosen or m not in self._excluded_nulls]
        return selected

    def start_analysis(self):
        """Run: check the data first and show what was found."""
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
        if self.output_folder and (Path(self.output_folder) / 'analysis_summary.tsv').exists() \
                and not ask_yes_no(self, TEXTS["overwrite_title"], TEXTS["overwrite_text"],
                                   yes=TEXTS["overwrite_yes"], no=TEXTS["overwrite_no"]):
            self.select_output_folder()
            return
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
        if report is not None:
            n_stops = sum(1 for i in report.issues if i.kind == 'stop_codon')
            self.stop_label.configure(text=TEXTS["status_stops_template"].format(n=n_stops))
        if error is not None:
            self.append_log(f"{error}", 'error')
        elif report is not None and not report.has_problems:
            self.append_log(TEXTS["preflight_ok"].format(n=len(report.genes)), 'ok')
        elif report is not None and report.has_problems:
            for line in report.format_text(get_language(), include_info=False).splitlines():
                self.append_log(line, 'warn')
            choice = PreflightDialog(self, report).show()
            if choice != 'continue':
                self.btn_run.configure(state="normal")
                return
            # "Continue anyway": codeml treats stop codons as missing data
            if any(i.kind == 'stop_codon' for i in report.issues):
                ignore_stops = True
        self._launch(selected, ignore_stops)

    def _launch(self, selected, ignore_stops: bool):
        self.stop_event.clear()
        self.pause_event.set()
        self._set_pause_button(False)
        self.progress_bar.set(0)
        self.progress_bar.configure(progress_color=self.COLORS['success'])
        self.progress_label.configure(text=TEXTS["progress_template"].format(done=0, total='…'))
        self.status_indicator.configure(text=TEXTS["status_running"], text_color=self.COLORS['success'])
        self.analysis_thread = threading.Thread(target=self._run_thread, args=(selected, ignore_stops),
                                                daemon=True)
        self.analysis_thread.start()

    def _timeout_seconds(self) -> int:
        """Minutes typed under Advanced settings, in seconds; 0 (automatic) when
        empty or invalid."""
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
        if summary.get('stopped'):
            self.progress_label.configure(text=TEXTS["progress_stopped"].format(
                ok=summary.get('ok', 0), total=total))
            self.progress_bar.configure(progress_color=self.COLORS['text_tertiary'])
        elif summary.get('failed'):
            self.progress_label.configure(text=TEXTS["progress_done_failed"].format(
                ok=summary.get('ok', 0), total=total, failed=summary['failed']))
            self.progress_bar.configure(progress_color=self.COLORS['warning'])
        else:
            self.progress_label.configure(text=TEXTS["progress_done"].format(done=done, total=total))
        if summary.get('failed'):
            def _item(g, r):
                why = ResultsViewerWindow._compact_reason(r).replace("\n", "\n    ")
                return f"• {g}\n    {why}"
            items = "\n".join(_item(g, r) for g, r in sorted(summary['failures'].items())[:30])
            if len(summary['failures']) > 30:
                items += f"\n… (+{len(summary['failures']) - 30})"
            show_message(self, TEXTS["failures_title"], TEXTS["failures_text"].format(
                failed=summary['failed'], total=total, items=items,
                path=Path(summary['output_folder']) / 'genes_status.tsv'), 'error')
        if summary.get('ok') and not summary.get('stopped'):
            self._open_results_viewer()

    @staticmethod
    def _clock(seconds: float) -> str:
        seconds = int(seconds)
        h, rest = divmod(seconds, 3600)
        m, s = divmod(rest, 60)
        return f"{h}:{m:02d}:{s:02d}" if h else f"{m}:{s:02d}"

    def _poll_stop_count(self):
        a = self.analysis_instance
        if a:
            cnt = getattr(a, 'current_stop_count', 0)
            self.stop_label.configure(text=TEXTS["status_stops_template"].format(n=cnt))
            if getattr(a, 'runs_total', 0) and not self.stop_event.is_set():
                self._show_live_progress(a)
        self.after(1000, self._poll_stop_count)

    def _show_live_progress(self, a):
        """Bar, percentage and time left from the backend's progress_snapshot (runs
        weighted by their expected time); label with the models running."""
        snap = a.progress_snapshot()
        running, frac, elapsed = snap['running'], snap['frac'], snap['elapsed']
        now = time.time()
        self.progress_bar.set(frac)
        if len(running) == 1:
            gene, (model, t0) = next(iter(running.items()))
            short = gene if len(gene) <= 24 else gene[:23] + "…"
            what = f"{model} ({short}, {self._clock(now - t0)})"
        elif running:
            by_model: dict = {}
            for model, _ in running.values():
                by_model[model] = by_model.get(model, 0) + 1
            what = ", ".join(TEXTS["progress_model_genes" if n > 1 else "progress_model_gene"].format(
                model=m, n=n) for m, n in by_model.items())
        else:
            what = "…"
        left = TEXTS["progress_left"].format(left=self._clock(snap['left'])) if snap['left'] else ""
        text = TEXTS["progress_running"].format(done=a.current_processed_genes,
                                                total=a.current_total_genes, pct=int(frac * 100),
                                                what=what, elapsed=self._clock(elapsed), left=left)
        self.progress_label.configure(text=text)


if __name__ == "__main__":
    app = App()
    app.mainloop()
