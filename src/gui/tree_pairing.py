"""Window that shows how the alignments were paired with the per-gene trees of a
folder, before the data check. Files are never renamed here: the user renames
them and checks again."""

import tkinter as tk
from pathlib import Path
from tkinter import ttk
from typing import Optional

import customtkinter as ctk

from backend.preflight import group_by_gene, list_alignment_files, pair_trees
from .gui_texts import TEXTS
from .summary_tab import ensure_tree_style
from .ui_helpers import FONT_SIZE, FONT_UI, PALETTE, RADIUS, SPACE, _button, _Modal, open_folder


class TreePairingDialog(_Modal):
    """Returns 'continue', or None when closed."""

    def __init__(self, parent, input_folder, tree_folder, general_tree: Optional[Path] = None):
        super().__init__(parent, TEXTS['pairing_title'], 1000, 580)
        self._input, self._folder, self._general = Path(input_folder), Path(tree_folder), general_tree
        self._summary = ctk.CTkLabel(self, text="", font=(FONT_UI, FONT_SIZE['md'], 'bold'), anchor='w',
                                     text_color=PALETTE['text_primary'])
        self._summary.pack(fill='x', padx=SPACE['xl'], pady=(SPACE['lg'], SPACE['xs']))
        ctk.CTkLabel(self, text=TEXTS['pairing_rules'], font=(FONT_UI, FONT_SIZE['sm']), anchor='w',
                     justify='left', wraplength=940, text_color=PALETTE['text_secondary']).pack(
            fill='x', padx=SPACE['xl'])

        box = ctk.CTkFrame(self, fg_color=PALETTE['bg_panel'], corner_radius=RADIUS['card'])
        box.pack(fill='both', expand=True, padx=SPACE['xl'], pady=SPACE['sm'])
        inner = tk.Frame(box, bg=PALETTE['bg_panel'])
        inner.pack(fill='both', expand=True, padx=SPACE['sm'], pady=SPACE['sm'])
        ensure_tree_style(self)
        self._tree = ttk.Treeview(inner, style='EP.Treeview', show='headings', selectmode='none',
                                  columns=('gene', 'tree', 'status'))
        sb = ctk.CTkScrollbar(inner, command=self._tree.yview)
        self._tree.configure(yscrollcommand=sb.set)
        sb.pack(side='right', fill='y')
        self._tree.pack(side='left', fill='both', expand=True)
        for cid, title, width in (('gene', TEXTS['pairing_col_gene'], 170), ('tree', TEXTS['pairing_col_tree'], 260),
                                  ('status', TEXTS['pairing_col_status'], 460)):
            self._tree.heading(cid, text=title, anchor='w')
            self._tree.column(cid, width=width, anchor='w', stretch=cid == 'status')
        self._tree.tag_configure('ok', foreground=PALETTE['success_fg'])
        self._tree.tag_configure('error', foreground=PALETTE['danger_fg'])
        self._tree.tag_configure('warn', foreground=PALETTE['warning_fg'])

        ctk.CTkLabel(self, text=TEXTS['pairing_branch_note'], font=(FONT_UI, FONT_SIZE['sm']), anchor='w',
                     justify='left', wraplength=940, text_color=PALETTE['text_tertiary']).pack(
            fill='x', padx=SPACE['xl'])
        row = ctk.CTkFrame(self, fg_color='transparent')
        row.pack(fill='x', padx=SPACE['xl'], pady=(SPACE['sm'], SPACE['lg']))
        _button(row, TEXTS['pairing_continue'], lambda: self._close('continue'),
                PALETTE['accent_fill']).pack(side='right')
        _button(row, TEXTS['pairing_check_again'], self._fill,
                PALETTE['neutral_fill']).pack(side='right', padx=(0, SPACE['sm']))
        _button(row, TEXTS['pairing_open_folder'], lambda: open_folder(self._folder, self),
                PALETTE['neutral_fill']).pack(side='right', padx=(0, SPACE['sm']))
        self._fill()

    def _fill(self):
        genes = sorted(group_by_gene(list_alignment_files(self._input))[0])
        pr = pair_trees(genes, self._input, self._folder)
        tree = self._tree
        tree.delete(*tree.get_children())
        general = self._general.name if self._general else None
        for gene in pr.missing:
            hint = pr.suggestions.get(gene)
            status = "✗ " + (TEXTS['pairing_did_you_mean'].format(name=hint.name) if hint
                             else TEXTS['pairing_no_tree'])
            if general:
                status += "  " + TEXTS['pairing_general_used'].format(tree=general)
            tree.insert('', 'end', values=(gene, '—', status), tags=('warn' if general else 'error',))
        for gene, files in pr.duplicates.items():
            tree.insert('', 'end', values=(gene, ", ".join(f.name for f in files),
                                           "✗ " + TEXTS['pairing_duplicate']), tags=('error',))
        for path in pr.orphans:
            tree.insert('', 'end', values=('—', path.name, "⚠ " + TEXTS['pairing_orphan']), tags=('warn',))
        for gene in sorted(pr.pairs):
            tool = pr.tools.get(gene)
            status = "✓ " + TEXTS['pairing_paired'] + (f" ({tool})" if tool else "")
            tree.insert('', 'end', values=(gene, pr.pairs[gene].name, status), tags=('ok',))
        self._summary.configure(text=TEXTS['pairing_summary'].format(
            paired=len(pr.pairs), total=len(genes), missing=len(pr.missing) + len(pr.duplicates),
            orphans=len(pr.orphans)))
