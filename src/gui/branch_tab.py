"""Branch tab of the results panel: for one gene, the tree with the labelled branches
in colour and their ω (Branch model), the Branch and Branch-site tests, and a table
of the branch groups."""

import re
import tkinter as tk
from pathlib import Path
from tkinter import ttk
from typing import Dict, List, Optional

import customtkinter as ctk
import numpy as np
import pandas as pd

from src.backend import lrt_stats

from . import charts
from .gui_texts import TEXTS
from .ui_helpers import (CURRENT_THEME, FONT_SIZE, FONT_UI, PALETTE, RADIUS, SPACE,
                         ask_save_file, show_message)
from .summary_tab import CHART_FORMATS, ensure_tree_style

# background, then #1, #2 … (Okabe-Ito, distinguishable with colour blindness)
LABEL_COLORS = ['#d55e00', '#0072b2', '#009e73', '#cc79a7', '#e69f00', '#56b4e9']


class Node:
    __slots__ = ('name', 'label', 'w', 'children', 'clade')

    def __init__(self):
        self.name, self.label, self.w, self.children, self.clade = '', None, None, [], False


def parse_newick(text: str) -> Optional[Node]:
    """Newick with codeml's extras: '#1'/'$1' branch labels and ' #0.61' ω values
    (from 'w ratios as labels'). Branch lengths are ignored."""
    s = text.strip()
    s = s[s.find('('):] if '(' in s else s
    s = s.rstrip().rstrip(';')
    pos = 0

    def node() -> Node:
        nonlocal pos
        n = Node()
        while pos < len(s) and s[pos].isspace():
            pos += 1
        if pos < len(s) and s[pos] == '(':
            pos += 1
            n.children.append(node())
            while pos < len(s) and s[pos] == ',':
                pos += 1
                n.children.append(node())
            while pos < len(s) and s[pos].isspace():
                pos += 1
            pos += 1                       # ')'
        start = pos
        while pos < len(s) and s[pos] not in ',)':
            pos += 1
        raw = s[start:pos].strip()
        raw = re.sub(r':[-\d.eE]+', '', raw)
        m = re.search(r'#\s*([\d.]+(?:[eE][-+]?\d+)?)', raw)
        if m and '.' in m.group(1):
            n.w = float(m.group(1))
            raw = raw[:m.start()] + raw[m.end():]
        lab = re.search(r'([#$])(\d+)\b', raw)
        if lab:
            n.label, n.clade = int(lab.group(2)), lab.group(1) == '$'
            raw = raw[:lab.start()] + raw[lab.end():]
        n.name = raw.strip().strip("'")
        return n
    try:
        root = node()
    except Exception:
        return None
    _spread_clade_labels(root, None)
    return root


def _spread_clade_labels(n: Node, inherited) -> None:
    """'$1' labels a whole clade in codeml ('#1' only the branch it sits on)."""
    if n.label is None and inherited is not None:
        n.label = inherited
    for c in n.children:
        _spread_clade_labels(c, n.label if n.clade else inherited)


def preorder(n: Node) -> List[Node]:
    out = [n]
    for c in n.children:
        out.extend(preorder(c))
    return out


class BranchTab:
    """Mixin for ResultsViewerWindow."""

    def _branch_genes(self) -> List[str]:
        have = set()
        for model in ('Branch', 'Branch-site'):
            col = f'{model}_lnL'
            if col in self.df.columns:
                have |= set(self.df.loc[self.df[col].notna(), 'Gene'])
        rank = getattr(self, '_gene_rank', {})
        qs = {g: v[2] for g, v in (self._pair_values('M0', 'Branch') or {}).items()}
        return sorted(have, key=lambda g: (qs.get(g, 2.0) if pd.notna(qs.get(g, 2.0)) else 2.0,
                                           rank.get(g, 0)))

    def _branch_tree(self, gene: str) -> Optional[Node]:
        """Tree with labels (from the tree codeml read) and ω per branch (Branch model)."""
        labelled = None
        for model in ('Branch', 'Branch-site'):
            rf = self._find_results_file(gene, model)
            if rf:
                tf = Path(rf).with_name(Path(rf).name.replace('_results.txt', '_tree.nwk'))
                if tf.exists():
                    labelled = parse_newick(tf.read_text(encoding='utf-8', errors='ignore'))
                    if labelled:
                        break
        rf = self._find_results_file(gene, 'Branch')
        omega = None
        if rf:
            text = Path(rf).read_text(encoding='utf-8', errors='ignore')
            m = re.search(r'w ratios as labels for TreeView:\s*\n(.+?;)', text)
            if m:
                omega = parse_newick(m.group(1))
        if labelled is None:
            return omega
        if omega is not None:
            a, b = preorder(labelled), preorder(omega)
            if len(a) == len(b):
                for x, y in zip(a, b):
                    x.w = y.w
        return labelled

    def _branchsite_info(self, gene: str) -> Dict:
        rf = self._find_results_file(gene, 'Branch-site')
        if not rf:
            return {}
        text = Path(rf).read_text(encoding='utf-8', errors='ignore')
        prop = re.search(r'^proportion\s+([\d.\s]+)$', text, re.M)
        fg = re.search(r'^foreground w\s+([\d.\s]+)$', text, re.M)
        out = {}
        if prop and fg:
            p, w = [float(x) for x in prop.group(1).split()], [float(x) for x in fg.group(1).split()]
            if len(p) >= 4 and len(w) >= 4:
                out['p2'] = p[2] + p[3]
                out['w2'] = w[2]
        counts = self._site_counts('Branch-site') if hasattr(self, '_site_counts') else {}
        out['sites'] = counts.get(gene, 0)
        return out

    def _create_tree_tab(self, parent):
        genes = self._branch_genes()
        if not genes:
            box = ctk.CTkFrame(parent, fg_color='transparent')
            box.pack(expand=True)
            ctk.CTkLabel(box, text=TEXTS["branch_no_data_title"], font=self._font('lg', 'bold'),
                         text_color=PALETTE['text_secondary']).pack(pady=(60, 6))
            ctk.CTkLabel(box, text=TEXTS["branch_no_data_hint"], font=self._font('sm'),
                         text_color=PALETTE['text_tertiary']).pack()
            return
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg as _Canvas
        from matplotlib.figure import Figure as _Figure

        bar = ctk.CTkFrame(parent, fg_color='transparent')
        bar.pack(fill='x', padx=SPACE['md'], pady=(SPACE['sm'], SPACE['xs']))
        ctk.CTkLabel(bar, text=TEXTS["branch_label_gene"], font=self._font('sm', 'bold'),
                     text_color=PALETTE['text_secondary']).pack(side='left', padx=(0, SPACE['sm']))
        combo = self._style_combo(ctk.CTkComboBox(bar, values=genes, width=320,
                                                  command=lambda g: self._branch_show(g)))
        combo.set(genes[0])
        combo.pack(side='left')
        self._ghost_button(bar, TEXTS["branch_export"], self._branch_export).pack(side='right')
        self._br_tests = ctk.CTkLabel(parent, text="", font=self._font('sm'), anchor='w', justify='left',
                                      text_color=PALETTE['text_secondary'], wraplength=1200)
        self._br_tests.pack(fill='x', padx=SPACE['md'], pady=(0, SPACE['xs']))

        body = ctk.CTkFrame(parent, fg_color='transparent')
        body.pack(fill='both', expand=True, padx=SPACE['md'], pady=(0, SPACE['sm']))
        side = ctk.CTkFrame(body, fg_color=PALETTE['bg_panel'], corner_radius=RADIUS['card'], width=360)
        side.pack(side='right', fill='y', padx=(SPACE['sm'], 0))
        side.pack_propagate(False)
        card = ctk.CTkFrame(body, fg_color=PALETTE['bg_panel'], corner_radius=RADIUS['card'])
        card.pack(side='left', fill='both', expand=True)

        c = charts.DARK if CURRENT_THEME['mode'] == 'dark' else charts.LIGHT
        self._br_colors = dict(c, bg=PALETTE['bg_panel'])
        fig = _Figure(figsize=(8, 4.5), facecolor=self._br_colors['bg'])
        self._br_canvas = _Canvas(fig, master=card)
        w = self._br_canvas.get_tk_widget()
        w.configure(bg=self._br_colors['bg'], highlightthickness=0)    # no white flash
        w.pack(fill='both', expand=True, padx=SPACE['sm'], pady=SPACE['sm'])

        ctk.CTkLabel(side, text=TEXTS["branch_groups"], font=self._font('sm', 'bold'),
                     text_color=PALETTE['text_primary'], anchor='w').pack(fill='x', padx=SPACE['md'],
                                                                         pady=(SPACE['md'], SPACE['xs']))
        inner = tk.Frame(side, bg=PALETTE['bg_panel'])
        inner.pack(fill='x', padx=SPACE['sm'])
        ensure_tree_style(self)
        tree = ttk.Treeview(inner, style='EP.Treeview', show='headings', height=6,
                            columns=('group', 'n', 'w'))
        for cid, title, width, anchor in (('group', TEXTS["branch_col_group"], 125, 'w'),
                                          ('n', TEXTS["branch_col_n"], 85, 'e'),
                                          ('w', TEXTS["branch_col_w"], 100, 'e')):
            tree.heading(cid, text=title, anchor=anchor)
            tree.column(cid, width=width, anchor=anchor, stretch=False)
        tree.tag_configure('odd', background=PALETTE['row_alt'])
        tree.pack(fill='x')
        self._br_table = tree
        self._br_bs = ctk.CTkLabel(side, text="", font=self._font('sm'), anchor='nw', justify='left',
                                   wraplength=330, text_color=PALETTE['text_secondary'])
        self._br_bs.pack(fill='x', padx=SPACE['md'], pady=(SPACE['md'], 0))
        self._branch_show(genes[0])

    def _branch_show(self, gene: str) -> None:
        self._br_gene = gene
        parts = []
        for pair in (('M0', 'Branch'), ('Branch-site_null', 'Branch-site')):
            v = (self._pair_values(*pair) or {}).get(gene)
            if v:
                lrt, p, q = v
                sig = pd.notna(q) and q < 0.05
                parts.append(f"{pair[1]} vs {pair[0].replace('_null', ' null')}: 2Δℓ {max(0.0, lrt):.2f}, "
                             f"p {lrt_stats.format_p(p)}, q {lrt_stats.format_p(q)}"
                             f" — {TEXTS['summary_verdict_sig'] if sig else TEXTS['summary_verdict_nonsig']}")
        self._br_tests.configure(text="   ·   ".join(parts))

        root = self._branch_tree(gene)
        fig = self._br_canvas.figure
        fig.clear()
        ax = fig.add_subplot(1, 1, 1)
        groups = self._draw_branch_tree(ax, root, self._br_colors) if root else {}
        fig.subplots_adjust(left=0.02, right=0.98, top=0.96, bottom=0.04)
        self._br_canvas.draw_idle()

        t = self._br_table
        t.delete(*t.get_children())
        for i, (lab, (n, w)) in enumerate(sorted(groups.items(), key=lambda kv: (kv[0] != 0, kv[0]))):
            name = TEXTS["branch_background"] if lab == 0 else f"#{lab}"
            t.insert('', 'end', values=(name, n, '' if w is None else f"{w:.4f}"),
                     tags=(('odd',) if i % 2 else ()))
        bs = self._branchsite_info(gene)
        if bs:
            self._br_bs.configure(text=TEXTS["branch_bs_info"].format(
                p2=f"{100 * bs['p2']:.1f}" if 'p2' in bs else '—',
                w2=f"{bs['w2']:.3g}" if 'w2' in bs else '—', n=bs.get('sites', 0)))
        else:
            self._br_bs.configure(text="")

    def _draw_branch_tree(self, ax, root: Node, c: Dict) -> Dict[int, tuple]:
        """Rectangular cladogram, tips aligned; background branches grey, labelled
        branches in the colour of their label with ω written on them. Returns
        {label: (n branches, ω)} (label 0 = background)."""
        tips = [n for n in preorder(root) if not n.children]
        y = {}
        for i, n in enumerate(tips):
            y[id(n)] = len(tips) - i

        def depth(n):
            return 0 if not n.children else 1 + max(depth(ch) for ch in n.children)
        total = depth(root)
        x = {}

        def place(n, level):
            if n.children:
                for ch in n.children:
                    place(ch, level + 1)
                y[id(n)] = sum(y[id(ch)] for ch in n.children) / len(n.children)
                x[id(n)] = level
            else:
                x[id(n)] = total
        place(root, 0)
        ax.set_facecolor(c['bg'])
        groups: Dict[int, list] = {}

        def color(lab):
            return c['ns'] if not lab else LABEL_COLORS[(lab - 1) % len(LABEL_COLORS)]

        def draw(n, parent=None):
            if parent is not None:
                lab = n.label or 0
                groups.setdefault(lab, []).append(n.w)
                col, lw = color(lab), (2.6 if lab else 1.4)
                px, py, nx, ny = x[id(parent)], y[id(parent)], x[id(n)], y[id(n)]
                ax.plot([px, px], [py, ny], color=c['ns'], lw=1.2, solid_capstyle='round', zorder=1)
                ax.plot([px, nx], [ny, ny], color=col, lw=lw, solid_capstyle='round', zorder=2 if lab else 1)
                if lab and n.w is not None:
                    ax.text((px + nx) / 2, ny + 0.18, f"#{lab}  ω {n.w:.3g}", color=col, fontsize=8,
                            ha='center', va='bottom', fontweight='bold')
            for ch in n.children:
                draw(ch, n)
        draw(root)
        for n in tips:
            ax.text(total + 0.12, y[id(n)], n.name, va='center', ha='left', fontsize=8.5, color=c['text'])
        ax.set_xlim(-0.2, total + 2.6)
        ax.set_ylim(0.3, len(tips) + 0.8)
        ax.axis('off')
        summary = {}
        for lab, ws in groups.items():
            vals = [w for w in ws if w is not None]
            summary[lab] = (len(ws), float(np.median(vals)) if vals else None)
        handles = []
        from matplotlib.lines import Line2D
        for lab in sorted(summary, key=lambda k: (k != 0, k)):
            n, w = summary[lab]
            name = TEXTS["branch_background"] if lab == 0 else f"#{lab}"
            handles.append(Line2D([0], [0], color=color(lab), lw=2.6 if lab else 1.4,
                                  label=f"{name}" + (f"   ω = {w:.3g}" if w is not None else "")))
        if handles:
            leg = ax.legend(handles=handles, loc='lower left', fontsize=8, frameon=True,
                            facecolor=c.get('box', c['bg']), edgecolor=c['box_edge'], labelcolor=c['text'])
            leg.get_frame().set_linewidth(0.6)
        return summary

    def _branch_export(self):
        gene = getattr(self, '_br_gene', None)
        if not gene:
            return
        path = ask_save_file(self, TEXTS["branch_export"], self.output_folder,
                             initialfile=f"{gene}_branches.pdf", defaultextension='.pdf',
                             filetypes=CHART_FORMATS)
        if not path:
            return
        try:
            from matplotlib.figure import Figure as _Figure
            root = self._branch_tree(gene)
            fig = _Figure(figsize=(7.1, 0.3 * max(10, len([n for n in preorder(root) if not n.children])) + 1),
                          facecolor='white')
            ax = fig.add_subplot(1, 1, 1)
            self._draw_branch_tree(ax, root, charts.LIGHT)
            ax.set_title(gene, fontsize=10, fontweight='bold', loc='left')
            ext = path.lower().rsplit('.', 1)[-1]
            kw = {'dpi': 600 if ext in ('tif', 'tiff') else 300, 'facecolor': 'white', 'bbox_inches': 'tight'}
            if ext in ('tif', 'tiff'):
                kw['pil_kwargs'] = {'compression': 'tiff_lzw'}
            fig.savefig(path, **kw)
            show_message(self, TEXTS["msg_success"], TEXTS["msg_exported_to"].format(path=path))
        except Exception as e:
            show_message(self, TEXTS["msg_error"], TEXTS["msg_export_err"].format(error=e), 'error')
