"""Summary tab of the results panel: one test at a time, with a chart (2Δℓ density
or ω per gene) and a table of genes that can be sorted, exported and opened in the
Positive Sites tab."""

import tkinter as tk
from tkinter import ttk
from typing import Dict, List, Optional, Tuple

import customtkinter as ctk
import numpy as np
import pandas as pd

from src.backend import lrt_stats

from . import charts
from .gui_texts import TEXTS
from .ui_helpers import (CURRENT_THEME, FONT_MONO, FONT_SIZE, FONT_UI, PALETTE, RADIUS, SPACE,
                         ask_save_file, mix, show_message)

# order of the test selector: the stricter site test first
TEST_ORDER = [('M8a', 'M8'), ('M1a', 'M2a'), ('M7', 'M8'),
              ('Branch-site_null', 'Branch-site'), ('M0', 'Branch'), ('M0', 'M1a')]
SITE_ALTS = ('M2a', 'M8')
POSITIVE_PAIRS = (('M1a', 'M2a'), ('M7', 'M8'), ('M8a', 'M8'))

# journal formats: vector PDF/SVG, 600 dpi TIFF, 300 dpi PNG; 180 mm wide (two columns)
CHART_FORMATS = [("PDF (vector)", "*.pdf"), ("TIFF, 600 dpi", "*.tif"),
                 ("PNG, 300 dpi", "*.png"), ("SVG (vector)", "*.svg")]
CHART_SIZE_IN = (7.1, 3.2)


def test_label(pair) -> str:
    null, alt = pair
    if pair == ('Branch-site_null', 'Branch-site'):
        return "Branch-site"
    if pair == ('M0', 'Branch'):
        return "Branch vs M0"
    return f"{alt} vs {null.replace('_null', ' null')}"


def gene_verdict(sig_by_pair: dict, rows: list) -> Tuple[str, str]:
    """Kind of conclusion for a gene from its positive-selection tests: 'neutral',
    'no_m8a', 'none', 'weak' or 'positive', and the sentence shown for it.
    rows: (test, significant, p text, q text, effect, n sites, ω, p₁)."""
    if sig_by_pair.get(('M7', 'M8')) and sig_by_pair.get(('M8a', 'M8')) is False \
            and not sig_by_pair.get(('M1a', 'M2a')):
        return 'neutral', TEXTS["conclusion_neutral"]
    if sig_by_pair.get(('M7', 'M8')) and ('M8a', 'M8') not in sig_by_pair \
            and not sig_by_pair.get(('M1a', 'M2a')):
        return 'no_m8a', TEXTS["conclusion_no_m8a"]
    sig_tests = [t for t, sig, *_ in rows if sig]
    if not sig_tests:
        return 'none', TEXTS["conclusion_none"]
    n = max((r[5] for r in rows if r[1] and r[5] is not None), default=None)
    if n == 0:
        w = max((r[6] for r in rows if r[1] and len(r) > 6 and pd.notna(r[6])), default=np.nan)
        p1 = max((r[7] for r in rows if r[1] and len(r) > 7 and pd.notna(r[7])), default=np.nan)
        detail = (TEXTS["conclusion_weak_class"].format(w=f"{w:.3g}", p1=f"{100 * p1:.1f}")
                  if pd.notna(w) and pd.notna(p1) else "")
        return 'weak', TEXTS["conclusion_weak"].format(tests=", ".join(sig_tests), detail=detail)
    sites = TEXTS["conclusion_sites"].format(n=n) if n is not None else ""
    return 'positive', TEXTS["conclusion_supported"].format(tests=", ".join(sig_tests), sites=sites)


def ensure_tree_style(widget) -> None:
    """ttk style 'EP.Treeview' in the program's theme (tables of the results panel)."""
    style = ttk.Style(widget)
    rowh = int(FONT_SIZE['sm'] * 2.1)
    style.configure('EP.Treeview', background=PALETTE['bg_panel'], fieldbackground=PALETTE['bg_panel'],
                    foreground=PALETTE['text_primary'], rowheight=rowh, borderwidth=0,
                    font=(FONT_MONO, -FONT_SIZE['sm']))
    style.configure('EP.Treeview.Heading', background=PALETTE['bg_elevated'],
                    foreground=PALETTE['text_secondary'], relief='flat', borderwidth=0,
                    font=(FONT_UI, -FONT_SIZE['xs'], 'bold'), padding=(6, 6))
    style.map('EP.Treeview.Heading', background=[('active', PALETTE['bg_elevated_hover'])])
    style.map('EP.Treeview', background=[('selected', PALETTE['accent_fill'])],
              foreground=[('selected', '#ffffff')])
    style.layout('EP.Treeview', [('Treeview.treearea', {'sticky': 'nswe'})])


class SummaryTab:
    """Mixin for ResultsViewerWindow."""

    # ── data ─────────────────────────────────────────────────────────

    def _available_tests(self) -> List[tuple]:
        known = [p for p in TEST_ORDER if self._pair_values(*p)]
        extra = [p for p in lrt_stats.PAIRS if p not in TEST_ORDER and self._pair_values(*p)]
        return known + extra

    def _site_counts(self, model: str) -> Dict[str, int]:
        """Sites with Pr(ω>1) ≥ 0.95 per gene for a model, from sites_BEB.tsv (one
        read for all genes); genes absent from it have none."""
        cache = getattr(self, '_site_count_cache', None)
        if cache is None:
            cache = self._site_count_cache = {}
            f = self.output_folder / 'sites_BEB.tsv'
            if f.exists():
                try:
                    t = pd.read_csv(f, sep='\t', dtype={'gene': str})
                    t = t[t['significance'].isin(['*', '**'])]
                    for (m, g), n in t.groupby(['model', 'gene']).size().items():
                        cache.setdefault(m, {})[g] = int(n)
                except Exception:
                    cache['_unreadable'] = True
        if cache.get('_unreadable') or not (self.output_folder / 'sites_BEB.tsv').exists():
            return {g: len(s) for g in self.df['Gene']
                    for s in [self._beb_sites(g, model)] if s is not None}
        return cache.get(model, {})

    def _positive_class(self, row, alt: str):
        w, p1 = row.get(f'{alt}_w_pos'), row.get(f'{alt}_p_pos')
        if (pd.isna(w) or pd.isna(p1)) and alt in SITE_ALTS:
            rf = self._find_results_file(row['Gene'], alt)
            if rf:
                from src.backend.sites_parser import SitesParser
                pc = SitesParser.extract_positive_class(rf) or {}
                w, p1 = pc.get('omega', np.nan), pc.get('p', np.nan)
        return w, p1

    def _gene_verdicts(self) -> Dict[str, Tuple[str, str]]:
        """{gene: (kind, sentence)} over the positive-selection tests that ran."""
        if getattr(self, '_verdict_cache', None) is not None:
            return self._verdict_cache
        tests = [p for p in POSITIVE_PAIRS if self._pair_values(*p)]
        counts = {alt: self._site_counts(alt) for alt in {a for _, a in tests}}
        out = {}
        for _, row in self.df.iterrows():
            gene = row['Gene']
            if row.get('status') == 'failed':
                out[gene] = ('failed', self._compact_reason(self._failure_reason(gene)))
                continue
            sig_by_pair, rows = {}, []
            for null, alt in tests:
                vals = self._pair_values(null, alt).get(gene)
                if not vals:
                    continue
                _, p, q = vals
                sig = bool(pd.notna(q) and q < 0.05)
                sig_by_pair[(null, alt)] = sig
                w, p1 = self._positive_class(row, alt)
                rows.append((test_label((null, alt)), sig, lrt_stats.format_p(p), lrt_stats.format_p(q),
                             '', counts[alt].get(gene, 0) if sig else None, w, p1))
            if rows:
                out[gene] = gene_verdict(sig_by_pair, rows)
        self._verdict_cache = out
        return out

    def _branch_tag_omegas(self) -> Dict[str, Dict[str, float]]:
        cache = getattr(self, '_branch_tag_cache', None)
        if cache is None:
            from src.backend.sites_parser import SitesParser
            cache = self._branch_tag_cache = {}
            for gene in self.df['Gene']:
                rf = self._find_results_file(gene, 'Branch')
                if rf:
                    try:
                        cache[gene] = SitesParser.extract_omega_by_tags(rf)
                    except Exception:
                        pass
        return cache

    # ── layout ───────────────────────────────────────────────────────

    def _create_summary_tab(self, parent):
        tests = self._available_tests()
        if not tests:
            ctk.CTkLabel(parent, text=TEXTS["summary_no_tests"], font=self._font('md'),
                         text_color=PALETTE['text_secondary']).pack(pady=40)
            return
        self._sum = {'test': tests[0], 'chart': 'lrt', 'sort': ('q', False)}
        self._summary_answer(parent)

        bar = ctk.CTkFrame(parent, fg_color='transparent')
        bar.pack(fill='x', padx=SPACE['md'], pady=(SPACE['sm'], 0))
        ctk.CTkLabel(bar, text=TEXTS["summary_test"], font=self._font('sm', 'bold'),
                     text_color=PALETTE['text_secondary']).pack(side='left', padx=(0, SPACE['sm']))
        by_label = {test_label(p): p for p in tests}
        if len(tests) <= 4:
            seg = ctk.CTkSegmentedButton(bar, values=list(by_label), font=self._font('sm', 'bold'), height=28,
                                         selected_color=PALETTE['accent_fill'],
                                         selected_hover_color=PALETTE['accent_fill'],
                                         unselected_color=PALETTE['bg_elevated'],
                                         unselected_hover_color=PALETTE['bg_elevated_hover'],
                                         text_color=PALETTE['text_primary'],
                                         command=lambda lab: self._summary_select(by_label[lab]))
        else:       # many tests: a list keeps room for the export buttons
            seg = self._style_combo(ctk.CTkComboBox(bar, values=list(by_label), width=200, state='readonly',
                                                    command=lambda lab: self._summary_select(by_label[lab])))
        seg.set(test_label(tests[0]))
        seg.pack(side='left')
        for text, cmd in ((TEXTS["summary_export_html"], self._export_html),
                          (TEXTS["summary_export_all"], self._export_excel),
                          (TEXTS["summary_export_table"], self._summary_export_table),
                          (TEXTS["summary_copy_genes"], self._summary_copy_genes)):
            self._ghost_button(bar, text, cmd).pack(side='right', padx=(SPACE['sm'], 0))
        hyp = ctk.CTkFrame(parent, fg_color='transparent')
        hyp.pack(fill='x', padx=SPACE['md'], pady=(SPACE['xs'], SPACE['xs']))
        self._sum_hyp = ctk.CTkLabel(hyp, text="", font=self._font('sm'), anchor='w', justify='left',
                                     text_color=PALETTE['text_secondary'], wraplength=1100)
        self._sum_hyp.pack(side='left')
        info = ctk.CTkLabel(hyp, text=" ? ", font=self._font('xs', 'bold'), corner_radius=8,
                            fg_color=PALETTE['bg_elevated'], text_color=PALETTE['text_primary'], cursor='hand2')
        info.pack(side='left', padx=(SPACE['sm'], 0))
        explain = ctk.CTkLabel(parent, text=TEXTS["summary_explain"], font=self._font('sm'), anchor='w',
                               justify='left', wraplength=1150, fg_color=PALETTE['bg_elevated'],
                               corner_radius=RADIUS['field'], text_color=PALETTE['text_primary'],
                               padx=SPACE['sm'], pady=SPACE['xs'])

        def _toggle(_e=None):
            if explain.winfo_ismapped():
                explain.pack_forget()
            else:
                explain.pack(fill='x', padx=SPACE['md'], pady=(0, SPACE['sm']), after=hyp)
        info.bind('<Button-1>', _toggle)
        explain.bind('<Button-1>', _toggle)

        self._summary_chart_card(parent)
        self._summary_table(parent)
        self._summary_select(tests[0])

    def _summary_answer_text(self):
        """(sentence, any_supported) over all positive-selection tests, from the
        same rules as the Evidence column; None when no such test ran."""
        verdicts = self._gene_verdicts()
        if not verdicts:
            return None
        by_kind: Dict[str, List[str]] = {}
        for gene, (kind, _) in verdicts.items():
            by_kind.setdefault(kind, []).append(gene)
        rank = getattr(self, '_gene_rank', {})
        for genes in by_kind.values():
            genes.sort(key=lambda g: (rank.get(g, len(rank)), g))
        counts = self._site_counts('M8') or self._site_counts('M2a')

        def names(genes, with_sites=False):
            """Up to three genes by name; more than that only as a count."""
            if len(genes) > 3:
                return ""
            return ": " + ", ".join(f"{g} ({TEXTS['answer_sites'].format(n=counts.get(g, 0))})"
                                    if with_sites else g for g in genes)

        pos = by_kind.get('positive', [])
        parts = [TEXTS['answer_positive'].format(n=len(pos), total=len(verdicts), genes=names(pos, True))
                 if pos else TEXTS['answer_none'].format(total=len(verdicts))]
        for kind in ('neutral', 'no_m8a', 'weak'):
            genes = by_kind.get(kind, [])
            if genes:
                parts.append(TEXTS[f'answer_{kind}'].format(n=len(genes), genes=names(genes)))
        return " ".join(parts), bool(pos)

    def _summary_answer(self, parent):
        found = self._summary_answer_text()
        if not found:
            return
        text, supported = found
        tone = PALETTE['success_fg'] if supported else PALETTE['text_secondary']
        box = ctk.CTkFrame(parent, fg_color=mix(PALETTE['bg_panel'], tone, 0.10), corner_radius=RADIUS['card'])
        box.pack(fill='x', padx=SPACE['md'], pady=(SPACE['sm'], 0))
        label = ctk.CTkLabel(box, text=text, font=self._font('sm', 'bold'), anchor='w', justify='left',
                             wraplength=1150, text_color=PALETTE['text_primary'])
        label.pack(fill='x', padx=SPACE['md'], pady=SPACE['xs'])
        box.bind('<Configure>', lambda e: label.configure(wraplength=max(300, e.width - 2 * SPACE['md'])), add='+')

    def _summary_chart_card(self, parent):
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg as _Canvas
        from matplotlib.figure import Figure as _Figure
        card = ctk.CTkFrame(parent, fg_color=PALETTE['bg_panel'], corner_radius=RADIUS['card'])
        card.pack(fill='x', padx=SPACE['md'], pady=(0, SPACE['sm']))
        top = ctk.CTkFrame(card, fg_color='transparent')
        top.pack(fill='x', padx=SPACE['md'], pady=(SPACE['sm'], 0))
        self._sum_chart_seg = ctk.CTkSegmentedButton(
            top, values=[TEXTS["chart_kind_lrt"], TEXTS["chart_kind_omega"]], height=24,
            font=self._font('xs', 'bold'), selected_color=PALETTE['accent_fill'],
            selected_hover_color=PALETTE['accent_fill'], unselected_color=PALETTE['bg_elevated'],
            unselected_hover_color=PALETTE['bg_elevated_hover'], text_color=PALETTE['text_primary'],
            command=lambda v: self._summary_draw('omega' if v == TEXTS["chart_kind_omega"] else 'lrt'))
        self._sum_chart_seg.set(TEXTS["chart_kind_lrt"])
        self._sum_chart_seg.pack(side='left')
        ctk.CTkLabel(top, text=TEXTS["chart_hint"], font=self._font('xs'),
                     text_color=PALETTE['text_tertiary']).pack(side='left', padx=SPACE['md'])
        self._ghost_button(top, TEXTS["chart_export"], self._summary_export_chart).pack(side='right')

        c = charts.DARK if CURRENT_THEME['mode'] == 'dark' else charts.LIGHT
        self._sum_colors = dict(c, bg=PALETTE['bg_panel'], box=PALETTE['bg_elevated'])
        fig = _Figure(figsize=(12, 1.8), facecolor=self._sum_colors['bg'])
        self._sum_canvas = _Canvas(fig, master=card)
        widget = self._sum_canvas.get_tk_widget()
        widget.configure(height=165, highlightthickness=0, bg=self._sum_colors['bg'])
        widget.pack(fill='x', padx=SPACE['sm'], pady=(0, SPACE['sm']))
        self._sum_hover = None

    def _ghost_button(self, parent, text, command):
        return ctk.CTkButton(parent, text=text, command=command, height=28, fg_color='transparent',
                             border_width=1, border_color=PALETTE['control_border'],
                             text_color=PALETTE['text_primary'], hover_color=PALETTE['bg_elevated'],
                             corner_radius=RADIUS['field'], font=self._font('sm', 'bold'))

    def _summary_table(self, parent):
        ensure_tree_style(self)
        # the detail line stays visible; the table takes what is left
        detail = ctk.CTkFrame(parent, fg_color='transparent')
        detail.pack(side='bottom', fill='x', padx=SPACE['md'], pady=(0, SPACE['sm']))
        box = ctk.CTkFrame(parent, fg_color=PALETTE['bg_panel'], corner_radius=RADIUS['card'])
        box.pack(fill='both', expand=True, padx=SPACE['md'], pady=(0, SPACE['xs']))
        inner = tk.Frame(box, bg=PALETTE['bg_panel'])
        inner.pack(fill='both', expand=True, padx=SPACE['sm'], pady=SPACE['sm'])
        tree = ttk.Treeview(inner, style='EP.Treeview', show='headings', selectmode='browse')
        sb = ctk.CTkScrollbar(inner, command=tree.yview)
        hb = ctk.CTkScrollbar(inner, command=tree.xview, orientation='horizontal', height=12)
        tree.configure(yscrollcommand=sb.set, xscrollcommand=hb.set)
        hb.pack(side='bottom', fill='x')
        sb.pack(side='right', fill='y')
        tree.pack(side='left', fill='both', expand=True)
        for tag, color in (('positive', PALETTE['success_fg']), ('weak', PALETTE['warning_fg']),
                           ('neutral', PALETTE['warning_fg']), ('no_m8a', PALETTE['warning_fg']),
                           ('failed', PALETTE['danger_fg']), ('sig', PALETTE['success_fg']),
                           ('none', PALETTE['text_secondary'])):
            tree.tag_configure(tag, foreground=color)
        tree.tag_configure('odd', background=PALETTE['row_alt'])
        tree.bind('<<TreeviewSelect>>', lambda e: self._summary_show_detail())
        tree.bind('<Double-1>', lambda e: self._summary_open_sites())
        self._sum_tree = tree
        self._summary_heading_hints(tree)

        self._sum_detail = ctk.CTkLabel(detail, text=TEXTS["summary_pick_gene"], font=self._font('sm'),
                                        anchor='w', justify='left', wraplength=900,
                                        text_color=PALETTE['text_tertiary'])
        self._sum_detail.pack(side='left', fill='x', expand=True)
        self._sum_sites_btn = self._ghost_button(detail, TEXTS["summary_open_sites"], self._summary_open_sites)


    def _summary_heading_hints(self, tree):
        tip = {'w': None, 'col': None}

        def hide(_e=None):
            if tip['w'] is not None:
                tip['w'].destroy()
            tip['w'] = tip['col'] = None

        def move(e):
            if tree.identify_region(e.x, e.y) != 'heading':
                hide()
                return
            col = tree.identify_column(e.x)
            key = tree.column(col, 'id') if col else None
            text = TEXTS["summary_col_hints"].get(key, '')
            if key == tip['col']:
                return
            hide()
            if not text:
                return
            w = tk.Toplevel(tree)
            w.overrideredirect(True)
            tk.Label(w, text=text, justify='left', wraplength=320, bg=PALETTE['bg_elevated'],
                     fg=PALETTE['text_primary'], font=(FONT_UI, -FONT_SIZE['sm']), padx=8, pady=5).pack()
            w.geometry(f"+{e.x_root + 12}+{e.y_root + 18}")
            tip['w'], tip['col'] = w, key
        tree.bind('<Motion>', move, add='+')
        tree.bind('<Leave>', hide, add='+')

    # ── behaviour ────────────────────────────────────────────────────

    def _summary_columns(self, pair) -> List[Tuple[str, str, int, str]]:
        """(id, title, width, anchor) of the table for a test."""
        null, alt = pair
        cols = [('gene', TEXTS["col_gene"], 132, 'w')]
        if pair in POSITIVE_PAIRS:
            cols.append(('conclusion', TEXTS["col_conclusion"], 195, 'w'))
        else:
            cols.append(('result', TEXTS["col_result"], 132, 'w'))
        cols += [('q', 'q (BH)', 74, 'e'), ('p', 'p', 74, 'e'), ('lrt', '2Δℓ', 66, 'e')]
        if pair == ('M0', 'Branch'):
            cols.append(('df', 'df', 40, 'e'))
        cols.append(('mean_w', TEXTS["col_mean_w"].format(model=alt), 92 + 8 * len(alt), 'e'))
        if alt in SITE_ALTS:
            cols += [('w_pos', TEXTS["col_w_pos"], 122, 'e'), ('p1', 'p₁', 51, 'e'),
                     ('sites', TEXTS["col_sites"], 95, 'e')]
        if pair == ('M0', 'Branch'):
            cols.append(('w_tags', TEXTS["col_w_tags"], 220, 'w'))
        for cid, model in (('lnl0', null.replace('_null', ' null')), ('lnl1', alt)):
            cols.append((cid, f"lnL {model}", max(88, 8 * len(model) + 36), 'e'))
        return cols

    def _summary_rows(self, pair) -> List[dict]:
        null, alt = pair
        vals = self._pair_values(null, alt) or {}
        verdicts = self._gene_verdicts() if pair in POSITIVE_PAIRS else {}
        counts = self._site_counts(alt) if alt in SITE_ALTS else {}
        tags = self._branch_tag_omegas() if pair == ('M0', 'Branch') else {}
        rows = []
        for _, row in self.df.iterrows():
            gene = row['Gene']
            failed = row.get('status') == 'failed'
            if gene not in vals and not failed:
                continue
            lrt, p, q = vals.get(gene, (np.nan, np.nan, np.nan))
            sig = bool(pd.notna(q) and q < 0.05)
            r = {'gene': gene, 'q_num': q if pd.notna(q) else 2.0, 'failed': failed, 'sig': sig,
                 'result': TEXTS["result_failed"] if failed else
                 (TEXTS["summary_verdict_sig"] if sig else TEXTS["summary_verdict_nonsig"]),
                 'q': '' if failed else lrt_stats.format_p(q), 'p': '' if failed else lrt_stats.format_p(p),
                 'lrt': '' if failed else f"{max(0.0, lrt):.2f}",
                 'lnl0': self._fmt_num(row.get(f'{null}_lnL'), 3),
                 'lnl1': self._fmt_num(row.get(f'{alt}_lnL'), 3),
                 'mean_w': self._fmt_num(row.get(f'{alt}_omega'), 3)}
            kind, sentence = verdicts.get(gene, ('failed' if failed else ('sig' if sig else 'none'), ''))
            r['kind'], r['sentence'] = kind, sentence
            r['conclusion'] = TEXTS["conclusion_short"].get(kind, '')
            if alt in SITE_ALTS:
                w, p1 = self._positive_class(row, alt)
                r['w_pos'], r['p1'] = self._fmt_num(w, 3), self._fmt_num(p1, 3)
                r['sites'] = str(counts.get(gene, 0)) if sig else ''
            if pair == ('M0', 'Branch'):
                # Branch uses the rooted labelled tree: one branch length more than M0
                np_b, np_0 = row.get('Branch_np'), row.get('M0_np')
                nt_b, nt_0 = row.get('Branch_ntime'), row.get('M0_ntime')
                if pd.notna(np_b) and pd.notna(np_0):
                    extra = (nt_b - nt_0) if pd.notna(nt_b) and pd.notna(nt_0) else 0
                    r['df'] = str(max(1, int(np_b - np_0 - extra)))
                else:
                    r['df'] = ''
                tg = tags.get(gene, {})
                r['w_tags'] = "  ".join(f"{'bg' if k == 'background' else k} {v:.3g}"
                                        for k, v in sorted(tg.items(), key=lambda kv: (kv[0] != 'background', kv[0])))
            if failed:
                r['sentence'] = self._compact_reason(self._failure_reason(gene)).replace("\n", "; ")
            notes = self._gene_notes(gene)
            r['warn'] = bool(notes)
            if notes:
                r['sentence'] = (r['sentence'] + "   ⚠ " + notes.replace(' | ', '; ')).strip()
            rows.append(r)
        return rows

    @staticmethod
    def _fmt_num(v, digits: int) -> str:
        try:
            return '' if v is None or pd.isna(v) else f"{float(v):.{digits}f}"
        except (TypeError, ValueError):
            return ''

    def _summary_select(self, pair):
        self._sum['test'] = pair
        info = lrt_stats.PAIRS.get(pair, {})
        df_txt = (f"df = {1 if info.get('boundary') else info.get('df')}" if info.get('df')
                  else TEXTS["summary_df_branch"])
        self._sum_hyp.configure(text=TEXTS["test_hypotheses"].get(pair, test_label(pair)) + "   ·   " + df_txt)
        has_omega = pair[1] in SITE_ALTS
        self._sum_chart_seg.configure(state='normal' if has_omega else 'disabled')
        if not has_omega:
            self._sum['chart'] = 'lrt'
            self._sum_chart_seg.set(TEXTS["chart_kind_lrt"])
        self._summary_fill_table()
        self._summary_draw(self._sum['chart'])

    def _summary_fill_table(self):
        tree, pair = self._sum_tree, self._sum['test']
        cols = self._summary_columns(pair)
        tree.delete(*tree.get_children())
        tree.configure(columns=[c[0] for c in cols])
        for cid, title, width, anchor in cols:
            tree.heading(cid, text=title, anchor=anchor, command=lambda c=cid: self._summary_sort(c))
            tree.column(cid, width=width, minwidth=width if cid == 'gene' else 40, anchor=anchor, stretch=cid == 'gene')
        rows = self._summary_rows(pair)
        key, rev = self._sum['sort']
        num = {'q': 'q_num'}.get(key)
        if num:
            rows.sort(key=lambda r: r[num], reverse=rev)
        elif key:
            def _k(r):
                v = r.get(key, '')
                try:
                    return (0, float(v))
                except (TypeError, ValueError):
                    return (1, str(v))
            rows.sort(key=_k, reverse=rev)
        rows.sort(key=lambda r: r['failed'])
        self._sum_rows = {}
        for i, r in enumerate(rows):
            kind = r['kind'] if pair in POSITIVE_PAIRS else ('failed' if r['failed'] else ('sig' if r['sig'] else 'none'))
            values = [r.get(c[0], '') for c in cols]
            if r['warn']:
                values[0] = "⚠ " + values[0]
            iid = tree.insert('', 'end', values=values,
                              tags=(kind,) + (('odd',) if i % 2 else ()))
            self._sum_rows[iid] = r
        hint = TEXTS["summary_pick_gene"]
        if any(r['warn'] for r in rows):
            hint += "   " + TEXTS["summary_warn_hint"]
        self._sum_detail.configure(text=hint, text_color=PALETTE['text_tertiary'])
        self._sum_sites_btn.pack_forget()

    def _summary_sort(self, col):
        key, rev = self._sum['sort']
        self._sum['sort'] = (col, not rev if key == col else False)
        self._summary_fill_table()

    def _summary_show_detail(self):
        sel = self._sum_tree.selection()
        if not sel:
            return
        r = self._sum_rows.get(sel[0])
        if not r:
            return
        color = {'positive': PALETTE['success_fg'], 'failed': PALETTE['danger_fg'],
                 'none': PALETTE['text_secondary'], 'sig': PALETTE['success_fg']}.get(r['kind'], PALETTE['warning_fg'])
        text = f"{r['gene']}:  {r['sentence'] or r['result']}"
        self._sum_detail.configure(text=text, text_color=color)
        if self._sum['test'][1] in SITE_ALTS + ('Branch-site',) and not r['failed']:
            self._sum_sites_btn.pack(side='right')
        else:
            self._sum_sites_btn.pack_forget()

    def _summary_open_sites(self):
        sel = self._sum_tree.selection()
        r = self._sum_rows.get(sel[0]) if sel else None
        if getattr(self, '_sites_goto', None) is None:
            self._ensure_tab(TEXTS["viewer_tab_sites"])
        goto = getattr(self, '_sites_goto', None)
        if r and goto and self._sum['test'][1] in SITE_ALTS + ('Branch-site',):
            goto(r['gene'], self._sum['test'][1])

    def _summary_points(self, pair):
        null, alt = pair
        info = lrt_stats.PAIRS.get(pair, {})
        vals = self._pair_values(null, alt) or {}
        df_ = 1 if info.get('boundary') else (info.get('df') or 1)
        test = charts.TestData(null, alt, df_)
        for gene, (lrt, p, q) in vals.items():
            test.points.append(charts.Point(gene, float(lrt), bool(pd.notna(q) and q < 0.05),
                                            f"{gene}  2Δℓ {max(0.0, lrt):.2f}  q {lrt_stats.format_p(q)}"))
        whole, positive = [], []
        if alt in SITE_ALTS:
            for _, row in self.df.iterrows():
                gene = row['Gene']
                if row.get('status') == 'failed' or gene not in vals:
                    continue
                q = vals[gene][2]
                sig = bool(pd.notna(q) and q < 0.05)
                mw = row.get(f'{alt}_omega')
                if pd.notna(mw):
                    whole.append(charts.Point(gene, float(mw), sig, f"{gene}  mean ω {float(mw):.3f} ({alt})"))
                w, p1 = self._positive_class(row, alt)
                if pd.notna(w):
                    positive.append(charts.Point(gene, float(w), sig, f"{gene}  ω {float(w):.2f}"
                                                 + (f"  p₁ {float(p1):.3f}" if pd.notna(p1) else "")
                                                 + f"  q {lrt_stats.format_p(q)}"))
        return test, whole, positive

    def _summary_draw(self, kind):
        self._sum['chart'] = kind
        fig, c = self._sum_canvas.figure, self._sum_colors
        if self._sum_hover is not None:
            self._sum_canvas.mpl_disconnect(self._sum_hover)
        fig.clear()
        ax = fig.add_subplot(1, 1, 1)
        test, whole, positive = self._summary_points(self._sum['test'])
        if kind == 'omega' and (whole or positive):
            pts = charts.draw_omega(ax, whole, positive, c, compact=True,
                                    whole_label=TEXTS["chart_mean_w"].format(model=self._sum['test'][1]),
                                    positive_label=TEXTS["chart_positive_class"])
        else:
            pts = charts.draw_lrt(ax, test, c, compact=True)
        fig.subplots_adjust(left=0.05, right=0.985, top=0.9, bottom=0.27)
        self._sum_hover = charts.enable_hover(fig, ax, pts, c)
        self._sum_canvas.draw_idle()

    def _summary_export_chart(self):
        pair = self._sum['test']
        kind = self._sum['chart']
        name = f"{'omega' if kind == 'omega' else 'LRT'}_{test_label(pair).replace(' ', '_')}"
        path = ask_save_file(self, TEXTS["chart_export"], self.output_folder, initialfile=f"{name}.pdf",
                             defaultextension='.pdf', filetypes=CHART_FORMATS)
        if not path:
            return
        try:
            from matplotlib.figure import Figure as _Figure
            c = charts.LIGHT
            fig = _Figure(figsize=CHART_SIZE_IN, facecolor=c['bg'])
            ax = fig.add_subplot(1, 1, 1)
            test, whole, positive = self._summary_points(pair)
            if kind == 'omega' and (whole or positive):
                charts.draw_omega(ax, whole, positive, c,
                                  whole_label=TEXTS["chart_mean_w"].format(model=pair[1]),
                                  positive_label=TEXTS["chart_positive_class"])
            else:
                charts.draw_lrt(ax, test, c)
            ext = path.lower().rsplit('.', 1)[-1]
            kw = {'dpi': 600 if ext in ('tif', 'tiff') else 300, 'facecolor': c['bg'], 'bbox_inches': 'tight'}
            if ext in ('tif', 'tiff'):
                kw['pil_kwargs'] = {'compression': 'tiff_lzw'}
            fig.savefig(path, **kw)
            show_message(self, TEXTS["msg_success"], TEXTS["msg_exported_to"].format(path=path))
        except Exception as e:
            show_message(self, TEXTS["msg_error"], TEXTS["msg_export_err"].format(error=e), 'error')

    def _summary_export_table(self):
        null, alt = self._sum['test']
        col = lrt_stats.lrt_column(null, alt)
        name = f"EasyPAML_{test_label((null, alt)).replace(' ', '_')}"
        path = ask_save_file(self, TEXTS["summary_export_table"], self.output_folder,
                             initialfile=f"{name}.csv", defaultextension='.csv',
                             filetypes=[("CSV", "*.csv"), ("Excel", "*.xlsx")])
        if not path:
            return
        try:
            out = self._build_export_df(col)
            verdicts = self._gene_verdicts() if (null, alt) in POSITIVE_PAIRS else {}
            if verdicts:
                out.insert(1, 'evidence_all_tests', [TEXTS["conclusion_short"].get(verdicts.get(g, ('',))[0], '')
                                              for g in out['Gene']])
            if path.lower().endswith('.xlsx'):
                out.to_excel(path, index=False, sheet_name=test_label((null, alt))[:31])
            else:
                out.to_csv(path, index=False, encoding='utf-8-sig')
            show_message(self, TEXTS["msg_success"], TEXTS["msg_exported_to"].format(path=path))
        except Exception as e:
            show_message(self, TEXTS["msg_error"], TEXTS["msg_export_err"].format(error=e), 'error')

    def _summary_copy_genes(self):
        """Genes significant in the selected test (positive selection only for site
        tests), one per line, for enrichment tools such as g:Profiler or PANTHER."""
        pair = self._sum['test']
        rows = self._summary_rows(pair)
        keep = [r['gene'] for r in rows if r['sig'] and not r['failed']
                and (pair not in POSITIVE_PAIRS or r['kind'] == 'positive')]
        self.clipboard_clear()
        self.clipboard_append("\n".join(keep))
        show_message(self, TEXTS["msg_success"], TEXTS["summary_genes_copied"].format(n=len(keep)))
