"""Charts of a run: density of the LRT statistic (2Δℓ) of each test and the ω of each
gene, drawn on matplotlib axes. The results panel draws a compact version with a
hover box; the export draws one panel per test for a figure."""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np
from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter, LogLocator, MaxNLocator, NullFormatter
from scipy import stats

MIN_GENES_FOR_DENSITY = 5

LIGHT = {'bg': '#ffffff', 'text': '#1f2430', 'muted': '#555d6b', 'axis': '#9aa1ad',
         'line': '#08306b', 'fill': '#e8edf4', 'reject': '#e69f00', 'crit': '#b2182b',
         'sig': '#d55e00', 'ns': '#8c8c8c', 'box': '#ffffff', 'box_edge': '#c9ced6'}
DARK = {'bg': '#16161c', 'text': '#eeeef2', 'muted': '#a2a2b6', 'axis': '#5c5c6e',
        'line': '#9ecae1', 'fill': '#23232e', 'reject': '#e69f00', 'crit': '#f4a582',
        'sig': '#ffb000', 'ns': '#7d7d90', 'box': '#20202a', 'box_edge': '#3a3a48'}


@dataclass
class Point:
    gene: str
    x: float
    significant: bool
    hover: str = ''


@dataclass
class TestData:
    null: str
    alt: str
    df: int
    points: List[Point] = field(default_factory=list)

    @property
    def title(self) -> str:
        return f"{self.alt} vs {self.null}"


def _density(values: np.ndarray, xs: np.ndarray, log: bool = False) -> Optional[np.ndarray]:
    v = np.log10(values) if log else values
    if v.size < MIN_GENES_FOR_DENSITY or np.ptp(v) == 0:
        return None
    try:
        kde = stats.gaussian_kde(v)
    except (np.linalg.LinAlgError, ValueError):
        return None
    return kde(np.log10(xs) if log else xs)


def _style(ax, c, compact: bool):
    ax.set_facecolor(c['bg'])
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(c['axis'])
    ax.tick_params(colors=c['muted'], labelsize=8 if compact else 8.5)
    ax.xaxis.label.set_color(c['muted'])
    ax.yaxis.label.set_color(c['muted'])
    ax.yaxis.set_major_locator(MaxNLocator(3 if compact else 4))


def _dots(ax, rows, c, compact):
    """Too few genes for a density: one labelled dot per gene, on its own line."""
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    n = len(rows)
    to_axes = (ax.transData + ax.transAxes.inverted()).transform
    for i, (x, gene, color) in enumerate(rows):
        y = (n - i) / (n + 1)
        ax.plot([x], [y], 'o', color=color, ms=6 if compact else 6.5, zorder=4)
        right = to_axes((x, y))[0] > 0.7     # label on the left near the right edge
        ax.annotate(gene, (x, y), xytext=(-6 if right else 6, 0), textcoords='offset points',
                    va='center', ha='right' if right else 'left',
                    fontsize=7.5 if compact else 8, color=c['text'])
    ax.set_ylim(0, 1.25)


def _rug(ax, xs, y, color, compact, strong=False):
    if len(xs):
        ax.plot(xs, np.full(len(xs), y), '|', color=color,
                ms=(7 if strong else 6) if compact else (6 if strong else 5),
                mew=1.2 if strong else 0.8, alpha=1 if strong else 0.75, zorder=3)


def draw_lrt(ax, test: TestData, c: Dict[str, str], compact: bool = False,
             label: Optional[str] = None) -> List[Point]:
    """Density of 2Δℓ with the rejection region beyond the χ² 0.05 critical value
    and one tick per gene (significant = q < 0.05). Returns the points drawn."""
    _style(ax, c, compact)
    x = np.array([max(0.0, p.x) for p in test.points], float)
    crit05 = stats.chi2.ppf(0.95, test.df)
    crit01 = stats.chi2.ppf(0.99, test.df)
    hi = max(crit01 * 1.6, float(x.max()) * 1.05 if x.size else 1.0)
    xs = np.linspace(0, hi, 600)
    ys = _density(x, xs)

    ax.set_xlim(0, hi)
    ns = [p for p in test.points if not p.significant]
    sig = [p for p in test.points if p.significant]
    if ys is not None:
        ymax = float(ys.max())
        ax.fill_between(xs, ys, color=c['fill'], zorder=1, lw=0)
        rej = xs >= crit05
        ax.fill_between(xs[rej], ys[rej], color=c['reject'], alpha=0.45, lw=0, zorder=2)
        ax.plot(xs, ys, color=c['line'], lw=1.3, zorder=4)
        y_rug = -ymax * 0.09
        ax.set_ylim(y_rug * 1.8, ymax * 1.32)
        ax.axhline(0, color=c['axis'], lw=0.6, zorder=0)
        _rug(ax, [max(0.0, p.x) for p in ns], y_rug, c['ns'], compact)
        _rug(ax, [max(0.0, p.x) for p in sig], y_rug, c['sig'], compact, strong=True)
        top = ymax * 1.18
    else:
        _dots(ax, [(max(0.0, p.x), p.gene, c['sig'] if p.significant else c['ns'])
                   for p in sorted(test.points, key=lambda p: -p.x)], c, compact)
        ax.axvspan(crit05, hi, color=c['reject'], alpha=0.12, lw=0, zorder=0)
        top = 1.13

    ax.axvline(crit05, color=c['crit'], lw=1.1, ls='--', zorder=5)
    ax.annotate(f"χ²₀.₀₅ = {crit05:.2f}", (crit05, top), xytext=(4, 0), textcoords='offset points',
                color=c['crit'], fontsize=7.5, ha='left', va='bottom')
    if not compact:
        ax.axvline(crit01, color=c['crit'], lw=0.9, ls=':', zorder=5)

    n_sig = len(sig)
    info = f"{len(test.points)} gene(s)\n{n_sig} significant (q < 0.05)"
    ax.text(0.98, 0.95, info, transform=ax.transAxes, ha='right', va='top',
            fontsize=7.5 if compact else 8, color=c['text'], linespacing=1.4,
            bbox=dict(boxstyle='round,pad=0.4', fc=c['box'], ec=c['box_edge'], lw=0.6))
    ax.set_title(test.title + ("" if compact else f"  (df = {test.df})"), fontsize=9.5,
                 fontweight='bold', color=c['text'], pad=6)
    ax.set_xlabel("LRT statistic, 2Δℓ", fontsize=8.5)
    if not compact and ys is not None:
        ax.set_ylabel("density", fontsize=8.5)
    if label:
        ax.text(-0.09, 1.04, label, transform=ax.transAxes, fontsize=12, fontweight='bold',
                color=c['text'], va='bottom', ha='left')
    for p in test.points:
        p.x = max(0.0, p.x)
    return test.points


def draw_omega(ax, whole: Sequence[Point], positive: Sequence[Point], c: Dict[str, str],
               compact: bool = False, label: Optional[str] = None,
               whole_label: str = "whole gene", positive_label: str = "positive class"):
    """ω of each gene on a log axis: the gene-wide ω (grey) and the ω of the class
    that may be under positive selection (coloured by q < 0.05), with ω = 1 marked.
    Returns the points drawn, for the hover box."""
    _style(ax, c, compact)
    w = np.array([p.x for p in whole if p.x and p.x > 0], float)
    pc = np.array([p.x for p in positive if p.x and p.x > 0], float)
    allv = np.concatenate([w, pc]) if (w.size or pc.size) else np.array([0.1, 10.0])
    lo = max(1e-3, min(float(allv.min()) / 1.6, 0.5))
    hi = max(float(allv.max()) * 1.6, 2.0)
    xs = np.logspace(np.log10(lo), np.log10(hi), 600)
    ax.set_xscale('log')

    ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.axvline(1.0, color=c['crit'], lw=1.0, ls='--', zorder=5)
    yw, yp = _density(w, xs, log=True), _density(pc, xs, log=True)
    if yw is None and yp is None:
        ax.set_xlim(lo, hi)
        rows = [(p.x, f"{p.gene} ({whole_label})", c['ns']) for p in whole if p.x and p.x > 0]
        rows += [(p.x, f"{p.gene} ({positive_label})", c['sig'] if p.significant else c['ns'])
                 for p in positive if p.x and p.x > 0]
        _dots(ax, rows, c, compact)
        ax.annotate("ω = 1 (neutral)", (1.0, 1.13), xytext=(4, 0), textcoords='offset points',
                    color=c['crit'], fontsize=7.5, ha='left', va='bottom')
        ax.set_title("ω per gene", fontsize=9.5, fontweight='bold', color=c['text'], pad=6)
        ax.set_xlabel("ω (dN/dS), log scale", fontsize=8.5)
        if label:
            ax.text(-0.09, 1.04, label, transform=ax.transAxes, fontsize=12, fontweight='bold',
                    color=c['text'], va='bottom', ha='left')
        return list(whole) + list(positive)
    peaks = [float(y.max()) for y in (yw, yp) if y is not None]
    ymax = max(peaks) if peaks else 1.0
    if yw is not None:
        ax.fill_between(xs, yw, color=c['fill'], lw=0, zorder=1)
        ax.plot(xs, yw, color=c['ns'], lw=1.1, zorder=3, label=whole_label)
    if yp is not None:
        ax.plot(xs, yp, color=c['sig'], lw=1.4, zorder=4, label=positive_label)
    if not peaks:
        ax.set_yticks([])

    y_w, y_p = -ymax * 0.08, -ymax * 0.19
    ax.set_ylim(y_p * 1.6, ymax * 1.3)
    ax.set_xlim(lo, hi)
    ax.axhline(0, color=c['axis'], lw=0.6, zorder=0)
    ax.annotate("ω = 1 (neutral)", (1.0, ymax * 1.16), xytext=(4, 0), textcoords='offset points',
                color=c['crit'], fontsize=7.5, ha='left', va='bottom')

    _rug(ax, [p.x for p in whole if p.x and p.x > 0], y_w, c['ns'], compact)
    _rug(ax, [p.x for p in positive if p.x and p.x > 0 and not p.significant], y_p, c['ns'], compact)
    _rug(ax, [p.x for p in positive if p.x and p.x > 0 and p.significant], y_p, c['sig'], compact,
         strong=True)
    if not compact:
        ax.text(lo, y_w, f" {whole_label}", fontsize=7, color=c['muted'], va='center', ha='left')
        ax.text(lo, y_p, f" {positive_label}", fontsize=7, color=c['muted'], va='center', ha='left')
    if peaks:
        leg = ax.legend(loc='upper right', fontsize=7.5, frameon=True, facecolor=c['box'],
                        edgecolor=c['box_edge'], labelcolor=c['text'])
        leg.get_frame().set_linewidth(0.6)
    ax.set_title("ω per gene", fontsize=9.5, fontweight='bold', color=c['text'], pad=6)
    ax.set_xlabel("ω (dN/dS), log scale", fontsize=8.5)
    if not compact:
        ax.set_ylabel("density", fontsize=8.5)
    if label:
        ax.text(-0.09, 1.04, label, transform=ax.transAxes, fontsize=12, fontweight='bold',
                color=c['text'], va='bottom', ha='left')
    return list(whole) + list(positive)


def enable_hover(fig: Figure, ax, points: Sequence[Point], c: Dict[str, str],
                 pixels: float = 6, max_lines: int = 10) -> int:
    """Show a box with the genes near the mouse along the x axis. Returns the
    callback id, for mpl_disconnect when the axes are redrawn."""
    box = ax.annotate("", xy=(0, 0), xycoords='axes fraction', xytext=(0, 0),
                      textcoords='offset points', fontsize=8, color=c['text'], zorder=20,
                      annotation_clip=False, linespacing=1.35,
                      bbox=dict(boxstyle='round,pad=0.5', fc=c['box'], ec=c['box_edge'], lw=0.8))
    box.set_visible(False)
    guide = ax.axvline(0, color=c['muted'], lw=0.6, alpha=0.6, zorder=6)
    guide.set_visible(False)
    pts = [p for p in points if p.x is not None and np.isfinite(p.x) and p.x > (0 if ax.get_xscale() == 'log' else -1)]

    def on_move(event):
        if event.inaxes is not ax or event.x is None:
            if box.get_visible():
                box.set_visible(False)
                guide.set_visible(False)
                fig.canvas.draw_idle()
            return
        to_px = ax.transData.transform
        near = sorted(((abs(to_px((p.x, 0))[0] - event.x), p) for p in pts),
                      key=lambda t: t[0])
        near = [p for d, p in near if d <= pixels]
        if not near:
            if box.get_visible():
                box.set_visible(False)
                guide.set_visible(False)
                fig.canvas.draw_idle()
            return
        lines = [p.hover or p.gene for p in near[:max_lines]]
        if len(near) > max_lines:
            lines.append(f"+{len(near) - max_lines} more")
        box.set_text("\n".join(lines))
        fx, fy = ax.transAxes.inverted().transform((event.x, event.y))
        box.xy = (fx, fy)
        box.set_position((-12, 10) if fx > 0.55 else (12, 10))
        box.set_horizontalalignment('right' if fx > 0.55 else 'left')
        box.set_verticalalignment('bottom' if fy < 0.5 else 'top')
        guide.set_xdata([near[0].x, near[0].x])
        box.set_visible(True)
        guide.set_visible(True)
        fig.canvas.draw_idle()

    return fig.canvas.mpl_connect('motion_notify_event', on_move)


def export_figure(path, tests: Sequence[TestData], whole: Sequence[Point],
                  positive: Sequence[Point]) -> None:
    """One panel per test and an ω panel, two columns, for a manuscript figure."""
    c = LIGHT
    has_omega = bool(whole or positive)
    cols = min(2, max(1, len(tests)))
    lrt_rows = (len(tests) + cols - 1) // cols
    rows = lrt_rows + (1 if has_omega else 0)
    fig = Figure(figsize=(5.3 * cols, 3.3 * rows), facecolor=c['bg'])
    grid = fig.add_gridspec(rows, cols, hspace=0.6, wspace=0.28)
    letters = "ABCDEFGHIJ"
    for i, test in enumerate(tests):
        draw_lrt(fig.add_subplot(grid[i // cols, i % cols]), test, c, label=letters[i])
    if has_omega:      # full width under the tests
        draw_omega(fig.add_subplot(grid[rows - 1, :]), whole, positive, c,
                   label=letters[len(tests)])
    dpi = 600 if str(path).lower().endswith(('.tif', '.tiff')) else 300
    fig.savefig(path, dpi=dpi, facecolor=c['bg'], bbox_inches='tight')
