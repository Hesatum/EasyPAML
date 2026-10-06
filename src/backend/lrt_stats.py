"""
LRT model pairs and p/q computation, shared by the backend
(LRT_results.txt, analysis_summary.tsv) and the results panel.

p-values use chi2.sf: 1 - chi2.cdf(x) rounds to 0.0 for large x, sf stays
accurate down to about 1e-300.
"""

from typing import Dict, List, Optional, Sequence

import numpy as np
from scipy import stats

# (null, alternative) -> df. 'boundary' marks tests whose null lies on the
# boundary of the parameter space (ω = 1 fixed): significance uses χ²₁, as the
# PAML manual recommends for branch-site, and the 50:50 χ²₀/χ²₁ mixture is
# reported for reference.
PAIRS = {
    # M0 has ω; M1a has p0 and ω0 (ω1 = 1 fixed): one more parameter. M1a
    # reduces to M0 at p0 = 1, on the boundary.
    ('M0', 'M1a'): {'df': 1, 'boundary': True},
    ('M1a', 'M2a'): {'df': 2, 'boundary': False},
    ('M7', 'M8'): {'df': 2, 'boundary': False},
    ('M8a', 'M8'): {'df': 1, 'boundary': True},
    ('M0', 'Branch'): {'df': None, 'boundary': False},   # df depende das marcas
    ('Branch-site_null', 'Branch-site'): {'df': 1, 'boundary': True},
}

# pairs that test site-wise positive selection (they have a BEB site table)
POSITIVE_SELECTION_PAIRS = [('M1a', 'M2a'), ('M7', 'M8'), ('M8a', 'M8')]


def pairs_for(models: Sequence[str]):
    present = set(models)
    return [p for p in PAIRS if p[0] in present and p[1] in present]


def lrt_column(null: str, alt: str) -> str:
    return f"lrt_{null}_vs_{alt}"


def q_column(null: str, alt: str) -> str:
    return f"q_{null}_vs_{alt}"


def p_column(null: str, alt: str) -> str:
    return f"p_{null}_vs_{alt}"


def p_value(lrt: float, df: int, boundary: bool = False) -> float:
    """LRT p-value; 1.0 when lrt <= 0 (the alternative did not improve)."""
    if lrt is None or not np.isfinite(lrt) or lrt <= 0:
        return 1.0
    return float(stats.chi2.sf(lrt, df=1 if boundary else df))


def p_value_mixture(lrt: float) -> float:
    """50:50 χ²₀/χ²₁ mixture (Self & Liang 1987), for reference only."""
    if lrt is None or not np.isfinite(lrt) or lrt <= 0:
        return 1.0
    return float(0.5 * stats.chi2.sf(lrt, df=1))


def bh_qvalues(pvals: Sequence[float]) -> List[float]:
    """Benjamini-Hochberg within one family (one model pair, all genes)."""
    if len(pvals) == 0:
        return []
    return [float(q) for q in stats.false_discovery_control(list(pvals), method='bh')]


def format_p(p: Optional[float]) -> str:
    """Readable notation: 4.7e-22, 0.031, 1."""
    if p is None or (isinstance(p, float) and not np.isfinite(p)):
        return "NA"
    if p == 0:
        return "< 1e-300"
    if p >= 0.001:
        return f"{p:.3g}"
    return f"{p:.2e}"


def format_p_unicode(p: Optional[float]) -> str:
    """4.7×10⁻²², for the summary line in the interface."""
    if p is None or (isinstance(p, float) and not np.isfinite(p)):
        return "NA"
    if p == 0:
        return "< 10⁻³⁰⁰"
    if p >= 0.001:
        return f"{p:.3g}"
    mant, exp = f"{p:.1e}".split('e')
    sup = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    return f"{mant}×10{str(int(exp)).translate(sup)}"
