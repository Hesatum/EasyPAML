"""
Pares de modelos do LRT e cálculo de p e q -- fonte única usada pelo
backend (LRT_results.txt, analysis_summary.tsv) e pelo painel de resultados,
para os dois nunca divergirem.

p-valores via chi2.sf (função de sobrevivência): 1 - chi2.cdf(x) vira 0.0
por arredondamento quando x é grande (o painel mostrava "p = 0"); sf dá o
valor certo até ~1e-300.
"""

from typing import Dict, List, Optional, Sequence

import numpy as np
from scipy import stats

# (nulo, alternativo) -> df e descrição. 'boundary' marca os testes cujo nulo
# fica na fronteira do espaço de parâmetros (ω = 1 fixo): para eles a
# significância usa χ²₁ puro, como o manual do PAML recomenda para o
# branch-site, e a mistura 50:50 χ²₀/χ²₁ é reportada só como referência.
PAIRS = {
    ('M0', 'M1a'): {'df': 2, 'boundary': False},
    ('M1a', 'M2a'): {'df': 2, 'boundary': False},
    ('M7', 'M8'): {'df': 2, 'boundary': False},
    ('M8a', 'M8'): {'df': 1, 'boundary': True},
    ('M0', 'Branch'): {'df': None, 'boundary': False},   # df depende das marcas
    ('Branch-site_null', 'Branch-site'): {'df': 1, 'boundary': True},
}

# Pares que testam seleção positiva por sítio (têm tabela de sítios BEB)
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
    """p do LRT. lrt <= 0 (alternativo não melhorou) -> 1.0."""
    if lrt is None or not np.isfinite(lrt) or lrt <= 0:
        return 1.0
    return float(stats.chi2.sf(lrt, df=1 if boundary else df))


def p_value_mixture(lrt: float) -> float:
    """Mistura 50:50 χ²₀/χ²₁ (Self & Liang 1987) -- só referência."""
    if lrt is None or not np.isfinite(lrt) or lrt <= 0:
        return 1.0
    return float(0.5 * stats.chi2.sf(lrt, df=1))


def bh_qvalues(pvals: Sequence[float]) -> List[float]:
    """Benjamini-Hochberg dentro da família (um par de modelos, todos os genes)."""
    if len(pvals) == 0:
        return []
    return [float(q) for q in stats.false_discovery_control(list(pvals), method='bh')]


def format_p(p: Optional[float]) -> str:
    """Notação científica legível: 4.7e-22, 0.031, 1."""
    if p is None or (isinstance(p, float) and not np.isfinite(p)):
        return "NA"
    if p == 0:
        return "< 1e-300"
    if p >= 0.001:
        return f"{p:.3g}"
    return f"{p:.2e}"


def format_p_unicode(p: Optional[float]) -> str:
    """4,7×10⁻²² -- para a frase-resumo da interface."""
    if p is None or (isinstance(p, float) and not np.isfinite(p)):
        return "NA"
    if p == 0:
        return "< 10⁻³⁰⁰"
    if p >= 0.001:
        return f"{p:.3g}"
    mant, exp = f"{p:.1e}".split('e')
    sup = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    return f"{mant}×10{str(int(exp)).translate(sup)}"
