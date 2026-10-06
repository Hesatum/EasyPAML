"""
Time limit of each codeml run, scaled to the model and gene size
(docs/timing_benchmark.md):

    time ≈ T_ref[model] × (taxa / 30)^2.73 × (codons / 500)^0.71

The limit is SLACK times this estimate and at least MIN_SECONDS. A stuck
codeml is caught by the idle check (idle_timeout); this limit only stops
optimizations that never finish.
"""

from typing import Optional

# seconds for 30 taxa × 500 codons (least-squares fit on the log scale)
T_REF = {
    'M0': 140,
    'M1a': 340,
    'M2a': 530,
    'M7': 1890,
    'M8': 1950,
    'M8a': 2400,
    'Branch': 160,
    'Branch-site': 890,
    'Branch-site_null': 810,
}
REF_TAXA = 30
REF_CODONS = 500
EXP_TAXA = 2.73
EXP_CODONS = 0.71
SLACK = 5.0
MIN_SECONDS = 1800


def estimate_seconds(model: str, n_taxa: int, n_codons: int) -> float:
    """Expected codeml time (fit median) for one gene and model."""
    ref = T_REF.get(model)
    if ref is None:
        ref = T_REF['Branch-site'] if str(model).startswith('Branch') else max(T_REF.values())
    taxa = max(int(n_taxa), 2)
    codons = max(int(n_codons), 1)
    return ref * (taxa / REF_TAXA) ** EXP_TAXA * (codons / REF_CODONS) ** EXP_CODONS


def auto_timeout(model: str, n_taxa: int, n_codons: int) -> int:
    return int(max(MIN_SECONDS, SLACK * estimate_seconds(model, n_taxa, n_codons)))


def user_timeout(value) -> int:
    """Seconds chosen by the user; 0 (automatic) when empty, None, 0 or invalid."""
    try:
        v = float(value) if value not in (None, '') else 0.0
    except (TypeError, ValueError):
        return 0
    return int(v) if v > 0 else 0


def resolve_timeout(user_value: Optional[float], model: str, n_taxa: int, n_codons: int) -> int:
    """The user value (> 0) wins; otherwise the automatic limit."""
    return user_timeout(user_value) or auto_timeout(model, n_taxa, n_codons)
