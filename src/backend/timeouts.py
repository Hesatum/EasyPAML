"""
Tempo limite de cada execução do codeml, proporcional ao tamanho do gene.

O tempo do codeml cresce com o número de táxons (forte) e de códons (fraco).
Medido em docs/benchmark_tempos.md (codeml 4.9j, dados simulados com o
evolverNSsites, 10/30/60 táxons × 150/500/1500 códons) e conferido em 25 genes
reais de exemplos_teste:

    tempo ≈ T_ref[modelo] × (táxons / 30)^2,73 × (códons / 500)^0,71

O limite é SLACK vezes essa estimativa, nunca menos que MIN_SECONDS. Quem
protege contra um codeml travado é a detecção de inatividade (idle_timeout);
o limite total só corta otimizações que não terminam.
"""

from typing import Optional

# Segundos para 30 táxons × 500 códons (ajuste de mínimos quadrados em log)
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
    """Tempo esperado do codeml (mediana do ajuste) para um gene × modelo."""
    ref = T_REF.get(model)
    if ref is None:
        ref = T_REF['Branch-site'] if str(model).startswith('Branch') else max(T_REF.values())
    taxa = max(int(n_taxa), 2)
    codons = max(int(n_codons), 1)
    return ref * (taxa / REF_TAXA) ** EXP_TAXA * (codons / REF_CODONS) ** EXP_CODONS


def auto_timeout(model: str, n_taxa: int, n_codons: int) -> int:
    return int(max(MIN_SECONDS, SLACK * estimate_seconds(model, n_taxa, n_codons)))


def user_timeout(value) -> int:
    """Segundos escolhidos pelo usuário; 0 quando vazio, None, 0 ou inválido (= automático)."""
    try:
        v = float(value) if value not in (None, '') else 0.0
    except (TypeError, ValueError):
        return 0
    return int(v) if v > 0 else 0


def resolve_timeout(user_value: Optional[float], model: str, n_taxa: int, n_codons: int) -> int:
    """Valor do usuário (> 0) tem prioridade; senão, o limite automático."""
    return user_timeout(user_value) or auto_timeout(model, n_taxa, n_codons)
