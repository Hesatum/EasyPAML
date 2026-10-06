"""Tempo limite automático por modelo × tamanho do gene (docs/benchmark_tempos.md)."""

import pytest

from backend.timeouts import (MIN_SECONDS, SLACK, auto_timeout, estimate_seconds,
                              resolve_timeout, user_timeout)

# (modelo, táxons, códons, segundos medidos) -- amostra de docs/benchmark_tempos.md
MEASURED = [
    ('M0', 60, 1500, 2563),
    ('M1a', 30, 500, 392),
    ('M2a', 60, 1500, 11347),
    ('M7', 10, 150, 98),
    ('M7', 60, 500, 12182),
    ('M8', 30, 500, 1991),
    ('M8', 60, 500, 14231),
    ('M8a', 60, 500, 15972),
    ('Branch', 60, 1500, 2378),
    ('Branch-site', 60, 150, 3674),
    ('Branch-site_null', 60, 1500, 13799),
]


@pytest.mark.parametrize('model,taxa,codons,seconds', MEASURED)
def test_limit_leaves_at_least_double_the_measured_time(model, taxa, codons, seconds):
    assert auto_timeout(model, taxa, codons) >= 2 * seconds


def test_limit_grows_with_taxa_and_codons():
    assert auto_timeout('M8', 60, 500) > auto_timeout('M8', 30, 500) > MIN_SECONDS
    assert auto_timeout('M8', 30, 1500) > auto_timeout('M8', 30, 500)


def test_small_genes_get_the_floor():
    assert auto_timeout('M0', 5, 100) == MIN_SECONDS


def test_slack_over_estimate():
    assert auto_timeout('M8', 60, 1500) == int(SLACK * estimate_seconds('M8', 60, 1500))


def test_unknown_model_uses_a_generous_reference():
    assert auto_timeout('M3', 30, 500) >= auto_timeout('M8a', 30, 500)
    assert auto_timeout('BranchSite_custom', 30, 500) == auto_timeout('Branch-site', 30, 500)


@pytest.mark.parametrize('value,expected', [(None, 0), ('', 0), (0, 0), (-5, 0), ('abc', 0),
                                            (120, 120), ('90', 90), (7.9, 7)])
def test_user_timeout_parsing(value, expected):
    assert user_timeout(value) == expected


def test_user_value_overrides_automatic():
    assert resolve_timeout(42, 'M8', 60, 1500) == 42
    assert resolve_timeout(None, 'M8', 60, 1500) == auto_timeout('M8', 60, 1500)
