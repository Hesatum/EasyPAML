"""methods_text.txt and the version commit."""
from src.backend.methods_text import build_methods_text
from src.backend.version import source_commit, version_string


def _text(models, sizes):
    return build_methods_text(
        version='0.3.0.dev0 (commit abc1234)', codeml_version='4.9j', models=models,
        ctl={'CodonFreq': 2, 'ncatG': 10, 'kappa': 2, 'cleandata': 1}, omega0=0.5,
        pruned=True, family_sizes=sizes, n_genes=3)


def test_methods_text_has_versions_tests_and_family_sizes():
    t = _text(['M7', 'M8', 'M8a'], {('M7', 'M8'): 3, ('M8a', 'M8'): 3})
    assert 'EasyPAML 0.3.0.dev0 (commit abc1234; https://github.com/Hesatum/EasyPAML)' in t and 'PAML 4.9j' in t
    assert 'F3x4 (CodonFreq = 2)' in t and 'ncatG = 10' in t and 'cleandata = 1' in t
    assert 'M8 vs M7 (df = 2; 3 gene(s))' in t
    assert 'M8 vs M8a (df = 1, χ²₁' in t
    assert 'Benjamini-Hochberg' in t and 'Bayes Empirical Bayes' in t


def test_methods_text_without_m8a_does_not_invent_it():
    t = _text(['M7', 'M8'], {('M7', 'M8'): 3})
    assert 'M8a' not in t


def test_methods_text_reports_masked_stops_and_excluded_sequences():
    t = build_methods_text(
        version='0.3.0', codeml_version='4.9j', models=['M7', 'M8'],
        ctl={'CodonFreq': 2, 'ncatG': 10, 'kappa': 2, 'cleandata': 1}, omega0=0.5,
        pruned=True, family_sizes={}, n_genes=2,
        masked_stops={'geneB': 1}, excluded_taxa={'geneA': ['Macaca_mulata']})
    assert '1 stop codon(s) in 1 gene(s) (geneB) were treated as missing data' in t
    assert 'Macaca_mulata from geneA' in t


def test_version_string_has_commit_in_a_clone():
    c = source_commit()
    assert c is None or len(c.split('-')[0]) == 40
    if c:
        assert f"commit {c[:7]}" in version_string()


def test_methods_text_mentions_warm_start_only_when_used():
    assert 'starting values' not in _text(['M0', 'M7', 'M8'], {})
    t = build_methods_text(
        version='0.3.0', codeml_version='4.9j', models=['M0', 'M7', 'M8'],
        ctl={'CodonFreq': 2, 'ncatG': 10, 'kappa': 2, 'cleandata': 1}, omega0=0.5,
        pruned=True, family_sizes={}, n_genes=2, warm_start=True, omega_starts=(0.2, 1.0, 2.5))
    assert 'estimated under M0 were used as starting values' in t and '0.2, 1 and 2.5' in t
