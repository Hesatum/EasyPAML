"""methods_text.txt e commit da versão (rodada 2 do teste de usabilidade: o
revisor não conseguia citar a versão exata nem dizer quantos genes entraram
na correção BH)."""
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


def test_version_string_has_commit_in_a_clone():
    c = source_commit()
    assert c is None or len(c.split('-')[0]) == 40
    if c:
        assert f"commit {c[:7]}" in version_string()


def test_interpretation_candidates_use_panel_q_and_m8a(tmp_path):
    """Aba Interpretação (rodada 2): usava max(LRT) com df = 2, ignorava o
    M8a×M8 e contava genes que falharam."""
    from src.backend.go_enrichment import rank_candidates
    tsv = tmp_path / 'analysis_summary.tsv'
    tsv.write_text(
        "Gene\tstatus\tp_M7_vs_M8\tq_M7_vs_M8\tp_M8a_vs_M8\tq_M8a_vs_M8\n"
        "gA\tok\t1e-10\t1e-9\t1e-8\t1e-7\n"        # candidato pelo M8a×M8
        "gB\tok\t1e-10\t1e-9\t0.4\t0.6\n"           # só M7×M8: sítios neutros, não é candidato
        "gC\tfailed\t1e-20\t1e-19\t1e-20\t1e-19\n"  # falhou: fora
    )
    ann = tmp_path / 'go.tsv'
    ann.write_text("gene_id_full\tgo_biological_process\tgo_cellular_component\tgo_molecular_function\n"
                   "gA\t\t\t\ngB\t\t\t\ngC\t\t\t\n")
    cand, _ = rank_candidates(tsv, ann)
    assert list(cand['Gene']) == ['gA']
    assert cand.iloc[0]['test'] == 'M8 vs M8a' and abs(cand.iloc[0]['q_value'] - 1e-7) < 1e-12
