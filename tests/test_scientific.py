"""Item 1 -- riscos científicos: CodonFreq, .ctl explícito, M8a, LRT."""
import math
import re
from pathlib import Path

from src.backend import lrt_stats
from src.backend.codeml_backend import CodemlBatchAnalysis
from src.backend.ctl_params import (CODONFREQ_OPTIONS, DEFAULT_CODONFREQ, codonfreq_label,
                                    parse_codonfreq_label, parse_ctl_text)

ROOT = Path(__file__).resolve().parent.parent


def test_codonfreq_codes_follow_paml():
    names = {v: n for v, n, _ in CODONFREQ_OPTIONS}
    assert names[2] == 'F3x4'
    assert names[7] == 'FMutSel'
    assert names[0] == 'Fequal' and names[1] == 'F1x4' and names[3] == 'F61'
    assert codonfreq_label(2) == '2 = F3x4'
    assert codonfreq_label(7) == '7 = FMutSel'
    assert parse_codonfreq_label('7 = FMutSel') == 7
    assert parse_codonfreq_label('2') == 2


def test_default_codonfreq_is_f3x4_everywhere():
    assert DEFAULT_CODONFREQ == 2
    for name, cfg in CodemlBatchAnalysis.MODEL_CONFIGS.items():
        assert cfg['CodonFreq'] == 2, name


def test_no_wrong_codonfreq_label_left_in_sources():
    """O rótulo '7=F3×4' (errado) não pode voltar em nenhum lugar."""
    bad = re.compile(r'7\s*=\s*F3\s*[x×]\s*4', re.IGNORECASE)
    for path in list((ROOT / 'src').rglob('*.py')) + [ROOT / 'README.md']:
        assert not bad.search(path.read_text(encoding='utf-8')), path


def _ctl(model, **cfg):
    app = CodemlBatchAnalysis()
    app.config = cfg
    text = app.generate_ctl_content('g.fasta', 'g.nwk', 'out.txt',
                                    CodemlBatchAnalysis.MODEL_CONFIGS[model], model_name=model)
    return parse_ctl_text(text)


def test_ctl_writes_every_relevant_parameter():
    ctl = _ctl('M8')
    for key in ('seqfile', 'treefile', 'outfile', 'seqtype', 'CodonFreq', 'estFreq',
                'model', 'NSsites', 'icode', 'fix_kappa', 'kappa', 'fix_omega', 'omega',
                'fix_alpha', 'alpha', 'ncatG', 'cleandata', 'fix_blength', 'method',
                'Small_Diff', 'getSE', 'RateAncestor', 'clock'):
        assert key in ctl, key
    assert ctl['ncatG'] == '10'
    assert ctl['CodonFreq'] == '2'
    assert ctl['fix_kappa'] == '0' and float(ctl['kappa']) == 2.0
    assert ctl['fix_blength'] == '0'
    assert ctl['method'] == '0'
    assert ctl['NSsites'] == '8' and ctl['fix_omega'] == '0'


def test_ctl_global_overrides():
    ctl = _ctl('M7', CodonFreq=7, ncatG=4)
    assert ctl['CodonFreq'] == '7' and ctl['ncatG'] == '4'


def test_m8a_is_m8_with_omega_fixed_at_one():
    ctl = _ctl('M8a')
    assert ctl['NSsites'] == '8'
    assert ctl['fix_omega'] == '1'
    assert float(ctl['omega']) == 1.0


def test_warm_start_uses_branch_lengths_as_initial_values_not_fixed():
    app = CodemlBatchAnalysis()
    app.config = {}
    text = app.generate_ctl_content('g.fasta', 'g.nwk', 'o.txt', CodemlBatchAnalysis.MODEL_CONFIGS['M8'],
                                    model_name='M8', kappa=2.5, fix_blength=1)
    ctl = parse_ctl_text(text)
    assert ctl['fix_blength'] == '1'          # 1 = initial (pamlDOC); 2 = fixed
    assert float(ctl['kappa']) == 2.5
    src = (ROOT / 'src' / 'backend' / 'codeml_backend.py').read_text(encoding='utf-8')
    assert 'fix_bl = 2' not in src


def test_m8_auto_adds_m7_and_m8a():
    models = CodemlBatchAnalysis.auto_complete_null_models(['M8'])
    assert set(models) == {'M8', 'M7', 'M8a'}


def test_lrt_pairs_and_df():
    assert lrt_stats.PAIRS[('M7', 'M8')]['df'] == 2
    assert lrt_stats.PAIRS[('M8a', 'M8')]['df'] == 1
    assert lrt_stats.PAIRS[('M8a', 'M8')]['boundary'] is True
    assert set(lrt_stats.pairs_for(['M7', 'M8', 'M8a'])) == {('M7', 'M8'), ('M8a', 'M8')}


def test_p_value_does_not_underflow_to_zero():
    p = lrt_stats.p_value(98.23, 2)
    assert 0 < p < 1e-20
    assert math.isclose(p, math.exp(-98.23 / 2), rel_tol=1e-9)   # chi2 df=2 tem forma fechada
    assert lrt_stats.p_value(-1.0, 2) == 1.0
    assert lrt_stats.format_p(p).startswith('4.') and 'e-22' in lrt_stats.format_p(p)
    assert lrt_stats.format_p_unicode(p) == '4.7×10⁻²²'


def test_boundary_pair_uses_chi2_1_and_reports_mixture():
    p = lrt_stats.p_value(3.0, 1, boundary=True)
    assert math.isclose(lrt_stats.p_value_mixture(3.0), p / 2)


def test_method_note_fixed():
    note = CodemlBatchAnalysis.LRT_METHOD_NOTE
    assert 'ntime=0' not in note
    assert 'M8a vs M8' in note


def test_m7_m8_references():
    info = CodemlBatchAnalysis.MODEL_INFO
    assert '2000' in info['M7']['references'] and '2005' not in info['M7']['references']
    assert '2000' in info['M8']['references']
    assert 'M8a' in info


def test_ctl_comment_never_glued_to_long_values():
    """Bug real: com nome longo o '*' do comentário colava no valor
    ('..._seq.fasta*') e o codeml não achava o arquivo."""
    app = CodemlBatchAnalysis()
    app.config = {}
    long = '25_PHOT2__phototropin2_chloroplast_avoidance_high_light_M8a_seq.fasta'
    text = app.generate_ctl_content(long, 't.nwk', 'o.txt', CodemlBatchAnalysis.MODEL_CONFIGS['M8'],
                                    model_name='M8')
    line = next(l for l in text.splitlines() if l.strip().startswith('seqfile'))
    assert f"{long}   *" in line
    assert parse_ctl_text(text)['seqfile'] == long


def test_regenerated_summary_keeps_tiny_p_values(tmp_path):
    """Regenerar o TSV não pode transformar p = 4e-23 em 0.000000."""
    import pandas as pd
    for model, lnl, np_ in (('M7', -4174.199719, 20), ('M8', -4122.628318, 22)):
        d = tmp_path / model
        d.mkdir()
        (d / f'g_{model}_results.txt').write_text(
            f"CODONML\nns =  10  ls = 300\nlnL(ntime: 17  np: {np_}):  {lnl}      +0.000000\n")
    CodemlBatchAnalysis.regenerate_summary_files(tmp_path)
    df = pd.read_csv(tmp_path / 'analysis_summary.tsv', sep='\t')
    p = float(df.loc[0, 'p_M7_vs_M8'])
    assert 1e-24 < p < 1e-22
    assert float(df.loc[0, 'q_M7_vs_M8']) > 0


def test_m8a_can_be_left_out():
    """Item 2 (2ª rodada): opção para não adicionar o M8a; o padrão continua com ele."""
    assert set(CodemlBatchAnalysis.auto_complete_null_models(['M8'], include_m8a=False)) == {'M8', 'M7'}
    assert set(CodemlBatchAnalysis.auto_complete_null_models(['M8'])) == {'M8', 'M7', 'M8a'}
    # M8a escolhido à mão continua mesmo com a opção desligada
    assert 'M8a' in CodemlBatchAnalysis.auto_complete_null_models(['M8', 'M8a'], include_m8a=False)
    cli = (ROOT / 'easypaml_cli.py').read_text(encoding='utf-8')
    assert '--no-m8a' in cli


def test_lrt_degrees_of_freedom():
    """Diferença de parâmetros livres entre os modelos aninhados (rodada 2 do
    teste de usabilidade: M0 vs M1a estava com df = 2; o certo é 1)."""
    from src.backend import lrt_stats
    assert lrt_stats.PAIRS[('M0', 'M1a')]['df'] == 1
    assert lrt_stats.PAIRS[('M1a', 'M2a')]['df'] == 2
    assert lrt_stats.PAIRS[('M7', 'M8')]['df'] == 2
    assert lrt_stats.PAIRS[('M8a', 'M8')]['df'] == 1
    assert lrt_stats.PAIRS[('Branch-site_null', 'Branch-site')]['df'] == 1


def test_model_help_in_both_languages():
    """Rodada 2: a ajuda dos modelos ficava em inglês com a janela em PT."""
    en, pt = CodemlBatchAnalysis.MODEL_INFO, CodemlBatchAnalysis.MODEL_INFO_PT
    assert set(en) == set(pt)
    for code in en:
        assert set(en[code]) == set(pt[code]), code
        assert pt[code]['references'] == en[code]['references'], code
    assert 'q < 0,05' in pt['M2a']['interpretation'] and 'M8a' in pt['M8']['interpretation']
