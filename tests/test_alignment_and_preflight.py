"""Alignment reading, site map and the data check."""
import json
from pathlib import Path

import pandas as pd
import pytest

from src.backend.alignment_io import (cleandata_kept_codons, find_stop_codons, parse_phylip,
                                      read_alignment)
from src.backend.preflight import run_preflight, suggest_name
from src.backend.site_map import attach_original_positions, to_original, write_sitemap

DATA = Path(__file__).resolve().parent / 'data'


# ── PHYLIP relaxado ──────────────────────────────────────────────────────────

def test_relaxed_phylip_long_names(tmp_path):
    p = tmp_path / 'g.phy'
    p.write_text(" 3 9\nHomo_sapiens  ATGAAACCC\nPan_troglodytes ATGAAGCCC\nGorilla_gorilla  ATGAAACCG\n")
    aln = read_alignment(p)
    assert aln.names == ['Homo_sapiens', 'Pan_troglodytes', 'Gorilla_gorilla']
    assert aln.is_aligned and aln.length == 9
    assert aln.phylip_variant == 'relaxed-sequential'


def test_strict_phylip_ten_columns():
    text = " 2 6\nHomo_sapieATGAAA\nPan_trogloATGAAG\n"
    names, seqs, variant = parse_phylip(text)
    assert names == ['Homo_sapie', 'Pan_troglo'] and variant == 'strict-sequential'


def test_relaxed_interleaved_phylip():
    text = " 2 12\nseqA ATGAAA\nseqB ATGAAG\nCCCGGG\nCCCGGA\n"
    names, seqs, variant = parse_phylip(text)
    assert seqs['seqA'] == 'ATGAAACCCGGG' and seqs['seqB'] == 'ATGAAGCCCGGA'
    assert variant == 'relaxed-interleaved'


def test_phylip_with_sequence_split_over_lines():
    text = " 2 12\nseqA ATGAAA\nCCCGGG\nseqB ATGAAG\nCCCGGA\n"
    names, seqs, variant = parse_phylip(text)
    assert seqs['seqB'] == 'ATGAAGCCCGGA'


# ── Stops e cleandata ────────────────────────────────────────────────────────

def test_stop_codon_position_in_problematic_data():
    aln = read_alignment(DATA / 'gene_problematic.fasta')
    stops = find_stop_codons(aln.names, aln.seqs)
    assert stops == [('Gorilla_gorilla', 150, 'TGA')]


def test_cleandata_kept_codons_rule():
    names = ['a', 'b', 'c']
    seqs = {'a': 'ATGAAACCCGGGTTT', 'b': 'ATG---CCCTGATTT', 'c': 'ATGAAANCCGGGTTT'}
    # codon 2: gap in b; codon 3: N in c; codon 4: stop in b -> removed
    assert cleandata_kept_codons(names, seqs) == [1, 5]


def test_sitemap_maps_codeml_positions_back(tmp_path):
    aln = read_alignment(DATA / 'gene_problematic.fasta')
    kept = cleandata_kept_codons(aln.names, aln.seqs)
    assert len(kept) == 299 and 150 not in kept
    results = tmp_path / 'g_M8_results.txt'
    results.write_text("ns =  10  ls = 299\n")
    write_sitemap(tmp_path / 'g_M8_sitemap.json', cleandata=1, n_codons=300,
                  kept_codons=kept, sequences=aln.names, codeml_sites=299)
    sm = json.loads((tmp_path / 'g_M8_sitemap.json').read_text())
    assert sm['verified'] is True
    # codeml site 220 is alignment codon 221 (column 150 was removed)
    assert to_original(220, sm) == 221
    assert to_original(65, sm) == 65
    df = attach_original_positions(pd.DataFrame({'position': [65, 220, 280]}), results)
    assert list(df['position_original']) == [65, 221, 281]
    assert df['position_mapped'].all()


def test_sitemap_missing_means_unmapped(tmp_path):
    df = attach_original_positions(pd.DataFrame({'position': [5]}), tmp_path / 'x_M8_results.txt')
    assert list(df['position_original']) == [5] and not df['position_mapped'].any()


# ── Preflight ────────────────────────────────────────────────────────────────

def _folder_with(tmp_path, files):
    d = tmp_path / 'in'
    d.mkdir()
    for name, src in files.items():
        (d / name).write_text(Path(src).read_text() if isinstance(src, Path) else src)
    return d


def test_preflight_problematic_data(tmp_path):
    d = _folder_with(tmp_path, {'meu_gene.fasta': DATA / 'gene_problematic.fasta'})
    rep = run_preflight(d, DATA / 'gene_example.nwk')
    kinds = {i.kind: i for i in rep.issues}
    stop = kinds['stop_codon']
    assert stop.data['sequence'] == 'Gorilla_gorilla'
    assert stop.data['codon_position'] == 150 and stop.data['nucleotide_position'] == 448
    name = kinds['name_not_in_tree']
    assert name.data['name'] == 'Macaca_mulata'
    assert name.data['suggestion'] == 'Macaca_mulatta'
    assert 'Macaca_mulatta' in kinds['tree_taxa_pruned'].data['names']
    assert rep.has_problems
    pt = rep.format_text('pt')
    assert 'códon 150' in pt and "você quis dizer 'Macaca_mulatta'" in pt and 'EXCLUÍDA' in pt
    en = rep.format_text('en')
    assert 'codon 150' in en and "did you mean 'Macaca_mulatta'" in en


def test_preflight_duplicates_and_length(tmp_path):
    good = (DATA / 'gene_example.fasta').read_text()
    phy_names = [l[1:].strip() for l in good.splitlines() if l.startswith('>')]
    d = _folder_with(tmp_path, {
        'g1.fasta': DATA / 'gene_example.fasta',
        'g1.phy': " 3 6\n" + "\n".join(f"{n}  ATGAAA" for n in phy_names[:3]) + "\n",
        'g2.fasta': ">Homo_sapiens\nATGAA\n>Pan_troglodytes\nATGAG\n>Gorilla_gorilla\nATGAC\n",
    })
    rep = run_preflight(d, DATA / 'gene_example.nwk')
    by = {(i.gene, i.kind) for i in rep.issues}
    assert ('g1', 'duplicate_gene') in by
    assert ('g2', 'not_multiple_of_3') in by
    assert rep.files['g1'].suffix == '.fasta'          # FASTA wins


def test_preflight_clean_data_has_no_problems(tmp_path):
    d = _folder_with(tmp_path, {'g.fasta': DATA / 'gene_example.fasta'})
    rep = run_preflight(d, DATA / 'gene_example.nwk')
    assert not rep.has_problems


def test_preflight_ignore_stop_codons_downgrades_to_info(tmp_path):
    d = _folder_with(tmp_path, {'g.fasta': DATA / 'gene_problematic.fasta'})
    rep = run_preflight(d, DATA / 'gene_example.nwk', ignore_stop_codons=True)
    assert [i.severity for i in rep.issues if i.kind == 'stop_codon'] == ['info']


def test_suggest_name():
    assert suggest_name('Macaca_mulata', ['Macaca_mulatta', 'Papio_anubis']) == 'Macaca_mulatta'
    assert suggest_name('Zebra', ['Macaca_mulatta']) is None


def test_tree_names_of_iqtree_and_raxml_pair_with_their_gene(tmp_path):
    from src.backend.preflight import pair_trees
    trees = tmp_path / 'trees'
    trees.mkdir()
    for name in ('PEPC1.fasta.treefile', 'RAxML_bestTree.PPCK1', 'almt2.nwk', 'CHS.raxml.bestTree',
                 'AQP_PIP14.nwk', 'notes.txt'):
        (trees / name).write_text('(a,b,c);')
    pr = pair_trees(['PEPC1', 'PPCK1', 'ALMT2', 'CHS', 'AQP_PIP1-4'], tree_folder=trees)
    assert {g: p.name for g, p in pr.pairs.items()} == {
        'PEPC1': 'PEPC1.fasta.treefile', 'PPCK1': 'RAxML_bestTree.PPCK1', 'ALMT2': 'almt2.nwk',
        'CHS': 'CHS.raxml.bestTree'}
    assert pr.tools['PEPC1'] == 'IQ-TREE' and pr.tools['PPCK1'] == 'RAxML' and pr.tools['ALMT2'] == ''
    assert pr.missing == ['AQP_PIP1-4']                         # a near name is not paired
    assert pr.suggestions['AQP_PIP1-4'].name == 'AQP_PIP14.nwk'
    assert [p.name for p in pr.orphans] == ['AQP_PIP14.nwk']


def test_two_trees_for_one_gene_are_reported_not_guessed(tmp_path):
    from src.backend.preflight import pair_trees
    (tmp_path / 'g1.nwk').write_text('(a,b,c);')
    (tmp_path / 'g1.fasta.treefile').write_text('(a,b,c);')
    pr = pair_trees(['g1'], tree_folder=tmp_path)
    assert 'g1' not in pr.pairs and len(pr.duplicates['g1']) == 2


def test_data_check_names_the_tree_to_rename(tmp_path):
    data = Path(__file__).resolve().parent / 'data'
    inp, trees = tmp_path / 'in', tmp_path / 'trees'
    inp.mkdir()
    trees.mkdir()
    (inp / 'AQP_PIP1-4.fasta').write_text((data / 'gene_example.fasta').read_text())
    (trees / 'AQP_PIP14.nwk').write_text((data / 'gene_example.nwk').read_text())
    report = run_preflight(inp, None, tree_folder=trees)
    msgs = {i.kind: i.message('en') for i in report.issues}
    assert "did you mean 'AQP_PIP14.nwk'? rename the file" in msgs['no_tree']
    assert 'AQP_PIP14.nwk' in msgs['tree_orphans']
    assert any(i.kind == 'no_tree' and i.severity == 'error' for i in report.issues)
