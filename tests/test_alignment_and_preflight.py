"""Itens 1c e 3 -- leitura de alinhamentos, mapa de sítios, validação."""
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
    aln = read_alignment(DATA / 'gene_problematico.fasta')
    stops = find_stop_codons(aln.names, aln.seqs)
    assert stops == [('Gorilla_gorilla', 150, 'TGA')]


def test_cleandata_kept_codons_rule():
    names = ['a', 'b', 'c']
    seqs = {'a': 'ATGAAACCCGGGTTT', 'b': 'ATG---CCCTGATTT', 'c': 'ATGAAANCCGGGTTT'}
    # códon 2: gap em b; códon 3: N em c; códon 4: stop em b -> saem
    assert cleandata_kept_codons(names, seqs) == [1, 5]


def test_sitemap_maps_codeml_positions_back(tmp_path):
    aln = read_alignment(DATA / 'gene_problematico.fasta')
    kept = cleandata_kept_codons(aln.names, aln.seqs)
    assert len(kept) == 299 and 150 not in kept
    results = tmp_path / 'g_M8_results.txt'
    results.write_text("ns =  10  ls = 299\n")
    write_sitemap(tmp_path / 'g_M8_sitemap.json', cleandata=1, n_codons=300,
                  kept_codons=kept, sequences=aln.names, codeml_sites=299)
    sm = json.loads((tmp_path / 'g_M8_sitemap.json').read_text())
    assert sm['verified'] is True
    # codeml diz 220 -> no alinhamento do usuário é 221 (coluna 150 removida antes)
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
    d = _folder_with(tmp_path, {'meu_gene.fasta': DATA / 'gene_problematico.fasta'})
    rep = run_preflight(d, DATA / 'gene_exemplo.nwk')
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
    good = (DATA / 'gene_exemplo.fasta').read_text()
    phy_names = [l[1:].strip() for l in good.splitlines() if l.startswith('>')]
    d = _folder_with(tmp_path, {
        'g1.fasta': DATA / 'gene_exemplo.fasta',
        'g1.phy': " 3 6\n" + "\n".join(f"{n}  ATGAAA" for n in phy_names[:3]) + "\n",
        'g2.fasta': ">Homo_sapiens\nATGAA\n>Pan_troglodytes\nATGAG\n>Gorilla_gorilla\nATGAC\n",
    })
    rep = run_preflight(d, DATA / 'gene_exemplo.nwk')
    by = {(i.gene, i.kind) for i in rep.issues}
    assert ('g1', 'duplicate_gene') in by
    assert ('g2', 'not_multiple_of_3') in by
    assert rep.files['g1'].suffix == '.fasta'          # FASTA tem preferência


def test_preflight_clean_data_has_no_problems(tmp_path):
    d = _folder_with(tmp_path, {'g.fasta': DATA / 'gene_exemplo.fasta'})
    rep = run_preflight(d, DATA / 'gene_exemplo.nwk')
    assert not rep.has_problems


def test_preflight_ignore_stop_codons_downgrades_to_info(tmp_path):
    d = _folder_with(tmp_path, {'g.fasta': DATA / 'gene_problematico.fasta'})
    rep = run_preflight(d, DATA / 'gene_exemplo.nwk', ignore_stop_codons=True)
    assert [i.severity for i in rep.issues if i.kind == 'stop_codon'] == ['info']


def test_suggest_name():
    assert suggest_name('Macaca_mulata', ['Macaca_mulatta', 'Papio_anubis']) == 'Macaca_mulatta'
    assert suggest_name('Zebra', ['Macaca_mulatta']) is None
