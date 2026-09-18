"""Enriquecimento de GO: candidatos (LRT significativo) vs background (todos
os genes testados). Fisher exato por termo -- scipy ja e dependencia, nada
novo. Le direto o TSV de anotacao no formato "descricao [GO:XXXXXXX]; ..."
por coluna de categoria (biological_process/cellular_component/molecular_function).
"""
import re
from pathlib import Path
from typing import Optional

import pandas as pd
from scipy import stats

GO_TERM_RE = re.compile(r'([^;]+?)\s*\[(GO:\d+)\]')
GO_COLUMNS = {
    'biological_process': 'go_biological_process',
    'cellular_component': 'go_cellular_component',
    'molecular_function': 'go_molecular_function',
}


def _parse_go_terms(cell) -> set:
    if not isinstance(cell, str) or not cell.strip():
        return set()
    return {m.group(2) for m in GO_TERM_RE.finditer(cell)}


def load_gene_to_go(annotation_file: Path, gene_id_col: str = 'gene_id_full') -> dict:
    """gene_id -> {go_id: description}. Junta as 3 categorias GO numa coisa so
    (biological_process/cellular_component/molecular_function) -- pra
    enriquecimento por termo, a categoria so importa como rotulo de exibicao."""
    df = pd.read_csv(annotation_file, sep='\t', dtype=str)
    gene_to_go = {}
    for _, row in df.iterrows():
        terms = {}
        for col in GO_COLUMNS.values():
            cell = row.get(col)
            if isinstance(cell, str):
                for m in GO_TERM_RE.finditer(cell):
                    terms[m.group(2)] = m.group(1).strip()
        gene_to_go[row[gene_id_col]] = terms
    return gene_to_go


def enrich(candidate_genes: set, background_genes: set, gene_to_go: dict,
           min_candidates: int = 2) -> pd.DataFrame:
    """Fisher exato (2x2: no termo / fora do termo  x  candidato / background)
    por termo GO. Retorna DataFrame ordenado por p-valor, uma linha por termo
    que aparece em pelo menos `min_candidates` genes candidatos (termos
    presentes numa unica amostra nao sao informativos e so inflam o teste)."""
    n_cand = len(candidate_genes)
    n_bg = len(background_genes)
    if n_cand == 0 or n_bg == 0:
        return pd.DataFrame(columns=['go_id', 'description', 'n_candidates', 'n_background', 'odds_ratio', 'p_value'])

    term_candidates: dict = {}
    term_desc: dict = {}
    for gene in candidate_genes & background_genes:
        for go_id, desc in gene_to_go.get(gene, {}).items():
            term_candidates.setdefault(go_id, set()).add(gene)
            term_desc[go_id] = desc

    term_background: dict = {}
    for gene in background_genes:
        for go_id in gene_to_go.get(gene, {}):
            term_background.setdefault(go_id, set()).add(gene)

    rows = []
    for go_id, cand_set in term_candidates.items():
        if len(cand_set) < min_candidates:
            continue
        bg_set = term_background.get(go_id, set())
        a = len(cand_set)                      # candidatos com o termo
        b = n_cand - a                         # candidatos sem o termo
        c = len(bg_set) - a                    # background (nao-candidato) com o termo
        d = (n_bg - n_cand) - c                # background (nao-candidato) sem o termo
        odds_ratio, p_value = stats.fisher_exact([[a, b], [max(c, 0), max(d, 0)]], alternative='greater')
        rows.append({
            'go_id': go_id, 'description': term_desc[go_id],
            'n_candidates': a, 'n_background': len(bg_set),
            'odds_ratio': odds_ratio, 'p_value': p_value,
        })
    return pd.DataFrame(rows).sort_values('p_value') if rows else pd.DataFrame(
        columns=['go_id', 'description', 'n_candidates', 'n_background', 'odds_ratio', 'p_value'])


def rank_candidates(lrt_summary_tsv: Path, annotation_file: Path,
                     lrt_columns=('lrt_M1a_vs_M2a', 'lrt_M7_vs_M8'), df: int = 2,
                     sig_threshold: float = 0.05) -> tuple:
    """Le o analysis_summary.tsv do EasyPAML + o TSV de anotacao GO, devolve
    (tabela de candidatos ranqueada por efeito, tabela de enriquecimento GO).
    Ponto de entrada unico que a aba da GUI chama -- toda a logica de
    verdade mora nas duas funcoes acima, testaveis sem Tkinter."""
    summary = pd.read_csv(lrt_summary_tsv, sep='\t')
    gene_to_go = load_gene_to_go(annotation_file)

    best_stat = pd.Series(0.0, index=summary.index)
    for col in lrt_columns:
        if col in summary.columns:
            best_stat = best_stat.combine(summary[col].fillna(0).clip(lower=0), max)
    summary['p_value'] = stats.chi2.sf(best_stat, df=df)
    summary['go_terms'] = summary['Gene'].map(
        lambda g: '; '.join(sorted(gene_to_go.get(g, {}).values())) or 'sem anotacao'
    )

    all_genes = set(summary['Gene'])
    sig_genes = set(summary.loc[summary['p_value'] < sig_threshold, 'Gene'])

    candidates = summary.loc[summary['Gene'].isin(sig_genes)].sort_values('p_value')
    go_table = enrich(sig_genes, all_genes, gene_to_go)
    return candidates, go_table


def _self_check():
    """python3 -m src.backend.go_enrichment -- roda com dado real do projeto,
    falha alto se a logica quebrar."""
    ann = Path('/home/user/Desktop/projetos/matheus_wgs/resultados_analises/13_te_contamination/go_annotation/final_5725_annotation.tsv')
    if not ann.exists():
        print('SKIP: arquivo de anotacao real nao encontrado neste ambiente.')
        return
    gene_to_go = load_gene_to_go(ann)
    assert len(gene_to_go) > 5000, f'esperava ~5725 genes, achou {len(gene_to_go)}'
    assert gene_to_go['C1_41851674_g1001']['GO:0015995'] == 'chlorophyll biosynthetic process'

    import random
    random.seed(0)
    all_genes = set(gene_to_go)
    candidates = set(random.sample(sorted(all_genes), 50))
    go_table = enrich(candidates, all_genes, gene_to_go)
    assert not go_table.empty, 'esperava pelo menos um termo GO com >=2 candidatos numa amostra de 50 genes reais'
    assert (go_table['p_value'] >= 0).all() and (go_table['p_value'] <= 1).all()
    print(f'OK: {len(gene_to_go)} genes anotados, {len(go_table)} termos GO testados na amostra de checagem.')


if __name__ == '__main__':
    _self_check()
