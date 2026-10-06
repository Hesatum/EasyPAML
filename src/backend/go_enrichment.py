"""GO enrichment: candidates (significant LRT) against all tested genes, with
Fisher's exact test per term. The annotation TSV has one column per GO
category with entries like "description [GO:XXXXXXX]; ...".
"""
import re
from pathlib import Path

import pandas as pd
from scipy import stats

GO_TERM_RE = re.compile(r'([^;]+?)\s*\[(GO:\d+)\]')
GO_COLUMNS = ('go_biological_process', 'go_cellular_component', 'go_molecular_function')


def load_gene_to_go(annotation_file: Path, gene_id_col: str = 'gene_id_full') -> dict:
    """gene_id -> {go_id: description}, merging the three GO categories."""
    df = pd.read_csv(annotation_file, sep='\t', dtype=str)
    gene_to_go = {}
    for _, row in df.iterrows():
        terms = {}
        for col in GO_COLUMNS:
            cell = row.get(col)
            if isinstance(cell, str):
                for m in GO_TERM_RE.finditer(cell):
                    terms[m.group(2)] = m.group(1).strip()
        gene_to_go[row[gene_id_col]] = terms
    return gene_to_go


_ENRICH_COLUMNS = ['go_id', 'description', 'n_candidates', 'n_background', 'odds_ratio', 'p_value', 'q_value']


def enrich(candidate_genes: set, background_genes: set, gene_to_go: dict,
           min_candidates: int = 2) -> pd.DataFrame:
    """Fisher's exact test (term / not term x candidate / background) for each GO
    term present in at least `min_candidates` candidates. Returns a DataFrame
    sorted by p-value; q_value is Benjamini-Hochberg across the tested terms."""
    n_cand = len(candidate_genes)
    n_bg = len(background_genes)
    if n_cand == 0 or n_bg == 0:
        return pd.DataFrame(columns=_ENRICH_COLUMNS)

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
        a = len(cand_set)                      # candidates with the term
        b = n_cand - a                         # candidates without the term
        c = len(bg_set) - a                    # background with the term
        d = (n_bg - n_cand) - c                # background without the term
        odds_ratio, p_value = stats.fisher_exact([[a, b], [max(c, 0), max(d, 0)]], alternative='greater')
        rows.append({
            'go_id': go_id, 'description': term_desc[go_id],
            'n_candidates': a, 'n_background': len(bg_set),
            'odds_ratio': odds_ratio, 'p_value': p_value,
        })
    if not rows:
        return pd.DataFrame(columns=_ENRICH_COLUMNS)
    table = pd.DataFrame(rows)
    table['q_value'] = stats.false_discovery_control(table['p_value'], method='bh')
    return table.sort_values('p_value')


# Tests used to pick candidates. M7 vs M8 counts only when M8a vs M8 was not
# run: M8 can beat M7 through neutral sites alone (see METHODS.md).
CANDIDATE_TESTS = (('M1a', 'M2a'), ('M8a', 'M8'), ('M7', 'M8'))


def rank_candidates(lrt_summary_tsv: Path, annotation_file: Path,
                    sig_threshold: float = 0.05) -> tuple:
    """Read analysis_summary.tsv and the GO annotation TSV; return (candidates,
    GO enrichment table).

    A candidate has q < sig_threshold in M1a vs M2a or M8a vs M8 (or M7 vs M8
    when M8a vs M8 was not run), with the same p and q as LRT_results.txt.
    Failed genes are left out. Columns: Gene, test, p_value, q_value, go_terms."""
    summary = pd.read_csv(lrt_summary_tsv, sep='\t')
    gene_to_go = load_gene_to_go(annotation_file)
    if 'status' in summary.columns:
        summary = summary[summary['status'].astype(str) != 'failed']

    tests = [t for t in CANDIDATE_TESTS if f"q_{t[0]}_vs_{t[1]}" in summary.columns]
    if ('M8a', 'M8') in tests:
        tests = [t for t in tests if t != ('M7', 'M8')]

    best = []
    for _, row in summary.iterrows():
        pick = None
        for null, alt in tests:
            q = pd.to_numeric(row.get(f"q_{null}_vs_{alt}"), errors='coerce')
            pv = pd.to_numeric(row.get(f"p_{null}_vs_{alt}"), errors='coerce')
            if pd.notna(q) and (pick is None or q < pick[2]):
                pick = (f"{alt} vs {null}", pv, q)
        best.append(pick or (None, float('nan'), float('nan')))
    summary = summary.assign(test=[b[0] for b in best], p_value=[b[1] for b in best],
                             q_value=[b[2] for b in best])
    summary['go_terms'] = summary['Gene'].map(
        lambda g: '; '.join(sorted(gene_to_go.get(g, {}).values())) or 'no annotation'
    )

    all_genes = set(summary['Gene'])
    sig_genes = set(summary.loc[summary['q_value'] < sig_threshold, 'Gene'])

    candidates = summary.loc[summary['Gene'].isin(sig_genes)].sort_values('q_value')
    go_table = enrich(sig_genes, all_genes, gene_to_go)
    return candidates, go_table
