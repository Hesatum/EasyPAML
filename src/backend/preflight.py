"""
Data check before codeml runs: stop codons, names missing from the tree, taxa
that pruning will remove, duplicate files of one gene, and lengths that are not
a multiple of 3, each with the sequence name and codon position. The window
shows the result in a dialog and the command line prints it.
"""

import difflib
import re
from dataclasses import dataclass, field
from io import StringIO
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from .alignment_io import (ALIGNMENT_SUFFIXES, AlignmentError, find_stop_codons,
                           read_alignment)

ERROR = 'error'      # the gene cannot run like this
WARNING = 'warning'  # it runs, but the user should know (may change the result)
INFO = 'info'        # information only


@dataclass
class Issue:
    gene: str
    kind: str
    severity: str
    data: Dict = field(default_factory=dict)

    def message(self, lang: str = 'en') -> str:
        return format_issue(self, lang)


@dataclass
class PreflightReport:
    genes: List[str]
    files: Dict[str, Path]
    issues: List[Issue]
    ignored_files: List[Path] = field(default_factory=list)

    @property
    def has_errors(self) -> bool:
        return any(i.severity == ERROR for i in self.issues)

    @property
    def has_problems(self) -> bool:
        return any(i.severity in (ERROR, WARNING) for i in self.issues)

    def by_gene(self) -> Dict[str, List[Issue]]:
        out: Dict[str, List[Issue]] = {}
        for i in self.issues:
            out.setdefault(i.gene, []).append(i)
        return out

    def format_text(self, lang: str = 'en', include_info: bool = True) -> str:
        lines = []
        for gene, issues in self.by_gene().items():
            shown = [i for i in issues if include_info or i.severity != INFO]
            if not shown:
                continue
            lines.append(gene if gene else ('(geral)' if lang == 'pt' else '(general)'))
            for i in shown:
                tag = {'error': 'ERRO' if lang == 'pt' else 'ERROR',
                       'warning': 'AVISO' if lang == 'pt' else 'WARNING',
                       'info': 'INFO'}[i.severity]
                lines.append(f"  [{tag}] {i.message(lang)}")
        return "\n".join(lines)


# ── Tree ──────────────────────────────────────────────────────────────────

_LABEL_RE = re.compile(r"\s*[#$]\d+$")


def read_tree_taxa(tree_file) -> List[str]:
    """Tip names of the tree, without #1/$1 labels or an 'N 1' header."""
    from Bio import Phylo
    raw = re.sub(r'\s+([#$]\d)', r'\1', Path(tree_file).read_text(encoding='utf-8', errors='replace'))
    lines = raw.splitlines()
    if lines and lines[0].strip() and lines[0].strip().split()[0].isdigit() \
            and not lines[0].strip().startswith('('):
        raw = "\n".join(lines[1:])
    tree = Phylo.read(StringIO(raw), 'newick')
    names = []
    for t in tree.get_terminals():
        if t.name:
            names.append(_LABEL_RE.sub('', t.name.strip().strip("'\"")))
    return names


def suggest_name(name: str, candidates: Sequence[str]) -> Optional[str]:
    """Nome mais parecido (ou None) -- ex.: 'Macaca_mulata' -> 'Macaca_mulatta'."""
    lower = {c.lower(): c for c in candidates}
    hit = difflib.get_close_matches(name.lower(), list(lower), n=1, cutoff=0.75)
    return lower[hit[0]] if hit else None


# ── Arquivos de alinhamento ─────────────────────────────────────────────────

def list_alignment_files(input_folder) -> List[Path]:
    folder = Path(input_folder)
    files = [p for p in folder.iterdir()
             if p.is_file() and p.suffix.lower() in ALIGNMENT_SUFFIXES]
    return sorted(files, key=lambda p: p.name.lower())


_FORMAT_PREFERENCE = {'.fasta': 0, '.fas': 1, '.fa': 2, '.fna': 3, '.phy': 4, '.phylip': 5}


def group_by_gene(files: Sequence[Path]):
    """{gene: chosen file}, [files ignored as duplicates].

    Two files with the same base name (gene.fasta and gene.phy) are one gene;
    only one is used, FASTA first.
    """
    groups: Dict[str, List[Path]] = {}
    for f in files:
        groups.setdefault(f.stem, []).append(f)
    chosen: Dict[str, Path] = {}
    ignored: List[Path] = []
    for gene, fs in groups.items():
        fs = sorted(fs, key=lambda p: _FORMAT_PREFERENCE.get(p.suffix.lower(), 9))
        chosen[gene] = fs[0]
        ignored.extend(fs[1:])
    return chosen, ignored


TREE_SUFFIXES = ('.nwk', '.tree', '.tre', '.newick', '.treefile', '.contree')
_RAXML_PREFIXES = ('raxml_besttree.', 'raxml_bipartitions.', 'raxml_result.')
_RAXML_SUFFIXES = ('.raxml.besttree', '.raxml.support', '.raxml.besttreecollapsed')
_ALN_SUFFIXES = ('.fasta', '.fas', '.fa', '.fna', '.phy', '.phylip', '.aln')


def is_tree_name(name: str) -> bool:
    low = name.lower()
    return low.endswith(TREE_SUFFIXES + _RAXML_SUFFIXES) or low.startswith(_RAXML_PREFIXES)


def tree_key(name: str):
    """(gene key, tool) of a tree file name: extensions and the prefixes and
    suffixes of IQ-TREE and RAxML removed, lower case. tool is 'IQ-TREE', 'RAxML'
    or ''."""
    low, tool = name.lower(), ''
    for pre in _RAXML_PREFIXES:
        if low.startswith(pre):
            low, tool = low[len(pre):], 'RAxML'
    changed = True
    while changed:
        changed = False
        for suf in _RAXML_SUFFIXES + TREE_SUFFIXES + _ALN_SUFFIXES:
            if low.endswith(suf) and len(low) > len(suf):
                if suf in _RAXML_SUFFIXES:
                    tool = 'RAxML'
                elif suf in ('.treefile', '.contree'):
                    tool = tool or 'IQ-TREE'
                low, changed = low[:-len(suf)], True
    return low, tool


@dataclass
class TreePairing:
    """Genes and per-gene tree files paired by name."""
    pairs: Dict[str, Path] = field(default_factory=dict)
    tools: Dict[str, str] = field(default_factory=dict)
    missing: List[str] = field(default_factory=list)
    suggestions: Dict[str, Path] = field(default_factory=dict)
    orphans: List[Path] = field(default_factory=list)
    duplicates: Dict[str, List[Path]] = field(default_factory=dict)


def pair_trees(genes: Sequence[str], alignment_folder=None, tree_folder=None) -> TreePairing:
    """Pair each gene with a tree file in the tree folder (or, without one, in the
    alignments folder). Names must match exactly after tree_key(); a near name is
    only offered as a suggestion."""
    out = TreePairing()
    by_key: Dict[str, List[Path]] = {}
    folders = [Path(tree_folder)] if tree_folder else []
    if alignment_folder:
        folders.append(Path(alignment_folder))
    seen = set()
    for folder in folders:
        if not folder.is_dir():
            continue
        for p in sorted(folder.iterdir()):
            if p.is_file() and is_tree_name(p.name) and p.resolve() not in seen:
                seen.add(p.resolve())
                by_key.setdefault(tree_key(p.name)[0], []).append(p)
    gene_keys = {g.lower(): g for g in genes}
    for key, files in by_key.items():
        gene = gene_keys.get(key)
        if gene is None:
            out.orphans.extend(files)
        elif len(files) > 1 and len({f.parent for f in files}) == 1:
            out.duplicates[gene] = files
        else:
            out.pairs[gene] = files[0]
            out.tools[gene] = tree_key(files[0].name)[1]
    orphan_keys = {tree_key(p.name)[0]: p for p in out.orphans}
    for gene in genes:
        if gene in out.pairs or gene in out.duplicates:
            continue
        out.missing.append(gene)
        hit = difflib.get_close_matches(gene.lower(), list(orphan_keys), n=1, cutoff=0.8)
        if hit:
            out.suggestions[gene] = orphan_keys[hit[0]]
    return out


def discover_per_gene_trees(alignment_folder, tree_folder=None, genes=None) -> Dict[str, Path]:
    """{gene: tree file} for genes with their own tree (see pair_trees). Without
    genes, those of the alignments folder."""
    if genes is None:
        genes = list(group_by_gene(list_alignment_files(alignment_folder))[0]) if alignment_folder else []
    return pair_trees(sorted(genes), alignment_folder, tree_folder).pairs


# ── Main check ──────────────────────────────────────────────────

def run_preflight(input_folder, tree_file, auto_prune: bool = True,
                  ignore_stop_codons: bool = False,
                  per_gene_trees: Optional[Dict[str, Path]] = None,
                  tree_folder=None) -> PreflightReport:
    """Check every alignment in the folder against its tree.

    per_gene_trees: {gene: tree file} for genes with their own tree; the
    others use tree_file. tree_folder: folder of per-gene trees (pair_trees).
    """
    files = list_alignment_files(input_folder)
    chosen, ignored = group_by_gene(files)
    issues: List[Issue] = []
    pairing = pair_trees(sorted(chosen), input_folder, tree_folder)
    if per_gene_trees is None:
        per_gene_trees = pairing.pairs
    for gene in chosen:
        if gene in per_gene_trees:
            continue
        twins = pairing.duplicates.get(gene)
        if twins:
            issues.append(Issue(gene, 'tree_duplicate', WARNING if tree_file else ERROR,
                                {'files': ", ".join(p.name for p in twins), 'general': bool(tree_file)}))
        elif tree_folder or not tree_file:
            hint = pairing.suggestions.get(gene)
            issues.append(Issue(gene, 'no_tree', WARNING if tree_file else ERROR,
                                {'suggestion': hint.name if hint else '', 'general': bool(tree_file)}))
    if tree_folder and pairing.orphans:
        issues.append(Issue('', 'tree_orphans', WARNING,
                            {'n': len(pairing.orphans), 'names': [p.name for p in pairing.orphans]}))
    if per_gene_trees:
        issues.append(Issue('', 'per_gene_trees', INFO, {'n': len(per_gene_trees)}))

    if not files:
        issues.append(Issue('', 'no_alignments', ERROR, {'folder': str(input_folder)}))

    for dup in ignored:
        issues.append(Issue(dup.stem, 'duplicate_gene', WARNING,
                            {'used': chosen[dup.stem].name, 'ignored': dup.name}))

    tree_cache: Dict[str, Optional[List[str]]] = {}

    def taxa_for(gene: str):
        tf = (per_gene_trees or {}).get(gene) or tree_file
        if not tf:
            return None
        key = str(tf)
        if key not in tree_cache:
            try:
                tree_cache[key] = read_tree_taxa(tf) if tf else None
            except Exception as exc:  # unreadable tree
                tree_cache[key] = None
                issues.append(Issue('', 'tree_unreadable', ERROR,
                                    {'tree': Path(tf).name, 'error': str(exc)}))
        return tree_cache[key]

    for gene, path in chosen.items():
        try:
            aln = read_alignment(path)
        except AlignmentError as exc:
            issues.append(Issue(gene, 'unreadable', ERROR, {'file': path.name, 'error': str(exc)}))
            continue

        for d in aln.duplicate_names:
            issues.append(Issue(gene, 'duplicate_name', ERROR, {'name': d}))

        if len(aln.names) < 3:
            issues.append(Issue(gene, 'too_few_sequences', ERROR, {'n': len(aln.names)}))
            continue

        if not aln.is_aligned:
            lens = sorted(set(aln.lengths))
            issues.append(Issue(gene, 'unaligned', ERROR, {'lengths': lens}))
            continue

        if aln.length % 3 != 0:
            issues.append(Issue(gene, 'not_multiple_of_3', ERROR,
                                {'length': aln.length, 'remainder': aln.length % 3}))

        for name, pos, codon in find_stop_codons(aln.names, aln.seqs):
            last = pos == aln.length // 3
            issues.append(Issue(gene, 'stop_codon',
                                INFO if ignore_stop_codons else WARNING,
                                {'sequence': name, 'codon_position': pos,
                                 'nucleotide_position': (pos - 1) * 3 + 1,
                                 'codon': codon, 'terminal': last}))

        tree_taxa = taxa_for(gene)
        if tree_taxa is None:
            continue
        aln_set, tree_set = set(aln.names), set(tree_taxa)
        for name in aln.names:
            if name not in tree_set:
                issues.append(Issue(
                    gene, 'name_not_in_tree', WARNING if auto_prune else ERROR,
                    {'name': name, 'suggestion': suggest_name(name, sorted(tree_set - aln_set))
                     or suggest_name(name, tree_taxa), 'auto_prune': auto_prune}))
        absent = [t for t in tree_taxa if t not in aln_set]
        if absent:
            issues.append(Issue(gene, 'tree_taxa_pruned', INFO if auto_prune else ERROR,
                                {'names': absent, 'auto_prune': auto_prune}))
        if len(aln_set & tree_set) < 3:
            issues.append(Issue(gene, 'too_few_shared_taxa', ERROR,
                                {'n': len(aln_set & tree_set)}))

    return PreflightReport(genes=list(chosen), files=chosen, issues=issues,
                           ignored_files=ignored)


# ── Mensagens (PT / EN) ─────────────────────────────────────────────────────

_MESSAGES = {
    'no_alignments': {
        'pt': "Nenhum alinhamento (.fasta, .fas, .phy, .phylip) encontrado em {folder}.",
        'en': "No alignment (.fasta, .fas, .phy, .phylip) found in {folder}.",
    },
    'duplicate_gene': {
        'pt': "{used} e {ignored} são o mesmo gene; será usado só {used}.",
        'en': "{used} and {ignored} are the same gene; only {used} will be used.",
    },
    'tree_unreadable': {
        'pt': "Não foi possível ler a árvore {tree}: {error}",
        'en': "Could not read the tree {tree}: {error}",
    },
    'unreadable': {
        'pt': "Não foi possível ler {file}: {error}",
        'en': "Could not read {file}: {error}",
    },
    'duplicate_name': {
        'pt': "O nome '{name}' aparece mais de uma vez no alinhamento.",
        'en': "The name '{name}' appears more than once in the alignment.",
    },
    'too_few_sequences': {
        'pt': "Só {n} sequência(s); são necessárias pelo menos 3.",
        'en': "Only {n} sequence(s); at least 3 are needed.",
    },
    'unaligned': {
        'pt': "As sequências têm comprimentos diferentes ({lengths}); o arquivo não está alinhado.",
        'en': "Sequences have different lengths ({lengths}); the file is not aligned.",
    },
    'not_multiple_of_3': {
        'pt': "O alinhamento tem {length} nucleotídeos, que não é múltiplo de 3 (sobram {remainder}).",
        'en': "The alignment has {length} nucleotides, which is not a multiple of 3 ({remainder} left over).",
    },
    'stop_codon': {
        'pt': "Stop codon {codon} na sequência {sequence}, códon {codon_position} "
              "(nucleotídeo {nucleotide_position}){terminal_pt}. O codeml trata essa coluna inteira como dado ausente. "
              "Se o códon está errado (erro de sequenciamento, montagem ou alinhamento), corrija o "
              "arquivo, por exemplo trocando esse códon por NNN nessa sequência. Se ele é real "
              "(por exemplo, um pseudogene), \"Continuar mesmo assim\" roda o gene com o códon "
              "mascarado, como \"Ignorar stop codons\"; o aviso fica em genes_status.tsv, no painel "
              "e em methods_text.txt.",
        'en': "Stop codon {codon} in sequence {sequence}, codon {codon_position} "
              "(nucleotide {nucleotide_position}){terminal_en}. codeml treats that whole column as missing data. "
              "If the codon is wrong (sequencing, assembly or alignment error), fix the file, for "
              "example by replacing that codon with NNN in that sequence. If it is real (e.g. a "
              "pseudogene), \"Continue anyway\" runs the gene with the codon masked, the same as "
              "\"Ignore stop codons\"; the warning is kept in genes_status.tsv, the panel and "
              "methods_text.txt.",
    },
    'name_not_in_tree': {
        'pt': "'{name}' está no alinhamento mas não na árvore{suggest_pt}.{prune_pt}",
        'en': "'{name}' is in the alignment but not in the tree{suggest_en}.{prune_en}",
    },
    'tree_taxa_pruned': {
        'pt': "{count} táxon(s) da árvore não estão neste gene{prune_pt}: {list}",
        'en': "{count} tree taxon/taxa are not in this gene{prune_en}: {list}",
    },
    'no_tree': {
        'pt': "Nenhuma árvore com o nome deste gene{suggest_pt}.{general_pt}",
        'en': "No tree named after this gene{suggest_en}.{general_en}",
    },
    'tree_duplicate': {
        'pt': "Mais de uma árvore com o nome deste gene ({files}); deixe só uma.{general_pt}",
        'en': "More than one tree named after this gene ({files}); keep only one.{general_en}",
    },
    'tree_orphans': {
        'pt': "{n} árvore(s) sem alinhamento com o mesmo nome: {list}",
        'en': "{n} tree(s) with no alignment of the same name: {list}",
    },
    'per_gene_trees': {
        'pt': "{n} gene(s) com árvore própria (GENE.nwk na pasta); para eles ela substitui o arquivo de árvore.",
        'en': "{n} gene(s) with their own tree (GENE.nwk in the folder); for them it replaces the tree file.",
    },
    'too_few_shared_taxa': {
        'pt': "Só {n} sequência(s) em comum entre alinhamento e árvore; são necessárias pelo menos 3.",
        'en': "Only {n} sequence(s) shared by alignment and tree; at least 3 are needed.",
    },
}


def format_issue(issue: Issue, lang: str = 'en') -> str:
    lang = 'pt' if lang == 'pt' else 'en'
    d = dict(issue.data)
    if issue.kind == 'stop_codon':
        d['terminal_pt'] = ", no fim da sequência" if d.get('terminal') else ""
        d['terminal_en'] = ", at the end of the sequence" if d.get('terminal') else ""
    if issue.kind == 'name_not_in_tree':
        s = d.get('suggestion')
        d['suggest_pt'] = f" (você quis dizer '{s}'?)" if s else ""
        d['suggest_en'] = f" (did you mean '{s}'?)" if s else ""
        d['prune_pt'] = (" Ela será EXCLUÍDA da análise." if d.get('auto_prune')
                         else " Sem a poda automática, o codeml vai falhar.")
        d['prune_en'] = (" It will be EXCLUDED from the analysis." if d.get('auto_prune')
                         else " Without automatic pruning, codeml will fail.")
    if issue.kind == 'tree_taxa_pruned':
        names = d.get('names', [])
        d['count'] = len(names)
        d['list'] = ", ".join(names[:8]) + (" …" if len(names) > 8 else "")
        d['prune_pt'] = " e serão podados da árvore" if d.get('auto_prune') else ""
        d['prune_en'] = " and will be pruned from the tree" if d.get('auto_prune') else ""
    if issue.kind in ('no_tree', 'tree_duplicate'):
        s = d.get('suggestion')
        d['suggest_pt'] = f" (você quis dizer '{s}'? renomeie o arquivo)" if s else ""
        d['suggest_en'] = f" (did you mean '{s}'? rename the file)" if s else ""
        d['general_pt'] = (" A árvore geral será usada." if d.get('general')
                           else " Escolha um arquivo de árvore geral ou renomeie a árvore como GENE.nwk.")
        d['general_en'] = (" The general tree will be used." if d.get('general')
                           else " Choose a general tree file or name the tree GENE.nwk.")
    if issue.kind == 'tree_orphans':
        names = d.get('names', [])
        d['list'] = ", ".join(names[:8]) + (" …" if len(names) > 8 else "")
    if issue.kind == 'unaligned':
        d['lengths'] = ", ".join(str(x) for x in d.get('lengths', []))
    template = _MESSAGES.get(issue.kind, {}).get(lang)
    if not template:
        return f"{issue.kind}: {issue.data}"
    try:
        return template.format(**d)
    except (KeyError, IndexError):
        return f"{issue.kind}: {issue.data}"
