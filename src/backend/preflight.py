"""
Validação ANTES de rodar o codeml.

Tudo que o teste de usabilidade mostrou que o programa só descobria durante
(ou depois de) a execução -- stop codons, nomes que não batem com a árvore,
táxons que a poda automática vai tirar, arquivos duplicados do mesmo gene,
comprimento que não é múltiplo de 3 -- é verificado aqui, com nome da
sequência e posição do códon, para a interface mostrar num diálogo e o CLI
imprimir antes de começar.

Não depende de Tk; testável isoladamente.
"""

import difflib
import re
from dataclasses import dataclass, field
from io import StringIO
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from .alignment_io import (ALIGNMENT_SUFFIXES, AlignmentError, find_stop_codons,
                           read_alignment)

ERROR = 'error'      # o gene não pode rodar assim
WARNING = 'warning'  # roda, mas o usuário precisa saber (pode mudar o resultado)
INFO = 'info'        # só informativo


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


# ── Árvore ──────────────────────────────────────────────────────────────────

_LABEL_RE = re.compile(r"[#$]\d+$")


def read_tree_taxa(tree_file) -> List[str]:
    """Nomes das folhas da árvore (sem marcas #1/$1 e sem cabeçalho 'N 1')."""
    from Bio import Phylo
    raw = Path(tree_file).read_text(encoding='utf-8', errors='replace')
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
    """{gene: arquivo escolhido}, [arquivos ignorados por duplicidade].

    Dois arquivos com o mesmo nome-base (gene.fasta + gene.phy) são o MESMO
    gene: usa-se um só (FASTA tem preferência) em vez de rodar duas vezes e
    sobrescrever resultados.
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


# ── Verificação principal ──────────────────────────────────────────────────

def run_preflight(input_folder, tree_file, auto_prune: bool = True,
                  ignore_stop_codons: bool = False,
                  per_gene_trees: Optional[Dict[str, Path]] = None) -> PreflightReport:
    """Verifica todos os alinhamentos da pasta contra a árvore.

    per_gene_trees: {gene: arquivo de árvore} para os genes que têm árvore
    própria (o resto usa tree_file).
    """
    files = list_alignment_files(input_folder)
    chosen, ignored = group_by_gene(files)
    issues: List[Issue] = []

    if not files:
        issues.append(Issue('', 'no_alignments', ERROR, {'folder': str(input_folder)}))

    for dup in ignored:
        issues.append(Issue(dup.stem, 'duplicate_gene', WARNING,
                            {'used': chosen[dup.stem].name, 'ignored': dup.name}))

    tree_cache: Dict[str, Optional[List[str]]] = {}

    def taxa_for(gene: str):
        tf = (per_gene_trees or {}).get(gene) or tree_file
        key = str(tf)
        if key not in tree_cache:
            try:
                tree_cache[key] = read_tree_taxa(tf) if tf else None
            except Exception as exc:  # árvore ilegível
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
              "(nucleotídeo {nucleotide_position}){terminal_pt}. O codeml trata essa coluna inteira como dado ausente.",
        'en': "Stop codon {codon} in sequence {sequence}, codon {codon_position} "
              "(nucleotide {nucleotide_position}){terminal_en}. codeml treats that whole column as missing data.",
    },
    'name_not_in_tree': {
        'pt': "'{name}' está no alinhamento mas não na árvore{suggest_pt}.{prune_pt}",
        'en': "'{name}' is in the alignment but not in the tree{suggest_en}.{prune_en}",
    },
    'tree_taxa_pruned': {
        'pt': "{count} táxon(s) da árvore não estão neste gene{prune_pt}: {list}",
        'en': "{count} tree taxon/taxa are not in this gene{prune_en}: {list}",
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
    if issue.kind == 'unaligned':
        d['lengths'] = ", ".join(str(x) for x in d.get('lengths', []))
    template = _MESSAGES.get(issue.kind, {}).get(lang)
    if not template:
        return f"{issue.kind}: {issue.data}"
    try:
        return template.format(**d)
    except (KeyError, IndexError):
        return f"{issue.kind}: {issue.data}"
