"""
Codon alignment reading (FASTA and PHYLIP) and codon utilities.

PHYLIP may be sequential or interleaved, with strict (10-column) or relaxed
(space-separated) names; the reading that gives every sequence the length in
the header wins. codeml always receives a FASTA copy, so every accepted
variant behaves the same.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

ALIGNMENT_SUFFIXES = ('.fas', '.fasta', '.fa', '.fna', '.phy', '.phylip')

STOP_CODONS = frozenset({'TAA', 'TAG', 'TGA'})
NUCLEOTIDES = frozenset('ACGT')

_WS = re.compile(r'\s+')


class AlignmentError(ValueError):
    """A file that cannot be read as an alignment."""


@dataclass
class Alignment:
    path: Path
    fmt: str                                # 'fasta' | 'phylip'
    names: List[str]
    seqs: Dict[str, str]                    # name -> sequence (upper case, no spaces)
    phylip_variant: Optional[str] = None    # 'relaxed-sequential', 'strict-interleaved', ...
    duplicate_names: List[str] = field(default_factory=list)

    @property
    def lengths(self) -> List[int]:
        return [len(self.seqs[n]) for n in self.names]

    @property
    def is_aligned(self) -> bool:
        return len(set(self.lengths)) <= 1

    @property
    def length(self) -> int:
        return self.lengths[0] if self.names else 0


def _clean(seq: str) -> str:
    return _WS.sub('', seq).upper()


def _is_phylip_header(line: str) -> bool:
    parts = line.split()
    return len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit()


def parse_fasta(text: str) -> Tuple[List[str], Dict[str, str], List[str]]:
    names: List[str] = []
    chunks: Dict[str, List[str]] = {}
    dups: List[str] = []
    current = None
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith('>'):
            header = line[1:].strip()
            name = header.split()[0] if header else f"seq{len(names) + 1}"
            if name in chunks:
                dups.append(name)
                name = f"{name}__dup{len(names) + 1}"
            names.append(name)
            chunks[name] = []
            current = name
        elif current is not None:
            chunks[current].append(line)
    return names, {n: _clean(''.join(chunks[n])) for n in names}, dups


def _phylip_sequential(lines: List[str], ntax: int, nchar: int, strict: bool):
    names, seqs = [], {}
    cur = None
    for line in lines:
        if cur is None or len(seqs[cur]) >= nchar:
            if len(names) == ntax:
                return None  # text left over: wrong reading
            if strict:
                name, rest = line[:10].strip(), line[10:]
            else:
                parts = line.strip().split(None, 1)
                name, rest = parts[0], (parts[1] if len(parts) > 1 else '')
            if not name or name in seqs:
                return None
            names.append(name)
            seqs[name] = _clean(rest)
            cur = name
        else:
            seqs[cur] += _clean(line)
    if len(names) != ntax or any(len(seqs[n]) != nchar for n in names):
        return None
    return names, seqs


def _phylip_interleaved(lines: List[str], ntax: int, nchar: int, strict: bool):
    if len(lines) < ntax:
        return None
    names, seqs = [], {}
    for line in lines[:ntax]:
        if strict:
            name, rest = line[:10].strip(), line[10:]
        else:
            parts = line.strip().split(None, 1)
            name, rest = parts[0], (parts[1] if len(parts) > 1 else '')
        if not name or name in seqs:
            return None
        names.append(name)
        seqs[name] = _clean(rest)
    for i, line in enumerate(lines[ntax:]):
        seqs[names[i % ntax]] += _clean(line)
    if any(len(seqs[n]) != nchar for n in names):
        return None
    return names, seqs


def parse_phylip(text: str) -> Tuple[List[str], Dict[str, str], str]:
    lines = [l.rstrip('\r') for l in text.splitlines() if l.strip()]
    if not lines or not _is_phylip_header(lines[0]):
        raise AlignmentError("missing PHYLIP header ('N  L')")
    ntax, nchar = (int(x) for x in lines[0].split()[:2])
    body = lines[1:]
    # try relaxed names first (codeml accepts them too), then strict 10-column names
    for variant, fn, strict in (
        ('relaxed-sequential', _phylip_sequential, False),
        ('strict-sequential', _phylip_sequential, True),
        ('relaxed-interleaved', _phylip_interleaved, False),
        ('strict-interleaved', _phylip_interleaved, True),
    ):
        parsed = fn(body, ntax, nchar, strict)
        if parsed:
            names, seqs = parsed
            return names, seqs, variant
    raise AlignmentError(
        f"PHYLIP não reconhecido: o cabeçalho declara {ntax} sequências de "
        f"{nchar} caracteres, mas nenhuma leitura (sequencial/intercalado, "
        f"nomes de 10 colunas ou separados por espaço) bate com isso"
    )


def read_alignment(path) -> Alignment:
    path = Path(path)
    try:
        text = path.read_text(encoding='utf-8', errors='replace')
    except OSError as exc:
        raise AlignmentError(f"não foi possível ler o arquivo: {exc}") from exc
    stripped = text.lstrip()
    if not stripped:
        raise AlignmentError("empty file")
    if stripped.startswith('>'):
        names, seqs, dups = parse_fasta(text)
        return Alignment(path, 'fasta', names, seqs, duplicate_names=dups)
    first = stripped.splitlines()[0]
    if _is_phylip_header(first):
        names, seqs, variant = parse_phylip(text)
        return Alignment(path, 'phylip', names, seqs, phylip_variant=variant)
    raise AlignmentError("unknown format (expected FASTA '>' or PHYLIP 'N  L')")


def to_fasta(names: List[str], seqs: Dict[str, str], width: int = 60) -> str:
    out = []
    for n in names:
        s = seqs[n]
        out.append(f">{n}")
        out.extend(s[i:i + width] for i in range(0, len(s), width))
    return "\n".join(out) + "\n"


def codons(seq: str) -> List[str]:
    return [seq[i:i + 3] for i in range(0, len(seq) - len(seq) % 3, 3)]


def find_stop_codons(names: List[str], seqs: Dict[str, str]) -> List[Tuple[str, int, str]]:
    """[(sequence, 1-based codon position, codon)] for every stop codon."""
    hits = []
    for n in names:
        for i, c in enumerate(codons(seqs[n].replace('U', 'T')), 1):
            if c in STOP_CODONS:
                hits.append((n, i, c))
    return hits


def cleandata_kept_codons(names: List[str], seqs: Dict[str, str]) -> List[int]:
    """1-based positions of the codons codeml keeps with cleandata = 1.

    codeml's rule (checked with 4.9j and 4.10): a codon column is removed if any
    sequence has something other than A/C/G/T there (gap, N, ?, ambiguity) or a
    stop codon. codeml numbers BEB sites over the remaining columns; this list
    maps them back to the alignment.
    """
    if not names:
        return []
    per_seq = [codons(seqs[n].replace('U', 'T')) for n in names]
    ncod = min(len(c) for c in per_seq)
    kept = []
    for i in range(ncod):
        ok = True
        for cs in per_seq:
            c = cs[i]
            if c in STOP_CODONS or any(ch not in NUCLEOTIDES for ch in c):
                ok = False
                break
        if ok:
            kept.append(i + 1)
    return kept
