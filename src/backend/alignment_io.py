"""
Leitura de alinhamentos de códons (FASTA e PHYLIP) e utilitários de códons.

PHYLIP: aceita sequencial e intercalado, com nomes "estritos" (10 colunas)
ou "relaxados" (nome separado da sequência por espaço, qualquer tamanho).
O formato é escolhido pela interpretação em que todas as sequências ficam
com o comprimento declarado no cabeçalho.

O resto do EasyPAML sempre passa ao codeml uma cópia em FASTA (ver
codeml_backend), então qualquer variante aceita aqui funciona igual.
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
    """Arquivo que não pode ser lido como alinhamento."""


@dataclass
class Alignment:
    path: Path
    fmt: str                                # 'fasta' | 'phylip'
    names: List[str]
    seqs: Dict[str, str]                    # nome -> sequência (maiúsculas, sem espaços)
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
                return None  # sobrou texto: interpretação errada
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
        raise AlignmentError("cabeçalho PHYLIP ('N  L') ausente")
    ntax, nchar = (int(x) for x in lines[0].split()[:2])
    body = lines[1:]
    # Ordem de tentativa: relaxado primeiro (é o que o codeml também aceita,
    # nomes separados por espaço), depois estrito de 10 colunas.
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
        raise AlignmentError("arquivo vazio")
    if stripped.startswith('>'):
        names, seqs, dups = parse_fasta(text)
        return Alignment(path, 'fasta', names, seqs, duplicate_names=dups)
    first = stripped.splitlines()[0]
    if _is_phylip_header(first):
        names, seqs, variant = parse_phylip(text)
        return Alignment(path, 'phylip', names, seqs, phylip_variant=variant)
    raise AlignmentError("formato não reconhecido (esperado FASTA '>' ou PHYLIP 'N  L')")


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
    """[(sequência, posição do códon 1-based, códon)] para cada stop codon."""
    hits = []
    for n in names:
        for i, c in enumerate(codons(seqs[n].replace('U', 'T')), 1):
            if c in STOP_CODONS:
                hits.append((n, i, c))
    return hits


def cleandata_kept_codons(names: List[str], seqs: Dict[str, str]) -> List[int]:
    """Posições (1-based) dos códons que o codeml mantém com cleandata = 1.

    Regra do codeml (verificada com 4.9j e 4.10.x): uma coluna de códon sai
    da análise se, em QUALQUER sequência, o códon tiver algo diferente de
    A/C/G/T (gap, N, ?, ambiguidade) ou for um stop codon (o codeml converte a
    coluna inteira em '???'). As posições BEB do codeml são contadas nas
    colunas que sobram; este mapa permite voltar à numeração do alinhamento
    do usuário.
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
