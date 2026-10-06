"""
codeml control file (.ctl) parameters.

Every parameter is written explicitly so results do not depend on codeml's
internal defaults. Without ncatG, for example, codeml (4.9j and 4.10.10) uses 4
beta categories when M7 or M8 runs alone and 10 when 'NSsites = 7 8' runs in
one .ctl, and lnL changes.

Defaults follow the reference codeml.ctl shipped with PAML, except CodonFreq
(F3x4 instead of F61) and ncatG = 10, the usual choice for M7/M8.
"""

from collections import OrderedDict
from typing import Dict, List, Tuple

# (value, short name, description), as coded in the codeml section of pamlDOC
CODONFREQ_OPTIONS: List[Tuple[int, str, str]] = [
    (0, "Fequal", "1/61 each"),
    (1, "F1x4", "F1x4"),
    (2, "F3x4", "F3x4"),
    (3, "F61", "codon table"),
    (4, "F1x4MG", "F1x4MG"),
    (5, "F3x4MG", "F3x4MG"),
    (6, "FMutSel0", "FMutSel0"),
    (7, "FMutSel", "FMutSel"),
]

DEFAULT_CODONFREQ = 2  # F3x4

CODONFREQ_NAMES: Dict[int, str] = {v: name for v, name, _ in CODONFREQ_OPTIONS}


def codonfreq_label(value) -> str:
    """'2 = F3x4', the label used in the interface and the log."""
    try:
        v = int(value)
    except (TypeError, ValueError):
        return str(value)
    name = CODONFREQ_NAMES.get(v)
    return f"{v} = {name}" if name else str(v)


def parse_codonfreq_label(label) -> int:
    """Inverse of codonfreq_label: '2 = F3x4' or '2' -> 2."""
    text = str(label).strip()
    head = text.split('=', 1)[0].strip()
    return int(head)


# .ctl parameters that do not depend on the model, in the order of the
# reference codeml.ctl
DEFAULT_CTL_PARAMS: "OrderedDict[str, object]" = OrderedDict([
    ('noisy', 1),
    ('verbose', 1),
    ('runmode', 0),
    ('seqtype', 1),
    ('CodonFreq', DEFAULT_CODONFREQ),
    ('estFreq', 0),
    ('ndata', 1),
    ('clock', 0),
    ('aaDist', 0),
    ('icode', 0),
    ('Mgene', 0),
    ('fix_kappa', 0),
    ('kappa', 2),
    ('fix_alpha', 1),
    ('alpha', 0),
    ('Malpha', 0),
    ('ncatG', 10),
    ('getSE', 0),
    ('RateAncestor', 0),
    ('Small_Diff', 0.5e-6),
    ('cleandata', 1),
    ('fix_blength', 0),
    ('method', 0),
])

_COMMENTS = {
    'seqfile': 'alignment (copy next to this .ctl)',
    'treefile': 'tree actually used (pruned/unrooted, next to this .ctl)',
    'outfile': 'raw codeml output',
    'noisy': '0-9: stdout detail',
    'verbose': '1: detailed output',
    'runmode': '0: user tree',
    'seqtype': '1: codons',
    'CodonFreq': '0:Fequal 1:F1x4 2:F3x4 3:F61 4:F1x4MG 5:F3x4MG 6:FMutSel0 7:FMutSel',
    'estFreq': '0: observed frequencies',
    'ndata': 'number of data sets',
    'clock': '0: no clock',
    'aaDist': '0: equal amino acid distances',
    'model': '0: one omega for all branches; 2: omega per branch group',
    'NSsites': '0:M0 1:M1a 2:M2a 7:M7 8:M8',
    'icode': '0: universal genetic code',
    'Mgene': '0',
    'fix_kappa': '0: kappa estimated',
    'kappa': 'initial kappa (ts/tv)',
    'fix_omega': '0: omega estimated; 1: omega fixed',
    'omega': 'initial omega (fixed value if fix_omega = 1)',
    'fix_alpha': '1: no gamma rate variation',
    'alpha': '0: no gamma',
    'Malpha': '0',
    'ncatG': 'categories of the beta distribution (M7/M8)',
    'getSE': '0: no standard errors',
    'RateAncestor': '0: no ancestral reconstruction',
    'Small_Diff': 'convergence criterion',
    'cleandata': '1: remove columns with gaps/ambiguities/stops',
    'fix_blength': '0: ignore branch lengths in tree; 1: use as initial values; 2: fixed',
    'method': '0: simultaneous optimisation',
}

# Ordem final das linhas no .ctl
_CTL_ORDER = [
    'seqfile', 'treefile', 'outfile',
    'noisy', 'verbose', 'runmode',
    'seqtype', 'CodonFreq', 'estFreq', 'ndata', 'clock', 'aaDist',
    'model', 'NSsites', 'icode', 'Mgene',
    'fix_kappa', 'kappa', 'fix_omega', 'omega',
    'fix_alpha', 'alpha', 'Malpha', 'ncatG',
    'getSE', 'RateAncestor', 'Small_Diff', 'cleandata', 'fix_blength', 'method',
]


def _fmt_value(key: str, value) -> str:
    if key == 'Small_Diff':
        return f"{float(value):.1e}".replace('e-0', 'e-')
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def build_ctl_text(params: Dict[str, object]) -> str:
    """The .ctl text, one parameter per line, aligned on '=' like the reference codeml.ctl."""
    lines = []
    keys = [k for k in _CTL_ORDER if k in params] + \
           [k for k in params if k not in _CTL_ORDER]
    for key in keys:
        val = _fmt_value(key, params[key])
        comment = _COMMENTS.get(key)
        line = f"{key:>13} = {val}"
        if comment:
            # at least 3 spaces before '*': glued to a value ("file.fasta*")
            # codeml reads the '*' as part of the file name
            line = line.ljust(max(44, len(line) + 3)) + f"* {comment}"
        lines.append(line)
    return "\n".join(lines) + "\n"


def parse_ctl_text(text: str) -> Dict[str, str]:
    """Read 'key = value' pairs from a .ctl, ignoring comments after '*'."""
    out: Dict[str, str] = {}
    for raw in text.splitlines():
        line = raw.split('*', 1)[0].strip()
        if '=' not in line:
            continue
        k, v = line.split('=', 1)
        out[k.strip()] = v.strip()
    return out
