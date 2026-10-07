"""
Methods paragraph built from what a run did: EasyPAML version and commit,
codeml version, models, .ctl parameters, each LRT with df and null
distribution, and the number of genes in each Benjamini-Hochberg correction.
Written to methods_text.txt in the output folder for the user to review.
"""

from typing import Dict, List, Optional, Sequence, Tuple

from . import lrt_stats
from .ctl_params import CODONFREQ_NAMES

_SITE_MODELS = ('M0', 'M1a', 'M2a', 'M7', 'M8', 'M8a')


def _join(items: Sequence[str]) -> str:
    items = list(items)
    if len(items) <= 1:
        return "".join(items)
    return ", ".join(items[:-1]) + " and " + items[-1]


def build_methods_text(*, version: str, codeml_version: Optional[str], models: List[str],
                       ctl: Dict[str, object], omega0: float, pruned: bool,
                       family_sizes: Dict[Tuple[str, str], int], n_genes: int,
                       beb: bool = True, masked_stops: Optional[Dict[str, int]] = None,
                       excluded_taxa: Optional[Dict[str, List[str]]] = None,
                       warm_start: bool = False, omega_starts: Sequence[float] = (),
                       n_own_trees: int = 0) -> str:
    cf = ctl.get('CodonFreq')
    cf_name = CODONFREQ_NAMES.get(int(cf), str(cf)) if cf is not None else '?'
    site = [m for m in _SITE_MODELS if m in models]
    other = [m.replace('_null', ' null') for m in models if m not in _SITE_MODELS]
    fitted = site + other

    url = "https://github.com/Hesatum/EasyPAML"
    cite = f"{version[:-1]}; {url})" if version.endswith(')') else f"{version} ({url})"
    s = [f"Selection was tested with EasyPAML {cite} "
         f"running codeml from PAML {codeml_version or 'v?'} (Yang 2007) on {n_genes} gene(s)."]
    s.append(f"Models {_join(fitted)} were fitted with codon frequencies {cf_name} "
             f"(CodonFreq = {cf})"
             + (f", the beta distribution discretised into {ctl.get('ncatG')} categories "
                f"(ncatG = {ctl.get('ncatG')})" if any(m in models for m in ('M7', 'M8', 'M8a')) else "")
             + f", κ and ω estimated from initial values {ctl.get('kappa')} and {omega0}, and "
             + ("alignment columns with gaps, ambiguous characters or stop codons removed "
                "(cleandata = 1)." if int(ctl.get('cleandata', 1)) == 1
                else "all alignment columns kept (cleandata = 0)."))
    if n_own_trees:
        s.append(("Each gene was analysed with its own tree" if n_own_trees >= n_genes else
                  f"{n_own_trees} of the {n_genes} gene(s) were analysed with their own tree and the "
                  "others with a single tree")
                 + " (the gene-tree pairs are listed in run_config.json).")
    tree = "Input trees were unrooted for the site models"
    if pruned:
        tree += " and pruned to the taxa present in each alignment"
    s.append(tree + ".")
    if warm_start:
        s.append("The κ and branch lengths estimated under M0 were used as starting values for the "
                 "other models (fix_blength = 1)"
                 + (f", each fitted from initial ω {_join([f'{w:g}' for w in omega_starts])} keeping "
                    "the highest log-likelihood" if omega_starts else "") + ".")
    if masked_stops:
        genes = sorted(masked_stops)
        n = sum(masked_stops.values())
        s.append(f"{n} stop codon(s) in {len(genes)} gene(s) ({_join(genes)}) were treated as "
                 "missing data.")
    if excluded_taxa:
        parts = [f"{', '.join(names)} from {gene}" for gene, names in sorted(excluded_taxa.items())]
        s.append("Sequences absent from the tree were excluded: " + "; ".join(parts) + ".")

    tests = []
    for (null, alt) in lrt_stats.pairs_for(models):
        info = lrt_stats.PAIRS[(null, alt)]
        df = info['df'] if info['df'] is not None else 'number of foreground branch groups'
        tests.append(f"{alt} vs {null} (df = {df}"
                     + (", χ²₁, conservative relative to the 50:50 χ²₀/χ²₁ mixture" if info['boundary'] else "")
                     + f"; {family_sizes.get((null, alt), 0)} gene(s))")
    if tests:
        s.append("Nested models were compared with likelihood ratio tests: " + "; ".join(tests) + ". "
                 "P-values were corrected for multiple testing with the Benjamini-Hochberg procedure "
                 "within each test, across the genes of the run (numbers above), and genes with "
                 "q < 0.05 were considered significant.")
    if beb and any(m in models for m in ('M2a', 'M8')):
        s.append("Sites under positive selection were identified with the Bayes Empirical Bayes "
                 "procedure (posterior probability ≥ 0.95).")
    return " ".join(s)
