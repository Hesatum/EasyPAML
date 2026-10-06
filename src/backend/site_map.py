"""
Map between codeml's site numbering and the user's alignment numbering.

With cleandata = 1 codeml removes columns (gaps, ambiguities, stop codons) and
numbers the BEB/NEB sites over the remaining columns, so codon 221 of the
alignment becomes site 220 if an earlier column was removed.

For each gene and model the backend writes MODEL/GENE_MODEL_sitemap.json with
the kept codons. This module reads it and adds 'position_original' to site
tables.
"""

import json
import re
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

SITEMAP_SUFFIX = "_sitemap.json"


def sitemap_path_for(results_file) -> Path:
    """MODELO/GENE_MODELO_results.txt -> MODELO/GENE_MODELO_sitemap.json"""
    p = Path(results_file)
    name = p.name
    if name.endswith("_results.txt"):
        name = name[: -len("_results.txt")]
    else:
        name = p.stem
    return p.with_name(name + SITEMAP_SUFFIX)


def write_sitemap(path, *, cleandata: int, n_codons: int, kept_codons: List[int],
                  sequences: List[str], codeml_sites: Optional[int] = None) -> Dict:
    verified = None if codeml_sites is None else (codeml_sites == len(kept_codons))
    data = {
        "description": "kept_codons[i] = codon position in the user's alignment of "
                       "codeml site i+1 (cleandata removes columns before numbering)",
        "cleandata": int(cleandata),
        "n_codons_in_alignment": int(n_codons),
        "n_codons_used_by_codeml": len(kept_codons),
        "codeml_reported_sites": codeml_sites,
        "verified": verified,
        "sequences_used": list(sequences),
        "kept_codons": list(kept_codons),
    }
    Path(path).write_text(json.dumps(data, indent=1), encoding="utf-8")
    return data


def read_sitemap(results_file) -> Optional[Dict]:
    p = sitemap_path_for(results_file)
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return None


def codeml_site_count(results_file) -> Optional[int]:
    """Number of sites (codons) codeml analysed: 'ls' in the 'ns = 10  ls = 299'
    header of its output. Used to check the map computed by EasyPAML."""
    try:
        text = Path(results_file).read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return None
    m = re.search(r"\bns\s*=\s*\d+\s+ls\s*=\s*(\d+)", text)
    return int(m.group(1)) if m else None


def to_original(position: int, sitemap: Optional[Dict]) -> Optional[int]:
    if not sitemap:
        return None
    kept = sitemap.get("kept_codons") or []
    if 1 <= position <= len(kept):
        return int(kept[position - 1])
    return None


def attach_original_positions(df: pd.DataFrame, results_file) -> pd.DataFrame:
    """Add 'position_original' (the user's alignment numbering).

    Without a sitemap, or if it failed its check, the column repeats 'position'
    and 'position_mapped' is False so the interface can warn."""
    if df is None or df.empty or "position" not in df.columns:
        return df
    sm = read_sitemap(results_file)
    df = df.copy()
    if sm and sm.get("verified") is not False:
        df["position_original"] = [to_original(int(p), sm) for p in df["position"]]
        df["position_mapped"] = True
    else:
        df["position_original"] = df["position"].astype(int)
        df["position_mapped"] = False
    return df
