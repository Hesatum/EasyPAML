"""
Mapa entre a numeração de sítios do codeml e a do alinhamento do usuário.

Com cleandata = 1 o codeml remove colunas (gaps, ambiguidades, stop codons)
e numera os sítios do BEB/NEB nas colunas que SOBRARAM. Sem correção, um
sítio que no alinhamento do usuário é o códon 221 aparece como 220 se uma
coluna anterior foi removida.

Para cada gene x modelo o backend grava, ao lado da saída bruta,
  MODELO/GENE_MODELO_sitemap.json
com a lista de códons mantidos. Este módulo lê esse arquivo e acrescenta a
coluna 'position_original' às tabelas de sítios.
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
    """Número de sítios (códons) que o codeml realmente analisou.

    O mlc traz, no cabeçalho, 'ns = 10  ls = 299' (ls = códons após a
    limpeza). Usado para verificar o mapa calculado pelo EasyPAML."""
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
    """Acrescenta 'position_original' (numeração do alinhamento do usuário).

    Sem sitemap (resultado antigo) ou com cleandata = 0, a coluna repete
    'position' e 'position_mapped' fica False para a interface avisar."""
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
