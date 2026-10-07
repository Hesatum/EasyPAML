"""Sites under selection, site classes and ω values from codeml output."""

import re
from pathlib import Path
from typing import Dict, List, Optional
import pandas as pd
import numpy as np


class SitesParser:
    """Readers for codeml output files."""
    
    @staticmethod
    def parse_sites_from_file(filepath: Path, method: str = "BEB") -> pd.DataFrame:
        """BEB or NEB site table: position, amino_acid, pr_w_gt_1, post_mean,
        post_se, omega_lower, omega_upper and significance (* or **)."""
        
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()

        # Normalizar line endings (CODEML no Windows pode gerar \r\n)
        content = content.replace('\r\n', '\n').replace('\r', '\n')

        # Pode ser "BEB analysis", "BEB) analysis" ou "Bayes Empirical Bayes (BEB) analysis"
        # \s* between header and data tolerates spacing differences between codeml versions
        if method == "BEB":
            pattern = (
                r"BEB\b.*?analysis"
                r".*?Positively selected sites"
                r".*?\(amino acids refer to[^\)]*\)"
                r".*?Pr\(w>1\)[^\n]*\n"                  # linha de header da tabela
                r"\s*(.*?)"
                r"(?:\n\s*\n|Time used:|$)"              # termina em linha vazia ou fim
            )
        else:  # NEB
            pattern = (
                r"NEB\b.*?analysis"
                r".*?Positively selected sites"
                r".*?\(amino acids refer to[^\)]*\)"
                r".*?Pr\(w>1\)[^\n]*\n"
                r"\s*(.*?)"
                r"(?:\n\s*\n|Bayes|Time used:|$)"
            )

        match = re.search(pattern, content, re.DOTALL | re.IGNORECASE)

        if match:
            sites_text = match.group(1).strip()
            # M2a/M8 site lines, with
            # media posterior de omega +- erro padrao):
            # "   159 R      0.990**       8.976 +- 1.573"
            site_pattern = r'\s*(\d+)\s+([A-Z])\s+([\d.]+)([\*]*)\s+([\d.]+)\s*\+-\s*([\d.]+)'

            sites = []
            for line in sites_text.split('\n'):
                if not line.strip():
                    continue
                m = re.search(site_pattern, line)
                if m:
                    position = int(m.group(1))
                    amino_acid = m.group(2)
                    pr_w_gt_1 = float(m.group(3))
                    significance = m.group(4)
                    post_mean = float(m.group(5))
                    post_se = float(m.group(6))
                    omega_lower = post_mean - post_se
                    omega_upper = post_mean + post_se
                    sites.append({
                        'position': position,
                        'amino_acid': amino_acid,
                        'pr_w_gt_1': pr_w_gt_1,
                        'post_mean': post_mean,
                        'post_se': post_se,
                        'omega_lower': max(0, omega_lower),
                        'omega_upper': omega_upper,
                        'significance': significance,
                        'is_significant_95': pr_w_gt_1 >= 0.95,
                        'is_significant_99': pr_w_gt_1 >= 0.99
                    })
            return pd.DataFrame(sites)

        # Branch-site (model = 2, NSsites = 2) uses a
        # cabecalho diferente ("Positive sites for foreground lineages
        # Prob(w>1):" header and each line has only position, amino acid and
        # probability, without the posterior mean of ω
        # (e.g. "   987 S 0.993**")
        section_header = "BEB" if method == "BEB" else "NEB"
        bs_pattern = (
            rf"{section_header}\b.*?analysis"
            r".*?Positive sites for foreground lineages Prob\(w>1\):\s*\n"
            r"(.*?)"
            r"(?:\n\s*\n|Time used:|$)"
        )
        bs_match = re.search(bs_pattern, content, re.DOTALL | re.IGNORECASE)
        if not bs_match:
            return pd.DataFrame()

        bs_site_pattern = r'\s*(\d+)\s+([A-Z])\s+([\d.]+)([\*]*)\s*$'
        sites = []
        for line in bs_match.group(1).strip().split('\n'):
            if not line.strip():
                continue
            m = re.search(bs_site_pattern, line)
            if m:
                position = int(m.group(1))
                pr_w_gt_1 = float(m.group(3))
                sites.append({
                    'position': position,
                    'amino_acid': m.group(2),
                    'pr_w_gt_1': pr_w_gt_1,
                    'post_mean': np.nan,   # not reported for branch-site
                    'post_se': np.nan,
                    'omega_lower': np.nan,
                    'omega_upper': np.nan,
                    'significance': m.group(4),
                    'is_significant_95': pr_w_gt_1 >= 0.95,
                    'is_significant_99': pr_w_gt_1 >= 0.99
                })
        return pd.DataFrame(sites)
    
    @staticmethod
    def extract_omega_global(filepath: Path) -> Optional[float]:
        """ω from the 'dN & dS for each branch' table, for models with one ω across the tree."""
        try:
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            
            # Encontrar a tabela "dN & dS for each branch"
            if 'dN & dS for each branch' not in content:
                return None
            
            lines = content.split('\n')
            omegas = []
            
            for i, line in enumerate(lines):
                if 'dN & dS for each branch' in line:
                    # Procurar primeira linha de dados (pula headers e linhas vazias)
                    for j in range(i + 1, min(i + 200, len(lines))):
                        data_line = lines[j].strip()
                        
                        # Pular linhas vazias e headers
                        if not data_line or 'branch' in data_line.lower():
                            continue
                        
                        # Parar se encontrar outro header ou fim de tabela
                        if any(x in data_line.lower() for x in ['tree length', 'mlc', 'model', '---', 'dS tree', 'dN tree']):
                            break
                        
                        # Tentar extrair dN/dS (5ª coluna, index 4)
                        try:
                            parts = data_line.split()
                            if len(parts) >= 5:
                                omega_str = parts[4]
                                omega = float(omega_str)
                                
                                if -10 <= omega <= 100:
                                    omegas.append(omega)
                        except (ValueError, IndexError):
                            pass
                    break
            
            if omegas:
                return float(np.median(omegas))

            return None
        except Exception:
            return None

    @staticmethod
    def extract_omega_by_branches(filepath: Path) -> Dict[str, float]:
        """{branch: dN/dS}."""
        try:
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            
            if 'dN & dS for each branch' not in content:
                return {}
            
            branch_omegas = {}
            lines = content.split('\n')
            
            for i, line in enumerate(lines):
                if 'dN & dS for each branch' in line:
                    # Procurar linhas de dados
                    for j in range(i + 1, min(i + 200, len(lines))):
                        data_line = lines[j].strip()
                        
                        if not data_line or 'branch' in data_line.lower():
                            continue
                        
                        # Parar se encontrar fim de tabela
                        if any(x in data_line.lower() for x in ['tree length', 'dS tree', 'dN tree', '---']):
                            break
                        
                        try:
                            parts = data_line.split()
                            if len(parts) >= 5:
                                branch = parts[0]  # Exemplo: "19..20", "26..13"
                                omega_str = parts[4]
                                omega = float(omega_str)
                                
                                if -10 <= omega <= 100:
                                    branch_omegas[branch] = omega
                        except (ValueError, IndexError):
                            pass
                    break
            
            return branch_omegas
        except Exception:
            return {}

    @staticmethod
    def extract_omega_values_from_model_params(filepath: Path) -> Optional[float]:
        """ω from the model parameter lines (e.g. "omega (w) for branches:")."""
        try:
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            
            omega_patterns = [
                r'omega \(w\) for branches:\s*([\d.]+)',
                r'w\s*=\s*([\d.]+)',
                r'dN/dS.*?=\s*([\d.]+)',
                r'w \(dN/dS\)\s*=\s*([\d.]+)',
                r'\bw\s+=\s*([\d.]+)',
            ]
            
            for pattern in omega_patterns:
                match = re.search(pattern, content, re.IGNORECASE)
                if match:
                    omega = float(match.group(1))
                    if -10 <= omega <= 100:
                        return omega
            
            return None
        except Exception:
            return None

    @staticmethod
    def extract_omega_robust(filepath: Path) -> Optional[float]:
        """ω from the branch table, else from the model parameters. None for the
        branch-site models, which have no single ω."""
        try:
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            
            if 'site class' in content and 'background w' in content and 'foreground w' in content:
                return None
        except Exception:
            pass

        omega = SitesParser.extract_omega_global(filepath)
        if omega is not None:
            return omega
        
        omega = SitesParser.extract_omega_values_from_model_params(filepath)
        if omega is not None:
            return omega
        
        return None
    
    @staticmethod
    def extract_site_classes(filepath: Path) -> Optional[Dict[str, List[float]]]:
        """Proportions and ω of each site class ('p:' and 'w:' after 'MLEs of
        dN/dS (w) for site classes'), or None."""
        try:
            text = Path(filepath).read_text(encoding='utf-8', errors='ignore')
        except OSError:
            return None
        matches = list(re.finditer(
            r'MLEs of dN/dS \(w\) for site classes[^\n]*\n\s*\n\s*p:\s*([^\n]+)\n\s*w:\s*([^\n]+)', text))
        m = matches[-1] if matches else None   # last block, in files with several models
        if not m:
            return None
        try:
            p = [float(x) for x in m.group(1).split()]
            w = [float(x) for x in m.group(2).split()]
        except ValueError:
            return None
        if not p or len(p) != len(w):
            return None
        return {'p': p, 'w': w}

    @staticmethod
    def extract_positive_class(filepath: Path) -> Optional[Dict[str, float]]:
        """ω and proportion (p₁) of the M2a/M8 class that may have ω > 1: the last
        class of 'MLEs of dN/dS (w) for site classes'. Unlike the mean ω, it is
        not diluted by the other classes."""
        classes = SitesParser.extract_site_classes(filepath)
        if not classes:
            return None
        return {'p': classes['p'][-1], 'omega': classes['w'][-1]}

    @staticmethod
    def extract_omega_by_tags(filepath: Path) -> Dict[str, float]:
        """ω per branch label, e.g. {'background': 0.42, '#1': 0.48} for Branch and
        {'background': 0.09, 'foreground': 1.81} for Branch-site. 999 is codeml's
        placeholder for a label without data and is kept."""
        try:
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            
            tag_omegas = {}
            
            # Exemplo: "background w     0.09233  1.00000  0.09233  1.00000"
            #          "foreground w     0.09233  1.00000  1.81018  1.81018"
            branchsite_bg_pattern = r'background\s+w\s+([\d.\s]+?)(?:\n|$)'
            branchsite_fg_pattern = r'foreground\s+w\s+([\d.\s]+?)(?:\n|$)'
            
            bg_match = re.search(branchsite_bg_pattern, content)
            fg_match = re.search(branchsite_fg_pattern, content)
            
            if bg_match and fg_match:
                bg_vals = bg_match.group(1).strip().split()
                fg_vals = fg_match.group(1).strip().split()
                
                if bg_vals and fg_vals:
                    try:
                        bg_omega = float(bg_vals[-1])
                        fg_omega = float(fg_vals[-1])

                        if -10 <= bg_omega <= 100:
                            tag_omegas['background'] = bg_omega
                        if -10 <= fg_omega <= 100:
                            tag_omegas['foreground'] = fg_omega

                        if tag_omegas:
                            return tag_omegas
                    except Exception:
                        pass
            
            # Linha: "w (dN/dS) for branches:  0.35052 0.06632 999.00000"
            multi_omega_pattern = r'w\s*\(dN/dS\)\s*for\s+branches?\s*:\s*([\d.\s]+?)(?:\n|$)'
            multi_match = re.search(multi_omega_pattern, content, re.IGNORECASE)
            
            if multi_match:
                values_str = multi_match.group(1).strip()
                omega_values = []
                for val_str in values_str.split():
                    try:
                        val = float(val_str)
                        if -10 <= val <= 100 or val == 999:
                            omega_values.append(val)
                    except Exception:
                        pass
                
                # Se encontrou 2 ou mais valores
                if len(omega_values) >= 2:
                    tag_omegas['background'] = omega_values[0]
                    
                    for i, omega in enumerate(omega_values[1:], 1):
                        tag_omegas[f'#{i}'] = omega
                    
                    return tag_omegas
                
                elif len(omega_values) == 1:
                    tag_omegas['background'] = omega_values[0]
                    return tag_omegas
            
            branches_dict = SitesParser.extract_omega_by_branches(filepath)
            if branches_dict:
                unique_omegas = {}
                for branch, omega in branches_dict.items():
                    omega_key = f"{omega:.6f}"
                    if omega_key not in unique_omegas:
                        unique_omegas[omega_key] = []
                    unique_omegas[omega_key].append(branch)
                
                if len(unique_omegas) > 1:
                    sorted_omegas = sorted([(float(k), v) for k, v in unique_omegas.items()], 
                                          key=lambda x: x[0])
                    
                    if sorted_omegas:
                        tag_omegas['background'] = sorted_omegas[0][0]
                        for i, (omega, branches_list) in enumerate(sorted_omegas[1:], 1):
                            tag_omegas[f'#{i}'] = omega
                    
                    return tag_omegas
                
                elif len(unique_omegas) == 1:
                    omega = float(list(unique_omegas.keys())[0])
                    tag_omegas['background'] = omega
                    return tag_omegas
            
            return tag_omegas
        
        except Exception as e:
            print(f"[WARN] could not read ω per label: {e}")
            return {}
    
    @staticmethod
    def filter_sites_by_pvalue(df: pd.DataFrame, p_threshold: float = 0.95) -> pd.DataFrame:
        """Sites with Pr(ω>1) >= p_threshold."""
        if df.empty:
            return df
        
        return df[df['pr_w_gt_1'] >= p_threshold].sort_values('pr_w_gt_1', ascending=False)
    
    @staticmethod
    def extract_branchsite_class_data(filepath: Path) -> Dict[str, dict]:
        """Branch-site classes (0, 1, 2a, 2b) with proportion, background ω and
        foreground ω, for example:

        {
            '0': {'prop': 0.76190, 'bg_w': 0.09233, 'fg_w': 0.09233},
            '1': {'prop': 0.22315, 'bg_w': 1.00000, 'fg_w': 1.00000},
            '2a': {'prop': 0.01156, 'bg_w': 0.09233, 'fg_w': 1.81018},
            '2b': {'prop': 0.00339, 'bg_w': 1.00000, 'fg_w': 1.81018},
        }"""
        try:
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            
            result = {}
            
            # Formato: "site class             0        1       2a       2b"
            # Segue: "proportion       0.76190  0.22315  0.01156  0.00339"
            
            lines = content.split('\n')
            
            for i, line in enumerate(lines):
                if 'site class' in line.lower():
                    if i + 1 < len(lines) and 'proportion' in lines[i + 1].lower():
                        # Extrair classes da linha de headers
                        # Exemplo: "site class             0        1       2a       2b"
                        class_line = line.strip()
                        # Remover o prefixo "site class" ou similar
                        class_parts = re.sub(r'site\s+class', '', class_line, flags=re.IGNORECASE).strip().split()
                        
                        prop_line = lines[i + 1].strip()
                        # Remover prefixo "proportion"
                        prop_parts = re.sub(r'proportion', '', prop_line, flags=re.IGNORECASE).strip().split()
                        
                        # Extrair background w
                        bg_line = None
                        fg_line = None
                        for j in range(i + 2, min(i + 20, len(lines))):
                            if 'background w' in lines[j].lower():
                                bg_line = lines[j]
                            if 'foreground w' in lines[j].lower():
                                fg_line = lines[j]
                        
                        if bg_line and fg_line:
                            # Extrair valores de background w
                            bg_vals = re.sub(r'background\s+w', '', bg_line, flags=re.IGNORECASE).strip().split()
                            # Extrair valores de foreground w
                            fg_vals = re.sub(r'foreground\s+w', '', fg_line, flags=re.IGNORECASE).strip().split()
                            
                            # Montar resultado
                            for idx, cls in enumerate(class_parts):
                                if idx < len(prop_parts) and idx < len(bg_vals) and idx < len(fg_vals):
                                    try:
                                        result[cls] = {
                                            'prop': float(prop_parts[idx]),
                                            'bg_w': float(bg_vals[idx]),
                                            'fg_w': float(fg_vals[idx])
                                        }
                                    except (ValueError, IndexError):
                                        pass
                        
                        return result
            
            return {}
        except Exception as e:
            print(f"[WARN] could not read branch-site classes in '{filepath.name}': {e}")
            return {}