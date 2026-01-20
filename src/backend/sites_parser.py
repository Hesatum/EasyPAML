"""
Parser para extrair sítios sob seleção positiva de arquivos CODEML
Extrai dados de M8, M2a, Branch e outros modelos
"""

import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pandas as pd
import numpy as np


class SitesParser:
    """Parser para sítios sob seleção de arquivos CODEML"""
    
    @staticmethod
    def parse_sites_from_file(filepath: Path, method: str = "BEB") -> pd.DataFrame:
        """
        Parse sítios sob seleção de um arquivo de resultados CODEML
        
        Parameters:
        -----------
        filepath : Path
            Caminho para o arquivo de resultados
        method : str
            "BEB" ou "NEB" - qual análise usar
            
        Returns:
        --------
        pd.DataFrame
            Dataframe com colunas:
            - position: posição do aminoácido
            - amino_acid: aminoácido de uma letra
            - pr_w_gt_1: Pr(w>1) - p-value
            - post_mean: post mean for w (estimativa de omega)
            - post_se: standard error (se for disponível)
            - omega_lower: post mean - SE (limite inferior)
            - omega_upper: post mean + SE (limite superior)
            - significance: * ou ** baseado em p-value
        """
        
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        
        # Encontrar a seção apropriada (BEB ou NEB)
        # Pode ser "BEB analysis" ou "BEB) analysis"
        if method == "BEB":
            pattern = r"BEB\)?.*?analysis.*?Positively selected sites.*?\n\s*\(amino acids refer to.*?\)\s*\n\s*Pr\(w>1\).*?\n\n(.*?)(?:\n\n|Time used:|$)"
        else:  # NEB
            pattern = r"NEB analysis.*?Positively selected sites.*?\n\s*\(amino acids refer to.*?\)\s*\n\s*Pr\(w>1\).*?\n\n(.*?)(?:\n\n|Bayes|Time used:|$)"
        
        match = re.search(pattern, content, re.DOTALL | re.IGNORECASE)
        
        if not match:
            return pd.DataFrame()
        
        sites_text = match.group(1).strip()
        
        # Parser linha por linha
        sites = []
        
        # Pattern para linhas com dados de sítios
        # Exemplos:
        # "   159 R      0.990**       8.976 +- 1.573"
        # "   233 P      0.989*        8.967 +- 1.596"
        site_pattern = r'\s*(\d+)\s+([A-Z])\s+([\d.]+)([\*]*)\s+([\d.]+)\s*\+-\s*([\d.]+)'
        
        for line in sites_text.split('\n'):
            if not line.strip():
                continue
            
            match = re.search(site_pattern, line)
            if match:
                position = int(match.group(1))
                amino_acid = match.group(2)
                pr_w_gt_1 = float(match.group(3))
                significance = match.group(4)  # * ou **
                post_mean = float(match.group(5))
                post_se = float(match.group(6))
                
                # Calcular limites de confiança
                omega_lower = post_mean - post_se
                omega_upper = post_mean + post_se
                
                sites.append({
                    'position': position,
                    'amino_acid': amino_acid,
                    'pr_w_gt_1': pr_w_gt_1,
                    'post_mean': post_mean,
                    'post_se': post_se,
                    'omega_lower': max(0, omega_lower),  # omega não pode ser negativo
                    'omega_upper': omega_upper,
                    'significance': significance,
                    'is_significant_95': pr_w_gt_1 >= 0.95,
                    'is_significant_99': pr_w_gt_1 >= 0.99
                })
        
        return pd.DataFrame(sites)
    
    @staticmethod
    def extract_omega_global(filepath: Path) -> Optional[float]:
        """
        Extrai omega global (dN/dS) da tabela 'dN & dS for each branch'
        Para modelos como M2a, M8, etc. onde omega é constante através da árvore
        """
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
                        # Padrão: branch t N S dN/dS dN dS N*dN S*dS
                        try:
                            parts = data_line.split()
                            if len(parts) >= 5:
                                # Encontrar coluna dN/dS (a 5ª coluna após divisão por espaço)
                                # Formato típico: "19..20    0.016  1148.8  351.2  0.4355  0.0040  0.0091  4.6  3.2"
                                omega_str = parts[4]
                                omega = float(omega_str)
                                
                                # Validar que é um valor razoável de omega
                                if -10 <= omega <= 100:
                                    omegas.append(omega)
                        except (ValueError, IndexError):
                            pass
                    break
            
            if omegas:
                # Retornar a mediana ou média dos omegas encontrados
                # Para modelos com omega constante, todos devem ser similares
                return float(np.median(omegas))
            
            return None
        except:
            return None
    
    @staticmethod
    def extract_omega_by_branches(filepath: Path) -> Dict[str, float]:
        """
        Extrai omega (dN/dS) para cada branch da árvore
        Útil para análise Branch-site model (branch model)
        
        Returns:
        --------
        Dict[str, float]
            Mapping de branch -> omega
        """
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
        except:
            return {}
    
    @staticmethod
    def extract_omega_values_from_model_params(filepath: Path) -> Optional[float]:
        """
        Extrai omega dos parâmetros do modelo (seção de model parameters)
        Útil para M2a, M8, M1a e outros modelos de sites
        Busca por padrões como "omega (w) for branches:"
        """
        try:
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            
            # Padrões comuns para omega em CODEML
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
        except:
            return None
    
    @staticmethod
    def extract_omega_robust(filepath: Path) -> Optional[float]:
        """
        Extrai omega usando múltiplas estratégias em ordem de preferência
        1. Tabela 'dN & dS for each branch'
        2. Parâmetros do modelo
        
        IMPORTANTE: Para modelos Branch-site e Branch-site_null, retorna None
        pois não há um único omega (múltiplas classes de sítios)
        
        Returns o primeiro valor válido encontrado, ou None se não aplicável
        """
        try:
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            
            # Branch-site Model (NSsites=2) com site classes - não há omega único
            # Detectar pela presença de "site class" + "background w" e "foreground w"
            if 'site class' in content and 'background w' in content and 'foreground w' in content:
                # Este é um Branch-site Model - não retorna omega único
                return None
        except:
            pass
        
        # Tentar primeiro método: tabela de branches
        omega = SitesParser.extract_omega_global(filepath)
        if omega is not None:
            return omega
        
        # Tentar segundo método: parâmetros do modelo
        omega = SitesParser.extract_omega_values_from_model_params(filepath)
        if omega is not None:
            return omega
        
        return None
    
    @staticmethod
    def extract_omega_by_tags(filepath: Path) -> Dict[str, float]:
        """
        Extrai omegas diferenciados por tags/labels de branches (Branch/BranchSite models)
        Para modelos onde diferentes branches têm omegas diferentes
        
        Returns:
        --------
        Dict[str, float]
            Mapping de tag/label -> omega
            Exemplos:
            - Branch: {'background': 0.42, '#1': 0.48, '#2': 0.28}
            - BranchSite: {'background': 0.09, 'foreground': 1.81}
            - Com placeholder: {'background': 0.42, '#1': 0.0663, '#2': 999.0}
        """
        try:
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            
            tag_omegas = {}
            
            # ═══ PRIMEIRO: BranchSite (padrão específico "background w" e "foreground w") ═══
            # Exemplo: "background w     0.09233  1.00000  0.09233  1.00000"
            #          "foreground w     0.09233  1.00000  1.81018  1.81018"
            # Pega o ÚLTIMO valor (que é o mais relevante para o modelo)
            branchsite_bg_pattern = r'background\s+w\s+([\d.\s]+?)(?:\n|$)'
            branchsite_fg_pattern = r'foreground\s+w\s+([\d.\s]+?)(?:\n|$)'
            
            bg_match = re.search(branchsite_bg_pattern, content)
            fg_match = re.search(branchsite_fg_pattern, content)
            
            if bg_match and fg_match:
                # Extrair último valor (mais relevante)
                bg_vals = bg_match.group(1).strip().split()
                fg_vals = fg_match.group(1).strip().split()
                
                if bg_vals and fg_vals:
                    try:
                        bg_omega = float(bg_vals[-1])  # Último valor
                        fg_omega = float(fg_vals[-1])  # Último valor
                        
                        if -10 <= bg_omega <= 100:
                            tag_omegas['background'] = bg_omega
                        if -10 <= fg_omega <= 100:
                            tag_omegas['foreground'] = fg_omega
                        
                        if tag_omegas:
                            return tag_omegas
                    except:
                        pass
            
            # ═══ SEGUNDO: Branch (múltiplos valores com placeholders) ═══
            # Linha: "w (dN/dS) for branches:  0.35052 0.06632 999.00000"
            # MANTÉM 999 para mostrar ao usuário que não há dados para aquela tag
            multi_omega_pattern = r'w\s*\(dN/dS\)\s*for\s+branches?\s*:\s*([\d.\s]+?)(?:\n|$)'
            multi_match = re.search(multi_omega_pattern, content, re.IGNORECASE)
            
            if multi_match:
                values_str = multi_match.group(1).strip()
                omega_values = []
                for val_str in values_str.split():
                    try:
                        val = float(val_str)
                        # ACEITA TODOS os valores, inclusive 999 (placeholder)
                        if -10 <= val <= 100 or val == 999:
                            omega_values.append(val)
                    except:
                        pass
                
                # Se encontrou 2 ou mais valores
                if len(omega_values) >= 2:
                    # Primeiro valor é sempre background
                    tag_omegas['background'] = omega_values[0]
                    
                    # Valores subsequentes são atribuídos a tags numéricas (#1, #2, etc.)
                    for i, omega in enumerate(omega_values[1:], 1):
                        tag_omegas[f'#{i}'] = omega
                    
                    return tag_omegas
                
                elif len(omega_values) == 1:
                    tag_omegas['background'] = omega_values[0]
                    return tag_omegas
            
            # ═══ FALLBACK: extrair valores únicos da tabela de branches ═══
            branches_dict = SitesParser.extract_omega_by_branches(filepath)
            if branches_dict:
                # Obter valores únicos
                unique_omegas = {}
                for branch, omega in branches_dict.items():
                    omega_key = f"{omega:.6f}"
                    if omega_key not in unique_omegas:
                        unique_omegas[omega_key] = []
                    unique_omegas[omega_key].append(branch)
                
                # Se há múltiplos valores únicos
                if len(unique_omegas) > 1:
                    sorted_omegas = sorted([(float(k), v) for k, v in unique_omegas.items()], 
                                          key=lambda x: x[0])
                    
                    # Primeiro valor (menor) é provavelmente background
                    if sorted_omegas:
                        tag_omegas['background'] = sorted_omegas[0][0]
                        # Outros são tags
                        for i, (omega, branches_list) in enumerate(sorted_omegas[1:], 1):
                            tag_omegas[f'#{i}'] = omega
                    
                    return tag_omegas
                
                # Se todos iguais, retornar como background
                elif len(unique_omegas) == 1:
                    omega = float(list(unique_omegas.keys())[0])
                    tag_omegas['background'] = omega
                    return tag_omegas
            
            return tag_omegas
        
        except Exception as e:
            print(f"⚠️ Erro ao extrair omegas por tags: {e}")
            return {}
    
    @staticmethod
    def parse_sites_from_results_folder(output_folder: Path, models: List[str] = None) -> Dict[str, Dict]:
        """
        Parse sítios de múltiplos modelos em uma pasta de resultados
        
        Parameters:
        -----------
        output_folder : Path
            Pasta com resultados CODEML
        models : List[str], optional
            Lista de modelos a buscar (ex: ['M8', 'M2a', 'Branch'])
            Se None, busca todos disponíveis
            
        Returns:
        --------
        Dict[str, Dict]
            Estrutura: {gene_name: {model: dataframe}}
        """
        
        if models is None:
            models = ['M8', 'M2a', 'Branch']
        
        results = {}
        
        # Procurar arquivos de resultados por gene
        for results_file in output_folder.glob("*_results.txt"):
            # Extrair nome do gene e modelo
            match = re.search(r'(.+?)_([A-Za-z0-9]+)_results\.txt', results_file.name)
            if not match:
                continue
            
            gene_name = match.group(1)
            model_name = match.group(2)
            
            if model_name not in models:
                continue
            
            if gene_name not in results:
                results[gene_name] = {}
            
            # Tentar BEB primeiro, depois NEB
            for method in ['BEB', 'NEB']:
                df = SitesParser.parse_sites_from_file(results_file, method=method)
                if not df.empty:
                    results[gene_name][f'{model_name}_{method}'] = df
                    break
        
        return results
    
    @staticmethod
    def filter_sites_by_pvalue(df: pd.DataFrame, p_threshold: float = 0.95) -> pd.DataFrame:
        """
        Filtrar sítios por p-value
        
        Parameters:
        -----------
        df : pd.DataFrame
            DataFrame de sítios
        p_threshold : float
            Limiar de p-value (default 0.95)
            
        Returns:
        --------
        pd.DataFrame
            Sítios filtrados
        """
        if df.empty:
            return df
        
        return df[df['pr_w_gt_1'] >= p_threshold].sort_values('pr_w_gt_1', ascending=False)
    
    @staticmethod
    def get_codons_for_sites(fasta_file: Path, positions: List[int], ref_seq_index: int = 0) -> Dict[int, str]:
        """
        Obter codons para posições específicas de um arquivo FASTA
        
        Parameters:
        -----------
        fasta_file : Path
            Arquivo FASTA com sequências
        positions : List[int]
            Posições de aminoácidos (1-indexed)
        ref_seq_index : int
            Índice da sequência de referência
            
        Returns:
        --------
        Dict[int, str]
            Mapping de posição -> codon
        """
        
        from Bio import SeqIO
        
        codons = {}
        sequences = list(SeqIO.parse(fasta_file, 'fasta'))
        
        if not sequences:
            return codons
        
        ref_seq = str(sequences[ref_seq_index].seq)
        
        for pos in positions:
            # Converter posição de aminoácido para codon (1-indexed)
            codon_start = (pos - 1) * 3
            codon_end = codon_start + 3
            
            if codon_end <= len(ref_seq):
                codons[pos] = ref_seq[codon_start:codon_end]
        
        return codons
    
    @staticmethod
    def enrich_sites_with_codons(df: pd.DataFrame, fasta_file: Path) -> pd.DataFrame:
        """
        Enriquecer DataFrame de sítios com codons
        
        Parameters:
        -----------
        df : pd.DataFrame
            DataFrame de sítios
        fasta_file : Path
            Arquivo FASTA
            
        Returns:
        --------
        pd.DataFrame
            DataFrame com coluna 'codon' adicionada
        """
        
        if df.empty or not fasta_file.exists():
            df['codon'] = ''
            return df
        
        try:
            codons = SitesParser.get_codons_for_sites(fasta_file, df['position'].tolist())
            df['codon'] = df['position'].map(lambda pos: codons.get(pos, 'N/A'))
        except Exception as e:
            print(f"Erro ao extrair codons: {e}")
            df['codon'] = ''
        
        return df    
    @staticmethod
    def extract_branchsite_class_data(filepath: Path) -> Dict[str, dict]:
        """
        Extrai dados estruturados de classes de sítios do modelo Branch-site
        
        Retorna dicionário com dados para cada classe (0, 1, 2a, 2b):
        - prop: proporção da classe
        - bg_w: ω background
        - fg_w: ω foreground
        
        Exemplo de retorno:
        {
            '0': {'prop': 0.76190, 'bg_w': 0.09233, 'fg_w': 0.09233},
            '1': {'prop': 0.22315, 'bg_w': 1.00000, 'fg_w': 1.00000},
            '2a': {'prop': 0.01156, 'bg_w': 0.09233, 'fg_w': 1.81018},
            '2b': {'prop': 0.00339, 'bg_w': 1.00000, 'fg_w': 1.81018},
        }
        """
        try:
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            
            result = {}
            
            # Pattern para extrair "site class" line com proporções
            # Formato: "site class             0        1       2a       2b"
            # Segue: "proportion       0.76190  0.22315  0.01156  0.00339"
            
            lines = content.split('\n')
            
            for i, line in enumerate(lines):
                if 'site class' in line.lower():
                    # Próxima linha deve ser "proportion"
                    if i + 1 < len(lines) and 'proportion' in lines[i + 1].lower():
                        # Extrair classes da linha de headers
                        # Exemplo: "site class             0        1       2a       2b"
                        class_line = line.strip()
                        # Remover o prefixo "site class" ou similar
                        class_parts = re.sub(r'site\s+class', '', class_line, flags=re.IGNORECASE).strip().split()
                        
                        # Extrair proporções
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
            return {}