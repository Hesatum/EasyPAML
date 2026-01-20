"""
Módulo para extrair e processar dados de branches (ramos) de análises CODEML.
Permite capturar omega e estatísticas para cada branch/clade da árvore filogenética.
"""

import re
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import json


class BranchExtractor:
    """Extrai informações detalhadas de branches a partir de arquivos CODEML"""
    
    @staticmethod
    def extract_branch_table(filepath: Path) -> pd.DataFrame:
        """
        Extrai a tabela completa "dN & dS for each branch" de um arquivo CODEML
        
        Parameters:
        -----------
        filepath : Path
            Caminho para arquivo de resultados CODEML
            
        Returns:
        --------
        pd.DataFrame
            Tabela com colunas: branch, t, N, S, dN/dS, dN, dS, N*dN, S*dS
        """
        
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        
        # Procurar a seção "dN & dS for each branch"
        pattern = r'dN & dS for each branch\s*\n\s*branch\s+t\s+N\s+S\s+dN/dS\s+dN\s+dS\s+N\*dN\s+S\*dS(.*?)(?:\n\n|$)'
        match = re.search(pattern, content, re.DOTALL)
        
        if not match:
            return pd.DataFrame()
        
        table_text = match.group(1)
        rows = []
        
        for line in table_text.strip().split('\n'):
            line = line.strip()
            if not line:
                continue
            
            # Parse cada linha: "  19..20      0.010   1148.7    351.3   0.4354   0.0026   0.0059    2.9    2.1"
            parts = line.split()
            if len(parts) < 9:
                continue
            
            try:
                row = {
                    'branch': parts[0],
                    't': float(parts[1]),
                    'N': float(parts[2]),
                    'S': float(parts[3]),
                    'dN_dS': float(parts[4]),
                    'dN': float(parts[5]),
                    'dS': float(parts[6]),
                    'N_dN': float(parts[7]),
                    'S_dS': float(parts[8])
                }
                rows.append(row)
            except (ValueError, IndexError):
                continue
        
        return pd.DataFrame(rows)
    
    @staticmethod
    def extract_tree_structure(filepath: Path) -> Optional[str]:
        """
        Extrai a string da árvore filogenética do arquivo CODEML
        
        Parameters:
        -----------
        filepath : Path
            Caminho para arquivo de resultados CODEML
            
        Returns:
        --------
        str
            String da árvore em formato Newick
        """
        
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()
        
        # Procurar por linhas com árvore (geralmente começam com "(" e terminam com ";")
        tree_pattern = r'^\(\(.*\);?$'
        
        for line in lines:
            line = line.strip()
            if re.match(tree_pattern, line):
                # Remover nomes de espécies longas, manter só estrutura
                return line
        
        return None
    
    @staticmethod
    def extract_branch_omega_map(filepath: Path) -> Dict[str, float]:
        """
        Cria um mapeamento de branch_id -> omega para a árvore
        
        Parameters:
        -----------
        filepath : Path
            Caminho para arquivo de resultados CODEML
            
        Returns:
        --------
        Dict[str, float]
            Mapeamento {branch_id: omega_value}
        """
        
        df = BranchExtractor.extract_branch_table(filepath)
        if df.empty:
            return {}
        
        # Retornar dicionário com branch IDs como chaves e dN/dS como valores
        return dict(zip(df['branch'], df['dN_dS']))
    
    @staticmethod
    def create_model_summary(results_folder: Path, model_name: str, 
                            gene_name: Optional[str] = None) -> pd.DataFrame:
        """
        Cria um sumário compilado de todos os branches para um modelo específico
        
        Parameters:
        -----------
        results_folder : Path
            Pasta com resultados CODEML
        model_name : str
            Nome do modelo (M8, M2a, Branch, M1a, M0, M7)
        gene_name : str, optional
            Se fornecido, processa apenas este gene
            
        Returns:
        --------
        pd.DataFrame
            Tabela compilada com todos os genes
        """
        
        model_dir = results_folder / model_name
        if not model_dir.exists():
            return pd.DataFrame()
        
        all_rows = []
        
        # Procurar arquivos de resultados
        for results_file in model_dir.glob("*_results.txt"):
            # Extrair gene name
            match = re.search(r'(.+?)_' + re.escape(model_name) + r'_results\.txt', results_file.name)
            if not match:
                continue
            
            file_gene_name = match.group(1)
            
            if gene_name and file_gene_name != gene_name:
                continue
            
            # Extrair tabela de branches
            df = BranchExtractor.extract_branch_table(results_file)
            if df.empty:
                continue
            
            # Adicionar nome do gene
            df.insert(0, 'Gene', file_gene_name)
            all_rows.append(df)
        
        if not all_rows:
            return pd.DataFrame()
        
        return pd.concat(all_rows, ignore_index=True)
    
    @staticmethod
    def save_model_summaries(results_folder: Path, output_folder: Optional[Path] = None):
        """
        Salva sumários separados para cada modelo
        
        Parameters:
        -----------
        results_folder : Path
            Pasta com resultados CODEML
        output_folder : Path, optional
            Pasta para salvar sumários. Se None, usa results_folder
        """
        
        if output_folder is None:
            output_folder = results_folder
        
        models = ['M8', 'M2a', 'Branch', 'M1a', 'M0', 'M7']
        
        for model in models:
            df = BranchExtractor.create_model_summary(results_folder, model)
            
            if df.empty:
                continue
            
            output_file = output_folder / f'{model}_branches_summary.tsv'
            df.to_csv(output_file, sep='\t', index=False)
            print(f"Salvo: {output_file}")
    
    @staticmethod
    def create_tree_json_with_omega(filepath: Path) -> Dict:
        """
        Cria estrutura JSON da árvore com valores de omega nos nós
        
        Parameters:
        -----------
        filepath : Path
            Caminho para arquivo de resultados CODEML
            
        Returns:
        --------
        Dict
            Estrutura hierárquica com omega em cada nó
        """
        
        tree_str = BranchExtractor.extract_tree_structure(filepath)
        branch_omega = BranchExtractor.extract_branch_omega_map(filepath)
        
        if not tree_str:
            return {}
        
        # Parse árvore Newick (simplificado)
        # Formato: ((A:0.1,B:0.2):0.05,C:0.3);
        
        def parse_newick_to_json(newick_str: str, branch_map: Dict[str, float]) -> Dict:
            """Converte Newick para JSON com omega"""
            newick_str = newick_str.strip().rstrip(';')
            
            stack = []
            current_node = None
            i = 0
            
            while i < len(newick_str):
                char = newick_str[i]
                
                if char == '(':
                    # Iniciar novo nó
                    current_node = {'children': [], 'omega': None, 'branch_id': None}
                    stack.append(current_node)
                    
                elif char == ',':
                    # Separador entre nós filhos
                    pass
                    
                elif char == ')':
                    # Fechar nó
                    node = stack.pop()
                    
                    # Extrair identificador do nó (e.g., "19..20:0.010")
                    j = i + 1
                    node_id = ""
                    while j < len(newick_str) and newick_str[j] not in '(),;':
                        node_id += newick_str[j]
                        j += 1
                    
                    if node_id:
                        # Tentar extrair branch ID (parte antes do :)
                        branch_id = node_id.split(':')[0]
                        if branch_id in branch_map:
                            node['omega'] = branch_map[branch_id]
                            node['branch_id'] = branch_id
                        node['label'] = node_id
                    
                    i = j - 1
                    
                    if stack:
                        stack[-1]['children'].append(node)
                    else:
                        current_node = node
                
                else:
                    # Folha (nome de espécie)
                    j = i
                    label = ""
                    while j < len(newick_str) and newick_str[j] not in '(),;:':
                        label += newick_str[j]
                        j += 1
                    
                    if label:
                        leaf = {'label': label, 'omega': None, 'children': []}
                        if stack:
                            stack[-1]['children'].append(leaf)
                        else:
                            current_node = leaf
                    
                    i = j - 1
                
                i += 1
            
            return current_node or {}
        
        return parse_newick_to_json(tree_str, branch_omega)
    
    @staticmethod
    def export_model_branches_json(results_folder: Path, model_name: str, 
                                   output_file: Optional[Path] = None) -> Dict:
        """
        Exporta branches com omega como JSON para cada gene
        
        Parameters:
        -----------
        results_folder : Path
            Pasta com resultados CODEML
        model_name : str
            Nome do modelo
        output_file : Path, optional
            Arquivo para salvar JSON
            
        Returns:
        --------
        Dict
            Dados compilados {gene: {branch_id: omega, ...}}
        """
        
        model_dir = results_folder / model_name
        data = {}
        
        for results_file in model_dir.glob("*_results.txt"):
            match = re.search(r'(.+?)_' + re.escape(model_name) + r'_results\.txt', results_file.name)
            if not match:
                continue
            
            gene_name = match.group(1)
            omega_map = BranchExtractor.extract_branch_omega_map(results_file)
            
            if omega_map:
                data[gene_name] = omega_map
        
        if output_file:
            with open(output_file, 'w') as f:
                json.dump(data, f, indent=2)
            print(f"Salvo: {output_file}")
        
        return data
