"""Per-branch statistics (t, N, S, dN/dS, dN, dS) from codeml output."""

import re
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import json


class BranchExtractor:
    """Branch tables from codeml output files."""
    
    @staticmethod
    def extract_branch_table(filepath: Path) -> pd.DataFrame:
        """The "dN & dS for each branch" table: branch, t, N, S, dN/dS, dN, dS, N*dN, S*dS."""
        
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        
        # The column header row is followed by a blank line before the data rows,
        # so consume that blank line in the fixed part, then capture the data.
        pattern = r'dN & dS for each branch\s*\n\s*branch\s+t\s+N\s+S\s+dN/dS\s+dN\s+dS\s+N\*dN\s+S\*dS\s*\n\s*\n(.*?)(?:\n\s*\n|\Z)'
        match = re.search(pattern, content, re.DOTALL)
        
        if not match:
            return pd.DataFrame()
        
        table_text = match.group(1)
        rows = []
        
        for line in table_text.strip().split('\n'):
            line = line.strip()
            if not line:
                continue
            
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
        """The Newick tree printed in a codeml output file."""
        
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()
        
        tree_pattern = r'^\(\(.*\);?$'
        
        for line in lines:
            line = line.strip()
            if re.match(tree_pattern, line):
                return line
        
        return None
    
    @staticmethod
    def extract_branch_omega_map(filepath: Path) -> Dict[str, float]:
        """{branch_id: dN/dS}."""
        
        df = BranchExtractor.extract_branch_table(filepath)
        if df.empty:
            return {}
        
        return dict(zip(df['branch'], df['dN_dS']))
    
    @staticmethod
    def create_model_summary(results_folder: Path, model_name: str, 
                            gene_name: Optional[str] = None) -> pd.DataFrame:
        """Branch tables of every gene (or only gene_name) for one model."""
        
        model_dir = results_folder / model_name
        if not model_dir.exists():
            return pd.DataFrame()
        
        all_rows = []
        
        # Procurar arquivos de resultados
        for results_file in model_dir.glob("*_results.txt"):
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
            
            df.insert(0, 'Gene', file_gene_name)
            all_rows.append(df)
        
        if not all_rows:
            return pd.DataFrame()
        
        return pd.concat(all_rows, ignore_index=True)
    
    @staticmethod
    def save_model_summaries(results_folder: Path, output_folder: Optional[Path] = None):
        """Write one branch summary per model (to results_folder by default)."""
        
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
    def export_model_branches_json(results_folder: Path, model_name: str, 
                                   output_file: Optional[Path] = None) -> Dict:
        """{gene: {branch_id: dN/dS}} for one model, optionally saved as JSON."""
        
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
