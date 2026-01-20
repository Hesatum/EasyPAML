"""
Report generator for structured analysis.
Creates separate tables for each model enabling detailed per-branch/clade analysis.
"""

import pandas as pd
from pathlib import Path
from typing import Dict, Optional
import json

try:
    from .sites_parser import SitesParser
    from .branch_extractor import BranchExtractor
except ImportError:
    from sites_parser import SitesParser
    from branch_extractor import BranchExtractor


class StructuredReportGenerator:
    """Generates structured reports separating analyses by model"""
    
    @staticmethod
    def generate_all_reports(results_folder: Path, 
                            output_folder: Optional[Path] = None) -> Dict[str, Path]:
        """
        Generate all structured reports
        
        Parameters:
        -----------
        results_folder : Path
            Root folder with CODEML results
        output_folder : Path, optional
            Output folder for reports. If None, uses results_folder/reports
            
        Returns:
        --------
        Dict[str, Path]
            Mapping {report_name: file_path}
        """
        
        if output_folder is None:
            output_folder = results_folder / 'reports'
            output_folder.mkdir(exist_ok=True)
        
        reports = {}
        models = ['M8', 'M2a', 'Branch', 'M1a', 'M0', 'M7']
        
        print(f"\n{'='*60}")
        print("GENERATING STRUCTURED REPORTS")
        print(f"{'='*60}\n")
        
        # 1. Generate branch summaries (only for Branch model)
        print("[1/4] Generating branch summaries...")
        # Only the 'Branch' model has per-branch statistics
        if (results_folder / 'Branch').exists():
            df = BranchExtractor.create_model_summary(results_folder, 'Branch')
            if not df.empty:
                file_path = output_folder / f'Branch_branches_summary.tsv'
                df.to_csv(file_path, sep='\t', index=False)
                reports[f'Branch_branches'] = file_path
                print(f"  [OK] Branch_branches_summary.tsv ({len(df)} rows)")
        
        # 2. Generate site tables per model
        print("\n[2/4] Generating positive sites tables...")
        sites_summary = StructuredReportGenerator._generate_sites_report(
            results_folder, output_folder, models
        )
        if sites_summary:
            reports.update(sites_summary)
        
        # 3. Generate structured LRT comparisons
        print("\n[3/4] Generating structured LRT comparisons...")
        lrt_report = StructuredReportGenerator._generate_lrt_report(
            results_folder, output_folder
        )
        if lrt_report:
            reports['lrt_comparisons'] = lrt_report
            print(f"  [OK] lrt_comparisons_detailed.tsv")
        
        # 4. Generate tree structures with omega (only for Branch model)
        print("\n[4/4] Generating tree structures with omega...")
        if (results_folder / 'Branch').exists():
            tree_data = BranchExtractor.export_model_branches_json(
                results_folder, 'Branch'
            )
            
            if tree_data:
                file_path = output_folder / f'Branch_tree_omega.json'
                with open(file_path, 'w') as f:
                    json.dump(tree_data, f, indent=2)
                reports[f'Branch_tree'] = file_path
                gene_count = len(tree_data)
                print(f"  [OK] Branch_tree_omega.json ({gene_count} genes)")
        
        # Final summary
        print(f"\n{'='*60}")
        print(f"[OK] Reports generated in: {output_folder}")
        print(f"  Total files: {len(reports)}")
        print(f"{'='*60}\n")
        
        return reports
    
    @staticmethod
    def _generate_sites_report(results_folder: Path, output_folder: Path, 
                              models: list) -> Dict[str, Path]:
        """Generate positively selected sites report per model"""
        
        reports = {}
        
        for model in models:
            model_dir = results_folder / model
            if not model_dir.exists():
                continue
            
            all_sites = []
            
            for results_file in model_dir.glob('*_results.txt'):
                gene_name = results_file.stem.rsplit('_', 2)[0]
                
                # Try BEB first, then NEB
                for method in ['BEB', 'NEB']:
                    df = SitesParser.parse_sites_from_file(results_file, method=method)
                    if not df.empty:
                        df.insert(0, 'Gene', gene_name)
                        df.insert(1, 'Model', model)
                        df.insert(2, 'Method', method)
                        all_sites.append(df)
                        break
            
            if all_sites:
                combined_df = pd.concat(all_sites, ignore_index=True)
                file_path = output_folder / f'{model}_positively_selected_sites.tsv'
                combined_df.to_csv(file_path, sep='\t', index=False)
                reports[f'{model}_sites'] = file_path
                site_count = len(combined_df)
                print(f"  [OK] {model}_positively_selected_sites.tsv ({site_count} sites)")
        
        return reports
    
    @staticmethod
    def _generate_lrt_report(results_folder: Path, output_folder: Path) -> Optional[Path]:
        """Generate structured LRT comparison report"""
        
        # Load main TSV for LRT comparisons
        tsv_path = results_folder / 'analysis_summary.tsv'
        if not tsv_path.exists():
            return None
        
        df = pd.read_csv(tsv_path, sep='\t')
        
        comparisons = []
        lrt_cols = [col for col in df.columns if col.startswith('lrt_')]
        
        for lrt_col in lrt_cols:
            comparison_name = lrt_col.replace('lrt_', '').replace('_', ' vs ')
            
            for idx, row in df.iterrows():
                gene = row['Gene']
                lrt_val = row[lrt_col]
                
                if pd.isna(lrt_val):
                    continue
                
                # Extract models
                parts = lrt_col.replace('lrt_', '').split('_vs_')
                if len(parts) == 2:
                    model1, model2 = parts[0].upper(), parts[1].upper()
                    
                    # Get omegas
                    col1 = f'{model1}_omega'
                    col2 = f'{model2}_omega'
                    
                    omega1 = row.get(col1, None)
                    omega2 = row.get(col2, None)
                    
                    comparisons.append({
                        'Gene': gene,
                        'Comparison': comparison_name,
                        'Model1': model1,
                        'Model2': model2,
                        'LRT_2deltalnL': lrt_val,
                        'omega_Model1': omega1,
                        'omega_Model2': omega2
                    })
        
        if not comparisons:
            return None
        
        comp_df = pd.DataFrame(comparisons)
        file_path = output_folder / 'lrt_comparisons_detailed.tsv'
        comp_df.to_csv(file_path, sep='\t', index=False)
        
        return file_path
    
    @staticmethod
    def _generate_tree_reports(results_folder: Path, output_folder: Path, 
                              models: list) -> Dict[str, Path]:
        """Generate tree structures with omega for each model"""
        
        reports = {}
        
        for model in models:
            model_dir = results_folder / model
            if not model_dir.exists():
                continue
            
            tree_data = BranchExtractor.export_model_branches_json(
                results_folder, model
            )
            
            if tree_data:
                file_path = output_folder / f'{model}_tree_omega.json'
                with open(file_path, 'w') as f:
                    json.dump(tree_data, f, indent=2)
                reports[f'{model}_tree'] = file_path
                gene_count = len(tree_data)
                print(f"  [OK] {model}_tree_omega.json ({gene_count} genes)")
        
        return reports
    
    @staticmethod
    def generate_model_comparison_table(results_folder: Path, 
                                       output_file: Optional[Path] = None) -> pd.DataFrame:
        """
        Generate model comparison table with per-gene statistics
        
        Parameters:
        -----------
        results_folder : Path
            Results folder
        output_file : Path, optional
            Output file path
            
        Returns:
        --------
        pd.DataFrame
            Comparison table
        """
        
        tsv_path = results_folder / 'analysis_summary.tsv'
        if not tsv_path.exists():
            return pd.DataFrame()
        
        df = pd.read_csv(tsv_path, sep='\t')
        
        # Build comparison table
        comparison_data = []
        
        for idx, row in df.iterrows():
            gene = row['Gene']
            
            comparison_data.append({
                'Gene': gene,
                'M0_omega': row.get('M0_omega', None),
                'M1a_omega': row.get('M1a_omega', None),
                'M2a_omega': row.get('M2a_omega', None),
                'M7_omega': row.get('M7_omega', None),
                'M8_omega': row.get('M8_omega', None),
                'Branch_omega': row.get('Branch_omega', None),
                'M0_vs_M1a': row.get('lrt_M0_vs_M1a', None),
                'M1a_vs_M2a': row.get('lrt_M1a_vs_M2a', None),
                'M7_vs_M8': row.get('lrt_M7_vs_M8', None),
                'M0_vs_Branch': row.get('lrt_M0_vs_Branch', None),
            })
        
        comp_df = pd.DataFrame(comparison_data)
        
        if output_file:
            comp_df.to_csv(output_file, sep='\t', index=False)
        
        return comp_df
