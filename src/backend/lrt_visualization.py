"""
LRT (Likelihood Ratio Test) Visualization Module
================================================

Creates publication-quality chi-square distribution plots for CODEML LRT comparisons.
Replicates the R visualization style using matplotlib.

Author: EasyPAML
License: MIT
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
import numpy as np
from typing import Dict, Tuple, Optional, List
import os


class LRTVisualizer:
    """Creates chi-square distribution plots for LRT comparisons."""
    
    # LRT comparison configurations
    LRT_COMPARISONS = {
        'M0_vs_M1a': {
            'title': 'A) M0 vs M1a',
            'xlabel': '2Δl',
            'df': 1,
            'null_model': 'M0',
            'alternative_model': 'M1a',
            'description': 'Site model: Neutral vs Nearly Neutral'
        },
        'M1a_vs_M2a': {
            'title': 'B) M1a vs M2a',
            'xlabel': '2Δl',
            'df': 2,
            'null_model': 'M1a',
            'alternative_model': 'M2a',
            'description': 'Site model: Nearly Neutral vs Positive Selection'
        },
        'M7_vs_M8': {
            'title': 'C) M7 vs M8',
            'xlabel': '2Δl',
            'df': 2,
            'null_model': 'M7',
            'alternative_model': 'M8',
            'description': 'Site model: Beta vs Beta+ω'
        }
    }
    
    # Critical values for alpha levels
    ALPHA_LEVELS = {
        0.05: {'name': 'α = 0.05', 'color': 'darkgray', 'label': 'χ²_{df,0.05}'},
        0.01: {'name': 'α = 0.01', 'color': 'brown', 'label': 'χ²_{df,0.01}'}
    }
    
    def __init__(self, figsize: Tuple[int, int] = (14, 4), dpi: int = 100):
        """Initialize the visualizer."""
        self.figsize = figsize
        self.dpi = dpi
        
    def plot_single_comparison(self, 
                              lrt_stat: float,
                              df: int,
                              comparison_name: str,
                              output_path: Optional[str] = None) -> Tuple[plt.Figure, plt.Axes]:
        """
        Create a single LRT comparison plot.
        
        Parameters:
        -----------
        lrt_stat : float
            The 2Δl (2 * log-likelihood difference) statistic
        df : int
            Degrees of freedom for the chi-square distribution
        comparison_name : str
            Name of comparison (e.g., 'M0_vs_M1a')
        output_path : str, optional
            Path to save the figure
            
        Returns:
        --------
        fig : matplotlib.figure.Figure
        ax : matplotlib.axes.Axes
        """
        config = self.LRT_COMPARISONS.get(comparison_name, {})
        
        fig, ax = plt.subplots(figsize=(6, 4), dpi=self.dpi)
        
        # Determine x-axis limits
        x_max = max(lrt_stat * 1.2, 20)
        x = np.linspace(0, x_max, 1000)
        y = stats.chi2.pdf(x, df)
        
        # Plot the chi-square distribution
        ax.plot(x, y, 'k-', linewidth=2)
        ax.fill_between(x, y, alpha=0.1, color='gray')
        
        # Calculate critical values
        crit_05 = stats.chi2.ppf(0.95, df)
        crit_01 = stats.chi2.ppf(0.99, df)
        
        # Add vertical lines for critical values and test statistic
        ax.axvline(crit_05, color='darkgray', linestyle='--', linewidth=1.5, label=f'χ²_{{{df},0.05}}={crit_05:.2f}')
        ax.axvline(crit_01, color='brown', linestyle='--', linewidth=1.5, label=f'χ²_{{{df},0.01}}={crit_01:.2f}')
        ax.axvline(lrt_stat, color='red', linewidth=2.5, label=f'2Δl = {lrt_stat:.2f}')
        
        # Calculate p-value
        p_value = 1 - stats.chi2.cdf(lrt_stat, df)
        
        # Add text annotations
        textstr = f"2Δl = {lrt_stat:.2f}\np-value = {p_value:.2e}" if p_value < 0.001 else f"2Δl = {lrt_stat:.2f}\np-value = {p_value:.4f}"
        ax.text(0.98, 0.97, textstr, 
               transform=ax.transAxes,
               verticalalignment='top',
               horizontalalignment='right',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
               fontsize=10)
        
        # Labels and title
        ax.set_xlabel('χ² value', fontsize=11, fontweight='bold')
        ax.set_ylabel('Density', fontsize=11, fontweight='bold')
        ax.set_title(config.get('title', comparison_name), fontsize=12, fontweight='bold')
        ax.legend(loc='upper right', fontsize=9)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if output_path:
            fig.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        
        return fig, ax
    
    def plot_multiple_comparisons(self,
                                 lrt_stats: Dict[str, float],
                                 output_path: Optional[str] = None) -> Tuple[plt.Figure, List[plt.Axes]]:
        """
        Create multiple LRT comparison plots (similar to R output).
        
        Parameters:
        -----------
        lrt_stats : dict
            Dictionary mapping comparison names to LRT statistics
            Example: {'M0_vs_M1a': 559.26, 'M1a_vs_M2a': 0.0, 'M7_vs_M8': 12.54}
        output_path : str, optional
            Path to save the figure
            
        Returns:
        --------
        fig : matplotlib.figure.Figure
        axes : list of matplotlib.axes.Axes
        """
        n_comparisons = len(lrt_stats)
        fig, axes = plt.subplots(1, n_comparisons, figsize=(5*n_comparisons, 4), dpi=self.dpi)
        
        if n_comparisons == 1:
            axes = [axes]
        
        for idx, (comparison_name, lrt_stat) in enumerate(lrt_stats.items()):
            ax = axes[idx]
            config = self.LRT_COMPARISONS.get(comparison_name, {})
            df = config.get('df', 1)
            
            # Determine x-axis limits
            x_max = max(lrt_stat * 1.2, 20)
            x = np.linspace(0, x_max, 1000)
            y = stats.chi2.pdf(x, df)
            
            # Plot chi-square distribution
            ax.plot(x, y, 'k-', linewidth=2)
            ax.fill_between(x, y, alpha=0.1, color='gray')
            
            # Calculate critical values
            crit_05 = stats.chi2.ppf(0.95, df)
            crit_01 = stats.chi2.ppf(0.99, df)
            
            # Add vertical lines
            ax.axvline(crit_05, color='darkgray', linestyle='--', linewidth=1.5)
            ax.axvline(crit_01, color='brown', linestyle='--', linewidth=1.5)
            ax.axvline(lrt_stat, color='red', linewidth=2.5)
            
            # Calculate p-value
            p_value = 1 - stats.chi2.cdf(lrt_stat, df)
            
            # Add annotations
            y_max = max(y)
            ax.text(crit_05 - 2, y_max * 0.9, f'χ²={crit_05:.2f}', 
                   fontsize=9, color='darkgray', ha='right')
            ax.text(crit_01 - 2, y_max * 0.8, f'χ²={crit_01:.2f}', 
                   fontsize=9, color='brown', ha='right')
            ax.text(lrt_stat, y_max * 0.95, f'2Δl={lrt_stat:.2f}', 
                   fontsize=9, color='red', fontweight='bold', ha='center',
                   bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))
            
            # P-value text
            p_str = f'{p_value:.2e}' if p_value < 0.001 else f'{p_value:.4f}'
            ax.text(0.98, 0.05, f'p-value = {p_str}', 
                   transform=ax.transAxes,
                   fontsize=10,
                   ha='right',
                   bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
            
            # Labels
            ax.set_xlabel(config.get('xlabel', '2Δl'), fontsize=10, fontweight='bold')
            ax.set_ylabel('Density', fontsize=10, fontweight='bold')
            ax.set_title(config.get('title', comparison_name), fontsize=11, fontweight='bold')
            ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if output_path:
            fig.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        
        return fig, axes
    
    def create_gene_lrt_report(self,
                              gene_name: str,
                              lrt_results: Dict[str, float],
                              output_dir: str) -> str:
        """
        Create a complete LRT report for a gene with all comparisons.
        
        Parameters:
        -----------
        gene_name : str
            Name of the gene
        lrt_results : dict
            Dictionary with LRT statistics for each comparison
        output_dir : str
            Directory to save the figure
            
        Returns:
        --------
        str : Path to the saved figure
        """
        # Filter only available comparisons
        available_comparisons = {k: v for k, v in lrt_results.items() 
                                if k in self.LRT_COMPARISONS}
        
        if not available_comparisons:
            return None
        
        fig, axes = self.plot_multiple_comparisons(available_comparisons)
        
        # Add overall title
        fig.suptitle(f'LRT Analysis - Gene: {gene_name}', 
                    fontsize=14, fontweight='bold', y=1.02)
        
        # Save figure
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f'LRT_{gene_name}.png')
        fig.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close(fig)
        
        return output_path


def extract_lrt_statistics_from_results(results_dict: Dict) -> Dict[str, float]:
    """
    Extract LRT statistics from results dictionary.
    
    Parameters:
    -----------
    results_dict : dict
        Dictionary containing LRT results (expected to have keys like 
        'lrt_M0_vs_M1a', 'lrt_M1a_vs_M2a', etc.)
        
    Returns:
    --------
    dict : Dictionary mapping comparison names to 2Δl values
    """
    lrt_stats = {}
    
    mapping = {
        'lrt_M0_vs_M1a': 'M0_vs_M1a',
        'lrt_M1a_vs_M2a': 'M1a_vs_M2a',
        'lrt_M7_vs_M8': 'M7_vs_M8'
    }
    
    for key, comparison_name in mapping.items():
        if key in results_dict and results_dict[key] is not None:
            try:
                lrt_stats[comparison_name] = float(results_dict[key])
            except (ValueError, TypeError):
                pass
    
    return lrt_stats


# Example usage:
if __name__ == '__main__':
    # Example LRT statistics (from the R script)
    example_lrt = {
        'M0_vs_M1a': 559.26,
        'M1a_vs_M2a': 0.0,
        'M7_vs_M8': 12.54
    }
    
    visualizer = LRTVisualizer()
    fig, axes = visualizer.plot_multiple_comparisons(example_lrt, 
                                                     output_path='lrt_comparison.png')
    print("✅ LRT comparison plots created successfully!")
