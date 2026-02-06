"""
Interactive Analysis of FR1 Lick Microstructure Results

Load and explore the processed data from batch_lick_microstructure_FR1.py

Author: Analysis for Logan's DREADD study
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats

# Set style
sns.set_style("whitegrid")
sns.set_context("notebook", font_scale=1.1)


class FR1ResultsAnalyzer:
    """
    Interactive analyzer for FR1 lick microstructure results.
    """
    
    def __init__(self, results_csv_path):
        """
        Load results from CSV file.
        
        Parameters:
        -----------
        results_csv_path : str or Path
            Path to the FR1_lick_microstructure_all_mice.csv file
        """
        self.df = pd.read_csv(results_csv_path)
        self.results_folder = Path(results_csv_path).parent
        
        print(f"Loaded {len(self.df)} sessions")
        print(f"Mice: {self.df['mouse_id'].nunique()}")
        print(f"Groups: {self.df['group'].unique()}")
        print(f"Days range: {self.df['day'].min()} to {self.df['day'].max()}")
    
    def compare_groups(self, metric, day=None, show_individuals=True):
        """
        Compare groups on a specific metric.
        
        Parameters:
        -----------
        metric : str
            Column name to analyze (e.g., 'total_licks', 'mean_burst_size')
        day : int, optional
            Specific day to analyze. If None, uses all days
        show_individuals : bool
            Whether to show individual data points
        """
        # Filter by day if specified
        if day is not None:
            data = self.df[self.df['day'] == day].copy()
            title_suffix = f" - Day {day}"
        else:
            data = self.df.copy()
            title_suffix = " - All Days"
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Box plot
        sns.boxplot(data=data, x='group', y=metric, ax=ax, palette='Set2')
        
        # Add individual points if requested
        if show_individuals:
            sns.stripplot(data=data, x='group', y=metric, ax=ax, 
                         color='black', alpha=0.4, size=4)
        
        # Statistical comparison
        groups = data['group'].unique()
        if len(groups) == 2:
            # T-test for two groups
            group1_data = data[data['group'] == groups[0]][metric].dropna()
            group2_data = data[data['group'] == groups[1]][metric].dropna()
            
            t_stat, p_val = stats.ttest_ind(group1_data, group2_data)
            
            ax.text(0.5, 0.95, f'p = {p_val:.4f}', 
                   transform=ax.transAxes, ha='center', va='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        ax.set_title(f'{metric.replace("_", " ").title()}{title_suffix}',
                    fontweight='bold', fontsize=14)
        ax.set_xlabel('Group', fontsize=12)
        ax.set_ylabel(metric.replace('_', ' ').title(), fontsize=12)
        
        plt.tight_layout()
        plt.show()
        
        # Print statistics
        print(f"\n{metric} {title_suffix}:")
        print(data.groupby('group')[metric].describe())
    
    def plot_learning_curves(self, metrics=None, separate_cohorts=False):
        """
        Plot learning curves across days for specified metrics.
        
        Parameters:
        -----------
        metrics : list of str, optional
            Metrics to plot. If None, uses default set
        separate_cohorts : bool
            Whether to separate cohorts in the plot
        """
        if metrics is None:
            metrics = ['total_licks', 'n_bursts', 'mean_burst_size', 'mean_intraburst_freq_hz']
        
        n_metrics = len(metrics)
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        axes = axes.flatten()
        
        for idx, metric in enumerate(metrics):
            ax = axes[idx] if idx < 4 else plt.subplot(2, 2, idx + 1)
            
            if separate_cohorts:
                # Plot each cohort separately
                for cohort in sorted(self.df['cohort'].unique()):
                    cohort_data = self.df[self.df['cohort'] == cohort]
                    for group in sorted(cohort_data['group'].unique()):
                        group_data = cohort_data[cohort_data['group'] == group]
                        group_mean = group_data.groupby('day')[metric].agg(['mean', 'sem'])
                        
                        label = f'{group}-{cohort}'
                        ax.plot(group_mean.index, group_mean['mean'], 
                               marker='o', label=label, linewidth=2, alpha=0.7)
                        ax.fill_between(group_mean.index,
                                      group_mean['mean'] - group_mean['sem'],
                                      group_mean['mean'] + group_mean['sem'],
                                      alpha=0.15)
            else:
                # Plot groups collapsed across cohorts
                for group in sorted(self.df['group'].unique()):
                    group_data = self.df[self.df['group'] == group]
                    group_mean = group_data.groupby('day')[metric].agg(['mean', 'sem'])
                    
                    ax.plot(group_mean.index, group_mean['mean'], 
                           marker='o', label=group, linewidth=2.5)
                    ax.fill_between(group_mean.index,
                                  group_mean['mean'] - group_mean['sem'],
                                  group_mean['mean'] + group_mean['sem'],
                                  alpha=0.2)
            
            ax.set_xlabel('Training Day', fontsize=11)
            ax.set_ylabel(metric.replace('_', ' ').title(), fontsize=11)
            ax.set_title(metric.replace('_', ' ').title(), fontweight='bold')
            ax.legend(fontsize=9)
            ax.grid(alpha=0.3)
        
        plt.suptitle('FR1 Training: Learning Curves', 
                    fontsize=16, fontweight='bold', y=1.00)
        plt.tight_layout()
        plt.show()
    
    def analyze_individual_mouse(self, mouse_id):
        """
        Analyze a single mouse's data across days.
        
        Parameters:
        -----------
        mouse_id : str
            Mouse identifier
        """
        mouse_data = self.df[self.df['mouse_id'] == mouse_id].sort_values('day')
        
        if len(mouse_data) == 0:
            print(f"No data found for mouse: {mouse_id}")
            return
        
        group = mouse_data['group'].iloc[0]
        cohort = mouse_data['cohort'].iloc[0]
        
        print(f"\nMouse: {mouse_id}")
        print(f"Group: {group}")
        print(f"Cohort: {cohort}")
        print(f"Days trained: {len(mouse_data)}")
        print(f"Day range: {mouse_data['day'].min()} to {mouse_data['day'].max()}")
        
        # Create multi-panel figure
        fig, axes = plt.subplots(2, 3, figsize=(15, 8))
        axes = axes.flatten()
        
        metrics = [
            ('total_licks', 'Total Licks'),
            ('lick_rate_per_s', 'Lick Rate (licks/s)'),
            ('n_bursts', 'Number of Bursts'),
            ('mean_burst_size', 'Mean Burst Size'),
            ('mean_intraburst_freq_hz', 'Intraburst Frequency (Hz)'),
            ('mean_ili_ms', 'Mean ILI (ms)')
        ]
        
        for idx, (metric, title) in enumerate(metrics):
            ax = axes[idx]
            ax.plot(mouse_data['day'], mouse_data[metric], 
                   marker='o', linewidth=2, markersize=8, color='steelblue')
            ax.set_xlabel('Training Day')
            ax.set_ylabel(title)
            ax.set_title(title, fontweight='bold')
            ax.grid(alpha=0.3)
        
        plt.suptitle(f'{mouse_id} ({group}, {cohort})', 
                    fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.show()
        
        # Print summary table
        print("\nSummary Statistics:")
        print(mouse_data[['day', 'total_licks', 'n_bursts', 'mean_burst_size', 
                         'mean_intraburst_freq_hz']].to_string(index=False))
    
    def get_last_day_analysis(self):
        """
        Analyze the last training day for each mouse.
        """
        # Get last day for each mouse
        last_days = self.df.groupby('mouse_id')['day'].max().reset_index()
        last_days.columns = ['mouse_id', 'last_day']
        
        # Merge to get data for last day only
        last_day_data = self.df.merge(last_days, on='mouse_id')
        last_day_data = last_day_data[last_day_data['day'] == last_day_data['last_day']]
        
        print("\nLast Day Performance by Group:")
        
        metrics = ['total_licks', 'n_bursts', 'mean_burst_size', 'mean_intraburst_freq_hz']
        
        summary = last_day_data.groupby('group')[metrics].agg(['mean', 'std', 'count'])
        print(summary)
        
        return last_day_data
    
    def compare_first_vs_last(self):
        """
        Compare first vs last training day for each mouse.
        """
        # Get first and last day for each mouse
        first_last = []
        
        for mouse_id in self.df['mouse_id'].unique():
            mouse_data = self.df[self.df['mouse_id'] == mouse_id].sort_values('day')
            
            first_day = mouse_data.iloc[0]
            last_day = mouse_data.iloc[-1]
            
            first_last.append({
                'mouse_id': mouse_id,
                'group': first_day['group'],
                'cohort': first_day['cohort'],
                'first_day': first_day['day'],
                'last_day': last_day['day'],
                'first_total_licks': first_day['total_licks'],
                'last_total_licks': last_day['total_licks'],
                'first_n_bursts': first_day['n_bursts'],
                'last_n_bursts': last_day['n_bursts'],
                'first_burst_size': first_day['mean_burst_size'],
                'last_burst_size': last_day['mean_burst_size']
            })
        
        comparison_df = pd.DataFrame(first_last)
        
        # Calculate changes
        comparison_df['change_licks'] = comparison_df['last_total_licks'] - comparison_df['first_total_licks']
        comparison_df['change_bursts'] = comparison_df['last_n_bursts'] - comparison_df['first_n_bursts']
        comparison_df['change_burst_size'] = comparison_df['last_burst_size'] - comparison_df['first_burst_size']
        
        # Visualize changes
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        
        changes = [
            ('change_licks', 'Change in Total Licks'),
            ('change_bursts', 'Change in Number of Bursts'),
            ('change_burst_size', 'Change in Burst Size')
        ]
        
        for idx, (metric, title) in enumerate(changes):
            ax = axes[idx]
            sns.boxplot(data=comparison_df, x='group', y=metric, ax=ax, palette='Set2')
            sns.stripplot(data=comparison_df, x='group', y=metric, ax=ax,
                         color='black', alpha=0.5, size=6)
            ax.axhline(y=0, color='red', linestyle='--', alpha=0.5)
            ax.set_title(title, fontweight='bold')
            ax.set_ylabel('Change (Last - First)')
        
        plt.suptitle('First vs Last Day Comparison', fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.show()
        
        return comparison_df
    
    def export_for_stats(self, output_path=None):
        """
        Export data in format ready for statistical analysis in R/SPSS/Prism.
        
        Parameters:
        -----------
        output_path : str or Path, optional
            Where to save the file. If None, saves to results folder
        """
        if output_path is None:
            output_path = self.results_folder / 'FR1_data_for_stats.csv'
        
        # Select key columns
        stats_df = self.df[[
            'mouse_id', 'group', 'cohort', 'day',
            'total_licks', 'lick_rate_per_s',
            'n_bursts', 'mean_burst_size', 'mean_burst_duration_ms',
            'mean_intraburst_freq_hz', 'mean_ili_ms',
            'mean_lick_duration_ms', 'percent_long_licks'
        ]].copy()
        
        stats_df.to_csv(output_path, index=False)
        print(f"Exported data for statistical analysis to:\n{output_path}")


# Quick analysis functions
def quick_load(results_folder):
    """
    Quickly load results from the FR1_Analysis_Results folder.
    
    Parameters:
    -----------
    results_folder : str or Path
        Path to FR1_Analysis_Results folder
        
    Returns:
    --------
    FR1ResultsAnalyzer instance
    """
    results_folder = Path(results_folder)
    csv_file = results_folder / 'FR1_lick_microstructure_all_mice.csv'
    
    if not csv_file.exists():
        print(f"Error: Could not find {csv_file}")
        return None
    
    return FR1ResultsAnalyzer(csv_file)


# Example usage
if __name__ == "__main__":
    # CHANGE THIS PATH
    results_path = r"C:\Users\lmman\iCloudDrive\Otis Lab\DS Habit Ensure Study\DREADD Cohort 1-3 Set Up for Code\FR1_Analysis_Results\FR1_lick_microstructure_all_mice.csv"
    
    # Load the data
    analyzer = FR1ResultsAnalyzer(results_path)
    
    # Example analyses
    print("\n" + "="*60)
    print("EXAMPLE ANALYSES")
    print("="*60)
    
    # 1. Plot learning curves
    print("\n1. Plotting learning curves...")
    analyzer.plot_learning_curves()
    
    # 2. Compare groups on total licks
    print("\n2. Comparing groups on total licks...")
    analyzer.compare_groups('total_licks')
    
    # 3. Analyze a specific mouse
    print("\n3. Analyzing individual mouse...")
    first_mouse = analyzer.df['mouse_id'].iloc[0]
    analyzer.analyze_individual_mouse(first_mouse)
    
    # 4. First vs last day comparison
    print("\n4. Comparing first vs last training day...")
    comparison = analyzer.compare_first_vs_last()
    
    # 5. Export for stats
    print("\n5. Exporting data for statistical analysis...")
    analyzer.export_for_stats()
    
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE!")
    print("="*60)
