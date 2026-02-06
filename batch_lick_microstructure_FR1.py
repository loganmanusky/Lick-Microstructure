"""
Batch Lick Microstructure Analysis for FR1 Data
Processes all mice across all days for FR1 training

Folder Structure Expected:
- Acq_FR1/
  - Mouse1_GROUP_C#_Acq/
    - Acq_FR1_D01/
      - Mouse1_FR1_D1_DATE_Lick Event Log.csv
      - Mouse1_FR1_D1_DATE_Lever Event Log.csv (if available)
    - Acq_FR1_D02/
      ...
  - Mouse2_GROUP_C#_Acq/
    ...

Author: Analysis for Logan's DREADD study
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Set plotting style
sns.set_style("whitegrid")
sns.set_context("notebook", font_scale=1.1)


class LickMicrostructureAnalyzer:
    """Analyzes lick microstructure from CSV files."""
    
    def __init__(self, ibi_threshold_ms=500, min_licks_per_burst=3, 
                 long_lick_threshold_ms=300):
        """
        Initialize analyzer with parameters.
        
        Parameters:
        -----------
        ibi_threshold_ms : float
            Interburst interval threshold in milliseconds (default: 500ms)
        min_licks_per_burst : int
            Minimum licks to constitute a burst (default: 3)
        long_lick_threshold_ms : float
            Threshold for identifying long licks (default: 300ms)
        """
        self.ibi_threshold = ibi_threshold_ms
        self.min_licks_per_burst = min_licks_per_burst
        self.long_lick_threshold = long_lick_threshold_ms
    
    def analyze_session(self, lick_csv_path, lever_csv_path=None):
        """
        Analyze a single session's lick microstructure.
        
        Parameters:
        -----------
        lick_csv_path : Path
            Path to lick event log CSV
        lever_csv_path : Path, optional
            Path to lever event log CSV (for FR1, may not be needed)
            
        Returns:
        --------
        dict with session metrics
        """
        # Load lick data
        lick_data = pd.read_csv(lick_csv_path)
        lick_onsets = lick_data['On Timestamp'].values
        lick_offsets = lick_data['Off Timestamp'].values
        
        # Convert to seconds for easier interpretation
        lick_onsets_s = lick_onsets / 1000.0
        lick_offsets_s = lick_offsets / 1000.0
        
        # Basic session info
        total_licks = len(lick_onsets)
        session_duration_s = lick_onsets_s[-1] if total_licks > 0 else 0
        
        # Calculate lick durations
        lick_durations = lick_offsets - lick_onsets
        mean_lick_duration = np.mean(lick_durations)
        median_lick_duration = np.median(lick_durations)
        
        # Identify long licks
        long_licks = lick_durations > self.long_lick_threshold
        n_long_licks = np.sum(long_licks)
        
        # Calculate interlick intervals (ILIs)
        if total_licks > 1:
            ilis = np.diff(lick_onsets)  # in milliseconds
            mean_ili = np.mean(ilis)
            median_ili = np.median(ilis)
        else:
            ilis = np.array([])
            mean_ili = np.nan
            median_ili = np.nan
        
        # Identify bursts
        bursts = self._identify_bursts(lick_onsets)
        n_bursts = len(bursts)
        
        # Analyze bursts
        if n_bursts > 0:
            burst_sizes = [len(b) for b in bursts]
            mean_burst_size = np.mean(burst_sizes)
            median_burst_size = np.median(burst_sizes)
            
            # Calculate intraburst frequency for each burst
            intraburst_freqs = []
            for burst in bursts:
                if len(burst) >= 2:
                    burst_ilis = np.diff(burst)
                    mean_burst_ili = np.mean(burst_ilis)
                    # Convert to Hz: 1000ms / mean_ILI_ms
                    freq_hz = 1000.0 / mean_burst_ili if mean_burst_ili > 0 else np.nan
                    intraburst_freqs.append(freq_hz)
            
            mean_intraburst_freq = np.mean(intraburst_freqs) if intraburst_freqs else np.nan
            
            # Calculate burst durations
            burst_durations = []
            for burst in bursts:
                if len(burst) >= 2:
                    duration = burst[-1] - burst[0]  # in ms
                    burst_durations.append(duration)
            
            mean_burst_duration = np.mean(burst_durations) if burst_durations else np.nan
        else:
            mean_burst_size = np.nan
            median_burst_size = np.nan
            mean_intraburst_freq = np.nan
            mean_burst_duration = np.nan
            burst_sizes = []
        
        # Calculate lick rate (licks per second)
        lick_rate = total_licks / session_duration_s if session_duration_s > 0 else 0
        
        # Compile results
        results = {
            # Session info
            'total_licks': total_licks,
            'session_duration_s': session_duration_s,
            'lick_rate_per_s': lick_rate,
            
            # Lick durations
            'mean_lick_duration_ms': mean_lick_duration,
            'median_lick_duration_ms': median_lick_duration,
            'n_long_licks': n_long_licks,
            'percent_long_licks': (n_long_licks / total_licks * 100) if total_licks > 0 else 0,
            
            # Interlick intervals
            'mean_ili_ms': mean_ili,
            'median_ili_ms': median_ili,
            
            # Burst metrics
            'n_bursts': n_bursts,
            'mean_burst_size': mean_burst_size,
            'median_burst_size': median_burst_size,
            'mean_burst_duration_ms': mean_burst_duration,
            'mean_intraburst_freq_hz': mean_intraburst_freq,
            
            # Raw data for further analysis
            '_lick_onsets': lick_onsets,
            '_lick_offsets': lick_offsets,
            '_bursts': bursts,
            '_ilis': ilis,
            '_lick_durations': lick_durations
        }
        
        return results
    
    def _identify_bursts(self, lick_times):
        """
        Identify lick bursts based on IBI threshold.
        
        Returns:
        --------
        list of arrays, each containing lick times for one burst
        """
        if len(lick_times) < self.min_licks_per_burst:
            return []
        
        # Calculate interlick intervals
        ilis = np.diff(lick_times)
        
        # Find burst boundaries (where ILI exceeds threshold)
        burst_breaks = np.where(ilis > self.ibi_threshold)[0] + 1
        
        # Split licks into bursts
        bursts = np.split(lick_times, burst_breaks)
        
        # Filter bursts by minimum size
        bursts = [b for b in bursts if len(b) >= self.min_licks_per_burst]
        
        return bursts


def parse_mouse_info(folder_name):
    """
    Extract mouse ID, group, and cohort from folder name.
    
    Example: 'Ensure1_GIG_C1_Acq' -> 
        mouse_id='Ensure1', group='GIG', cohort='C1'
    """
    parts = folder_name.split('_')
    
    if len(parts) >= 3:
        mouse_id = parts[0]
        group = parts[1]
        cohort = parts[2]
    else:
        mouse_id = folder_name
        group = 'Unknown'
        cohort = 'Unknown'
    
    return mouse_id, group, cohort


def find_lick_csv(day_folder):
    """Find the lick event log CSV in a day folder."""
    lick_files = list(day_folder.glob("*Lick Event Log.csv")) + \
                 list(day_folder.glob("*Lick_Event_Log.csv"))
    
    if lick_files:
        return lick_files[0]
    return None


def batch_process_fr1(base_folder, output_folder=None):
    """
    Batch process all FR1 data.
    
    Parameters:
    -----------
    base_folder : str or Path
        Path to Acq_FR1 folder containing all mouse folders
    output_folder : str or Path, optional
        Where to save results (default: creates 'FR1_Analysis_Results' in base_folder)
        
    Returns:
    --------
    DataFrame with all results
    """
    base_folder = Path(base_folder)
    
    if output_folder is None:
        output_folder = base_folder.parent / 'FR1_Analysis_Results'
    else:
        output_folder = Path(output_folder)
    
    output_folder.mkdir(exist_ok=True, parents=True)
    
    print(f"Processing FR1 data from: {base_folder}")
    print(f"Output will be saved to: {output_folder}\n")
    
    # Initialize analyzer
    analyzer = LickMicrostructureAnalyzer()
    
    # Find all mouse folders
    mouse_folders = [f for f in base_folder.iterdir() if f.is_dir()]
    
    print(f"Found {len(mouse_folders)} mouse folders\n")
    
    # Storage for results
    all_results = []
    
    # Process each mouse
    for mouse_folder in sorted(mouse_folders):
        mouse_id, group, cohort = parse_mouse_info(mouse_folder.name)
        
        print(f"Processing: {mouse_id} (Group: {group}, Cohort: {cohort})")
        
        # Find all day folders
        day_folders = [f for f in mouse_folder.iterdir() 
                      if f.is_dir() and 'Acq_FR1_D' in f.name]
        
        print(f"  Found {len(day_folders)} training days")
        
        # Process each day
        for day_folder in sorted(day_folders):
            # Extract day number
            day_name = day_folder.name
            day_num = int(day_name.split('_D')[1][:2])
            
            # Find lick CSV
            lick_csv = find_lick_csv(day_folder)
            
            if lick_csv is None:
                print(f"    Day {day_num}: No lick CSV found")
                continue
            
            try:
                # Analyze this session
                results = analyzer.analyze_session(lick_csv)
                
                # Remove raw data arrays (too large for DataFrame)
                raw_data = {k: results.pop(k) for k in list(results.keys()) 
                           if k.startswith('_')}
                
                # Add metadata
                results['mouse_id'] = mouse_id
                results['group'] = group
                results['cohort'] = cohort
                results['day'] = day_num
                results['date'] = lick_csv.stem.split('_')[-4] if len(lick_csv.stem.split('_')) > 4 else 'Unknown'
                results['file_path'] = str(lick_csv)
                
                all_results.append(results)
                
                print(f"    Day {day_num}: ✓ ({results['total_licks']} licks, "
                      f"{results['n_bursts']} bursts)")
                
            except Exception as e:
                print(f"    Day {day_num}: ✗ Error - {str(e)}")
        
        print()  # Blank line between mice
    
    # Convert to DataFrame
    df_results = pd.DataFrame(all_results)
    
    # Save results
    output_csv = output_folder / 'FR1_lick_microstructure_all_mice.csv'
    df_results.to_csv(output_csv, index=False)
    print(f"\n✓ Saved all results to: {output_csv}")
    
    # Create summary statistics by group
    summary_stats = create_summary_statistics(df_results, output_folder)
    
    # Create visualizations
    create_visualizations(df_results, output_folder)
    
    return df_results


def create_summary_statistics(df, output_folder):
    """Create summary statistics grouped by group and day."""
    
    print("\nCreating summary statistics...")
    
    # Define metrics to summarize
    metrics = [
        'total_licks', 'lick_rate_per_s', 'mean_ili_ms', 
        'n_bursts', 'mean_burst_size', 'mean_intraburst_freq_hz',
        'mean_lick_duration_ms', 'percent_long_licks'
    ]
    
    # Group by group and day
    summary = df.groupby(['group', 'day'])[metrics].agg(['mean', 'std', 'count'])
    summary.to_csv(output_folder / 'FR1_summary_by_group_day.csv')
    
    # Overall group summary (across all days)
    group_summary = df.groupby('group')[metrics].agg(['mean', 'std', 'count'])
    group_summary.to_csv(output_folder / 'FR1_summary_by_group_overall.csv')
    
    # Individual mouse summary (across days)
    mouse_summary = df.groupby(['mouse_id', 'group'])[metrics].agg(['mean', 'std', 'count'])
    mouse_summary.to_csv(output_folder / 'FR1_summary_by_mouse.csv')
    
    print(f"  ✓ Saved summary statistics")
    
    return summary


def create_visualizations(df, output_folder):
    """Create key visualizations for the data."""
    
    print("\nCreating visualizations...")
    
    viz_folder = output_folder / 'visualizations'
    viz_folder.mkdir(exist_ok=True)
    
    # 1. Total licks across days by group
    fig, ax = plt.subplots(figsize=(10, 6))
    for group in df['group'].unique():
        group_data = df[df['group'] == group]
        group_mean = group_data.groupby('day')['total_licks'].agg(['mean', 'sem'])
        ax.plot(group_mean.index, group_mean['mean'], marker='o', label=group, linewidth=2)
        ax.fill_between(group_mean.index, 
                        group_mean['mean'] - group_mean['sem'],
                        group_mean['mean'] + group_mean['sem'],
                        alpha=0.2)
    ax.set_xlabel('Training Day', fontsize=12)
    ax.set_ylabel('Total Licks', fontsize=12)
    ax.set_title('FR1 Training: Total Licks Across Days', fontsize=14, fontweight='bold')
    ax.legend(title='Group')
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(viz_folder / 'total_licks_across_days.png', dpi=300)
    plt.close()
    
    # 2. Number of bursts across days
    fig, ax = plt.subplots(figsize=(10, 6))
    for group in df['group'].unique():
        group_data = df[df['group'] == group]
        group_mean = group_data.groupby('day')['n_bursts'].agg(['mean', 'sem'])
        ax.plot(group_mean.index, group_mean['mean'], marker='o', label=group, linewidth=2)
        ax.fill_between(group_mean.index,
                        group_mean['mean'] - group_mean['sem'],
                        group_mean['mean'] + group_mean['sem'],
                        alpha=0.2)
    ax.set_xlabel('Training Day', fontsize=12)
    ax.set_ylabel('Number of Bursts', fontsize=12)
    ax.set_title('FR1 Training: Bursts Across Days', fontsize=14, fontweight='bold')
    ax.legend(title='Group')
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(viz_folder / 'bursts_across_days.png', dpi=300)
    plt.close()
    
    # 3. Mean burst size across days
    fig, ax = plt.subplots(figsize=(10, 6))
    for group in df['group'].unique():
        group_data = df[df['group'] == group]
        group_mean = group_data.groupby('day')['mean_burst_size'].agg(['mean', 'sem'])
        ax.plot(group_mean.index, group_mean['mean'], marker='o', label=group, linewidth=2)
        ax.fill_between(group_mean.index,
                        group_mean['mean'] - group_mean['sem'],
                        group_mean['mean'] + group_mean['sem'],
                        alpha=0.2)
    ax.set_xlabel('Training Day', fontsize=12)
    ax.set_ylabel('Mean Burst Size (licks)', fontsize=12)
    ax.set_title('FR1 Training: Burst Size Across Days', fontsize=14, fontweight='bold')
    ax.legend(title='Group')
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(viz_folder / 'burst_size_across_days.png', dpi=300)
    plt.close()
    
    # 4. Intraburst frequency across days
    fig, ax = plt.subplots(figsize=(10, 6))
    for group in df['group'].unique():
        group_data = df[df['group'] == group]
        group_mean = group_data.groupby('day')['mean_intraburst_freq_hz'].agg(['mean', 'sem'])
        ax.plot(group_mean.index, group_mean['mean'], marker='o', label=group, linewidth=2)
        ax.fill_between(group_mean.index,
                        group_mean['mean'] - group_mean['sem'],
                        group_mean['mean'] + group_mean['sem'],
                        alpha=0.2)
    ax.axhline(y=6.5, color='gray', linestyle='--', alpha=0.5, label='Typical range (6-8 Hz)')
    ax.axhline(y=8, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Training Day', fontsize=12)
    ax.set_ylabel('Intraburst Frequency (Hz)', fontsize=12)
    ax.set_title('FR1 Training: Lick Frequency Within Bursts', fontsize=14, fontweight='bold')
    ax.legend(title='Group')
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(viz_folder / 'intraburst_frequency_across_days.png', dpi=300)
    plt.close()
    
    # 5. Individual mouse trajectories - Total Licks
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    
    for idx, group in enumerate(sorted(df['group'].unique())):
        ax = axes[idx]
        group_data = df[df['group'] == group]
        
        for mouse in group_data['mouse_id'].unique():
            mouse_data = group_data[group_data['mouse_id'] == mouse]
            ax.plot(mouse_data['day'], mouse_data['total_licks'], 
                   marker='o', alpha=0.7, label=mouse)
        
        ax.set_xlabel('Training Day')
        ax.set_ylabel('Total Licks')
        ax.set_title(f'Group: {group}', fontweight='bold')
        ax.legend(fontsize=8, ncol=2)
        ax.grid(alpha=0.3)
    
    plt.suptitle('Individual Mouse Trajectories - Total Licks', 
                fontsize=16, fontweight='bold', y=1.00)
    plt.tight_layout()
    plt.savefig(viz_folder / 'individual_trajectories_licks.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved visualizations to {viz_folder}")


def main():
    """
    Main function - modify this with your paths!
    """
    
    # CHANGE THIS TO YOUR PATH
    base_folder = r"C:\Users\lmman\iCloudDrive\Otis Lab\DS Habit Ensure Study\DREADD Cohort 1-3 Set Up for Code\Acq_FR1"
    
    # Run batch processing
    df_results = batch_process_fr1(base_folder)
    
    print("\n" + "="*60)
    print("PROCESSING COMPLETE!")
    print("="*60)
    print(f"\nTotal sessions processed: {len(df_results)}")
    print(f"Unique mice: {df_results['mouse_id'].nunique()}")
    print(f"Groups: {', '.join(df_results['group'].unique())}")
    print(f"Cohorts: {', '.join(df_results['cohort'].unique())}")
    print("\nCheck the 'FR1_Analysis_Results' folder for:")
    print("  - CSV files with all data")
    print("  - Summary statistics")
    print("  - Visualization plots")
    
    return df_results


if __name__ == "__main__":
    df = main()
