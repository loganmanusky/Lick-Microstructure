# FR1 Lick Microstructure Analysis Guide

Analysis pipeline for lick microstructure data from FR1 operant training sessions.

## Overview

This pipeline processes lick event logs from multiple mice across multiple training days and extracts key microstructure metrics including:
- Total licks and lick rate
- Burst analysis (number, size, duration)
- Intraburst lick frequency
- Interlick intervals
- Lick duration analysis

## File Structure Expected

```
Acq_FR1/
├── Ensure1_GIG_C1_Acq/
│   ├── Acq_FR1_D01/
│   │   └── Ensure1_FR1_D1_2024-07-22_Lick Event Log.csv
│   ├── Acq_FR1_D02/
│   │   └── Ensure1_FR1_D2_2024-07-23_Lick Event Log.csv
│   └── ...
├── Ensure2_YFPG_C1_Acq/
│   └── ...
└── ...
```

## Step 1: Batch Process All Data

Run `batch_lick_microstructure_FR1.py` to process all mice and all days:

```python
python batch_lick_microstructure_FR1.py
```

**Before running:**
1. Open the script
2. Update line 545 with your path:
   ```python
   base_folder = r"YOUR_PATH_HERE\Acq_FR1"
   ```

**What it does:**
- Finds all mouse folders in Acq_FR1
- Processes each training day for each mouse
- Extracts lick microstructure metrics
- Saves results to `FR1_Analysis_Results/` folder

**Output files:**
- `FR1_lick_microstructure_all_mice.csv` - All data in one file
- `FR1_summary_by_group_day.csv` - Grouped by experimental group and day
- `FR1_summary_by_group_overall.csv` - Overall group statistics
- `FR1_summary_by_mouse.csv` - Individual mouse summaries
- `visualizations/` folder with plots

## Step 2: Interactive Analysis

Use `analyze_FR1_results.py` to explore the processed data:

```python
from analyze_FR1_results import FR1ResultsAnalyzer

# Load your results
analyzer = FR1ResultsAnalyzer('FR1_Analysis_Results/FR1_lick_microstructure_all_mice.csv')

# Plot learning curves
analyzer.plot_learning_curves()

# Compare groups on a specific metric
analyzer.compare_groups('total_licks')
analyzer.compare_groups('mean_burst_size', day=5)  # Specific day

# Analyze individual mouse
analyzer.analyze_individual_mouse('Ensure1')

# Compare first vs last day
comparison = analyzer.compare_first_vs_last()

# Get last day performance
last_day = analyzer.get_last_day_analysis()

# Export for statistical software
analyzer.export_for_stats()
```

## Key Metrics Explained

### Session-Level Metrics:
- **total_licks**: Total number of licks in the session
- **session_duration_s**: Length of session in seconds
- **lick_rate_per_s**: Average licks per second

### Burst Metrics:
- **n_bursts**: Number of burst episodes
- **mean_burst_size**: Average licks per burst
- **mean_burst_duration_ms**: Average duration of bursts
- **mean_intraburst_freq_hz**: Lick frequency within bursts (should be ~6-8 Hz)

### Lick Characteristics:
- **mean_ili_ms**: Mean interlick interval in milliseconds
- **mean_lick_duration_ms**: Average tongue contact duration
- **percent_long_licks**: % of licks exceeding 300ms (may indicate sipper issues)

## Analysis Parameters

You can adjust these in the scripts:

- **ibi_threshold_ms** (default: 500): Gap between licks defining burst boundaries
- **min_licks_per_burst** (default: 3): Minimum licks to count as a burst
- **long_lick_threshold_ms** (default: 300): Duration threshold for "long licks"

## Experimental Groups

Your study has 4 groups:
- **GIG**: Gi-DREADD, Goal-directed condition
- **YFPG**: YFP control, Goal-directed condition  
- **GIH**: Gi-DREADD, Habit condition
- **YFPH**: YFP control, Habit condition

Across 3 cohorts (C1, C2, C3)

## Common Analyses

### 1. Did groups acquire FR1 similarly?
```python
analyzer.plot_learning_curves(['total_licks', 'n_bursts'])
analyzer.compare_groups('total_licks', day=None)  # All days
```

### 2. Were there group differences on last day?
```python
last_day = analyzer.get_last_day_analysis()
analyzer.compare_groups('mean_burst_size', day=last_day['day'].max())
```

### 3. Did groups show similar learning trajectories?
```python
comparison = analyzer.compare_first_vs_last()
print(comparison.groupby('group')[['change_licks', 'change_bursts']].describe())
```

### 4. Individual mouse check
```python
# Check if any mice are outliers
for mouse in analyzer.df['mouse_id'].unique():
    analyzer.analyze_individual_mouse(mouse)
```

## Statistical Analysis

Export data for your stats software:

```python
analyzer.export_for_stats()
```

This creates `FR1_data_for_stats.csv` with clean format for:
- SPSS
- GraphPad Prism
- R
- Python (statsmodels, pingouin)

**Recommended statistical approach:**
- Mixed effects model with mouse as random effect
- Fixed effects: Group, Day, Group×Day interaction
- Can also analyze specific days with ANOVA/t-tests

## Troubleshooting

**"No lick CSV found"**
- Check that files are named correctly with "Lick Event Log" or "Lick_Event_Log"
- Verify folder structure matches expected format

**Missing days**
- This is normal! Mice may have different numbers of training days
- The script handles this automatically

**Plots look strange**
- Check for outlier mice with `analyze_individual_mouse()`
- Verify group assignments parsed correctly from folder names
- Consider adjusting IBI threshold if burst detection seems off

## Next Steps

After FR1 analysis, you can adapt these scripts for:
1. VI15 data
2. VI30 data
3. Combined analysis across training stages

The same basic structure will work - just change the folder paths!

## Questions?

Key things to check:
1. Are intraburst frequencies in the 6-8 Hz range? (indicates normal licking motor pattern)
2. Are there many long licks? (>10% suggests sipper placement issues)
3. Do groups show similar acquisition? (no baseline differences before experimental manipulation)

---

**Author:** Analysis pipeline for Logan's DREADD habit study
**Date:** February 2025
