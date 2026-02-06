# Lick Microstructure Web Analyzer

A custom web application for analyzing lick microstructure in rodent behavioral experiments. Built with Plotly Dash and tailored for Otis Lab data format.

## Features

### 📊 **Comprehensive Analysis**
- **Session Overview**: Total licks, lick rate, temporal patterns
- **Burst Analysis**: Number of bursts, mean burst size, burst duration
- **Microstructure Metrics**: Intraburst frequency, interlick intervals
- **Statistical Modeling**: Weibull distribution fitting for burst probability
- **Lick Characteristics**: Duration analysis, long lick detection

### 🎨 **Interactive Visualizations**
- Real-time session histogram (standard or cumulative)
- Intraburst frequency distribution
- Burst size histogram
- Weibull probability plot
- Lick duration distribution

### ⚙️ **Adjustable Parameters**
- **IBI Threshold** (100-2000ms): Define burst boundaries
- **Minimum Burst Size** (1-10 licks): Filter small lick clusters
- **Bin Size** (10-300s): Session histogram resolution
- **Cumulative Option**: Toggle cumulative lick display

### 📥 **Data Export**
- Download complete analysis results as CSV
- Includes all metrics and statistical parameters

## Installation

### Prerequisites
- Python 3.8 or newer
- pip package manager

### Setup

1. **Install dependencies:**
```bash
pip install -r requirements_webapp.txt
```

2. **Run the application:**
```bash
python lick_microstructure_app.py
```

3. **Access the app:**
Open your web browser and navigate to:
```
http://127.0.0.1:8050
```

## Usage

### 1. Upload Data

Drag and drop or click to select your CSV file. The file must contain:
- **Required**: `On Timestamp` column (lick onset times in milliseconds)
- **Optional**: `Off Timestamp` column (lick offset times for duration analysis)

**Supported Format:**
```csv
Unnamed: 0,On Timestamp,Off Timestamp
0,92342,92542
1,92829,92916
2,92955,92994
...
```

### 2. Adjust Parameters

Use the sliders to customize your analysis:

**Interburst Interval (IBI)** 
- Default: 500ms
- Range: 100-2000ms
- Purpose: Defines the minimum gap between licks that separates bursts

**Minimum Licks per Burst**
- Default: 3 licks
- Range: 1-10 licks
- Purpose: Filters out very small lick clusters

**Histogram Bin Size**
- Default: 60s
- Range: 10-300s
- Purpose: Time resolution for session overview

### 3. View Results

The app displays:

**Statistics Panel:**
- Session metrics (total licks, duration, rate)
- Burst metrics (count, size, frequency)
- Lick characteristics (ILI, duration, long licks)

**Interactive Graphs:**
- **Session Overview**: Temporal pattern of licking
- **Intraburst Frequency**: ILI distribution within bursts
- **Burst Size Distribution**: Histogram of licks per burst
- **Weibull Fit**: Burst probability model (requires ≥10 bursts)
- **Lick Durations**: Contact time distribution

### 4. Download Results

Click "📥 Download Results (CSV)" to export:
- All calculated metrics
- Weibull parameters
- Summary statistics

## Data Format

### Your Current Format (Otis Lab)
Your data comes from operant boxes with this structure:

```
Acq_FR1/
├── MouseID_Group_Cohort_Acq/
│   ├── Acq_FR1_D01/
│   │   └── MouseID_FR1_D1_DATE_Lick Event Log.csv
│   └── Acq_FR1_D02/
│       └── MouseID_FR1_D2_DATE_Lick Event Log.csv
```

**CSV Structure:**
- Column 1: Index (can be ignored)
- Column 2: `On Timestamp` - lick onset in milliseconds from session start
- Column 3: `Off Timestamp` - lick offset in milliseconds from session start

### Preprocessing for Batch Analysis

If you want to analyze multiple files:
1. Use the `batch_lick_microstructure_FR1.py` script (provided separately)
2. Or upload individual files to this web app one at a time

## Understanding the Metrics

### Session-Level Metrics
- **Total Licks**: Count of all lick events
- **Session Duration**: Time from first to last lick
- **Lick Rate**: Licks per second (total licks / duration)

### Burst Metrics
- **Number of Bursts**: Count of lick clusters separated by IBI threshold
- **Mean Burst Size**: Average licks per burst
- **Mean Burst Duration**: Average time span of bursts (ms)
- **Intraburst Frequency**: Lick rate within bursts (typically 6-8 Hz in rodents)

### Microstructure
- **Mean ILI**: Average inter-lick interval across all licks
- **Weibull α (alpha)**: Scale parameter of burst size distribution
- **Weibull β (beta)**: Shape parameter of burst size distribution
- **Weibull r²**: Goodness of fit for Weibull model

### Lick Characteristics
- **Mean Lick Duration**: Average tongue contact time (ms)
- **Long Licks**: Licks exceeding 300ms threshold
- **% Long Licks**: Percentage of licks that are "long"

> **Note**: Long licks (>300ms) may indicate suboptimal sipper placement or fluid bridges

## Interpreting Results

### Typical Values (Rodents)
- **Intraburst Frequency**: 6-8 Hz (normal licking rhythm)
- **Lick Duration**: 50-200ms
- **Long Licks**: <10% (higher suggests sipper issues)

### Common Patterns
- **More bursts, smaller size**: May indicate lower palatability or satiation
- **Fewer bursts, larger size**: May indicate higher palatability or motivation
- **Changes in intraburst frequency**: Typically stable; changes suggest motor issues

### Weibull Analysis
- Requires ≥10 bursts for reliable fitting
- **β > 1**: Burst size increases over session (learning/motivation)
- **β < 1**: Burst size decreases over session (satiation)
- **β ≈ 1**: Exponential distribution (random burst termination)

## Comparison to lickcalc

This app is inspired by [lickcalc](https://lickcalc.uit.no) but customized for your data:

**Similarities:**
- Same core microstructure metrics
- Weibull distribution fitting
- Interactive parameter adjustment
- Real-time visualization updates

**Differences:**
- **Simplified file format**: Works directly with your CSV format
- **No Med Associates parsing**: Designed for your specific data structure
- **Streamlined interface**: Focused on essential metrics
- **Custom for operant tasks**: Tailored for FR1/VI schedule analysis

## Troubleshooting

**App won't start:**
- Check that all dependencies are installed
- Verify Python version is 3.8+
- Try: `pip install --upgrade -r requirements_webapp.txt`

**File upload fails:**
- Verify CSV has `On Timestamp` column
- Check that timestamps are numeric (milliseconds)
- Ensure no missing/corrupt data in CSV

**No bursts detected:**
- Decrease IBI threshold (try 300-400ms)
- Decrease minimum burst size (try 2 licks)
- Check that you have enough licks (need at least min_burst_size × 2)

**Weibull plot shows "Need ≥10 bursts":**
- This is normal for short sessions or very restrictive parameters
- Adjust IBI threshold and min burst size to detect more bursts
- Weibull fitting requires sufficient data for reliable statistics

**Graphs look strange:**
- Check for outliers in your data (extreme timestamps)
- Verify timestamps are in milliseconds (not seconds)
- Try different bin sizes for session histogram

## Advanced Usage

### Deploying Online

To deploy this app online (e.g., on Heroku, PythonAnywhere):

1. **Create `Procfile`:**
```
web: gunicorn lick_microstructure_app:server
```

2. **Add to requirements:**
```
gunicorn==21.2.0
```

3. **Follow your hosting platform's deployment guide**

### Customization

The app is built modularly. You can customize:

**Analysis Parameters:**
Edit default values in the sliders (lines 500-600)

**Visual Styling:**
Change Bootstrap theme in `app = dash.Dash(...)` line

**Additional Metrics:**
Add new calculations to the `LickAnalyzer` class

**Export Format:**
Modify the `download_results` callback to include additional data

## Citations

If you use this tool in your research, please cite:

**Original lickcalc:**
```
Volcko KL & McCutcheon JE. lickcalc: Easy analysis of lick microstructure 
in experiments of rodent ingestive behaviour. (2025)
https://github.com/uit-no/lickcalc
```

**Key References:**
- Davis & Smith (1992) - Foundational microstructure work
- Johnson (2018) - Comprehensive microstructure review
- Naneix et al. (2020) - Applications to ingestive behavior

## Support

For issues or questions:
1. Check the troubleshooting section above
2. Review the lickcalc documentation: [https://lickcalc.uit.no/help](https://lickcalc.uit.no/help)
3. Consult with your lab colleagues or PI

## File Structure

```
lick-microstructure-app/
├── lick_microstructure_app.py    # Main application
├── requirements_webapp.txt        # Python dependencies
└── README_webapp.md              # This file
```

## License

GPL-3.0 License (following lickcalc)

---

**Built for:** Logan's DREADD habit study at MUSC
**Lab:** Dr. James Otis' Laboratory
**Date:** February 2025
