"""
Lick Microstructure Analysis Web Application
Custom version for Otis Lab data format

Based on lickcalc architecture but adapted for CSV format from operant boxes
Author: Custom build for Logan's DREADD study
"""

import dash
from dash import dcc, html, Input, Output, State, callback_context
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

import pandas as pd
import numpy as np
from scipy import stats
import base64
import io
from pathlib import Path
import json

# Initialize Dash app with Bootstrap theme
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.title = "Lick Microstructure Analyzer"

# Server for deployment
server = app.server


# ===========================
# ANALYSIS FUNCTIONS
# ===========================

class LickAnalyzer:
    """Core analysis functions for lick microstructure."""
    
    def __init__(self, lick_onsets, lick_offsets=None, ibi_threshold=500, 
                 min_burst_size=3, long_lick_threshold=300):
        """
        Initialize analyzer.
        
        Parameters:
        -----------
        lick_onsets : array
            Lick onset timestamps in milliseconds
        lick_offsets : array, optional
            Lick offset timestamps in milliseconds
        ibi_threshold : float
            Interburst interval threshold in ms
        min_burst_size : int
            Minimum licks per burst
        long_lick_threshold : float
            Threshold for long licks in ms
        """
        self.lick_onsets = np.array(lick_onsets)
        self.lick_offsets = np.array(lick_offsets) if lick_offsets is not None else None
        self.ibi_threshold = ibi_threshold
        self.min_burst_size = min_burst_size
        self.long_lick_threshold = long_lick_threshold
        
        # Calculate basic metrics
        self.total_licks = len(self.lick_onsets)
        self.session_duration_s = self.lick_onsets[-1] / 1000.0 if self.total_licks > 0 else 0
        
        # Calculate derived metrics
        self._calculate_ilis()
        self._identify_bursts()
        self._analyze_lick_durations()
        
    def _calculate_ilis(self):
        """Calculate interlick intervals."""
        if self.total_licks > 1:
            self.ilis = np.diff(self.lick_onsets)
            self.mean_ili = np.mean(self.ilis)
            self.median_ili = np.median(self.ilis)
        else:
            self.ilis = np.array([])
            self.mean_ili = np.nan
            self.median_ili = np.nan
    
    def _identify_bursts(self):
        """Identify lick bursts."""
        if self.total_licks < self.min_burst_size:
            self.bursts = []
            self.n_bursts = 0
            return
        
        # Find burst boundaries
        burst_breaks = np.where(self.ilis > self.ibi_threshold)[0] + 1
        bursts = np.split(self.lick_onsets, burst_breaks)
        
        # Filter by minimum size
        self.bursts = [b for b in bursts if len(b) >= self.min_burst_size]
        self.n_bursts = len(self.bursts)
        
        # Calculate burst metrics
        if self.n_bursts > 0:
            self.burst_sizes = np.array([len(b) for b in self.bursts])
            self.mean_burst_size = np.mean(self.burst_sizes)
            
            # Calculate intraburst frequency for each burst
            self.intraburst_freqs = []
            self.burst_durations = []
            
            for burst in self.bursts:
                if len(burst) >= 2:
                    burst_ilis = np.diff(burst)
                    mean_burst_ili = np.mean(burst_ilis)
                    freq_hz = 1000.0 / mean_burst_ili if mean_burst_ili > 0 else np.nan
                    self.intraburst_freqs.append(freq_hz)
                    
                    duration = burst[-1] - burst[0]
                    self.burst_durations.append(duration)
            
            self.mean_intraburst_freq = np.mean([f for f in self.intraburst_freqs if not np.isnan(f)])
            self.mean_burst_duration = np.mean(self.burst_durations)
        else:
            self.burst_sizes = np.array([])
            self.mean_burst_size = np.nan
            self.mean_intraburst_freq = np.nan
            self.mean_burst_duration = np.nan
    
    def _analyze_lick_durations(self):
        """Analyze lick durations if offset data available."""
        if self.lick_offsets is not None and len(self.lick_offsets) == len(self.lick_onsets):
            self.lick_durations = self.lick_offsets - self.lick_onsets
            self.mean_lick_duration = np.mean(self.lick_durations)
            self.median_lick_duration = np.median(self.lick_durations)
            
            # Long licks
            long_licks = self.lick_durations > self.long_lick_threshold
            self.n_long_licks = np.sum(long_licks)
            self.percent_long_licks = (self.n_long_licks / self.total_licks * 100) if self.total_licks > 0 else 0
        else:
            self.lick_durations = None
            self.mean_lick_duration = np.nan
            self.median_lick_duration = np.nan
            self.n_long_licks = 0
            self.percent_long_licks = 0
    
    def fit_weibull(self):
        """Fit Weibull distribution to burst sizes (if enough bursts)."""
        if self.n_bursts < 10:  # Minimum bursts for reliable fit
            return None, None, None
        
        try:
            # Calculate survival probabilities
            sorted_sizes = np.sort(self.burst_sizes)
            n = len(sorted_sizes)
            survival_prob = 1 - (np.arange(1, n + 1) - 0.5) / n
            
            # Log-log transform for linear fit
            log_sizes = np.log(sorted_sizes)
            log_log_surv = np.log(-np.log(survival_prob))
            
            # Linear regression
            valid = np.isfinite(log_log_surv)
            if np.sum(valid) < 3:
                return None, None, None
            
            slope, intercept, r_value, _, _ = stats.linregress(
                log_sizes[valid], log_log_surv[valid]
            )
            
            # Weibull parameters
            beta = slope  # Shape parameter
            alpha = np.exp(-intercept / slope)  # Scale parameter
            r_squared = r_value ** 2
            
            return alpha, beta, r_squared
        except:
            return None, None, None
    
    def get_session_histogram(self, bin_size_s=60, cumulative=False):
        """
        Get session histogram data.
        
        Parameters:
        -----------
        bin_size_s : float
            Bin size in seconds
        cumulative : bool
            Whether to return cumulative counts
        """
        if self.total_licks == 0:
            return np.array([]), np.array([])
        
        # Convert to seconds
        lick_times_s = self.lick_onsets / 1000.0
        
        # Create bins
        max_time = np.ceil(self.session_duration_s)
        bins = np.arange(0, max_time + bin_size_s, bin_size_s)
        
        # Histogram
        counts, _ = np.histogram(lick_times_s, bins=bins)
        
        if cumulative:
            counts = np.cumsum(counts)
        
        # Bin centers for plotting
        bin_centers = bins[:-1] + bin_size_s / 2
        
        return bin_centers, counts


# ===========================
# FILE PARSING
# ===========================

def parse_csv_file(contents, filename):
    """Parse uploaded CSV file."""
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    
    try:
        # Read CSV
        df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
        
        # Check for required columns
        if 'On Timestamp' not in df.columns:
            return None, "Error: 'On Timestamp' column not found"
        
        lick_onsets = df['On Timestamp'].values
        
        # Check for offset column
        lick_offsets = df['Off Timestamp'].values if 'Off Timestamp' in df.columns else None
        
        return (lick_onsets, lick_offsets, filename), None
    
    except Exception as e:
        return None, f"Error parsing file: {str(e)}"


# ===========================
# LAYOUT COMPONENTS
# ===========================

def create_header():
    """Create app header."""
    return dbc.Row([
        dbc.Col([
            html.H1("🔬 Lick Microstructure Analyzer", className="text-primary mb-2"),
            html.P("Comprehensive analysis of licking behavior in rodent experiments", 
                   className="text-muted"),
            html.Hr()
        ])
    ])


def create_upload_section():
    """Create file upload section."""
    return dbc.Card([
        dbc.CardHeader(html.H4("📁 Upload Data", className="mb-0")),
        dbc.CardBody([
            dcc.Upload(
                id='upload-data',
                children=html.Div([
                    html.I(className="fas fa-cloud-upload-alt fa-3x mb-3"),
                    html.Br(),
                    'Drag and Drop or ',
                    html.A('Select File', style={'color': '#007bff', 'cursor': 'pointer'}),
                    html.Br(),
                    html.Small('Accepts CSV files with "On Timestamp" and "Off Timestamp" columns',
                              className="text-muted")
                ]),
                style={
                    'width': '100%',
                    'height': '150px',
                    'lineHeight': '150px',
                    'borderWidth': '2px',
                    'borderStyle': 'dashed',
                    'borderRadius': '10px',
                    'textAlign': 'center',
                    'backgroundColor': '#f8f9fa'
                },
                multiple=False
            ),
            html.Div(id='upload-status', className="mt-3")
        ])
    ], className="mb-4")


def create_parameter_controls():
    """Create parameter control sliders."""
    return dbc.Card([
        dbc.CardHeader(html.H4("⚙️ Analysis Parameters", className="mb-0")),
        dbc.CardBody([
            # IBI Threshold
            html.Label("Interburst Interval (IBI) Threshold", className="fw-bold"),
            html.P("Minimum gap between licks defining burst boundaries", 
                   className="text-muted small"),
            dcc.Slider(
                id='ibi-slider',
                min=100,
                max=2000,
                step=50,
                value=500,
                marks={i: f'{i}ms' for i in range(100, 2001, 300)},
                tooltip={"placement": "bottom", "always_visible": True}
            ),
            html.Hr(),
            
            # Min Burst Size
            html.Label("Minimum Licks per Burst", className="fw-bold mt-3"),
            html.P("Minimum number of licks to constitute a burst", 
                   className="text-muted small"),
            dcc.Slider(
                id='min-burst-slider',
                min=1,
                max=10,
                step=1,
                value=3,
                marks={i: str(i) for i in range(1, 11)},
                tooltip={"placement": "bottom", "always_visible": True}
            ),
            html.Hr(),
            
            # Bin Size
            html.Label("Session Histogram Bin Size", className="fw-bold mt-3"),
            html.P("Time bin size for session overview", 
                   className="text-muted small"),
            dcc.Slider(
                id='bin-size-slider',
                min=10,
                max=300,
                step=10,
                value=60,
                marks={i: f'{i}s' for i in range(0, 301, 60)},
                tooltip={"placement": "bottom", "always_visible": True}
            ),
            html.Hr(),
            
            # Cumulative toggle
            dbc.Checklist(
                options=[{"label": " Show cumulative licks", "value": 1}],
                value=[],
                id="cumulative-toggle",
                className="mt-3"
            )
        ])
    ], className="mb-4")


def create_results_section():
    """Create results display section."""
    return dbc.Card([
        dbc.CardHeader(html.H4("📊 Analysis Results", className="mb-0")),
        dbc.CardBody([
            dbc.Row([
                dbc.Col([
                    html.H5("Session Overview", className="text-primary"),
                    html.Div(id='session-stats')
                ], width=4),
                dbc.Col([
                    html.H5("Burst Metrics", className="text-primary"),
                    html.Div(id='burst-stats')
                ], width=4),
                dbc.Col([
                    html.H5("Lick Characteristics", className="text-primary"),
                    html.Div(id='lick-stats')
                ], width=4)
            ])
        ])
    ], className="mb-4")


# ===========================
# APP LAYOUT
# ===========================

app.layout = dbc.Container([
    # Store for data
    dcc.Store(id='stored-data'),
    
    # Header
    create_header(),
    
    # Upload section
    create_upload_section(),
    
    # Parameter controls
    create_parameter_controls(),
    
    # Results
    create_results_section(),
    
    # Graphs
    dbc.Row([
        dbc.Col([
            dbc.Card([
                dbc.CardHeader(html.H5("Session Overview")),
                dbc.CardBody([
                    dcc.Graph(id='session-histogram')
                ])
            ])
        ], width=12)
    ], className="mb-4"),
    
    dbc.Row([
        dbc.Col([
            dbc.Card([
                dbc.CardHeader(html.H5("Intraburst Frequency")),
                dbc.CardBody([
                    dcc.Graph(id='ili-histogram')
                ])
            ])
        ], width=6),
        dbc.Col([
            dbc.Card([
                dbc.CardHeader(html.H5("Burst Size Distribution")),
                dbc.CardBody([
                    dcc.Graph(id='burst-histogram')
                ])
            ])
        ], width=6)
    ], className="mb-4"),
    
    dbc.Row([
        dbc.Col([
            dbc.Card([
                dbc.CardHeader(html.H5("Burst Probability (Weibull Fit)")),
                dbc.CardBody([
                    dcc.Graph(id='weibull-plot')
                ])
            ])
        ], width=6),
        dbc.Col([
            dbc.Card([
                dbc.CardHeader(html.H5("Lick Durations")),
                dbc.CardBody([
                    dcc.Graph(id='duration-histogram')
                ])
            ])
        ], width=6)
    ], className="mb-4"),
    
    # Download button
    dbc.Row([
        dbc.Col([
            dbc.Button("📥 Download Results (CSV)", id="download-button", color="primary", size="lg"),
            dcc.Download(id="download-data")
        ], className="text-center")
    ])
    
], fluid=True, className="p-4")


# ===========================
# CALLBACKS
# ===========================

@app.callback(
    [Output('stored-data', 'data'),
     Output('upload-status', 'children')],
    Input('upload-data', 'contents'),
    State('upload-data', 'filename')
)
def upload_file(contents, filename):
    """Handle file upload."""
    if contents is None:
        raise PreventUpdate
    
    result, error = parse_csv_file(contents, filename)
    
    if error:
        return None, dbc.Alert(error, color="danger")
    
    lick_onsets, lick_offsets, fname = result
    
    # Store data
    data = {
        'lick_onsets': lick_onsets.tolist(),
        'lick_offsets': lick_offsets.tolist() if lick_offsets is not None else None,
        'filename': fname
    }
    
    status = dbc.Alert([
        html.I(className="fas fa-check-circle me-2"),
        f"Successfully loaded: {fname}",
        html.Br(),
        html.Small(f"{len(lick_onsets)} licks detected")
    ], color="success")
    
    return data, status


@app.callback(
    [Output('session-stats', 'children'),
     Output('burst-stats', 'children'),
     Output('lick-stats', 'children'),
     Output('session-histogram', 'figure'),
     Output('ili-histogram', 'figure'),
     Output('burst-histogram', 'figure'),
     Output('weibull-plot', 'figure'),
     Output('duration-histogram', 'figure')],
    [Input('stored-data', 'data'),
     Input('ibi-slider', 'value'),
     Input('min-burst-slider', 'value'),
     Input('bin-size-slider', 'value'),
     Input('cumulative-toggle', 'value')]
)
def update_analysis(stored_data, ibi_threshold, min_burst_size, bin_size, cumulative):
    """Update all analyses and plots."""
    if stored_data is None:
        # Return empty figures
        empty_fig = go.Figure()
        empty_fig.update_layout(
            title="Upload data to begin analysis",
            xaxis_title="",
            yaxis_title="",
            template="plotly_white"
        )
        return (
            html.P("Upload data to see results"), 
            html.P("Upload data to see results"),
            html.P("Upload data to see results"),
            empty_fig, empty_fig, empty_fig, empty_fig, empty_fig
        )
    
    # Load data
    lick_onsets = np.array(stored_data['lick_onsets'])
    lick_offsets = np.array(stored_data['lick_offsets']) if stored_data['lick_offsets'] else None
    
    # Analyze
    analyzer = LickAnalyzer(
        lick_onsets, lick_offsets,
        ibi_threshold=ibi_threshold,
        min_burst_size=min_burst_size
    )
    
    # Session stats
    session_stats = html.Div([
        html.P([html.Strong("Total Licks: "), f"{analyzer.total_licks}"]),
        html.P([html.Strong("Session Duration: "), f"{analyzer.session_duration_s:.1f} s"]),
        html.P([html.Strong("Lick Rate: "), f"{analyzer.total_licks/analyzer.session_duration_s:.2f} licks/s"])
    ])
    
    # Burst stats
    alpha, beta, r2 = analyzer.fit_weibull()
    weibull_str = f"α={alpha:.2f}, β={beta:.2f}, r²={r2:.3f}" if alpha else "N/A (need >10 bursts)"
    
    burst_stats = html.Div([
        html.P([html.Strong("Number of Bursts: "), f"{analyzer.n_bursts}"]),
        html.P([html.Strong("Mean Burst Size: "), f"{analyzer.mean_burst_size:.1f}" if not np.isnan(analyzer.mean_burst_size) else "N/A"]),
        html.P([html.Strong("Intraburst Freq: "), f"{analyzer.mean_intraburst_freq:.2f} Hz" if not np.isnan(analyzer.mean_intraburst_freq) else "N/A"]),
        html.P([html.Strong("Weibull: "), weibull_str], className="small")
    ])
    
    # Lick stats
    lick_stats = html.Div([
        html.P([html.Strong("Mean ILI: "), f"{analyzer.mean_ili:.1f} ms" if not np.isnan(analyzer.mean_ili) else "N/A"]),
        html.P([html.Strong("Mean Duration: "), f"{analyzer.mean_lick_duration:.1f} ms" if not np.isnan(analyzer.mean_lick_duration) else "N/A"]),
        html.P([html.Strong("Long Licks: "), f"{analyzer.n_long_licks} ({analyzer.percent_long_licks:.1f}%)"])
    ])
    
    # === FIGURES ===
    
    # 1. Session histogram
    bin_centers, counts = analyzer.get_session_histogram(bin_size_s=bin_size, cumulative=(1 in cumulative))
    fig_session = go.Figure()
    fig_session.add_trace(go.Bar(
        x=bin_centers,
        y=counts,
        marker_color='steelblue'
    ))
    fig_session.update_layout(
        title="Licks Over Time" + (" (Cumulative)" if 1 in cumulative else ""),
        xaxis_title="Time (s)",
        yaxis_title="Cumulative Licks" if 1 in cumulative else "Licks per Bin",
        template="plotly_white",
        hovermode='x unified'
    )
    
    # 2. ILI histogram (intraburst only)
    fig_ili = go.Figure()
    if analyzer.n_bursts > 0:
        # Get ILIs within bursts only
        intraburst_ilis = []
        for burst in analyzer.bursts:
            if len(burst) >= 2:
                intraburst_ilis.extend(np.diff(burst))
        
        if intraburst_ilis:
            fig_ili.add_trace(go.Histogram(
                x=intraburst_ilis,
                nbinsx=50,
                marker_color='coral'
            ))
            fig_ili.update_layout(
                title=f"Intraburst ILI Distribution (Mean: {analyzer.mean_intraburst_freq:.1f} Hz)",
                xaxis_title="Interlick Interval (ms)",
                yaxis_title="Count",
                template="plotly_white"
            )
    else:
        fig_ili.update_layout(title="No bursts detected", template="plotly_white")
    
    # 3. Burst size histogram
    fig_burst = go.Figure()
    if analyzer.n_bursts > 0:
        fig_burst.add_trace(go.Histogram(
            x=analyzer.burst_sizes,
            nbinsx=min(30, int(np.max(analyzer.burst_sizes))),
            marker_color='mediumseagreen'
        ))
        fig_burst.update_layout(
            title=f"Burst Size Distribution (Mean: {analyzer.mean_burst_size:.1f} licks)",
            xaxis_title="Licks per Burst",
            yaxis_title="Count",
            template="plotly_white"
        )
    else:
        fig_burst.update_layout(title="No bursts detected", template="plotly_white")
    
    # 4. Weibull plot
    fig_weibull = go.Figure()
    if alpha is not None:
        sorted_sizes = np.sort(analyzer.burst_sizes)
        n = len(sorted_sizes)
        survival_prob = 1 - (np.arange(1, n + 1) - 0.5) / n
        
        fig_weibull.add_trace(go.Scatter(
            x=sorted_sizes,
            y=survival_prob,
            mode='markers',
            name='Data',
            marker=dict(color='steelblue', size=6)
        ))
        
        # Fitted line
        x_fit = np.linspace(min(sorted_sizes), max(sorted_sizes), 100)
        y_fit = np.exp(-(x_fit / alpha) ** beta)
        fig_weibull.add_trace(go.Scatter(
            x=x_fit,
            y=y_fit,
            mode='lines',
            name=f'Weibull Fit (r²={r2:.3f})',
            line=dict(color='red', width=2)
        ))
        
        fig_weibull.update_layout(
            title=f"Burst Probability (Weibull: α={alpha:.2f}, β={beta:.2f})",
            xaxis_title="Burst Size (licks)",
            yaxis_title="P(Burst ≥ size)",
            template="plotly_white",
            yaxis_type="log"
        )
    else:
        fig_weibull.update_layout(
            title="Need ≥10 bursts for Weibull fit",
            template="plotly_white"
        )
    
    # 5. Duration histogram
    fig_duration = go.Figure()
    if analyzer.lick_durations is not None:
        fig_duration.add_trace(go.Histogram(
            x=analyzer.lick_durations,
            nbinsx=50,
            marker_color='mediumpurple'
        ))
        fig_duration.add_vline(
            x=300, line_dash="dash", line_color="red",
            annotation_text="Long lick threshold"
        )
        fig_duration.update_layout(
            title=f"Lick Duration Distribution (Mean: {analyzer.mean_lick_duration:.1f} ms)",
            xaxis_title="Duration (ms)",
            yaxis_title="Count",
            template="plotly_white"
        )
    else:
        fig_duration.update_layout(
            title="No offset data available",
            template="plotly_white"
        )
    
    return (session_stats, burst_stats, lick_stats,
            fig_session, fig_ili, fig_burst, fig_weibull, fig_duration)


@app.callback(
    Output("download-data", "data"),
    Input("download-button", "n_clicks"),
    [State('stored-data', 'data'),
     State('ibi-slider', 'value'),
     State('min-burst-slider', 'value')],
    prevent_initial_call=True
)
def download_results(n_clicks, stored_data, ibi_threshold, min_burst_size):
    """Download results as CSV."""
    if stored_data is None:
        raise PreventUpdate
    
    # Load data
    lick_onsets = np.array(stored_data['lick_onsets'])
    lick_offsets = np.array(stored_data['lick_offsets']) if stored_data['lick_offsets'] else None
    
    # Analyze
    analyzer = LickAnalyzer(
        lick_onsets, lick_offsets,
        ibi_threshold=ibi_threshold,
        min_burst_size=min_burst_size
    )
    
    # Prepare results
    alpha, beta, r2 = analyzer.fit_weibull()
    
    results = {
        'Metric': [
            'Total Licks',
            'Session Duration (s)',
            'Lick Rate (licks/s)',
            'Number of Bursts',
            'Mean Burst Size',
            'Mean Burst Duration (ms)',
            'Mean Intraburst Frequency (Hz)',
            'Mean ILI (ms)',
            'Mean Lick Duration (ms)',
            'Long Licks (n)',
            'Long Licks (%)',
            'Weibull Alpha',
            'Weibull Beta',
            'Weibull R²'
        ],
        'Value': [
            analyzer.total_licks,
            f"{analyzer.session_duration_s:.2f}",
            f"{analyzer.total_licks/analyzer.session_duration_s:.2f}",
            analyzer.n_bursts,
            f"{analyzer.mean_burst_size:.2f}" if not np.isnan(analyzer.mean_burst_size) else "N/A",
            f"{analyzer.mean_burst_duration:.2f}" if not np.isnan(analyzer.mean_burst_duration) else "N/A",
            f"{analyzer.mean_intraburst_freq:.2f}" if not np.isnan(analyzer.mean_intraburst_freq) else "N/A",
            f"{analyzer.mean_ili:.2f}" if not np.isnan(analyzer.mean_ili) else "N/A",
            f"{analyzer.mean_lick_duration:.2f}" if not np.isnan(analyzer.mean_lick_duration) else "N/A",
            analyzer.n_long_licks,
            f"{analyzer.percent_long_licks:.2f}",
            f"{alpha:.3f}" if alpha else "N/A",
            f"{beta:.3f}" if beta else "N/A",
            f"{r2:.3f}" if r2 else "N/A"
        ]
    }
    
    df = pd.DataFrame(results)
    
    return dcc.send_data_frame(df.to_csv, f"lick_analysis_{stored_data['filename'].replace('.csv', '')}_results.csv", index=False)


# ===========================
# RUN APP
# ===========================

if __name__ == '__main__':
    app.run(debug=True, host='127.0.0.1', port=8050)
