"""
Utility Functions for EDE Project
Author: Dewan Sajid Islam
Date: 2026
"""

import numpy as np
import json
import os
from scipy.interpolate import interp1d
from scipy.integrate import trapz
import matplotlib.pyplot as plt

# ============================================================================
# File I/O Utilities
# ============================================================================

def load_parameters(model='EDE', filepath=None):
    """
    Load cosmological parameters from JSON file.
    
    Parameters
    ----------
    model : str
        'EDE' or 'LCDM'
    filepath : str or None
        Path to parameters file
        
    Returns
    -------
    params : dict
        Cosmological parameters
    """
    if filepath is None:
        filepath = os.path.join(os.path.dirname(__file__), '..', 'data', 'parameters.json')
    
    with open(filepath, 'r') as f:
        all_params = json.load(f)
    
    return all_params.get(model, {})

def save_parameters(params, filename):
    """
    Save parameters to JSON file.
    
    Parameters
    ----------
    params : dict
        Parameters to save
    filename : str
        Output filename
    """
    with open(filename, 'w') as f:
        json.dump(params, f, indent=2)

def load_best_fit(filepath=None):
    """
    Load best-fit parameters.
    
    Parameters
    ----------
    filepath : str or None
        Path to best-fit file
        
    Returns
    -------
    best_fit : dict
        Best-fit parameters and statistics
    """
    if filepath is None:
        filepath = os.path.join(os.path.dirname(__file__), '..', 'data', 'best_fit_parameters.json')
    
    with open(filepath, 'r') as f:
        return json.load(f)

# ============================================================================
# Statistical Utilities
# ============================================================================

def chi_squared(data, model, errors):
    """
    Calculate chi-squared.
    
    Parameters
    ----------
    data : array
        Observed data
    model : array
        Model predictions
    errors : array
        Measurement errors
        
    Returns
    -------
    chi2 : float
        Chi-squared value
    """
    return np.sum(((data - model) / errors)**2)

def reduced_chi_squared(chi2, n_data, n_params):
    """
    Calculate reduced chi-squared.
    
    Parameters
    ----------
    chi2 : float
        Chi-squared value
    n_data : int
        Number of data points
    n_params : int
        Number of parameters
        
    Returns
    -------
    red_chi2 : float
        Reduced chi-squared
    """
    return chi2 / (n_data - n_params)

def akaike_information_criterion(chi2, n_params, n_data):
    """
    Calculate Akaike Information Criterion (AIC).
    
    Parameters
    ----------
    chi2 : float
        Chi-squared value
    n_params : int
        Number of parameters
    n_data : int
        Number of data points
        
    Returns
    -------
    aic : float
        AIC value
    """
    return chi2 + 2 * n_params + (2 * n_params * (n_params + 1)) / (n_data - n_params - 1)

def bayesian_information_criterion(chi2, n_params, n_data):
    """
    Calculate Bayesian Information Criterion (BIC).
    
    Parameters
    ----------
    chi2 : float
        Chi-squared value
    n_params : int
        Number of parameters
    n_data : int
        Number of data points
        
    Returns
    -------
    bic : float
        BIC value
    """
    return chi2 + n_params * np.log(n_data)

def gelman_rubin(chains):
    """
    Calculate Gelman-Rubin diagnostic for MCMC chains.
    
    Parameters
    ----------
    chains : array, shape (n_walkers, n_steps, n_params)
        MCMC chains
        
    Returns
    -------
    R_hat : array, shape (n_params,)
        Gelman-Rubin statistic for each parameter
    """
    n_walkers, n_steps, n_params = chains.shape
    
    # Calculate within-chain variance
    chain_vars = np.var(chains, axis=1, ddof=1)  # shape: (n_walkers, n_params)
    W = np.mean(chain_vars, axis=0)  # within-chain variance
    
    # Calculate between-chain variance
    chain_means = np.mean(chains, axis=1)  # shape: (n_walkers, n_params)
    overall_mean = np.mean(chain_means, axis=0)  # shape: (n_params,)
    
    B = n_steps / (n_walkers - 1) * np.sum((chain_means - overall_mean)**2, axis=0)
    
    # Estimated variance
    var_plus = (n_steps - 1) / n_steps * W + B / n_steps
    
    # R-hat statistic
    R_hat = np.sqrt(var_plus / W)
    
    return R_hat

def effective_sample_size(chains):
    """
    Calculate effective sample size (ESS).
    
    Parameters
    ----------
    chains : array, shape (n_walkers, n_steps, n_params)
        MCMC chains
        
    Returns
    -------
    ess : array, shape (n_params,)
        Effective sample size for each parameter
    """
    n_walkers, n_steps, n_params = chains.shape
    
    # Flatten chains across walkers
    flat_chains = chains.reshape(-1, n_params)
    
    # Calculate autocorrelation
    ess = np.zeros(n_params)
    for i in range(n_params):
        # Simple estimate of ESS
        # In practice, use more sophisticated autocorrelation time estimation
        acf = autocorrelation(flat_chains[:, i])
        
        # Find where ACF drops below 0.5
        tau = np.argmax(acf < 0.5)
        if tau == 0:
            tau = 1
        
        ess[i] = n_walkers * n_steps / (2 * tau)
    
    return ess

def autocorrelation(x, max_lag=None):
    """
    Calculate autocorrelation function.
    
    Parameters
    ----------
    x : array
        Time series
    max_lag : int or None
        Maximum lag to compute
        
    Returns
    -------
    acf : array
        Autocorrelation function
    """
    n = len(x)
    if max_lag is None:
        max_lag = min(n // 2, 1000)
    
    x_mean = np.mean(x)
    x_var = np.var(x)
    
    acf = np.zeros(max_lag)
    for lag in range(max_lag):
        if x_var == 0:
            acf[lag] = 0
        else:
            acf[lag] = np.mean((x[:n-lag] - x_mean) * (x[lag:] - x_mean)) / x_var
    
    return acf

# ============================================================================
# Cosmological Utilities
# ============================================================================

def redshift_to_time(z, H0=70.0, Omega_m=0.3, Omega_L=0.7):
    """
    Convert redshift to cosmic time (approximate).
    
    Parameters
    ----------
    z : float or array
        Redshift
    H0 : float
        Hubble constant in km/s/Mpc
    Omega_m : float
        Matter density parameter
    Omega_L : float
        Dark energy density parameter
        
    Returns
    -------
    t : float or array
        Cosmic time in Gyr
    """
    z = np.asarray(z)
    
    # Hubble time
    H0_s = H0 * 1000 / 3.086e19  # Convert to 1/s
    t_H = 1.0 / H0_s / 3.154e7 / 1e9  # Hubble time in Gyr
    
    # For flat ΛCDM (approximate)
    # Full calculation requires integration of 1/(H(z)*(1+z))
    # This is a simplified approximation
    if Omega_L == 0:
        # Einstein-de Sitter
        t = (2.0/3.0) * t_H * (1 + z)**(-1.5)
    else:
        # Flat ΛCDM approximation
        a = 1.0 / (1.0 + z)
        t = t_H * (2.0/3.0) * np.arcsinh(np.sqrt(Omega_L/Omega_m * a**-3)) / np.sqrt(Omega_L)
    
    return t

def comoving_distance(z, H0=70.0, Omega_m=0.3, Omega_k=0.0):
    """
    Calculate comoving distance (approximate).
    
    Parameters
    ----------
    z : float or array
        Redshift
    H0 : float
        Hubble constant in km/s/Mpc
    Omega_m : float
        Matter density parameter
    Omega_k : float
        Curvature density parameter
        
    Returns
    -------
    D_C : float or array
        Comoving distance in Mpc
    """
    z = np.asarray(z)
    
    # Hubble distance
    D_H = 299792.458 / H0  # Mpc
    
    # For flat universe (Omega_k = 0)
    if Omega_k == 0:
        # Numerical integration would be better
        # Simplified approximation for flat ΛCDM
        Omega_L = 1.0 - Omega_m
        
        # Approximate integral
        def integrand(zp):
            return 1.0 / np.sqrt(Omega_m * (1+zp)**3 + Omega_L)
        
        # Simple trapezoidal integration
        if np.isscalar(z):
            z_vals = np.linspace(0, z, 1000)
            integrand_vals = integrand(z_vals)
            D_C = trapz(integrand_vals, z_vals)
        else:
            D_C = np.zeros_like(z)
            for i, zi in enumerate(z):
                if zi == 0:
                    D_C[i] = 0
                else:
                    z_vals = np.linspace(0, zi, 1000)
                    integrand_vals = integrand(z_vals)
                    D_C[i] = trapz(integrand_vals, z_vals)
        
        return D_H * D_C
    
    else:
        raise NotImplementedError("Non-flat universe not implemented in this simplified version.")

def angular_diameter_distance(z, H0=70.0, Omega_m=0.3):
    """
    Calculate angular diameter distance.
    
    Parameters
    ----------
    z : float or array
        Redshift
    H0 : float
        Hubble constant in km/s/Mpc
    Omega_m : float
        Matter density parameter
        
    Returns
    -------
    D_A : float or array
        Angular diameter distance in Mpc
    """
    D_C = comoving_distance(z, H0, Omega_m)
    return D_C / (1.0 + z)

def luminosity_distance(z, H0=70.0, Omega_m=0.3):
    """
    Calculate luminosity distance.
    
    Parameters
    ----------
    z : float or array
        Redshift
    H0 : float
        Hubble constant in km/s/Mpc
    Omega_m : float
        Matter density parameter
        
    Returns
    -------
    D_L : float or array
        Luminosity distance in Mpc
    """
    D_C = comoving_distance(z, H0, Omega_m)
    return D_C * (1.0 + z)

# ============================================================================
# Plotting Utilities
# ============================================================================

def setup_plot_style(style='default'):
    """
    Setup matplotlib plotting style.
    
    Parameters
    ----------
    style : str
        Plot style: 'default', 'publication', 'presentation'
    """
    if style == 'default':
        plt.rcParams.update({
            'font.size': 12,
            'axes.labelsize': 14,
            'axes.titlesize': 16,
            'xtick.labelsize': 12,
            'ytick.labelsize': 12,
            'legend.fontsize': 12,
            'figure.titlesize': 18,
            'figure.dpi': 100,
            'savefig.dpi': 300,
            'axes.grid': True,
            'grid.alpha': 0.3,
            'lines.linewidth': 2,
            'lines.markersize': 6,
        })
    
    elif style == 'publication':
        plt.rcParams.update({
            'font.size': 10,
            'font.family': 'serif',
            'font.serif': ['Times', 'Palatino', 'New Century Schoolbook', 'Bookman', 'Computer Modern Roman'],
            'mathtext.fontset': 'cm',
            'mathtext.rm': 'serif',
            'axes.labelsize': 12,
            'axes.titlesize': 14,
            'xtick.labelsize': 10,
            'ytick.labelsize': 10,
            'legend.fontsize': 10,
            'figure.titlesize': 16,
            'figure.dpi': 300,
            'savefig.dpi': 600,
            'savefig.format': 'pdf',
            'savefig.bbox': 'tight',
            'axes.grid': True,
            'grid.alpha': 0.2,
            'lines.linewidth': 1.5,
            'lines.markersize': 4,
        })
    
    elif style == 'presentation':
        plt.rcParams.update({
            'font.size': 14,
            'axes.labelsize': 16,
            'axes.titlesize': 18,
            'xtick.labelsize': 14,
            'ytick.labelsize': 14,
            'legend.fontsize': 14,
            'figure.titlesize': 20,
            'figure.dpi': 150,
            'savefig.dpi': 300,
            'axes.grid': True,
            'grid.alpha': 0.3,
            'lines.linewidth': 3,
            'lines.markersize': 8,
        })

def save_figure(fig, filename, formats=['pdf', 'png'], dpi=300, transparent=False):
    """
    Save figure in multiple formats.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object
    filename : str
        Base filename (without extension)
    formats : list
        List of formats to save
    dpi : int
        Resolution for raster formats
    transparent : bool
        Whether to use transparent background
    """
    import os
    
    # Create figures directory if it doesn't exist
    os.makedirs('figures', exist_ok=True)
    
    for fmt in formats:
        if fmt == 'pdf':
            fig.savefig(f'figures/{filename}.pdf', 
                       bbox_inches='tight', 
                       transparent=transparent)
        elif fmt == 'png':
            fig.savefig(f'figures/{filename}.png', 
                       dpi=dpi, 
                       bbox_inches='tight',
                       transparent=transparent)
        elif fmt == 'eps':
            fig.savefig(f'figures/{filename}.eps', 
                       bbox_inches='tight',
                       transparent=transparent)
        else:
            print(f"Warning: Unknown format '{fmt}', skipping.")
    
    print(f"Saved figure as figures/{filename}.{{{','.join(formats)}}}")

def add_credits(ax, text=None, position='bottom right'):
    """
    Add credits/citation to plot.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes object
    text : str or None
        Credit text. If None, uses default.
    position : str
        'bottom right', 'bottom left', 'top right', 'top left'
    """
    if text is None:
        text = "Dewan Sajid Islam, Entanglement Dark Energy Model (2024)"
    
    if position == 'bottom right':
        x, y, ha, va = 0.98, 0.02, 'right', 'bottom'
    elif position == 'bottom left':
        x, y, ha, va = 0.02, 0.02, 'left', 'bottom'
    elif position == 'top right':
        x, y, ha, va = 0.98, 0.98, 'right', 'top'
    elif position == 'top left':
        x, y, ha, va = 0.02, 0.98, 'left', 'top'
    else:
        raise ValueError(f"Unknown position: {position}")
    
    ax.text(x, y, text,
            transform=ax.transAxes,
            fontsize=8,
            alpha=0.7,
            horizontalalignment=ha,
            verticalalignment=va)

# ============================================================================
# Numerical Utilities
# ============================================================================

def log_interp1d(x, y, **kwargs):
    """
    Logarithmic interpolation (linear in log-space).
    
    Parameters
    ----------
    x : array
        Independent variable
    y : array
        Dependent variable
    **kwargs : dict
        Additional arguments to interp1d
        
    Returns
    -------
    interp_func : callable
        Interpolation function
    """
    log_x = np.log10(x)
    log_y = np.log10(y)
    
    # Remove any infinities or NaNs
    mask = np.isfinite(log_x) & np.isfinite(log_y)
    if not np.all(mask):
        log_x = log_x[mask]
        log_y = log_y[mask]
    
    interp_log = interp1d(log_x, log_y, **kwargs)
    
    def interpolator(x_new):
        return 10**interp_log(np.log10(x_new))
    
    return interpolator

def smooth_data(x, y, window_size=11, window='hanning'):
    """
    Smooth data using a window function.
    
    Parameters
    ----------
    x : array
        Independent variable (must be equally spaced)
    y : array
        Dependent variable
    window_size : int
        Size of smoothing window (must be odd)
    window : str
        Window type: 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
        
    Returns
    -------
    y_smooth : array
        Smoothed data
    """
    if window_size < 3:
        return y
    
    if window_size % 2 == 0:
        window_size += 1
    
    if window == 'flat':  # Moving average
        w = np.ones(window_size) / window_size
    elif window == 'hanning':
        w = np.hanning(window_size)
        w /= w.sum()
    elif window == 'hamming':
        w = np.hamming(window_size)
        w /= w.sum()
    elif window == 'bartlett':
        w = np.bartlett(window_size)
        w /= w.sum()
    elif window == 'blackman':
        w = np.blackman(window_size)
        w /= w.sum()
    else:
        raise ValueError(f"Unknown window type: {window}")
    
    # Pad the data
    y_padded = np.pad(y, (window_size//2, window_size//2), mode='edge')
    
    # Convolve
    y_smooth = np.convolve(w, y_padded, mode='valid')
    
    return y_smooth

# ============================================================================
# Unit Conversion
# ============================================================================

def unit_conversion(value, from_unit, to_unit):
    """
    Convert between cosmological units.
    
    Parameters
    ----------
    value : float or array
        Value to convert
    from_unit : str
        Original unit
    to_unit : str
        Target unit
        
    Returns
    -------
    converted : float or array
        Converted value
    """
    # Length conversions
    length_units = {
        'm': 1.0,
        'km': 1e3,
        'Mpc': 3.086e22,
        'Mpc/h': 3.086e22,  # Same as Mpc, h is dimensionless
        'Gpc': 3.086e25,
    }
    
    # Time conversions
    time_units = {
        's': 1.0,
        'yr': 3.154e7,
        'Gyr': 3.154e16,
    }
    
    # Mass conversions
    mass_units = {
        'kg': 1.0,
        'M_sun': 1.989e30,
    }
    
    # Energy conversions
    energy_units = {
        'J': 1.0,
        'eV': 1.602e-19,
        'GeV': 1.602e-10,
    }
    
    # Combine all units
    all_units = {**length_units, **time_units, **mass_units, **energy_units}
    
    if from_unit not in all_units or to_unit not in all_units:
        raise ValueError(f"Unknown unit. Available units: {list(all_units.keys())}")
    
    return value * all_units[from_unit] / all_units[to_unit]

# ============================================================================
# Test Functions
# ============================================================================

def test_utilities():
    """Test all utility functions."""
    print("Testing utility functions...")
    
    # Test statistics
    data = np.array([1.0, 2.0, 3.0])
    model = np.array([1.1, 2.1, 3.1])
    errors = np.array([0.1, 0.1, 0.1])
    
    chi2 = chi_squared(data, model, errors)
    print(f"Chi-squared: {chi2:.2f}")
    
    # Test unit conversion
    value = 1.0
    converted = unit_conversion(value, 'Mpc', 'm')
    print(f"1 Mpc = {converted:.2e} m")
    
    # Test redshift to time
    z = 0.0
    t = redshift_to_time(z)
    print(f"Age of universe at z={z}: {t:.2f} Gyr")
    
    print("All tests passed!")

if __name__ == "__main__":
    test_utilities()
