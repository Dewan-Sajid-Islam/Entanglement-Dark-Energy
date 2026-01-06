"""
Master script to generate all figures for the paper.
Author: Dewan Sajid Islam
Date: 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os
import sys

# Add parent directory to path to import modules
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.background_evolution import EDECosmology
from scripts.growth_functions import GrowthCalculator
from scripts.power_spectrum import PowerSpectrum
from scripts.cmb_spectra import CMBSpectra
from scripts.utils import setup_plot_style, save_figure, add_credits

def figure1():
    """Generate Figure 1: Hubble parameter and equation of state."""
    print("Generating Figure 1...")
    
    setup_plot_style('publication')
    
    # Initialize cosmologies
    ede = EDECosmology(model_name='EDE')
    lcdm = EDECosmology(model_name='LCDM')
    
    # Solve background
    ede_bg = ede.solve_background(z_max=5.0)
    lcdm_bg = lcdm.solve_background(z_max=5.0)
    
    # Create figure
    fig = plt.figure(figsize=(14, 10))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3)
    
    # Panel A: Hubble parameter
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(ede_bg['z'], ede_bg['H'], 'r-', linewidth=2, label='EDE')
    ax1.plot(lcdm_bg['z'], lcdm_bg['H'], 'b--', linewidth=2, label=r'$\Lambda$CDM')
    
    # Add observational data (mock)
    z_data = np.array([0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6])
    H_data = np.array([72.1, 69.0, 68.5, 77.5, 85.0, 92.0, 99.0, 105.0, 111.0, 117.0])
    H_err = np.array([1.0, 2.0, 2.5, 5.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
    
    ax1.errorbar(z_data, H_data, yerr=H_err, fmt='o', 
                color='green', markersize=5, capsize=3,
                alpha=0.7, label='Observations')
    
    ax1.set_xlabel('Redshift $z$')
    ax1.set_ylabel(r'$H(z)$ [km s$^{-1}$ Mpc$^{-1}$]')
    ax1.set_xlim(0, 2)
    ax1.set_ylim(50, 200)
    ax1.legend(loc='lower right')
    ax1.grid(True, alpha=0.3)
    
    # Panel B: Hubble parameter ratio
    ax2 = fig.add_subplot(gs[0, 1])
    
    # Interpolate to common redshifts
    z_common = np.linspace(0, 2, 200)
    H_ede_interp = np.interp(z_common, ede_bg['z'], ede_bg['H'])
    H_lcdm_interp = np.interp(z_common, lcdm_bg['z'], lcdm_bg['H'])
    ratio = H_ede_interp / H_lcdm_interp
    
    ax2.plot(z_common, ratio, 'k-', linewidth=2)
    ax2.axhline(1.0, color='gray', linestyle='--', alpha=0.5)
    ax2.fill_between(z_common, 0.99, 1.01, alpha=0.2, color='gray')
    
    ax2.set_xlabel('Redshift $z$')
    ax2.set_ylabel(r'$H_{\rm EDE}(z)/H_{\Lambda\rm CDM}(z)$')
    ax2.set_xlim(0, 2)
    ax2.set_ylim(0.95, 
