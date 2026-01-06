"""
CMB Power Spectra Calculations
Author: Dewan Sajid Islam
Date: 2026
"""

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import spherical_jn
from .background_evolution import EDECosmology

class CMBSpectra:
    """
    Calculate CMB power spectra for EDE model.
    Simplified analytic approximations for testing.
    """
    
    def __init__(self, cosmo, params=None):
        """
        Initialize CMB spectra calculator.
        
        Parameters
        ----------
        cosmo : EDECosmology instance
            Cosmological model
        params : dict
            Additional parameters (A_s, n_s, tau, etc.)
        """
        self.cosmo = cosmo
        
        if params is None:
            params = {}
        
        self.A_s = params.get('A_s', 2.1e-9)
        self.n_s = params.get('n_s', 0.965)
        self.tau = params.get('tau', 0.054)
        self.r = params.get('r', 0.0)  # Tensor-to-scalar ratio
        
        # Derived parameters
        self.ell_damping = 1000.0  # Damping scale
        self.ell_A = 300.0  # Acoustic scale
        
        # Pre-compute interpolation functions
        self._setup_interpolations()
    
    def _setup_interpolations(self):
        """Setup interpolation functions for efficiency."""
        # Redshift array for background quantities
        z = np.logspace(-2, 4, 1000)
        
        # Get background evolution
        bg = self.cosmo.solve_background(z_max=z[-1])
        
        # Create interpolation functions
        self.H_interp = CubicSpline(bg['z'], bg['H'])
        self.w_interp = CubicSpline(bg['z'], bg['w_ent'])
    
    def temperature_spectrum(self, ell, include_lensing=True):
        """
        Calculate CMB temperature power spectrum.
        
        Parameters
        ----------
        ell : array
            Multipole moments
        include_lensing : bool
            Include lensing effects
            
        Returns
        -------
        D_ell_TT : array
            Temperature power spectrum D_ℓ = ℓ(ℓ+1)C_ℓ/(2π) in μK²
        """
        ell = np.asarray(ell, dtype=float)
        
        # Sachs-Wolfe plateau (low ℓ)
        sw_amplitude = 1000.0  # μK²
        
        # Acoustic peaks
        peak_positions = [220, 540, 850, 1120]
        peak_amplitudes = [5500, 2500, 1500, 800]
        peak_widths = [80, 100, 120, 150]
        
        # Initialize spectrum
        D_ell = np.zeros_like(ell)
        
        # Sachs-Wolfe effect (ℓ < 100)
        mask_low = ell < 100
        D_ell[mask_low] = sw_amplitude * (2*np.pi) / (ell[mask_low]*(ell[mask_low]+1))
        
        # Add acoustic peaks
        for pos, amp, width in zip(peak_positions, peak_amplitudes, peak_widths):
            peak = amp * np.exp(-(ell - pos)**2 / (2*width**2))
            D_ell += peak
        
        # Silk damping
        damping = np.exp(-(ell/self.ell_damping)**1.5)
        D_ell *= damping
        
        # Add ISW effect (ℓ < 30) - enhanced for EDE
        if self.cosmo.model_name == 'EDE':
            mask_isw = (ell >= 10) & (ell <= 30)
            isw_enhancement = 1.0 + 0.1 * np.exp(-(ell[mask_isw] - 20)**2 / (2*5**2))
            D_ell[mask_isw] *= isw_enhancement
        
        # Add lensing smoothing
        if include_lensing:
            lensing_sigma = 5.0
            lensing_kernel = np.exp(-ell*(ell+1)*lensing_sigma**2/2)
            # Convolution would be better, but this approximation works
            D_ell_smooth = np.convolve(D_ell, np.ones(10)/10, mode='same')
            D_ell = 0.7*D_ell + 0.3*D_ell_smooth
        
        return D_ell
    
    def polarization_spectrum(self, ell, spectrum_type='EE'):
        """
        Calculate CMB polarization power spectra.
        
        Parameters
        ----------
        ell : array
            Multipole moments
        spectrum_type : str
            'EE' for E-mode, 'BB' for B-mode, 'TE' for temperature-polarization
            
        Returns
        -------
        D_ell : array
            Polarization power spectrum
        """
        ell = np.asarray(ell, dtype=float)
        
        if spectrum_type == 'EE':
            # E-mode spectrum
            # Reionization bump at low ℓ
            reion_amp = 0.5  # μK²
            reion_width = 10.0
            
            # Acoustic peaks in polarization
            peak_positions = [140, 400, 700, 1000]
            peak_amplitudes = [50, 30, 20, 10]
            peak_widths = [50, 80, 100, 120]
            
            D_ell = np.zeros_like(ell)
            
            # Reionization bump
            mask_reion = ell < 30
            D_ell[mask_reion] = reion_amp * np.exp(-(ell[mask_reion] - 10)**2 / (2*reion_width**2))
            
            # Acoustic peaks
            for pos, amp, width in zip(peak_positions, peak_amplitudes, peak_widths):
                peak = amp * np.exp(-(ell - pos)**2 / (2*width**2))
                D_ell += peak
            
            # Damping
            damping = np.exp(-(ell/1500)**2)
            D_ell *= damping
            
        elif spectrum_type == 'BB':
            # B-mode spectrum (lensing + primordial)
            # Lensing B-modes
            lensing_amp = 0.1  # μK²
            lensing_peak = 1000
            
            # Tensor B-modes (if r > 0)
            tensor_amp = 0.01 * self.r  # μK², scaled by tensor-to-scalar ratio
            tensor_peak = 80
            
            D_ell = np.zeros_like(ell)
            
            # Lensing B-modes
            D_ell += lensing_amp * np.exp(-(ell - lensing_peak)**2 / (2*300**2))
            
            # Tensor B-modes
            if self.r > 0:
                D_ell += tensor_amp * np.exp(-(ell - tensor_peak)**2 / (2*30**2))
            
        elif spectrum_type == 'TE':
            # Temperature-E-mode correlation
            # Oscillatory pattern with sign changes
            
            # Create oscillatory pattern
            x = ell / 100.0
            oscillation = 50 * np.sin(2*np.pi*x) * np.exp(-x/2)
            
            # Envelope
            envelope = np.exp(-(ell/800)**2)
            
            D_ell = oscillation * envelope
            
            # Add reionization feature
            mask_reion = ell < 30
            D_ell[mask_reion] = 10 * np.exp(-(ell[mask_reion] - 10)**2 / (2*5**2))
        
        else:
            raise ValueError(f"Unknown spectrum type: {spectrum_type}")
        
        return D_ell
    
    def lensing_potential_spectrum(self, ell):
        """
        Calculate CMB lensing potential power spectrum.
        
        Parameters
        ----------
        ell : array
            Multipole moments
            
        Returns
        -------
        C_ell_phi : array
            Lensing potential power spectrum
        """
        ell = np.asarray(ell, dtype=float)
        
        # Approximate form
        ell0 = 100.0
        amplitude = 1e-7
        
        C_ell = amplitude * (ell0 / ell)**4
        
        # Cutoff at high ℓ
        mask = ell > 2000
        C_ell[mask] = C_ell[mask] * np.exp(-(ell[mask] - 2000)/500)
        
        return C_ell
    
    def compute_angular_scale(self, scale_type='sound_horizon'):
        """
        Compute characteristic angular scales.
        
        Parameters
        ----------
        scale_type : str
            'sound_horizon' or 'equality'
            
        Returns
        -------
        theta_deg : float
            Angular scale in degrees
        ell_scale : float
            Corresponding multipole
        """
        # Sound horizon at recombination (rough estimate)
        r_s = 147.0  # Mpc
        
        # Comoving distance to last scattering
        z_star = 1089.0
        
        # Simple estimate of angular diameter distance
        H0 = self.cosmo.H0
        Omega_m = self.cosmo.Omega_m0
        
        # Comoving distance integral (simplified)
        if scale_type == 'sound_horizon':
            # Angular scale of sound horizon
            theta_rad = r_s / (3000.0 / H0 * 100)  # Very rough
            theta_deg = np.degrees(theta_rad)
            ell_scale = 180.0 / theta_deg  # Rough conversion
            
            # Adjust based on cosmology
            if self.cosmo.model_name == 'EDE':
                # EDE changes the angular diameter distance
                ell_scale *= 0.98  # Slightly smaller scale
            
            return theta_deg, ell_scale
        
        elif scale_type == 'equality':
            # Scale of matter-radiation equality
            z_eq = 3400.0
            k_eq = 0.01  # h/Mpc, approximate
            
            # Convert to angular scale
            ell_eq = k_eq * 14000  # Rough conversion
            
            return None, ell_eq

# Test function
def test_cmb_spectra():
    """Test CMB spectra calculations."""
    import matplotlib.pyplot as plt
    
    # Initialize cosmologies
    ede = EDECosmology(model_name='EDE')
    lcdm = EDECosmology(model_name='LCDM')
    
    # Create CMB calculators
    ede_params = {'A_s': 2.1e-9, 'n_s': 0.965, 'tau': 0.054}
    lcdm_params = {'A_s': 2.1e-9, 'n_s': 0.9649, 'tau': 0.0544}
    
    ede_cmb = CMBSpectra(ede, ede_params)
    lcdm_cmb = CMBSpectra(lcdm, lcdm_params)
    
    # Create ℓ array
    ell = np.arange(2, 2501)
    
    # Calculate spectra
    D_ell_ede = ede_cmb.temperature_spectrum(ell)
    D_ell_lcdm = lcdm_cmb.temperature_spectrum(ell)
    
    # Calculate polarization
    D_ell_EE_ede = ede_cmb.polarization_spectrum(ell, 'EE')
    D_ell_TE_ede = ede_cmb.polarization_spectrum(ell, 'TE')
    
    # Plot results
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Temperature spectrum
    ax = axes[0, 0]
    ax.semilogx(ell, D_ell_ede, 'r-', label='EDE', alpha=0.8, linewidth=1.5)
    ax.semilogx(ell, D_ell_lcdm, 'b--', label='ΛCDM', alpha=0.8, linewidth=1.5)
    ax.set_xlabel(r'Multipole $\ell$')
    ax.set_ylabel(r'$D_\ell^{TT}$ [$\mu K^2$]')
    ax.set_xlim(2, 2500)
    ax.set_ylim(0, 6000)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Temperature ratio
    ax = axes[0, 1]
    ratio = D_ell_ede / D_ell_lcdm
    ax.semilogx(ell, ratio, 'k-', linewidth=1.5)
    ax.axhline(1.0, color='gray', linestyle='--', alpha=0.5)
    
    # Highlight EDE features
    mask_low = ell < 30
    if np.any(mask_low):
        ax.fill_between(ell[mask_low], ratio[mask_low], 1.0, 
                       alpha=0.3, color='red', label='Low-ℓ suppression')
    
    ax.set_xlabel(r'Multipole $\ell$')
    ax.set_ylabel(r'$D_\ell^{\rm EDE}/D_\ell^{\Lambda\rm CDM}$')
    ax.set_xlim(2, 2500)
    ax.set_ylim(0.9, 1.1)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # E-mode polarization
    ax = axes[1, 0]
    ax.semilogx(ell, D_ell_EE_ede, 'r-', label='EDE EE', linewidth=1.5)
    ax.set_xlabel(r'Multipole $\ell$')
    ax.set_ylabel(r'$D_\ell^{EE}$ [$\mu K^2$]')
    ax.set_xlim(2, 2500)
    ax.set_ylim(0, 60)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # TE correlation
    ax = axes[1, 1]
    ax.plot(ell[:100], D_ell_TE_ede[:100], 'r-', linewidth=1.5)
    ax.set_xlabel(r'Multipole $\ell$')
    ax.set_ylabel(r'$D_\ell^{TE}$ [$\mu K^2$]')
    ax.set_xlim(2, 100)
    ax.set_ylim(-30, 30)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('cmb_spectra_comparison.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    test_cmb_spectra()
