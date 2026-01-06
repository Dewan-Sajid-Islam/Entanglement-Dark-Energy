"""
Matter Power Spectrum Calculations
Author: Dewan Sajid Islam
Date: 2026
"""

import numpy as np
from scipy.interpolate import CubicSpline
from .background_evolution import EDECosmology
from .growth_functions import GrowthCalculator

class PowerSpectrum:
    """
    Calculate matter power spectrum for EDE model.
    """
    
    def __init__(self, cosmo, A_s=2.1e-9, n_s=0.965):
        """
        Initialize power spectrum calculator.
        
        Parameters
        ----------
        cosmo : EDECosmology instance
            Cosmological model
        A_s : float
            Primordial amplitude
        n_s : float
            Spectral index
        """
        self.cosmo = cosmo
        self.A_s = A_s
        self.n_s = n_s
        self.growth_calc = GrowthCalculator(cosmo)
        
        # Eisenstein & Hu transfer function parameters
        self.theta_cmb = 2.7255 / 2.7  # T_CMB/2.7K
        
    def eisenstein_hu_transfer(self, k, z=0.0):
        """
        Eisenstein & Hu (1998) transfer function.
        
        Parameters
        ----------
        k : array
            Wavenumbers in h/Mpc
        z : float
            Redshift
            
        Returns
        -------
        T(k) : array
            Transfer function
        """
        # Convert k to 1/Mpc
        k_Mpc = k * self.cosmo.h  # Now in 1/Mpc
        
        # Shape parameter (Eisenstein & Hu 1998)
        Omega_m = self.cosmo.Omega_m0
        Omega_b = self.cosmo.Omega_b0
        h = self.cosmo.h
        
        # Effective shape parameter
        theta = self.theta_cmb
        Omega_m_h2 = Omega_m * h**2
        Omega_b_h2 = Omega_b * h**2
        
        # Equation 26 from EH98
        Gamma_eff = Omega_m_h2 * (np.sqrt(0.5/theta) / theta) * \
                   (1.0 + (0.5/theta)**0.7 * Omega_b_h2/Omega_m_h2)
        
        # Sound horizon scale
        k_eq = 0.0746 * Omega_m_h2 * theta**(-2)  # Equation 3
        s = 44.5 * np.log(9.83/Omega_m_h2) / np.sqrt(1.0 + 10.0 * Omega_b_h2**0.75)  # Equation 6
        
        # Transfer function (BBKS inspired but with EH corrections)
        q = k_Mpc / (Gamma_eff * h)
        
        # Baryon suppression (Equation 17)
        alpha_Gamma = 1.0 - 0.328 * np.log(431.0 * Omega_m_h2) * Omega_b_h2/Omega_m_h2 + \
                      0.38 * np.log(22.3 * Omega_m_h2) * (Omega_b_h2/Omega_m_h2)**2
        
        # Main transfer function
        T = np.log(1.0 + 2.34*q) / (2.34*q)
        T *= (1.0 + 3.89*q + (16.1*q)**2 + (5.46*q)**3 + (6.71*q)**4)**(-0.25)
        
        # Baryon corrections
        if Omega_b_h2 > 0:
            # Silk damping scale (Equation 11)
            k_silk = 1.6 * (Omega_b_h2)**0.52 * (Omega_m_h2)**0.73 * \
                     (1.0 + (10.4*Omega_m_h2)**(-0.95))
            
            # Baryon suppression (Equation 20)
            f_baryon = 1.0 / (1.0 + (k_Mpc*s/5.4)**4)
            T_baryon = T / (1.0 + (k_Mpc*s/5.2)**2) + \
                       f_baryon / (1.0 + (k_Mpc*s/5.4)**2)
            
            T = alpha_Gamma * T + (1.0 - alpha_Gamma) * T_baryon
        
        # Apply growth factor for redshift
        D = self.growth_calc.get_growth_factor(z)
        
        return T * D
    
    def primordial_power(self, k):
        """
        Primordial power spectrum.
        
        Parameters
        ----------
        k : array
            Wavenumbers in h/Mpc
            
        Returns
        -------
        P_prim(k) : array
            Primordial power spectrum
        """
        k_Mpc = k * self.cosmo.h  # Convert to 1/Mpc
        k_pivot = 0.05  # Mpc^-1, Planck pivot scale
        
        # Power-law spectrum
        P_prim = self.A_s * (k_Mpc / k_pivot)**(self.n_s - 1.0)
        
        return P_prim
    
    def linear_power_spectrum(self, k, z=0.0):
        """
        Linear matter power spectrum.
        
        Parameters
        ----------
        k : array
            Wavenumbers in h/Mpc
        z : float
            Redshift
            
        Returns
        -------
        P(k,z) : array
            Power spectrum in (Mpc/h)^3
        """
        # Primordial spectrum
        P_prim = self.primordial_power(k)
        
        # Transfer function
        T = self.eisenstein_hu_transfer(k, z)
        
        # Power spectrum
        k_Mpc = k * self.cosmo.h  # 1/Mpc
        
        # Normalization (Equation 29 from EH98)
        delta_H = 1.94e-5 * self.cosmo.Omega_m0**(-0.785 - 0.05*np.log(self.cosmo.Omega_m0)) * \
                  np.exp(-0.95*(self.n_s - 1) - 0.169*(self.n_s - 1)**2)
        
        # Final power spectrum
        P_k = (2*np.pi**2 / k_Mpc**3) * delta_H**2 * (k_Mpc/0.002)**self.n_s * T**2
        
        # Convert to (Mpc/h)^3
        P_k = P_k / self.cosmo.h**3
        
        return P_k
    
    def sigma_R(self, R, z=0.0):
        """
        Calculate σ_R, the RMS fluctuation in spheres of radius R.
        
        Parameters
        ----------
        R : float
            Radius in Mpc/h
        z : float
            Redshift
            
        Returns
        -------
        sigma_R : float
            RMS fluctuation
        """
        # Create k array
        k_min = 1e-4
        k_max = 10.0
        k = np.logspace(np.log10(k_min), np.log10(k_max), 1000)
        
        # Power spectrum
        P_k = self.linear_power_spectrum(k, z)
        
        # Window function (top-hat in Fourier space)
        x = k * R
        W = 3.0 * (np.sin(x) - x*np.cos(x)) / x**3
        
        # Integral for sigma^2
        integrand = k**2 * P_k * W**2 / (2*np.pi**2)
        sigma2 = np.trapz(integrand, k)
        
        return np.sqrt(sigma2)
    
    def sigma8(self, z=0.0):
        """
        Calculate σ_8 at redshift z.
        
        Parameters
        ----------
        z : float
            Redshift
            
        Returns
        -------
        sigma8 : float
            σ_8 value
        """
        return self.sigma_R(8.0, z)
    
    def nonlinear_power(self, k, z=0.0, method='halofit'):
        """
        Nonlinear power spectrum using HALOFIT approximation.
        
        Parameters
        ----------
        k : array
            Wavenumbers in h/Mpc
        z : float
            Redshift
        method : str
            'halofit' for Takahashi et al. (2012) implementation
            
        Returns
        -------
        P_nl(k,z) : array
            Nonlinear power spectrum
        """
        # For now, return linear spectrum
        # Full HALOFIT implementation is complex
        return self.linear_power_spectrum(k, z)

# Test function
def test_power_spectrum():
    """Test power spectrum calculations."""
    import matplotlib.pyplot as plt
    
    # Initialize cosmology
    ede = EDECosmology(model_name='EDE')
    lcdm = EDECosmology(model_name='LCDM')
    
    # Create power spectrum calculators
    ede_ps = PowerSpectrum(ede)
    lcdm_ps = PowerSpectrum(lcdm)
    
    # Create k array
    k = np.logspace(-4, 1, 500)  # h/Mpc
    
    # Calculate power spectra
    Pk_ede = ede_ps.linear_power_spectrum(k, z=0.0)
    Pk_lcdm = lcdm_ps.linear_power_spectrum(k, z=0.0)
    
    # Calculate sigma8
    sigma8_ede = ede_ps.sigma8(z=0.0)
    sigma8_lcdm = lcdm_ps.sigma8(z=0.0)
    
    print(f"EDE σ_8 = {sigma8_ede:.3f}")
    print(f"ΛCDM σ_8 = {sigma8_lcdm:.3f}")
    
    # Plot results
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Power spectra
    ax = axes[0]
    ax.loglog(k, Pk_ede, 'r-', label='EDE', linewidth=2)
    ax.loglog(k, Pk_lcdm, 'b--', label='ΛCDM', linewidth=2)
    ax.set_xlabel(r'$k$ [$h$ Mpc$^{-1}$]')
    ax.set_ylabel(r'$P(k)$ [$(h^{-1}$ Mpc)$^3$]')
    ax.set_xlim(1e-4, 10)
    ax.set_ylim(1e-4, 1e5)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Ratio
    ax = axes[1]
    ratio = Pk_ede / Pk_lcdm
    ax.semilogx(k, ratio, 'k-', linewidth=2)
    ax.axhline(1.0, color='gray', linestyle='--', alpha=0.5)
    
    # Highlight suppression region
    mask = k < 0.01
    if np.any(mask):
        ax.fill_between(k[mask], ratio[mask], 1.0, 
                       alpha=0.3, color='red', label='EDE suppression')
    
    ax.set_xlabel(r'$k$ [$h$ Mpc$^{-1}$]')
    ax.set_ylabel(r'$P_{\rm EDE}(k)/P_{\Lambda\rm CDM}(k)$')
    ax.set_xlim(1e-4, 10)
    ax.set_ylim(0.9, 1.1)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('power_spectrum_comparison.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    test_power_spectrum()
