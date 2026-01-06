"""
Entanglement Dark Energy Background Evolution Solver
Author: Dewan Sajid Islam
Date: 2024
"""

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import root
from typing import Tuple, Callable

class EDECosmology:
    """
    Solves the background evolution equations for the Entanglement Dark Energy model.
    """
    
    def __init__(self, H0: float = 72.1, Omega_m0: float = 0.294, 
                 Omega_r0: float = 8.24e-5, Omega_k0: float = 0.0):
        """
        Initialize cosmology with given parameters.
        
        Parameters
        ----------
        H0 : float
            Hubble constant in km/s/Mpc
        Omega_m0 : float
            Present matter density parameter
        Omega_r0 : float
            Present radiation density parameter
        Omega_k0 : float
            Present curvature density parameter (0 for flat)
        """
        self.H0 = H0
        self.h = H0 / 100.0
        self.Omega_m0 = Omega_m0
        self.Omega_r0 = Omega_r0
        self.Omega_k0 = Omega_k0
        
        # Constants
        self.c = 2.99792458e5  # km/s
        self.G = 4.30091e-9    # Mpc (km/s)^2 / M_sun
        self.Mpc_to_km = 3.08567758e19  # 1 Mpc in km
        
        # Derived
        self.H0_s = H0 / self.Mpc_to_km  # H0 in 1/s
        self.rho_crit0 = 3 * self.H0_s**2 / (8 * np.pi * self.G)  # Critical density
        
    def H_LCDM(self, z: float) -> float:
        """ΛCDM Hubble parameter."""
        a = 1.0 / (1.0 + z)
        return self.H0 * np.sqrt(self.Omega_m0 * a**-3 + self.Omega_r0 * a**-4 
                                + self.Omega_k0 * a**-2 + (1 - self.Omega_m0 - self.Omega_r0 - self.Omega_k0))
    
    def future_horizon_integral(self, z: float, H_func: Callable) -> float:
        """Calculate future event horizon R_h(z)."""
        def integrand(zp):
            return 1.0 / H_func(zp)
        
        # Adaptive integration
        result, _ = quad(integrand, z, np.inf, epsabs=1e-10, epsrel=1e-10, 
                        limit=1000, points=[z, z+0.1, z+1.0, z+10.0])
        return result / (1.0 + z)
    
    def rho_ent(self, z: float, H: float, Rh: float, dRh_dz: float) -> float:
        """Entanglement dark energy density (Eq. 8)."""
        Hz = H * self.Mpc_to_km  # Convert to 1/s for density calculation
        Rh_mpc = Rh  # Already in Mpc
        
        # Time derivative: dRh/dt = -H(1+z)dRh/dz
        dRh_dt = -Hz * (1 + z) * dRh_dz
        
        term1 = 1.0 / (Rh_mpc**2)
        term2 = dRh_dt / (Hz * Rh_mpc**3)
        
        rho = (3.0 / (8.0 * np.pi * self.G)) * (term1 - term2)
        return rho
    
    def w_ent(self, z: float, H: float, Rh: float, dRh_dz: float, d2Rh_dz2: float) -> float:
        """Equation of state parameter (Eq. 11)."""
        Hz = H * self.Mpc_to_km
        Rh_mpc = Rh
        
        # First and second derivatives
        dRh_dt = -Hz * (1 + z) * dRh_dz
        d2Rh_dt2 = Hz**2 * (1 + z) * ((1 + z) * d2Rh_dz2 + dRh_dz)
        
        HRh = Hz * Rh_mpc
        dHRh_dt = dRh_dt + Rh_mpc * Hz
        
        # Avoid division by zero
        if abs(HRh - 1) < 1e-10:
            return -1.0
        
        w = -1.0 + (2.0/3.0) * (dRh_dt / (Hz * Rh_mpc)) * \
            (1.0 - 1.0/HRh + d2Rh_dt2/(2.0 * Hz * dRh_dt) - dRh_dt/(2.0 * Hz * Rh_mpc)) / \
            (1.0 - dRh_dt/(Hz * Rh_mpc))
        
        # Physical bounds
        return np.clip(w, -1.5, 0.0)
    
    def solve_background(self, z_max: float = 10.0, n_points: int = 1000) -> dict:
        """
        Solve the coupled background equations for EDE.
        
        Returns
        -------
        dict with arrays for z, H, Rh, Omega_ent, w_ent
        """
        # Log-spaced redshift array
        z = np.logspace(-4, np.log10(z_max), n_points)
        z = np.sort(np.concatenate([z, np.linspace(0, 2, 500)]))  # Denser at low-z
        
        # Initial guess using LCDM
        H_guess = self.H_LCDM(z)
        Rh_guess = np.array([self.future_horizon_integral(zi, self.H_LCDM) for zi in z])
        
        # Iterative solution
        for iteration in range(20):
            H_new = np.zeros_like(z)
            
            # Interpolate for derivative calculations
            from scipy.interpolate import CubicSpline
            H_interp = CubicSpline(z, H_guess)
            Rh_interp = CubicSpline(z, Rh_guess)
            
            for i, zi in enumerate(z):
                # Calculate derivatives
                if i == 0:
                    dRh_dz = (Rh_guess[1] - Rh_guess[0]) / (z[1] - z[0])
                    d2Rh_dz2 = 0.0
                elif i == len(z) - 1:
                    dRh_dz = (Rh_guess[-1] - Rh_guess[-2]) / (z[-1] - z[-2])
                    d2Rh_dz2 = 0.0
                else:
                    dRh_dz = (Rh_guess[i+1] - Rh_guess[i-1]) / (z[i+1] - z[i-1])
                    d2Rh_dz2 = (Rh_guess[i+1] - 2*Rh_guess[i] + Rh_guess[i-1]) / ((z[i+1] - z[i])**2)
                
                # Calculate entanglement density
                rho_ent_val = self.rho_ent(zi, H_guess[i], Rh_guess[i], dRh_dz)
                
                # Calculate total density
                a = 1.0 / (1.0 + zi)
                rho_m = self.Omega_m0 * a**-3 * self.rho_crit0
                rho_r = self.Omega_r0 * a**-4 * self.rho_crit0
                rho_tot = rho_m + rho_r + rho_ent_val
                
                # New Hubble parameter
                H_new[i] = np.sqrt(8 * np.pi * self.G * rho_tot / 3) * self.Mpc_to_km
            
            # Update future horizon with new H(z)
            def H_func_new(zp):
                return np.interp(zp, z, H_new)
            
            Rh_new = np.array([self.future_horizon_integral(zi, H_func_new) for zi in z])
            
            # Check convergence
            delta_H = np.max(np.abs(H_new - H_guess) / H_guess)
            delta_Rh = np.max(np.abs(Rh_new - Rh_guess) / np.abs(Rh_guess))
            
            if delta_H < 1e-6 and delta_Rh < 1e-6:
                print(f"Converged after {iteration+1} iterations")
                H_guess, Rh_guess = H_new, Rh_new
                break
            
            # Damped update
            H_guess = 0.7 * H_guess + 0.3 * H_new
            Rh_guess = 0.7 * Rh_guess + 0.3 * Rh_new
        
        # Calculate derived quantities
        Omega_ent = np.zeros_like(z)
        w = np.zeros_like(z)
        
        H_interp = CubicSpline(z, H_guess)
        Rh_interp = CubicSpline(z, Rh_guess)
        
        for i, zi in enumerate(z):
            if i == 0:
                dRh_dz = (Rh_guess[1] - Rh_guess[0]) / (z[1] - z[0])
                d2Rh_dz2 = (Rh_guess[2] - 2*Rh_guess[1] + Rh_guess[0]) / ((z[1] - z[0])**2)
            elif i == len(z) - 1:
                dRh_dz = (Rh_guess[-1] - Rh_guess[-2]) / (z[-1] - z[-2])
                d2Rh_dz2 = (Rh_guess[-1] - 2*Rh_guess[-2] + Rh_guess[-3]) / ((z[-1] - z[-2])**2)
            else:
                dRh_dz = (Rh_guess[i+1] - Rh_guess[i-1]) / (z[i+1] - z[i-1])
                d2Rh_dz2 = (Rh_guess[i+1] - 2*Rh_guess[i] + Rh_guess[i-1]) / ((z[i+1] - z[i])**2)
            
            rho_ent_val = self.rho_ent(zi, H_guess[i], Rh_guess[i], dRh_dz)
            H_s = H_guess[i] / self.Mpc_to_km
            rho_crit = 3 * H_s**2 / (8 * np.pi * self.G)
            
            Omega_ent[i] = rho_ent_val / rho_crit
            w[i] = self.w_ent(zi, H_guess[i], Rh_guess[i], dRh_dz, d2Rh_dz2)
        
        return {
            'z': z,
            'H': H_guess,
            'Rh': Rh_guess,
            'Omega_ent': Omega_ent,
            'w_ent': w,
            'Omega_m': self.Omega_m0 * (1 + z)**3 / (H_guess/self.H0)**2,
            'Omega_r': self.Omega_r0 * (1 + z)**4 / (H_guess/self.H0)**2
          }
