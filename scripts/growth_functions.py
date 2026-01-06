"""
Growth Functions for Structure Formation
Author: Dewan Sajid Islam
Date: 2026
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from .background_evolution import EDECosmology

class GrowthCalculator:
    """
    Calculate linear growth of structure for EDE model.
    """
    
    def __init__(self, cosmo):
        """
        Initialize with cosmology.
        
        Parameters
        ----------
        cosmo : EDECosmology instance
            Cosmological model
        """
        self.cosmo = cosmo
        self.growth_data = None
        
    def solve_growth_equation(self, z_max=10.0, n_points=1000):
        """
        Solve the linear growth equation.
        
        The equation is:
        d²δ/da² + (3/a + dlnH/da) dδ/da - (3Ω_m/(2a²)) δ = 0
        
        Parameters
        ----------
        z_max : float
            Maximum redshift
        n_points : int
            Number of points
            
        Returns
        -------
        dict with z, D, f
        """
        print("Solving growth equation...")
        
        # Create scale factor array (from early to late)
        a_min = 1e-3
        a_max = 1.0
        a = np.logspace(np.log10(a_min), np.log10(a_max), n_points)
        z = 1.0/a - 1.0
        
        # Get background evolution
        if self.cosmo.H_background is None:
            self.cosmo.solve_background(z_max=z_max)
        
        # Interpolate H(z) and Ω_m(z)
        H_interp = CubicSpline(self.cosmo.z_background, self.cosmo.H_background)
        
        def omega_m(z_val):
            """Matter density parameter at redshift z."""
            H_z = H_interp(z_val) * 1000 / self.cosmo.Mpc_to_m  # Convert to 1/s
            a_val = 1.0 / (1.0 + z_val)
            rho_m = self.cosmo.Omega_m0 * a_val**-3 * self.cosmo.rho_crit0
            rho_crit = 3 * H_z**2 / (8 * np.pi * self.cosmo.G)
            return rho_m / rho_crit
        
        # Growth ODE in terms of scale factor
        def growth_ode(a_val, y):
            """ODE for growth factor D and growth rate f."""
            D, f = y
            z_val = 1.0/a_val - 1.0
            
            # Hubble parameter and derivative
            H = H_interp(z_val) * 1000 / self.cosmo.Mpc_to_m
            
            # Numerical derivative dH/da
            da = a_val * 1e-6
            H_plus = H_interp(1.0/(a_val + da) - 1.0) * 1000 / self.cosmo.Mpc_to_m
            H_minus = H_interp(1.0/(a_val - da) - 1.0) * 1000 / self.cosmo.Mpc_to_m
            dH_da = (H_plus - H_minus) / (2 * da)
            
            # dlnH/da
            dlnH_dlna = a_val * dH_da / H
            
            # Ω_m at this redshift
            Omega_m = omega_m(z_val)
            
            # Growth rate derivative
            df_dlna = f**2 + f * (0.5 - 1.5 * Omega_m) + 1.5 * Omega_m
            
            # Convert to derivatives with respect to a
            dD_da = f * D / a_val
            df_da = df_dlna * f / a_val
            
            return [dD_da, df_da]
        
        # Initial conditions during matter domination (a << 1)
        a_start = a_min
        D_start = a_start  # D ∝ a during matter domination
        f_start = 1.0  # f = dlnD/dlna ≈ 1 during matter domination
        
        # Solve ODE
        sol = solve_ivp(
            growth_ode,
            [a_start, a_max],
            [D_start, f_start],
            t_eval=a,
            method='RK45',
            rtol=1e-8,
            atol=1e-10
        )
        
        if not sol.success:
            print(f"Growth ODE solver failed: {sol.message}")
            return None
        
        # Extract solution
        D = sol.y[0]
        f = sol.y[1]
        
        # Normalize D(z=0) = 1
        D = D / D[-1]
        
        # Store results
        self.growth_data = {
            'z': z,
            'a': a,
            'D': D,
            'f': f,
            'sigma8': None  # Will be computed separately
        }
        
        return self.growth_data
    
    def compute_sigma8(self, sigma8_0=0.83):
        """
        Compute σ8(z) from growth factor.
        
        Parameters
        ----------
        sigma8_0 : float
            σ8 at z=0
            
        Returns
        -------
        sigma8_z : array
            σ8(z)
        """
        if self.growth_data is None:
            self.solve_growth_equation()
        
        # σ8(z) = σ8(0) * D(z)
        sigma8_z = sigma8_0 * self.growth_data['D']
        
        self.growth_data['sigma8'] = sigma8_z
        
        return sigma8_z
    
    def get_growth_factor(self, z):
        """Get growth factor D(z)."""
        if self.growth_data is None:
            self.solve_growth_equation()
        
        return np.interp(z, self.growth_data['z'], self.growth_data['D'])
    
    def get_growth_rate(self, z):
        """Get growth rate f(z)."""
        if self.growth_data is None:
            self.solve_growth_equation()
        
        return np.interp(z, self.growth_data['z'], self.growth_data['f'])
    
    def get_fsigma8(self, z, sigma8_0=0.83):
        """Get fσ8(z)."""
        if self.growth_data is None:
            self.solve_growth_equation()
        
        if self.growth_data['sigma8'] is None:
            self.compute_sigma8(sigma8_0)
        
        f = np.interp(z, self.growth_data['z'], self.growth_data['f'])
        sigma8 = np.interp(z, self.growth_data['z'], self.growth_data['sigma8'])
        
        return f * sigma8

# Test function
def compare_growth_models():
    """Compare growth in EDE vs LCDM."""
    import matplotlib.pyplot as plt
    
    # Initialize cosmologies
    ede = EDECosmology(model_name='EDE')
    lcdm = EDECosmology(model_name='LCDM')
    
    # Create growth calculators
    ede_growth = GrowthCalculator(ede)
    lcdm_growth = GrowthCalculator(lcdm)
    
    # Solve growth equations
    ede_data = ede_growth.solve_growth_equation(z_max=5.0)
    lcdm_data = lcdm_growth.solve_growth_equation(z_max=5.0)
    
    # Compute σ8
    ede_sigma8 = ede_growth.compute_sigma8(sigma8_0=0.829)
    lcdm_sigma8 = lcdm_growth.compute_sigma8(sigma8_0=0.811)
    
    # Plot results
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Growth factor
    ax = axes[0]
    ax.plot(ede_data['z'], ede_data['D'], 'r-', label='EDE', linewidth=2)
    ax.plot(lcdm_data['z'], lcdm_data['D'], 'b--', label='ΛCDM', linewidth=2)
    ax.set_xlabel('Redshift z')
    ax.set_ylabel('Growth Factor D(z)')
    ax.set_xlim(0, 5)
    ax.set_ylim(0, 1.1)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Growth rate
    ax = axes[1]
    ax.plot(ede_data['z'], ede_data['f'], 'r-', label='EDE', linewidth=2)
    ax.plot(lcdm_data['z'], lcdm_data['f'], 'b--', label='ΛCDM', linewidth=2)
    ax.set_xlabel('Redshift z')
    ax.set_ylabel('Growth Rate f(z)')
    ax.set_xlim(0, 5)
    ax.set_ylim(0.4, 1.1)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # fσ8
    ax = axes[2]
    z = ede_data['z']
    fsigma8_ede = ede_growth.get_fsigma8(z, sigma8_0=0.829)
    fsigma8_lcdm = lcdm_growth.get_fsigma8(z, sigma8_0=0.811)
    
    ax.plot(z, fsigma8_ede, 'r-', label='EDE', linewidth=2)
    ax.plot(z, fsigma8_lcdm, 'b--', label='ΛCDM', linewidth=2)
    ax.set_xlabel('Redshift z')
    ax.set_ylabel(r'$f\sigma_8(z)$')
    ax.set_xlim(0, 2)
    ax.set_ylim(0.3, 0.6)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('growth_comparison.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    compare_growth_models()
