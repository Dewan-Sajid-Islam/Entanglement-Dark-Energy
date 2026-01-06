"""
Full MCMC Analysis for Entanglement Dark Energy Model
Author: Dewan Sajid Islam
Date: 2024
"""

import numpy as np
import emcee
import corner
import matplotlib.pyplot as plt
from scipy.stats import multivariate_normal
from scipy.optimize import minimize
import h5py
import json
from datetime import datetime
from background_evolution import EDECosmology

class EDELikelihood:
    """
    Likelihood calculation for EDE model including CMB, SN, BAO, and H0 data.
    """
    
    def __init__(self):
        # Load observational data
        self.load_data()
        
        # Priors
        self.prior_bounds = {
            'H0': (60.0, 80.0),
            'Omega_m': (0.1, 0.5),
            'Omega_b': (0.02, 0.06),
            'sigma8': (0.6, 0.9),
            'n_s': (0.92, 1.01),
            'tau': (0.01, 0.2),
            'w0': (-1.2, -0.8),
            'wa': (-1.0, 1.0)
        }
        
    def load_data(self):
        """Load observational data sets."""
        # Planck CMB likelihood (simplified - in practice use full Planck likelihood)
        self.cmb_data = {
            'ell': np.arange(2, 2501),
            'Cl_tt': None,  # Would load from Planck data
            'Cl_ee': None,
            'Cl_te': None,
            'covariance': None
        }
        
        # Pantheon+ SN data (simplified)
        self.sn_data = {
            'z': np.array([0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0]),
            'mu': np.array([33.0, 34.0, 35.5, 37.0, 38.5, 39.5, 40.5, 41.2, 42.5, 43.8]),
            'mu_err': np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2])
        }
        
        # BAO data
        self.bao_data = {
            'z': [0.106, 0.15, 0.38, 0.51, 0.61],
            'DV_over_rs': [2.98, 4.47, 8.42, 10.75, 12.85],
            'DV_over_rs_err': [0.13, 0.16, 0.14, 0.13, 0.13]
        }
        
        # H0 measurement
        self.h0_data = {
            'value': 73.04,
            'error': 1.04
        }
        
    def compute_distance_modulus(self, z, cosmo):
        """Compute distance modulus for SN data."""
        from scipy.integrate import quad
        
        def integrand(zp):
            return 1.0 / cosmo.H(zp)
        
        mu = np.zeros_like(z)
        for i, zi in enumerate(z):
            dc, _ = quad(integrand, 0, zi)
            dl = dc * (1 + zi)
            mu[i] = 5 * np.log10(dl * cosmo.Mpc_to_km) + 25  # dl in Mpc
        return mu
    
    def compute_bao_observable(self, z, cosmo):
        """Compute BAO observable DV/rs."""
        from scipy.integrate import quad
        
        def integrand(zp):
            return 1.0 / cosmo.H(zp)
        
        DV_over_rs = np.zeros_like(z)
        for i, zi in enumerate(z):
            # Sound horizon at drag epoch (approximate)
            rs = 147.0  # Mpc, approximate value
            
            # Comoving angular diameter distance
            dc, _ = quad(integrand, 0, zi)
            da = dc / (1 + zi)
            
            # DV = [cz * DA^2 / H(z)]^(1/3)
            DV = ((cosmo.c * zi * da**2) / cosmo.H(zi))**(1/3)
            DV_over_rs[i] = DV / rs
        
        return DV_over_rs
    
    def log_prior(self, theta):
        """Log prior probability."""
        H0, Omega_m, Omega_b, sigma8, n_s, tau, w0, wa = theta
        
        # Check bounds
        for i, (key, val) in enumerate(self.prior_bounds.items()):
            if not (val[0] <= theta[i] <= val[1]):
                return -np.inf
        
        # Additional physical priors
        if Omega_m + Omega_b > 0.99:
            return -np.inf
        
        # Flat priors within bounds
        return 0.0
    
    def log_likelihood(self, theta):
        """Log likelihood for given parameters."""
        # Priors
        lp = self.log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        
        H0, Omega_m, Omega_b, sigma8, n_s, tau, w0, wa = theta
        
        # Initialize cosmology
        cosmo = EDECosmology(H0=H0, Omega_m0=Omega_m, Omega_r0=2.47e-5/H0**2)
        
        try:
            # Solve background
            bg = cosmo.solve_background()
            
            # Calculate observables
            chi2_total = 0.0
            
            # 1. H0 likelihood (Gaussian)
            chi2_h0 = ((H0 - self.h0_data['value']) / self.h0_data['error'])**2
            chi2_total += chi2_h0
            
            # 2. SN likelihood
            mu_pred = self.compute_distance_modulus(self.sn_data['z'], cosmo)
            chi2_sn = np.sum(((mu_pred - self.sn_data['mu']) / self.sn_data['mu_err'])**2)
            chi2_total += chi2_sn
            
            # 3. BAO likelihood
            dv_pred = self.compute_bao_observable(np.array(self.bao_data['z']), cosmo)
            chi2_bao = np.sum(((dv_pred - self.bao_data['DV_over_rs']) / self.bao_data['DV_over_rs_err'])**2)
            chi2_total += chi2_bao
            
            # 4. CMB likelihood (simplified - would use full Planck likelihood)
            # Here we approximate with distance to last scattering and acoustic scale
            z_star = 1089.0  # Redshift of last scattering
            mu_star = self.compute_distance_modulus([z_star], cosmo)[0]
            
            # Approximate CMB chi2 using acoustic scale constraints
            # In practice, use full Planck likelihood code
            chi2_cmb = ((mu_star - 43.0) / 0.3)**2  # Simplified
            
            chi2_total += chi2_cmb
            
            return -0.5 * chi2_total
            
        except Exception as e:
            print(f"Error in likelihood calculation: {e}")
            return -np.inf
    
    def log_probability(self, theta):
        """Total log probability (prior + likelihood)."""
        lp = self.log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.log_likelihood(theta)

class EDEMCMC:
    """
    MCMC sampler for EDE parameter estimation.
    """
    
    def __init__(self, n_walkers=32, n_dim=8):
        self.n_walkers = n_walkers
        self.n_dim = n_dim
        self.likelihood = EDELikelihood()
        
        # Initial guess (LCDM values)
        self.initial_guess = [67.36, 0.3156, 0.0493, 0.811, 0.9649, 0.054, -1.0, 0.0]
        
        # Initialize walkers around guess
        self.initial_pos = self.initial_guess + 1e-3 * np.random.randn(n_walkers, n_dim)
        
        # Sampler
        self.sampler = emcee.EnsembleSampler(
            n_walkers, n_dim, self.likelihood.log_probability,
            moves=emcee.moves.StretchMove(a=2.0)
        )
    
    def run(self, n_steps=5000, burn_in=1000, progress=True):
        """Run MCMC sampling."""
        print(f"Starting MCMC with {self.n_walkers} walkers for {n_steps} steps")
        
        # Run burn-in
        if burn_in > 0:
            print("Running burn-in...")
            pos, prob, state = self.sampler.run_mcmc(
                self.initial_pos, burn_in, progress=progress, skip_initial_state_check=True
            )
            self.sampler.reset()
            initial_pos_burned = pos
        else:
            initial_pos_burned = self.initial_pos
        
        # Main run
        print("Running main MCMC...")
        self.sampler.run_mcmc(initial_pos_burned, n_steps, progress=progress)
        
        # Calculate acceptance fraction
        af = self.sampler.acceptance_fraction
        print(f"Acceptance fraction: {np.mean(af):.3f} ± {np.std(af):.3f}")
        
        # Effective sample size
        tau = self.sampler.get_autocorr_time(quiet=True)
        print(f"Autocorrelation time: {tau}")
        
        return self.sampler
    
    def save_results(self, filename="ede_chains.h5"):
        """Save MCMC chains to HDF5 file."""
        with h5py.File(f"data/chains/{filename}", "w") as f:
            # Save chains
            f.create_dataset("chains", data=self.sampler.get_chain())
            f.create_dataset("log_prob", data=self.sampler.get_log_prob())
            
            # Save metadata
            f.attrs["n_walkers"] = self.n_walkers
            f.attrs["n_dim"] = self.n_dim
            f.attrs["n_steps"] = self.sampler.iteration
            f.attrs["date"] = datetime.now().isoformat()
            
            # Parameter names
            param_names = ['H0', 'Omega_m', 'Omega_b', 'sigma8', 'n_s', 'tau', 'w0', 'wa']
            f.create_dataset("param_names", data=np.array(param_names, dtype='S'))
    
    def analyze_chains(self):
        """Analyze MCMC chains and compute statistics."""
        chains = self.sampler.get_chain(discard=500, thin=50, flat=True)
        
        param_names = ['H0', 'Ω_m', 'Ω_b', 'σ_8', 'n_s', 'τ', 'w_0', 'w_a']
        
        # Compute statistics
        means = np.mean(chains, axis=0)
        medians = np.median(chains, axis=0)
        stds = np.std(chains, axis=0)
        
        # 68% credible intervals
        lower = np.percentile(chains, 16, axis=0)
        upper = np.percentile(chains, 84, axis=0)
        
        print("\n" + "="*60)
        print("MCMC RESULTS (68% CREDIBLE INTERVALS)")
        print("="*60)
        for i, name in enumerate(param_names):
            print(f"{name:6s}: {medians[i]:.4f} +{upper[i]-medians[i]:.4f} -{medians[i]-lower[i]:.4f}")
        
        # Create corner plot
        fig = corner.corner(
            chains,
            labels=param_names,
            quantiles=[0.16, 0.5, 0.84],
            show_titles=True,
            title_kwargs={"fontsize": 12},
            label_kwargs={"fontsize": 14}
        )
        
        plt.savefig("figures/corner_plot.pdf", dpi=300, bbox_inches='tight')
        plt.savefig("figures/corner_plot.png", dpi=300, bbox_inches='tight')
        
        # Save results to JSON
        results = {
            'parameters': param_names,
            'means': means.tolist(),
            'medians': medians.tolist(),
            'stds': stds.tolist(),
            'lower_68': lower.tolist(),
            'upper_68': upper.tolist()
        }
        
        with open('data/best_fit_parameters.json', 'w') as f:
            json.dump(results, f, indent=2)
        
        return results

# Main execution
if __name__ == "__main__":
    # Initialize and run MCMC
    mcmc = EDEMCMC(n_walkers=32, n_dim=8)
    
    # Run for demonstration (use fewer steps for quick run)
    sampler = mcmc.run(n_steps=2000, burn_in=500, progress=True)
    
    # Analyze and save results
    results = mcmc.analyze_chains()
    mcmc.save_results()
    
    print("\nMCMC analysis complete! Results saved to data/")
