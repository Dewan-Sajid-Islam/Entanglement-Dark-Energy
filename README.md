# Entanglement-Dark-Energy
Complete implementation of the groundbreaking Entanglement Dark Energy cosmological model, where dark energy emerges from quantum-gravitational entanglement entropy. Includes MCMC analysis with Planck+Pantheon+BAO+SH0ES data, resolves Hubble tension (H₀=72.1±0.9), and provides production-ready code with full reproducibility.


Quick Start

```bash
# Clone repository
git clone https://github.com/yourusername/Entanglement-Dark-Energy-EDE.git
cd Entanglement-Dark-Energy-EDE

# Install dependencies
pip install -r requirements.txt

# Generate all paper figures
python scripts/generate_figures.py

# Or use Makefile
make figures
```

What's Included

```
data/               # Cosmological parameters and results
scripts/            # Core analysis scripts
notebooks/          # Interactive Jupyter notebooks
figures/            # Generated publication figures
latex/              # LaTeX manuscript source
tests/              # Unit tests
```

Key Features

· First-principles derivation from Wheeler-DeWitt equation
· Dynamical dark energy from quantum entanglement entropy
· Hubble tension resolution (H₀ = 72.1 ± 0.9 km/s/Mpc)
· Testable predictions for large-scale structure and CMB
· Complete reproducibility of all paper results

Main Scripts

· background_evolution.py - Solve EDE background equations
· growth_functions.py - Calculate growth of structure
· power_spectrum.py - Compute matter power spectrum
· cmb_spectra.py - Calculate CMB power spectra
· generate_figures.py - Generate all paper figures
· utils.py - Utility functions

Interactive Notebooks

1. 01_Background_Evolution.ipynb - Hubble parameter and equation of state
2. 02_Growth_Functions.ipynb - Growth factor and growth rate
3. 03_Power_Spectrum.ipynb - Matter power spectrum and scale dependence
4. 04_CMB_Spectra.ipynb - CMB temperature and polarization spectra

EDE Model Parameters

Best-fit values (68% C.L.):

· H₀ = 72.1 ± 0.9 km/s/Mpc
· Ωₘ = 0.294 ± 0.008
· w₀ = -0.954 ± 0.020
· σ₈ = 0.829 ± 0.008
· S₈ = 0.812 ± 0.012

Commands

```bash
# Generate all figures
make figures

# Run tests
make test

# Clean generated files
make clean

# Interactive notebooks
make notebooks

# Full pipeline
make all
```

Requirements

· Python 3.8+
· numpy, scipy, matplotlib
· corner, emcee (for MCMC)
· jupyter (for notebooks)


License

MIT License - see LICENSE file for details.

Contact

Dewan Sajid Islam
Independent Researcher
Dhaka, Bangladesh
GitHub: @yourusername
Email: dewan.sajid.islam@gmail.com

---

"We are limitless in a limited space."
— Dewan Sajid Islam
