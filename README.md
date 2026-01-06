🌌 Entanglement Dark Energy (EDE) Model

<div align="center">

"We are limitless in a limited space."
— Dewan Sajid Islam

</div>

🚀 Quick Start

```bash
# Clone the repository
git clone https://github.com/dewansajidislam/Entanglement-Dark-Energy-EDE.git
cd Entanglement-Dark-Energy-EDE

# Install dependencies
pip install -r requirements.txt

# Generate all figures from the paper
make figures

# Or run everything at once
make all
```

📊 What This Repository Contains

Directory Contents Description
data/ parameters.json, best_fit_parameters.json Cosmological parameters, MCMC results, observational data
scripts/ background_evolution.py, growth_functions.py, power_spectrum.py, cmb_spectra.py Core computational engine for EDE model
notebooks/ 01_Background_Evolution.ipynb, 02_Growth_Functions.ipynb, 03_Power_Spectrum.ipynb, 04_CMB_Spectra.ipynb Interactive Jupyter notebooks for exploration
figures/ figure1.pdf, figure2.pdf, figure3.pdf, figure4.pdf Publication-ready figures (generated)
latex/ main.tex, references.bib Complete LaTeX manuscript source
tests/ test_background.py, test_growth.py Unit tests for validation

🧮 Key Results: EDE vs ΛCDM

Cosmological Parameters (68% CL)

Parameter EDE Model ΛCDM Model Units Improvement
Hubble Constant H₀ = 72.1 ± 0.9 H₀ = 67.36 ± 0.54 km/s/Mpc ✅ Hubble tension resolved
Matter Density Ωₘ = 0.294 ± 0.008 Ωₘ = 0.3156 ± 0.0074 — —
Baryon Density Ω_b = 0.048 ± 0.001 Ω_b = 0.0493 ± 0.0002 — —
Fluctuation Amplitude σ₈ = 0.829 ± 0.008 σ₈ = 0.811 ± 0.006 — —
Structure Parameter S₈ = 0.812 ± 0.012 S₈ = 0.832 ± 0.013 — ✅ Reduces S₈ tension
Dark Energy EoS w₀ = -0.954 ± 0.020 w = -1 (fixed) — ✅ Dynamical dark energy
Spectral Index n_s = 0.965 ± 0.004 n_s = 0.9649 ± 0.0042 — Consistent
Optical Depth τ = 0.054 ± 0.007 τ = 0.0544 ± 0.0073 — Consistent

Derived Quantities

Quantity EDE Value ΛCDM Value Units Physical Meaning
Age of Universe 13.65 ± 0.05 Gyr 13.80 ± 0.02 Gyr Gyr Time since Big Bang
Sound Horizon r_s = 147.0 ± 0.3 r_s = 147.0 ± 0.3 Mpc BAO scale at drag epoch
Angular Scale θ* = 0.010413 ± 0.000005 θ* = 0.010413 ± 0.000005 rad CMB acoustic scale
Redshift of Equality z_eq = 3400 ± 50 z_eq = 3402 ± 26 — Matter-radiation equality
Redshift of Reionization z_reion = 7.7 ± 0.8 z_reion = 7.7 ± 0.8 — Epoch of reionization

🔬 The EDE Model: Core Equations

1. Entanglement Energy Density

```
ρ_ent = (3/(8πG)) * [1/R_h² - (Ṙ_h)/(H R_h³)]
```

Where R_h is the future event horizon

2. Equation of State

```
w_ent = -1 + (2/3)*(Ṙ_h/(H R_h))*[1 - 1/(H R_h) + (̈R_h)/(2H Ṙ_h) - Ṙ_h/(2H R_h)]/[1 - Ṙ_h/(H R_h)]
```

3. Future Event Horizon

```
R_h(t) = a(t) ∫_t^∞ dt'/a(t')
```

🛠️ How to Use This Repository

Command Line Interface

Command Action Output
make install Install all dependencies Installed packages
make test Run unit tests Test results
make figures Generate all paper figures PDF/PNG in figures/
make paper Compile LaTeX manuscript PDF in latex/
make clean Remove generated files Clean workspace
make all Install, test, and generate Complete pipeline

Python API Examples

```python
# Import the EDE model
from scripts.background_evolution import EDECosmology

# Initialize EDE cosmology
ede = EDECosmology(model_name='EDE')
lcdm = EDECosmology(model_name='LCDM')

# Solve background evolution
bg_ede = ede.solve_background(z_max=5.0)
bg_lcdm = lcdm.solve_background(z_max=5.0)

# Get Hubble parameter and equation of state
H_ede = bg_ede['H']           # H(z) for EDE
w_ede = bg_ede['w_ent']       # w(z) for EDE
H_lcdm = bg_lcdm['H']         # H(z) for ΛCDM

print(f"Hubble tension reduction: {(ede.H0/lcdm.H0 - 1)*100:.1f}%")
print(f"EDE w(z=0) = {w_ede[0]:.3f} (dynamical!)")
```

📈 Testable Predictions

1. Scale-Dependent Matter Clustering

Scale (k) EDE Prediction ΛCDM Prediction Observable
k < 0.001 h/Mpc Suppressed by 2-3% No suppression CMB low-ℓ, galaxy surveys
0.001 < k < 0.01 h/Mpc Transition region Scale-invariant Future LSS surveys
k > 0.01 h/Mpc Identical to ΛCDM Scale-invariant Current constraints

2. CMB Anomalies

Feature EDE Signature ΛCDM Expectation Detection Method
Low-ℓ TT Power Suppressed quadrupole Standard Sachs-Wolfe Planck, future CMB
TE Correlation (ℓ<30) Enhanced signal at ℓ≈25 Standard correlation CMB polarization
ISW Effect Unique correlation pattern Standard ISW CMB×LSS cross-correlation

3. Growth of Structure

Redshift fσ₈(EDE) fσ₈(ΛCDM) Enhancement
z = 0.0 0.471 0.462 +2.0%
z = 0.5 0.478 0.470 +1.7%
z = 1.0 0.485 0.479 +1.3%
z = 2.0 0.492 0.488 +0.8%

🎯 Key Features of the EDE Model

Feature EDE Implementation Advantage
Origin of Dark Energy Quantum entanglement entropy First-principles, no new fields
Hubble Tension Natural resolution H₀ = 72.1 km/s/Mpc No fine-tuning required
S₈ Tension Reduced to 0.812 ± 0.012 Better agreement with weak lensing
Dynamical w(z) Derived from entanglement evolution Testable with future surveys
Early Universe Negligible at high-z Preserves BBN and recombination
Theoretical Basis Wheeler-DeWitt + Holography Connects QG to cosmology

📁 Complete File Reference

Essential Files for Reproduction

File Purpose Usage
scripts/background_evolution.py Solve EDE background equations python -c "from scripts.background_evolution import EDECosmology; ede=EDECosmology(); bg=ede.solve_background()"
scripts/growth_functions.py Calculate growth factors python -c "from scripts.growth_functions import GrowthCalculator; from scripts.background_evolution import EDECosmology; ede=EDECosmology(); growth=GrowthCalculator(ede); data=growth.solve_growth_equation()"
scripts/power_spectrum.py Compute matter power spectrum python -c "from scripts.power_spectrum import PowerSpectrum; from scripts.background_evolution import EDECosmology; ede=EDECosmology(); ps=PowerSpectrum(ede); k=[0.01,0.1,1.0]; Pk=ps.linear_power_spectrum(k); print(Pk)"
scripts/cmb_spectra.py Calculate CMB spectra python -c "from scripts.cmb_spectra import CMBSpectra; from scripts.background_evolution import EDECosmology; ede=EDECosmology(); cmb=CMBSpectra(ede); ell=range(2,2501); D_ell=cmb.temperature_spectrum(ell)"
scripts/generate_figures.py Generate all paper figures python scripts/generate_figures.py
data/parameters.json Cosmological parameters Load with json.load(open('data/parameters.json'))

Notebook Workflow

1. Start with: notebooks/01_Background_Evolution.ipynb
2. Then explore: notebooks/02_Growth_Functions.ipynb
3. Continue to: notebooks/03_Power_Spectrum.ipynb
4. Finish with: notebooks/04_CMB_Spectra.ipynb

🧪 Validation & Testing

Test Command Expected Result
Background evolution python -m pytest tests/test_background.py -v All tests pass
Growth functions python -c "from scripts.growth_functions import GrowthCalculator; from scripts.background_evolution import EDECosmology; test_growth()" D(z=0)=1.0, f(z=0)≈0.5
Power spectrum python -c "from scripts.power_spectrum import PowerSpectrum; from scripts.background_evolution import EDECosmology; test_power_spectrum()" P(k) positive, σ₈≈0.8
Figure generation python scripts/generate_figures.py 4 figures + corner plot

🌟 Why EDE is Different

Aspect Traditional Models Entanglement Dark Energy
Theoretical Basis Ad-hoc scalar fields, modified gravity Emergent from quantum gravity
Hubble Tension Requires early dark energy, extra relics Natural outcome of formalism
Fine-tuning Cosmological constant: 10¹²⁰ fine-tuning No fine-tuning: scale set by horizon
Predictions Generic w(z) evolution, few distinctive signals Specific scale-dependent clustering, CMB signatures
Connection to QG None Direct: entanglement entropy → dark energy

🔮 Future Extensions

Extension Status Potential Impact
Full Boltzmann code integration Planned Accurate CMB and LSS predictions
Bayesian evidence calculation In progress Quantitative model comparison
Perturbation theory for EDE Planned Small-scale clustering predictions
Connection to string theory Future work Fundamental quantum gravity basis
Observational forecasts (Euclid, Roman) Planned Testable predictions for future surveys

📞 Contact & Collaboration

Dewan Sajid Islam
Independent Cosmology Researcher
Dhaka, Bangladesh

· Email: johnnykitty06@gmail.com

📜 License

This work is released under the MIT License - see the LICENSE file for details.

🙏 Acknowledgments

· The Planck Collaboration for CMB data
· Pantheon+ team for supernova observations
· BOSS/eBOSS collaborations for BAO measurements
· SHOES team for local H₀ measurements
· Developers of CAMB, CLASS, MontePython
· The open-source scientific Python community

---

<div align="center">

"We are limitless in a limited space."

The universe's acceleration emerges from its quantum information content.

</div>

---

🚨 Troubleshooting

Issue Solution
ModuleNotFoundError Run pip install -r requirements.txt
Memory errors on Android/Pydroid Reduce grid points in scripts (change n_points parameters)
Matplotlib font warnings Ignore - figures will still generate correctly
Slow execution Use make quick-test for basic checks
Notebooks not loading Ensure Jupyter is installed: pip install jupyter

📊 Performance Benchmarks

Operation Typical Runtime Memory Usage
Background solution (z_max=5) 2-3 seconds ~100 MB
Growth functions calculation 1-2 seconds ~50 MB
Power spectrum (500 k-points) 0.5 seconds ~30 MB
CMB spectra (ℓ=2-2500) 0.3 seconds ~20 MB
Full figure generation 10-15 seconds ~200 MB

---

<div align="center">

⚡ Ready to Explore?

```bash
git clone https://github.com/dewansajidislam/Entanglement-Dark-Energy-EDE.git
cd Entanglement-Dark-Energy-EDE
make all
```

Discover how quantum entanglement shapes our accelerating universe.

</div>
