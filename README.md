Entanglement Dark Energy (EDE) Model

<div align="center">

"We are limitless in a limited space."
— Dewan Sajid Islam

</div>

---

QUICK START

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

---

REPOSITORY CONTENTS

data/ - Cosmological parameters, MCMC results, observational data

· parameters.json - EDE and ΛCDM parameters
· best_fit_parameters.json - MCMC best-fit results
· chains/ - MCMC chain files (generated)

scripts/ - Core computational engine

· background_evolution.py - Solve EDE background equations
· growth_functions.py - Calculate growth factors
· power_spectrum.py - Compute matter power spectrum
· cmb_spectra.py - Calculate CMB spectra
· utils.py - Utility functions
· generate_figures.py - Generate all paper figures

notebooks/ - Interactive Jupyter notebooks

· 01_Background_Evolution.ipynb - Hubble parameter and w(z)
· 02_Growth_Functions.ipynb - Growth factor and fσ₈(z)
· 03_Power_Spectrum.ipynb - Matter power spectrum
· 04_CMB_Spectra.ipynb - CMB temperature and polarization

figures/ - Generated publication figures

· figure1.pdf/.png - Hubble parameter and equation of state
· figure2.pdf/.png - Growth functions
· figure3.pdf/.png - Matter power spectrum
· figure4.pdf/.png - CMB spectra
· corner_plot.pdf/.png - Parameter constraints

latex/ - Complete LaTeX manuscript

· main.tex - Manuscript source
· references.bib - Bibliography

tests/ - Unit tests for validation

· test_background.py - Background evolution tests
· test_growth.py - Growth function tests

---

KEY RESULTS: EDE vs ΛCDM

Cosmological Parameters (68% Confidence Level)

Hubble Constant

· EDE Model: H₀ = 72.1 ± 0.9 km/s/Mpc
· ΛCDM Model: H₀ = 67.36 ± 0.54 km/s/Mpc
· Improvement: Hubble tension naturally resolved

Matter Density

· EDE Model: Ωₘ = 0.294 ± 0.008
· ΛCDM Model: Ωₘ = 0.3156 ± 0.0074

Structure Parameters

· EDE Model: σ₈ = 0.829 ± 0.008, S₈ = 0.812 ± 0.012
· ΛCDM Model: σ₈ = 0.811 ± 0.006, S₈ = 0.832 ± 0.013
· Improvement: Reduces S₈ tension with weak lensing

Dark Energy Equation of State

· EDE Model: w₀ = -0.954 ± 0.020 (dynamical)
· ΛCDM Model: w = -1 (fixed, cosmological constant)

Primordial Parameters

· Spectral Index: n_s = 0.965 ± 0.004 (both models)
· Optical Depth: τ = 0.054 ± 0.007 (both models)

Derived Quantities

Age of Universe

· EDE: 13.65 ± 0.05 Gyr
· ΛCDM: 13.80 ± 0.02 Gyr

Sound Horizon Scale

· Both models: r_s = 147.0 ± 0.3 Mpc

Angular Sound Horizon

· Both models: θ* = 0.010413 ± 0.000005 rad

Redshift of Matter-Radiation Equality

· EDE: z_eq = 3400 ± 50
· ΛCDM: z_eq = 3402 ± 26

---

THE EDE MODEL: CORE EQUATIONS

1. Entanglement Energy Density

```
ρ_ent = (3/(8πG)) * [1/R_h² - (Ṙ_h)/(H R_h³)]
```

Where R_h is the future event horizon

2. Equation of State

```
w_ent = -1 + (2/3)*(Ṙ_h/(H R_h))*[1 - 1/(H R_h) + (R̈_h)/(2H Ṙ_h) - Ṙ_h/(2H R_h)]/[1 - Ṙ_h/(H R_h)]
```

3. Future Event Horizon

```
R_h(t) = a(t) ∫_t^∞ dt'/a(t')
```

---

HOW TO USE THIS REPOSITORY

Command Line Interface

Installation

```bash
make install        # Install all dependencies
```

Testing

```bash
make test           # Run unit tests
make quick-test     # Quick functionality test
```

Figure Generation

```bash
make figures        # Generate all paper figures
```

Manuscript Compilation

```bash
make paper          # Compile LaTeX manuscript
```

Maintenance

```bash
make clean          # Remove generated files
make all            # Install, test, and generate figures
```

Python API Examples

Basic Usage

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

Growth Functions

```python
from scripts.growth_functions import GrowthCalculator

ede_growth = GrowthCalculator(ede)
growth_data = ede_growth.solve_growth_equation(z_max=5.0)
fsigma8 = ede_growth.get_fsigma8(growth_data['z'], sigma8_0=0.829)
```

Power Spectrum

```python
from scripts.power_spectrum import PowerSpectrum

ede_ps = PowerSpectrum(ede)
k = np.logspace(-4, 1, 500)  # h/Mpc
Pk = ede_ps.linear_power_spectrum(k, z=0.0)
sigma8 = ede_ps.sigma8(z=0.0)
```

---

TESTABLE PREDICTIONS

Scale-Dependent Matter Clustering

Large Scales (k < 0.001 h/Mpc)

· Prediction: Suppressed by 2-3% compared to ΛCDM
· Observable: CMB low-ℓ, future galaxy surveys
· Cause: Entanglement correction to primordial power spectrum

Intermediate Scales (0.001 < k < 0.01 h/Mpc)

· Prediction: Transition region
· Observable: Large-scale structure surveys
· Cause: Scale-dependent evolution of entanglement entropy

Small Scales (k > 0.01 h/Mpc)

· Prediction: Identical to ΛCDM
· Observable: Current clustering constraints
· Cause: Standard CDM behavior at small scales

CMB Anomalies

Low-ℓ Temperature Power

· Prediction: Suppressed quadrupole
· Detection Method: Planck, future CMB missions
· Significance: Explains observed low quadrupole

TE Correlation (ℓ < 30)

· Prediction: Enhanced signal at ℓ ≈ 25
· Detection Method: CMB polarization measurements
· Significance: Unique EDE signature

ISW Effect

· Prediction: Distinct correlation pattern
· Detection Method: CMB × LSS cross-correlation
· Significance: Testable with current surveys

Growth of Structure

Redshift Evolution of fσ₈(z)

· z = 0.0: fσ₈ = 0.471 (EDE) vs 0.462 (ΛCDM) → +2.0% enhancement
· z = 0.5: fσ₈ = 0.478 (EDE) vs 0.470 (ΛCDM) → +1.7% enhancement
· z = 1.0: fσ₈ = 0.485 (EDE) vs 0.479 (ΛCDM) → +1.3% enhancement
· z = 2.0: fσ₈ = 0.492 (EDE) vs 0.488 (ΛCDM) → +0.8% enhancement

---

KEY FEATURES OF THE EDE MODEL

Origin of Dark Energy

· EDE Implementation: Quantum entanglement entropy
· Advantage: First-principles derivation, no new fields required

Hubble Tension Resolution

· EDE Implementation: Natural prediction H₀ = 72.1 ± 0.9 km/s/Mpc
· Advantage: No fine-tuning, emerges from formalism

S₈ Tension Reduction

· EDE Implementation: S₈ = 0.812 ± 0.012
· Advantage: Better agreement with weak lensing surveys

Dynamical Equation of State

· EDE Implementation: Derived from entanglement evolution
· Advantage: Testable prediction for future surveys

Early Universe Compatibility

· EDE Implementation: Negligible at high redshift
· Advantage: Preserves BBN and recombination successes

Theoretical Basis

· EDE Implementation: Wheeler-DeWitt + Holographic principle
· Advantage: Connects quantum gravity to observational cosmology

---

COMPLETE FILE REFERENCE

Essential Files for Reproduction

Background Evolution

· File: scripts/background_evolution.py
· Purpose: Solve EDE background equations
· Usage: python -c "from scripts.background_evolution import EDECosmology; ede=EDECosmology(); bg=ede.solve_background()"

Growth Functions

· File: scripts/growth_functions.py
· Purpose: Calculate growth factors and rates
· Usage: python -c "from scripts.growth_functions import GrowthCalculator; from scripts.background_evolution import EDECosmology; ede=EDECosmology(); growth=GrowthCalculator(ede); data=growth.solve_growth_equation()"

Power Spectrum

· File: scripts/power_spectrum.py
· Purpose: Compute matter power spectrum
· Usage: python -c "from scripts.power_spectrum import PowerSpectrum; from scripts.background_evolution import EDECosmology; ede=EDECosmology(); ps=PowerSpectrum(ede); k=[0.01,0.1,1.0]; Pk=ps.linear_power_spectrum(k); print(Pk)"

CMB Spectra

· File: scripts/cmb_spectra.py
· Purpose: Calculate CMB temperature and polarization spectra
· Usage: python -c "from scripts.cmb_spectra import CMBSpectra; from scripts.background_evolution import EDECosmology; ede=EDECosmology(); cmb=CMBSpectra(ede); ell=range(2,2501); D_ell=cmb.temperature_spectrum(ell)"

Figure Generation

· File: scripts/generate_figures.py
· Purpose: Generate all paper figures
· Usage: python scripts/generate_figures.py

Parameter Files

· File: data/parameters.json
· Purpose: Cosmological parameters for both models
· Usage: import json; params = json.load(open('data/parameters.json'))

Notebook Workflow

1. Start with: notebooks/01_Background_Evolution.ipynb - Understand background evolution
2. Then explore: notebooks/02_Growth_Functions.ipynb - Study growth of structure
3. Continue to: notebooks/03_Power_Spectrum.ipynb - Analyze matter clustering
4. Finish with: notebooks/04_CMB_Spectra.ipynb - Examine CMB predictions

---

VALIDATION & TESTING

Background Evolution Tests

· Command: python -m pytest tests/test_background.py -v
· Expected: All tests pass
· Checks: Parameter initialization, background solution, physical bounds

Growth Function Validation

· Command: python -c "from scripts.growth_functions import GrowthCalculator; from scripts.background_evolution import EDECosmology; test_growth()"
· Expected: D(z=0)=1.0, f(z=0)≈0.5
· Checks: Growth normalization, rate calculations

Power Spectrum Checks

· Command: python -c "from scripts.power_spectrum import PowerSpectrum; from scripts.background_evolution import EDECosmology; test_power_spectrum()"
· Expected: P(k) positive for all k, σ₈≈0.8
· Checks: Power spectrum positivity, σ₈ calculation

Figure Generation

· Command: python scripts/generate_figures.py
· Expected: 4 main figures + corner plot generated
· Checks: All figures created without errors

---

WHY EDE IS DIFFERENT

Theoretical Basis

· Traditional Models: Ad-hoc scalar fields, modified gravity
· Entanglement Dark Energy: Emergent from quantum gravity principles

Hubble Tension Approach

· Traditional Models: Require early dark energy, extra relics
· Entanglement Dark Energy: Natural outcome of the formalism

Fine-tuning Problem

· Traditional Models: Cosmological constant requires 10¹²⁰ fine-tuning
· Entanglement Dark Energy: No fine-tuning - scale set by cosmic horizon

Predictive Power

· Traditional Models: Generic w(z) evolution, few distinctive signals
· Entanglement Dark Energy: Specific scale-dependent clustering, unique CMB signatures

Connection to Quantum Gravity

· Traditional Models: None
· Entanglement Dark Energy: Direct - entanglement entropy generates dark energy

---

FUTURE EXTENSIONS

Full Boltzmann Code Integration

· Status: Planned
· Impact: Accurate CMB and LSS predictions

Bayesian Evidence Calculation

· Status: In progress
· Impact: Quantitative model comparison with ΛCDM

Perturbation Theory for EDE

· Status: Planned
· Impact: Small-scale clustering predictions

Connection to String Theory

· Status: Future work
· Impact: Fundamental quantum gravity basis

Observational Forecasts

· Status: Planned (Euclid, Roman, CMB-S4)
· Impact: Testable predictions for future surveys

---

CONTACT

Dewan Sajid Islam
Independent Cosmology Researcher
Dhaka, Bangladesh

· Email: johnnykitty06@gmail.com 
· OrcID: https://orcid.org/0009-0008-8784-8480 
· DOI:  https://doi.org/10.5281/zenodo.18166694

---

LICENSE

This work is released under the MIT License - see the LICENSE file for details.

---

ACKNOWLEDGMENTS

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

Discover how quantum entanglement shapes our accelerating universe.

</div>

---

TROUBLESHOOTING

ModuleNotFoundError

· Solution: Run pip install -r requirements.txt

Memory Errors on Android/Pydroid

· Solution: Reduce grid points in scripts (change n_points parameters)

Matplotlib Font Warnings

· Solution: Ignore - figures will still generate correctly

Slow Execution

· Solution: Use make quick-test for basic checks

Notebooks Not Loading

· Solution: Ensure Jupyter is installed: pip install jupyter

---

PERFORMANCE BENCHMARKS

Background Solution (z_max=5)

· Runtime: 2-3 seconds
· Memory: ~100 MB

Growth Functions Calculation

· Runtime: 1-2 seconds
· Memory: ~50 MB

Power Spectrum (500 k-points)

· Runtime: 0.5 seconds
· Memory: ~30 MB

CMB Spectra (ℓ=2-2500)

· Runtime: 0.3 seconds
· Memory: ~20 MB

Full Figure Generation

· Runtime: 10-15 seconds
· Memory: ~200 MB
