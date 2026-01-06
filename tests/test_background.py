"""
Test background evolution calculations.
"""

import numpy as np
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.background_evolution import EDECosmology

def test_initialization():
    """Test cosmology initialization."""
    ede = EDECosmology(model_name='EDE')
    lcdm = EDECosmology(model_name='LCDM')
    
    assert ede.H0 == 72.1
    assert lcdm.H0 == 67.36
    print("✓ Initialization test passed")

def test_background_solution():
    """Test background solution."""
    ede = EDECosmology(model_name='EDE')
    bg = ede.solve_background(z_max=2.0)
    
    # Check shapes
    assert len(bg['z']) > 0
    assert len(bg['H']) == len(bg['z'])
    assert len(bg['w_ent']) == len(bg['z'])
    
    # Check physical bounds
    assert np.all(bg['H'] > 0)
    assert np.all(bg['w_ent'] <= 0)  # w should be negative
    
    print("✓ Background solution test passed")

def test_w_ent_evolution():
    """Test equation of state evolution."""
    ede = EDECosmology(model_name='EDE')
    bg = ede.solve_background(z_max=2.0)
    
    # w(z=0) should be around -0.95
    w0 = bg['w_ent'][0]
    assert -1.0 < w0 < -0.9
    
    print(f"✓ w(z=0) = {w0:.3f} (expected ~ -0.954)")

def run_all_tests():
    """Run all tests."""
    print("Running background evolution tests...")
    print("-" * 40)
    
    test_initialization()
    test_background_solution()
    test_w_ent_evolution()
    
    print("-" * 40)
    print("All tests passed! ✓")

if __name__ == "__main__":
    run_all_tests()
