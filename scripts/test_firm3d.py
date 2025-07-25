#!/usr/bin/env python3
"""Test firm3d import and basic functionality."""

import sys
import os

# Add firm3d to path
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "codes", "firm3d", "src"))

try:
    print("Testing firm3d imports...")
    from simsopt.field.boozermagneticfield import BoozerRadialInterpolant
    print("✓ BoozerRadialInterpolant imported")
    
    # Try loading the boozmn file
    boozmn_file = "booz_xform/boozmn_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"
    print(f"\nTesting BoozerRadialInterpolant with {boozmn_file}")
    
    try:
        bri = BoozerRadialInterpolant(boozmn_file, order=3, no_K=True)
        print("✓ BoozerRadialInterpolant created successfully")
        print(f"  nfp = {bri.nfp}")
        print(f"  psi0 = {bri.psi0}")
    except Exception as e:
        print(f"✗ Error creating BoozerRadialInterpolant: {e}")
        
except ImportError as e:
    print(f"✗ Import error: {e}")
    print("\nTrying to understand the issue...")
    print(f"Python path includes:")
    for p in sys.path:
        if "firm3d" in p:
            print(f"  {p}")