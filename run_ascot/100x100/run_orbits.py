#!/usr/bin/env python3
"""
Run ASCOT5 orbit tracing with same parameters as SIMPLE.

This script:
1. Loads ASCOT5 HDF5 file with markers and field
2. Sets simulation parameters matching SIMPLE runtime
3. Traces orbits and records trajectory data
4. Saves results to orbits.h5 for extraction

Same runtime as SIMPLE: 1e-3 seconds (10001 macrosteps at ~1e-4 s per step)
"""
import numpy as np
from pathlib import Path
import unyt
from a5py import Ascot

# Configuration
H5_FILE = "ascot.h5"
RESULTS_FILE = "orbits.h5"

# Simulation parameters (matched to SIMPLE)
SIM_TIME = 1e-3 * unyt.s          # Total simulation time (1 ms, same as SIMPLE)
ORBITWRITE_INTERVAL = 1e-7 * unyt.s  # Write orbit every 0.1 microseconds
N_ORBIT_POINTS = 10001            # Should match SIMPLE's ntimstep


def main():
    print("=" * 60)
    print("ASCOT5 Orbit Tracing")
    print("=" * 60)

    # Check HDF5 file
    h5_path = Path(H5_FILE)
    if not h5_path.exists():
        raise FileNotFoundError(f"HDF5 file not found: {H5_FILE}")

    print(f"\n1. Opening ASCOT5 HDF5 file: {H5_FILE}")
    a5 = Ascot(H5_FILE)

    print(f"\n2. Setting simulation options...")
    print(f"   Simulation time: {SIM_TIME}")
    print(f"   Orbit write interval: {ORBITWRITE_INTERVAL}")
    print(f"   Expected orbit points: {N_ORBIT_POINTS}")

    # Initialize options for guiding center simulation with orbit recording
    opt = {
        'SIM_MODE': 2,                    # Guiding center mode
        'ENABLE_ORBITWRITE': 1,           # Enable orbit diagnostics
        'ORBITWRITE_MODE': 1,             # Time interval mode
        'ORBITWRITE_NPOINT': N_ORBIT_POINTS,  # Max points per marker
        'ORBITWRITE_INTERVAL': ORBITWRITE_INTERVAL,  # Write interval
        'ENDCOND_SIMTIMELIM': 1,          # Stop at sim time limit
        'ENDCOND_MAXMILEAGE': SIM_TIME,   # Max simulation time
    }

    try:
        a5.simulation_initoptions(**opt)
        print("   Options initialized")
    except Exception as e:
        print(f"   Options initialization error: {e}")
        print("   Attempting with minimal options...")
        opt_minimal = {
            'SIM_MODE': 2,
            'ENABLE_ORBITWRITE': 1,
            'ENDCOND_SIMTIMELIM': 1,
            'ENDCOND_MAXMILEAGE': SIM_TIME,
        }
        a5.simulation_initoptions(**opt_minimal)

    print(f"\n3. Initializing simulation inputs...")
    a5.simulation_initinputs(bfield=True, plasma=False, wall=False,
                             efield=False, neutral=False)
    print("   Inputs initialized")

    print(f"\n4. Reading markers from HDF5...")
    a5.simulation_initmarkers()
    print("   Markers initialized")

    print(f"\n5. Running orbit tracing simulation...")
    print(f"   This may take several minutes...")
    run = a5.simulation_run(printsummary=True)
    print(f"   Simulation completed!")

    print(f"\n6. Processing results...")
    print(f"   Results saved to: orbits.h5")

    # Export results to separate file
    a5.simulation_exportdata(RESULTS_FILE)
    print(f"   Results exported to: {RESULTS_FILE}")

    # Free resources
    print(f"\n7. Cleaning up...")
    a5.simulation_free(inputs=True, markers=True, diagnostics=True)

    print("\n" + "=" * 60)
    print("Orbit tracing complete!")
    print("=" * 60)
    print(f"\nResults saved to: {RESULTS_FILE}")
    print(f"Use 'make extract' to extract orbit data")


if __name__ == "__main__":
    main()
