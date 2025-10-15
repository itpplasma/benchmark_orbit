#!/usr/bin/env python3
"""
Run ASCOT5 orbit tracing with same parameters as SIMPLE.

This script:
1. Loads ASCOT5 HDF5 file with markers and field
2. Creates required simulation inputs (plasma, E-field, etc.)
3. Sets simulation parameters matching SIMPLE runtime
4. Runs ascot5_main executable to trace orbits
5. Results are stored in the HDF5 file for extraction

Same runtime as SIMPLE: 1e-3 seconds (1 ms)
"""
import numpy as np
from pathlib import Path
import unyt
import subprocess
import os
from a5py import Ascot
from a5py.ascot5io.options import Opt

# Configuration
H5_FILE = "ascot.h5"

# Simulation parameters (matched to SIMPLE)
SIM_TIME = 1e-3  # seconds (1 ms, same as SIMPLE)


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

    # Check for required inputs
    if not hasattr(a5.data, 'bfield') or a5.data.bfield.active is None:
        raise ValueError("No magnetic field found. Run 'make setup' first.")
    if not hasattr(a5.data, 'marker_gc') or a5.data.marker_gc.active is None:
        raise ValueError("No markers found. Run 'make setup_mark' first.")

    print(f"   Using field: {a5.data.bfield.active.get_desc()}")
    print(f"   Using markers: {a5.data.marker_gc.active.get_desc()}")

    # Create missing inputs
    print(f"\n2. Creating simulation inputs...")

    # Plasma (dummy - not used in GC simulation without collisions)
    a5.data.create_input("plasma flat")

    # Electric field (zero)
    a5.data.create_input("E_TC", exyz=np.array([0,0,0]))

    # Dummy inputs (not used but required)
    a5.data.create_input("N0_1D")
    a5.data.create_input("Boozer")
    a5.data.create_input("MHD_STAT")
    a5.data.create_input("asigma_loc")

    print(f"   Created: plasma, E-field, neutrals, Boozer, MHD, asigma")

    # Set up simulation options
    print(f"\n3. Configuring simulation options...")
    print(f"   Simulation time: {SIM_TIME} s (same as SIMPLE)")

    opt = Opt.get_default()
    opt.update({
        # Guiding-center mode
        "SIM_MODE": 2,
        "ENABLE_ADAPTIVE": 1,

        # End conditions
        "ENDCOND_SIMTIMELIM": 1,
        "ENDCOND_MAX_MILEAGE": SIM_TIME,  # 1 ms physical time (same as SIMPLE)
        "ENDCOND_CPUTIMELIM": 1,
        "ENDCOND_MAX_CPUTIME": 3600.0,  # 1 hour CPU time max
        "ENDCOND_WALLHIT": 0,  # Don't stop on wall hit
        "ENDCOND_RHOLIM": 1,  # Enable rho limits
        "ENDCOND_MAX_RHO": 1.5,  # Stop well outside LCFS

        # Physics
        "ENABLE_ORBIT_FOLLOWING": 1,
        "ENABLE_COULOMB_COLLISIONS": 0,  # No collisions (like SIMPLE)
        "ENABLE_ATOMIC": 0,  # No atomic reactions

        # Orbit recording
        "ENABLE_ORBITWRITE": 1,
        "ORBITWRITE_MODE": 1,  # Time interval mode
        "ORBITWRITE_INTERVAL": 1e-7,  # 0.1 microsecond intervals
        "ORBITWRITE_NPOINT": 20000,  # Max points per marker
    })
    a5.data.create_input("opt", **opt, desc="GC_ORBIT_TRACE")

    print(f"   Options configured")

    # Find ascot5_main executable
    print(f"\n4. Finding ascot5_main executable...")
    ascot5_paths = [
        os.path.expandvars("$CODE/ascot5/build/ascot5_main"),
        "/home/ert/code/ascot5/build/ascot5_main",
        os.path.expandvars("$CODE/ascot5/src/ascot5_main"),
        "/home/ert/code/ascot5/src/ascot5_main",
        "../../../build/ascot5_main",
        "../../../src/ascot5_main"
    ]

    ascot5_main = None
    for path in ascot5_paths:
        if Path(path).exists():
            ascot5_main = path
            break

    if ascot5_main is None:
        raise FileNotFoundError(
            "Could not find ascot5_main executable.\n"
            "Please build ASCOT5 first:\n"
            "  cd /home/ert/code/ascot5\n"
            "  mkdir -p build && cd build\n"
            "  cmake .. && make -j"
        )

    print(f"   Using: {ascot5_main}")

    # Run simulation
    print(f"\n5. Running orbit tracing simulation...")
    print(f"   This will take several minutes for 10,000 markers...")
    print(f"   (Output is live from ascot5_main)")

    result = subprocess.run(
        [ascot5_main, f"--d=ORBIT_TRACE"],
        cwd=Path.cwd(),
        text=True
    )

    if result.returncode != 0:
        print(f"\nSimulation failed with return code {result.returncode}")
        return

    print(f"\n6. Simulation completed successfully!")
    print(f"   Results saved in: {H5_FILE}")

    print("\n" + "=" * 60)
    print("Orbit tracing complete!")
    print("=" * 60)
    print(f"\nUse 'make extract' to extract orbit data to NetCDF")


if __name__ == "__main__":
    main()
