#!/usr/bin/env python3
"""
Setup ASCOT5 markers from SIMPLE initial conditions.

This script:
1. Reads SIMPLE start.dat (VMEC coordinates: s, theta, phi, vnorm, vthnorm)
2. Converts VMEC coordinates to ASCOT cylindrical using input_rhotheta2rz
3. Creates guiding center markers for ASCOT5 HDF5 file
4. Uses same initial conditions as SIMPLE for comparison

Note: VMEC s = s_tor (toroidal flux), and ASCOT rho = sqrt(s_tor).
"""
import numpy as np
from pathlib import Path
import unyt
from a5py import Ascot
from a5py.ascot5io.marker import Marker

# Configuration
START_FILE = "../start.dat"  # From SIMPLE setup
H5_FILE = "ascot.h5"

# Physics parameters (alpha particles)
SPECIES = "alpha"
ENERGY = 3.5e6 * unyt.eV


def read_simple_start(filename):
    """Read SIMPLE start.dat file with VMEC initial conditions."""
    data = np.loadtxt(filename)
    # Columns: s, theta, phi, vnorm, vthnorm
    return {
        's': data[:, 0],
        'theta': data[:, 1],
        'phi': data[:, 2],
        'vnorm': data[:, 3],
        'vthnorm': data[:, 4],
    }


def main():
    print("=" * 60)
    print("ASCOT5 Marker Setup from SIMPLE Initial Conditions")
    print("=" * 60)

    # Check input files
    start_path = Path(START_FILE)
    if not start_path.exists():
        raise FileNotFoundError(f"START file not found: {START_FILE}")

    h5_path = Path(H5_FILE)
    if not h5_path.exists():
        raise FileNotFoundError(f"HDF5 file not found: {H5_FILE}")

    print(f"\n1. Reading SIMPLE start file: {START_FILE}")
    vmec_coords = read_simple_start(START_FILE)
    n_markers = len(vmec_coords['s'])
    print(f"   Found {n_markers} markers")
    print(f"   s range: [{vmec_coords['s'].min():.3f}, {vmec_coords['s'].max():.3f}]")
    print(f"   theta range: [{vmec_coords['theta'].min():.3f}, {vmec_coords['theta'].max():.3f}]")
    print(f"   phi range: [{vmec_coords['phi'].min():.3f}, {vmec_coords['phi'].max():.3f}]")

    print(f"\n2. Opening ASCOT5 HDF5 file: {H5_FILE}")
    a5 = Ascot(H5_FILE)

    print(f"\n3. Initializing magnetic field for coordinate conversion...")
    a5.input_init(bfield=True)

    print(f"\n4. Converting VMEC (s, theta, phi) to cylindrical (R, Z)...")
    print(f"   Note: ASCOT rho = sqrt(s_tor), using {n_markers} markers")

    # Convert s_tor to rho: rho = sqrt(s)
    rho_coords = np.sqrt(vmec_coords['s'])
    theta_coords = vmec_coords['theta'] * unyt.rad
    phi_coords = (vmec_coords['phi'] * unyt.rad).to('deg')

    # Use ASCOT's coordinate conversion (handles VMEC->cylindrical properly)
    r_coords, z_coords = a5.input_rhotheta2rz(
        rho_coords,
        theta_coords,
        phi_coords,
        0*unyt.s  # Time (irrelevant for static equilibrium)
    )

    print(f"   R range: [{r_coords.min():.2f}, {r_coords.max():.2f}]")
    print(f"   Z range: [{z_coords.min():.2f}, {z_coords.max():.2f}]")

    # Free field (no longer needed)
    a5.input_free()

    print(f"\n5. Creating marker dictionary...")
    # Generate markers using ASCOT's template
    mrk = Marker.generate("gc", n=n_markers, species=SPECIES)

    # Set coordinates
    mrk["r"][:] = r_coords
    mrk["z"][:] = z_coords
    mrk["phi"][:] = phi_coords

    # Set energy (constant for all markers in SIMPLE case)
    mrk["energy"][:] = ENERGY

    # Set pitch angle from SIMPLE's vthnorm
    mrk["pitch"][:] = vmec_coords['vthnorm']

    # Set random gyroangle (not specified in SIMPLE)
    mrk["zeta"][:] = 2*np.pi * np.random.rand(n_markers) * unyt.rad

    # Set marker weights (total power / number of markers)
    POWER = 10e6 * unyt.W  # 10 MW
    mrk["weight"][:] = (POWER / (n_markers * mrk["energy"])).to("particles/s")

    print(f"\n6. Writing {n_markers} markers to HDF5...")
    a5.data.create_input(
        "gc",
        **mrk,
        desc=f"GC markers from SIMPLE: {n_markers} particles on s={vmec_coords['s'].mean():.2f} surface"
    )

    print("\n" + "=" * 60)
    print("Marker setup complete!")
    print("=" * 60)
    print(f"\nMarkers ready: {n_markers} GC alpha particles")
    print(f"Energy: {ENERGY.to('MeV')}")
    print(f"Pitch range: [{mrk['pitch'].min():.3f}, {mrk['pitch'].max():.3f}]")
    print(f"\nUse 'make run' to trace orbits")


if __name__ == "__main__":
    main()
