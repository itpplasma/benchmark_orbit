#!/usr/bin/env python3
"""
Setup ASCOT5 markers from SIMPLE initial conditions.

This script:
1. Reads SIMPLE start.dat (VMEC coordinates: s, theta, phi, vnorm, vthnorm)
2. Converts VMEC coordinates to ASCOT cylindrical coordinates (r, phi, z)
3. Creates guiding center markers for ASCOT5 HDF5 file
4. Uses same initial conditions as SIMPLE for comparison

VMEC to cylindrical conversion uses VMEC equilibrium file.
"""
import numpy as np
from pathlib import Path
import xarray as xr
import unyt
from a5py import Ascot

# Configuration
START_FILE = "../start.dat"  # From SIMPLE setup
VMEC_FILE = "wout.nc"
H5_FILE = "ascot.h5"

# Physics parameters
MASS_ALPHA = 4.0 * unyt.amu  # Alpha particle mass
CHARGE_ALPHA = 2 * unyt.e    # Alpha particle charge
ANUM_ALPHA = 4
ZNUM_ALPHA = 2


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


def vmec_to_cylindrical(vmec_file, s, theta, phi):
    """
    Convert VMEC flux coordinates to cylindrical (R, phi, Z).

    Uses Fourier expansions from VMEC output.
    s: normalized poloidal flux (0 at axis, 1 at LCFS)
    theta: poloidal angle [rad]
    phi: toroidal angle [rad]
    """
    # Load VMEC file
    vmec_ds = xr.open_dataset(vmec_file)

    # Get Fourier coefficients
    rmnc = vmec_ds['rmnc'].values  # (ns, mn)
    zmns = vmec_ds['zmns'].values  # (ns, mn)
    xm = vmec_ds['xm'].values      # (mn,)
    xn = vmec_ds['xn'].values      # (mn,)

    # Get radial grid
    s_vmec = vmec_ds['s_full'].values if 's_full' in vmec_ds else np.linspace(0, 1, len(rmnc))

    # Interpolate to requested s values (simple linear for now)
    r_coeffs = np.interp(s, s_vmec, rmnc, axis=0)
    z_coeffs = np.interp(s, s_vmec, zmns, axis=0)

    # Evaluate Fourier series: R = sum(rmnc * cos(m*theta - n*phi))
    n_particles = len(s)
    r_cyl = np.zeros(n_particles)
    z_cyl = np.zeros(n_particles)

    for i in range(len(xm)):
        m = int(xm[i])
        n = int(xn[i])
        angle = m * theta - n * phi

        r_cyl += r_coeffs[:, i] * np.cos(angle)
        z_cyl += z_coeffs[:, i] * np.sin(angle)

    vmec_ds.close()

    return r_cyl, z_cyl


def vmec_energy_from_vnorm(vnorm, mass):
    """
    Convert SIMPLE velocity normalization to ASCOT energy.

    SIMPLE uses normalized velocity vnorm where 1.0 ≈ thermal velocity.
    Convert to kinetic energy for ASCOT.
    """
    # For quasi-isotropic distribution, assume thermal energy ~3 keV
    e_thermal = 3.5e3 * unyt.eV
    energy = 0.5 * mass * (vnorm**2 * (2 * e_thermal / mass))
    return energy.to('eV')


def main():
    print("=" * 60)
    print("ASCOT5 Marker Setup from SIMPLE Initial Conditions")
    print("=" * 60)

    # Check input files
    start_path = Path(START_FILE)
    if not start_path.exists():
        raise FileNotFoundError(f"START file not found: {START_FILE}")

    vmec_path = Path(VMEC_FILE)
    if not vmec_path.exists():
        raise FileNotFoundError(f"VMEC file not found: {VMEC_FILE}")

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

    print(f"\n2. Converting VMEC to cylindrical coordinates...")
    print(f"   Using VMEC file: {VMEC_FILE}")
    r_cyl, z_cyl = vmec_to_cylindrical(
        VMEC_FILE,
        vmec_coords['s'],
        vmec_coords['theta'],
        vmec_coords['phi']
    )
    print(f"   R range: [{r_cyl.min():.2f}, {r_cyl.max():.2f}] m")
    print(f"   Z range: [{z_cyl.min():.2f}, {z_cyl.max():.2f}] m")

    print(f"\n3. Computing marker energies and pitch angles...")
    energy = vmec_energy_from_vnorm(vmec_coords['vnorm'], MASS_ALPHA)
    energy = np.atleast_1d(energy)

    # Compute pitch angle from vthnorm
    # In SIMPLE: vthnorm = v_parallel / v_total
    pitch = vmec_coords['vthnorm']  # cos(alpha) where alpha is pitch angle
    print(f"   Energy range: [{energy.min():.1f}, {energy.max():.1f}]")
    print(f"   Pitch range: [{pitch.min():.3f}, {pitch.max():.3f}]")

    print(f"\n4. Opening ASCOT5 HDF5 file: {H5_FILE}")
    a5 = Ascot(H5_FILE)

    print(f"\n5. Writing markers to HDF5...")
    # Create marker data dictionary
    markers = {
        'n': n_markers,
        'ids': np.arange(1, n_markers + 1, dtype=np.int64),
        'r': r_cyl * unyt.m,
        'phi': vmec_coords['phi'] * unyt.rad,
        'z': z_cyl * unyt.m,
        'energy': energy,
        'pitch': pitch,
        'zeta': np.zeros(n_markers) * unyt.rad,
        'mass': MASS_ALPHA * np.ones(n_markers),
        'charge': CHARGE_ALPHA * np.ones(n_markers, dtype=np.int16),
        'anum': ANUM_ALPHA * np.ones(n_markers, dtype=np.int16),
        'znum': ZNUM_ALPHA * np.ones(n_markers, dtype=np.int16),
        'weight': 1e18 * np.ones(n_markers) * unyt.particles / unyt.s,
        'time': np.zeros(n_markers) * unyt.s,
    }

    # Write markers to HDF5
    a5.data.create_input(
        'gc',  # Guiding center markers
        **markers,
        activate=True,
        desc=f"GC markers from SIMPLE: {n_markers} particles on s={vmec_coords['s'].mean():.2f} surface"
    )

    print(f"\n   Markers written to: marker_gc (QID: {a5.data.marker_gc.active.get_qid()})")

    print("\n" + "=" * 60)
    print("Marker setup complete!")
    print("=" * 60)
    print(f"\nMarkers ready: {n_markers} GC particles")
    print(f"Use 'make run' to trace orbits")


if __name__ == "__main__":
    main()
