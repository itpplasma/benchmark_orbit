#!/usr/bin/env python3
"""
Bundle fort.90XXX orbit output files into a single NetCDF file.

Data columns from SIMPLE orbit output:
  1: time (t)
  2-4: reference frame coordinates (s, theta, phi)
  5: momentum magnitude (p_abs)
  6: parallel velocity fraction (v_par/v)
"""

import glob
import numpy as np
import xarray as xr
import sys
from pathlib import Path


def parse_orbit_file(filepath):
    """Parse a single fort.90XXX orbit file."""
    try:
        data = np.loadtxt(filepath)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        return data
    except Exception as e:
        print(f"Warning: Failed to parse {filepath}: {e}", file=sys.stderr)
        return None


def main():
    # Find all fort.90XXX files
    fort_files = sorted(glob.glob("fort.9[0-9][0-9][0-9][0-9]"))
    if not fort_files:
        print("Error: No fort.9XXXX files found in current directory", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(fort_files)} orbit files")

    # Extract particle IDs from filenames
    orbit_data = {}
    max_timesteps = 0

    for filepath in fort_files:
        # Extract particle ID from fort.90XXX (e.g., "fort.90001" -> 1)
        # Remove "fort.9" prefix and convert to int (removes leading zeros)
        particle_id = int(filepath.split('.')[-1][1:])

        data = parse_orbit_file(filepath)
        if data is not None and len(data) > 0:
            orbit_data[particle_id] = data
            max_timesteps = max(max_timesteps, len(data))
            print(f"  Particle {particle_id}: {len(data)} timesteps")

    if not orbit_data:
        print("Error: No valid data found in orbit files", file=sys.stderr)
        sys.exit(1)

    # Create arrays padded to max_timesteps
    n_particles = len(orbit_data)
    particle_ids = sorted(orbit_data.keys())

    # Initialize arrays with NaN (will be masked)
    time_data = np.full((n_particles, max_timesteps), np.nan)
    s_data = np.full((n_particles, max_timesteps), np.nan)
    theta_data = np.full((n_particles, max_timesteps), np.nan)
    phi_data = np.full((n_particles, max_timesteps), np.nan)
    p_abs_data = np.full((n_particles, max_timesteps), np.nan)
    v_par_data = np.full((n_particles, max_timesteps), np.nan)

    # Fill arrays
    for i, pid in enumerate(particle_ids):
        data = orbit_data[pid]
        n_steps = len(data)
        time_data[i, :n_steps] = data[:, 0]
        s_data[i, :n_steps] = data[:, 1]
        theta_data[i, :n_steps] = data[:, 2]
        phi_data[i, :n_steps] = data[:, 3]
        p_abs_data[i, :n_steps] = data[:, 4]
        v_par_data[i, :n_steps] = data[:, 5]

    # Create xarray Dataset
    ds = xr.Dataset(
        {
            "time": (["particle", "timestep"], time_data),
            "s": (["particle", "timestep"], s_data),
            "theta": (["particle", "timestep"], theta_data),
            "phi": (["particle", "timestep"], phi_data),
            "p_abs": (["particle", "timestep"], p_abs_data),
            "v_par": (["particle", "timestep"], v_par_data),
        },
        coords={
            "particle": particle_ids,
            "timestep": np.arange(max_timesteps),
        },
        attrs={
            "description": "SIMPLE orbit tracer output (macrostep mode)",
            "source": "fort.90XXX files",
            "n_particles": n_particles,
            "n_timesteps": max_timesteps,
        },
    )

    # Save to NetCDF
    output_file = "orbits.nc"
    ds.to_netcdf(output_file)
    print(f"\nSuccessfully created {output_file}")
    print(f"  Particles: {n_particles}")
    print(f"  Max timesteps: {max_timesteps}")
    print(f"  File size: {Path(output_file).stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
