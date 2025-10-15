#!/usr/bin/env python3
"""
Extract orbit trajectories from ASCOT5 results and save as NetCDF.

This script:
1. Reads orbits.h5 from ASCOT5 simulation
2. Extracts trajectory data for each marker
3. Saves to orbits_ascot.nc in xarray format for easy plotting

Output format matches SIMPLE orbit output for comparison.
"""
import numpy as np
from pathlib import Path
import xarray as xr
from a5py import Ascot

# Configuration
RESULTS_FILE = "orbits.h5"
OUTPUT_NC = "orbits_ascot.nc"


def main():
    print("=" * 60)
    print("Extracting ASCOT5 Orbit Trajectories")
    print("=" * 60)

    # Check results file
    results_path = Path(RESULTS_FILE)
    if not results_path.exists():
        raise FileNotFoundError(f"Results file not found: {RESULTS_FILE}")

    print(f"\n1. Opening ASCOT5 results: {RESULTS_FILE}")
    a5 = Ascot(RESULTS_FILE)

    # Get simulation data
    print(f"\n2. Accessing orbit data...")
    try:
        inistate = a5.results.orbits.initial
        endstate = a5.results.orbits.final
        marker_ids = a5.results.orbits.list()

        print(f"   Found {len(marker_ids)} markers")
    except Exception as e:
        print(f"   Error accessing orbit data: {e}")
        print(f"   This may indicate the simulation didn't complete properly")
        raise

    print(f"\n3. Extracting trajectory data...")

    # Collect all orbit data
    orbit_data = {}
    max_timesteps = 0

    for marker_id in marker_ids:
        try:
            # Extract trajectory quantities
            time = a5.results.orbits.get(
                inistate, endstate, 'mileage', ids=marker_id
            )
            r = a5.results.orbits.get(
                inistate, endstate, 'r', ids=marker_id
            )
            phi = a5.results.orbits.get(
                inistate, endstate, 'phi', ids=marker_id
            )
            z = a5.results.orbits.get(
                inistate, endstate, 'z', ids=marker_id
            )

            if time is not None and len(time) > 0:
                orbit_data[marker_id] = {
                    'time': time,
                    'r': r,
                    'phi': phi,
                    'z': z,
                }
                max_timesteps = max(max_timesteps, len(time))
                if len(marker_ids) <= 10 or marker_id % 100 == 0:
                    print(f"   Marker {marker_id}: {len(time)} timesteps")
        except Exception as e:
            if len(marker_ids) <= 10:
                print(f"   Marker {marker_id}: No data - {e}")

    if not orbit_data:
        raise RuntimeError("No orbit data extracted")

    print(f"   Total markers with data: {len(orbit_data)}")
    print(f"   Max timesteps: {max_timesteps}")

    print(f"\n4. Creating NetCDF dataset...")

    # Create arrays padded to max_timesteps
    n_particles = len(orbit_data)
    marker_list = sorted(orbit_data.keys())

    # Initialize arrays with NaN
    time_data = np.full((n_particles, max_timesteps), np.nan)
    r_data = np.full((n_particles, max_timesteps), np.nan)
    phi_data = np.full((n_particles, max_timesteps), np.nan)
    z_data = np.full((n_particles, max_timesteps), np.nan)

    # Fill arrays
    for i, marker_id in enumerate(marker_list):
        data = orbit_data[marker_id]
        n_steps = len(data['time'])
        time_data[i, :n_steps] = data['time']
        r_data[i, :n_steps] = data['r']
        phi_data[i, :n_steps] = data['phi']
        z_data[i, :n_steps] = data['z']

    # Create xarray Dataset
    ds = xr.Dataset(
        {
            "time": (["particle", "timestep"], time_data),
            "r": (["particle", "timestep"], r_data),
            "phi": (["particle", "timestep"], phi_data),
            "z": (["particle", "timestep"], z_data),
        },
        coords={
            "particle": marker_list,
            "timestep": np.arange(max_timesteps),
        },
        attrs={
            "description": "ASCOT5 orbit tracer output (GC mode)",
            "source": "orbits.h5",
            "n_particles": n_particles,
            "n_timesteps": max_timesteps,
            "coordinate_system": "cylindrical (R, phi, Z)",
        },
    )

    # Save to NetCDF
    print(f"\n5. Saving to: {OUTPUT_NC}")
    ds.to_netcdf(OUTPUT_NC)
    print(f"   NetCDF saved successfully!")
    print(f"   File size: {Path(OUTPUT_NC).stat().st_size / 1e6:.1f} MB")

    print("\n" + "=" * 60)
    print("Extraction complete!")
    print("=" * 60)
    print(f"\nOrbit data ready: {OUTPUT_NC}")
    print(f"Particles: {n_particles}")
    print(f"Timesteps: {max_timesteps}")
    print(f"Use 'make plot_orbit ORBIT=<id>' to visualize")


if __name__ == "__main__":
    main()
