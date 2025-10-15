#!/usr/bin/env python3
"""
Plot ASCOT5 orbit trajectory from orbits_ascot.nc.

Compares with SIMPLE orbits in cylindrical coordinates.

Usage: python plot_orbit_ascot.py <particle_id>

Creates four subplots showing orbit in R-Z plane and coordinate evolution.
"""
import sys
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np


def main():
    if len(sys.argv) != 2:
        print("Usage: python plot_orbit_ascot.py <particle_id>")
        sys.exit(1)

    try:
        particle_id = int(sys.argv[1])
    except ValueError:
        print(f"Error: particle_id must be an integer, got '{sys.argv[1]}'")
        sys.exit(1)

    # Load NetCDF file
    try:
        ds = xr.open_dataset("orbits_ascot.nc")
    except FileNotFoundError:
        print("Error: orbits_ascot.nc not found. Run 'make run' first.")
        sys.exit(1)

    # Check if particle exists
    if particle_id not in ds.particle.values:
        print(f"Error: Particle {particle_id} not found in orbits_ascot.nc")
        print(f"Available particles: {sorted(ds.particle.values)}")
        sys.exit(1)

    # Extract orbit data
    orbit = ds.sel(particle=particle_id)
    time = orbit["time"].values
    r = orbit["r"].values
    phi = orbit["phi"].values
    z = orbit["z"].values

    # Remove NaN values
    valid = ~np.isnan(time)
    time = time[valid]
    r = r[valid]
    phi = phi[valid]
    z = z[valid]

    if len(time) == 0:
        print(f"Error: No valid data for particle {particle_id}")
        sys.exit(1)

    # Create figure with 4 subplots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f"ASCOT5 Orbit for Marker {particle_id} (GC mode)", fontsize=14, fontweight="bold")

    # Subplot 1: R-Z plane (poloidal cross-section)
    ax = axes[0, 0]
    ax.plot(r, z, "b,", markersize=2)
    ax.scatter(r[0], z[0], color="green", s=100, marker="o", label="Start", zorder=5)
    ax.scatter(r[-1], z[-1], color="red", s=100, marker="x", label="End", zorder=5)
    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    ax.set_title("Poloidal Cross-Section (R-Z)")
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_aspect('equal')

    # Subplot 2: R vs time
    ax = axes[0, 1]
    ax.plot(time, r, "r,", markersize=2)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("R [m]")
    ax.set_title("Major Radius Evolution")
    ax.grid(True, alpha=0.3)

    # Subplot 3: Z vs time
    ax = axes[1, 0]
    ax.plot(time, z, "g,", markersize=2)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Z [m]")
    ax.set_title("Vertical Position Evolution")
    ax.grid(True, alpha=0.3)

    # Subplot 4: phi vs time
    ax = axes[1, 1]
    ax.plot(time, phi, "m,", markersize=2)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\phi$ [rad]")
    ax.set_title("Toroidal Angle Evolution")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save as PNG
    png_filename = f"orbit_ascot_{particle_id}.png"
    plt.savefig(png_filename, dpi=150, bbox_inches="tight")
    print(f"Saved plot to {png_filename}")

    # Show plot
    plt.show()

    ds.close()


if __name__ == "__main__":
    main()
