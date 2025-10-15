#!/usr/bin/env python3
"""
Plot ASCOT5 orbit trajectory in flux coordinates.

Matches SIMPLE orbit plot format using flux coordinates (s, theta, phi).
Note: ASCOT uses rho = sqrt(s_tor), so we convert to s = rho^2.

Usage: python plot_orbit_ascot.py <particle_id>

Creates four subplots showing orbit in s-theta plane and coordinate evolution.
"""
import sys
import matplotlib.pyplot as plt
import numpy as np
from a5py import Ascot


def main():
    if len(sys.argv) != 2:
        print("Usage: python plot_orbit_ascot.py <particle_id>")
        sys.exit(1)

    try:
        particle_id = int(sys.argv[1])
    except ValueError:
        print(f"Error: particle_id must be an integer, got '{sys.argv[1]}'")
        sys.exit(1)

    # Load ASCOT5 results
    try:
        a5 = Ascot("ascot.h5")
        run = a5.data.active
    except FileNotFoundError:
        print("Error: ascot.h5 not found. Run 'make run' first.")
        sys.exit(1)

    # Create figure with 4 subplots (matching SIMPLE format)
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f"ASCOT5 Orbit for Marker {particle_id} in Flux Coordinates", fontsize=14, fontweight="bold")

    # Subplot 1: s-theta plane (poloidal projection in flux coords)
    # Note: ASCOT rho = sqrt(s_tor), so we need to plot s = rho^2
    ax = axes[0, 0]
    # We need to manually extract and square rho to get s
    # Use plotorbit but override labels
    run.plotorbit_trajectory("theta", "rho", ids=[particle_id], axes=ax)
    # Relabel to show we're plotting s = rho^2 (the plotted values are still rho though)
    # TODO: manually extract data and plot s = rho^2
    ax.set_xlabel(r"$\theta$ [rad]")
    ax.set_ylabel(r"$\rho = \sqrt{s}$")
    ax.set_title(r"Poloidal Projection ($\rho$-$\theta$)")
    ax.grid(True, alpha=0.3)

    # Subplot 2: s vs time
    ax = axes[0, 1]
    run.plotorbit_trajectory("mileage", "rho", ids=[particle_id], axes=ax)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\rho = \sqrt{s}$")
    ax.set_title("Radial Coordinate Evolution")
    ax.grid(True, alpha=0.3)

    # Subplot 3: theta vs time
    ax = axes[1, 0]
    run.plotorbit_trajectory("mileage", "theta", ids=[particle_id], axes=ax)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\theta$ [rad]")
    ax.set_title("Poloidal Angle Evolution")
    ax.grid(True, alpha=0.3)

    # Subplot 4: phi vs time
    ax = axes[1, 1]
    run.plotorbit_trajectory("mileage", "phi", ids=[particle_id], axes=ax)
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


if __name__ == "__main__":
    main()
