#!/usr/bin/env python3
"""
Plot FIRM3D orbit trajectory in flux coordinates.

Usage: python plot_orbit_firm3d.py <trajectory_file>

Creates four subplots showing orbit in s-theta plane and coordinate evolution.
Matches SIMPLE/ASCOT5 plot format for comparison.
"""
import sys
import matplotlib.pyplot as plt
import numpy as np


def main():
    if len(sys.argv) != 2:
        print("Usage: python plot_orbit_firm3d.py <trajectory_file>")
        print("\nExample:")
        print("  python plot_orbit_firm3d.py trajectory_data_vmec_tol_1e-06_resolution_48_tmax_0.001_trapped.txt")
        sys.exit(1)

    trajectory_file = sys.argv[1]

    # Load FIRM3D trajectory data
    try:
        data = np.loadtxt(trajectory_file)
    except FileNotFoundError:
        print(f"Error: File not found: {trajectory_file}")
        sys.exit(1)

    # Extract columns: time, s, theta, phi, [extra]
    time = data[:, 0]
    s = data[:, 1]
    theta = data[:, 2]
    phi = data[:, 3]

    # Create figure with 4 subplots (matching SIMPLE/ASCOT format)
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Extract file description for title
    import os
    basename = os.path.basename(trajectory_file)
    is_vmec = 'vmec' in basename
    coord_type = "VMEC" if is_vmec else "Boozer"

    fig.suptitle(f"FIRM3D Orbit in {coord_type} Coordinates\n{basename}",
                 fontsize=14, fontweight="bold")

    # Subplot 1: s-theta plane (poloidal projection)
    ax = axes[0, 0]
    ax.plot(theta, s, "b,", markersize=2)
    ax.scatter(theta[0], s[0], color="green", s=100, marker="o", label="Start", zorder=5)
    ax.scatter(theta[-1], s[-1], color="red", s=100, marker="x", label="End", zorder=5)
    ax.set_xlabel(r"$\theta$ [rad]")
    ax.set_ylabel(r"$s$")
    ax.set_title(r"Poloidal Projection ($s$-$\theta$)")
    ax.grid(True, alpha=0.3)
    ax.legend()

    # Subplot 2: s vs time
    ax = axes[0, 1]
    ax.plot(time, s, "r,", markersize=2)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$s$")
    ax.set_title("Flux Surface Evolution")
    ax.grid(True, alpha=0.3)

    # Subplot 3: theta vs time
    ax = axes[1, 0]
    ax.plot(time, theta, "g,", markersize=2)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\theta$ [rad]")
    ax.set_title("Poloidal Angle Evolution")
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
    png_filename = trajectory_file.replace('.txt', '.png')
    plt.savefig(png_filename, dpi=150, bbox_inches="tight")
    print(f"Saved plot to {png_filename}")

    # Show plot
    plt.show()


if __name__ == "__main__":
    main()
