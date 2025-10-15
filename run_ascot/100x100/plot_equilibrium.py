#!/usr/bin/env python3
"""
Plot VMEC equilibrium from ASCOT5 HDF5 file.

This script:
1. Loads an existing ASCOT5 HDF5 file with VMEC magnetic field
2. Plots the magnetic field and flux surfaces in R-Z planes
3. Displays and exports the plot to PNG

Requires: ascot.h5 (run setup_ascot.py first)
"""
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from a5py import Ascot

# Configuration
OUTPUT_H5 = "ascot.h5"
OUTPUT_PNG = "equilibrium.png"


def main():
    print("=" * 60)
    print("ASCOT5 VMEC Equilibrium Plotting")
    print("=" * 60)

    # Check if HDF5 file exists
    h5_path = Path(OUTPUT_H5)
    if not h5_path.exists():
        raise FileNotFoundError(f"HDF5 file not found: {OUTPUT_H5}")
        print(f"Run 'make setup' first to generate {OUTPUT_H5}")

    print(f"\n1. Opening HDF5 file: {OUTPUT_H5}")
    a5 = Ascot(OUTPUT_H5)

    # Initialize the field for evaluation
    print(f"\n2. Initializing field for plotting...")
    a5.input_init(bfield=True)

    # Read field data to get extents
    field_data = a5.data.bfield.active.read()
    rmin = field_data["b_rmin"].item() if hasattr(field_data["b_rmin"], 'item') else field_data["b_rmin"]
    rmax = field_data["b_rmax"].item() if hasattr(field_data["b_rmax"], 'item') else field_data["b_rmax"]
    zmin = field_data["b_zmin"].item() if hasattr(field_data["b_zmin"], 'item') else field_data["b_zmin"]
    zmax = field_data["b_zmax"].item() if hasattr(field_data["b_zmax"], 'item') else field_data["b_zmax"]

    print(f"   Field extent: R=[{rmin:.2f}, {rmax:.2f}] m, Z=[{zmin:.2f}, {zmax:.2f}] m")

    # Create R and Z grids for plotting
    r_grid = np.linspace(rmin, rmax, 200)
    z_grid = np.linspace(zmin, zmax, 200)

    print(f"\n3. Generating plots...")

    # Determine field period from VMEC data
    # For stellarators, one field period = 360°/nfp
    import unyt
    nfp = field_data.get("nfp", 5)  # Default to 5 if not found
    if isinstance(nfp, np.ndarray):
        nfp = int(nfp.item())
    else:
        nfp = int(nfp)

    field_period = 360.0 / nfp
    print(f"   Field period: {field_period:.1f}° (nfp={nfp})")

    # Plot at 4 toroidal angles spanning one field period
    phi_angles = np.linspace(0, field_period * 0.75, 4) * unyt.deg

    # Create figure with subplots (4 rows = 4 phi angles, 2 columns = |B| and rho)
    fig, axes = plt.subplots(4, 2, figsize=(12, 16))
    fig.suptitle(
        f"VMEC Equilibrium\n(One field period: 0° to {field_period:.1f}°)",
        fontsize=13, fontweight='bold'
    )

    for i, phi in enumerate(phi_angles):
        phi_deg = phi.to('deg').v
        print(f"   - Plotting at phi={phi_deg:.1f}°...")

        # Plot |B| (magnetic field magnitude)
        a5.input_plotrz(r_grid, z_grid, "bnorm", phi=phi, axes=axes[i, 0])
        axes[i, 0].set_title(f"|B| at phi={phi_deg:.1f}° [T]", fontsize=11, fontweight='bold')
        axes[i, 0].set_xlabel("R [m]")
        axes[i, 0].set_ylabel("Z [m]")
        axes[i, 0].set_aspect('equal')

        # Plot rho (normalized flux coordinate) with contours
        a5.input_plotrz(r_grid, z_grid, "rho", phi=phi, axes=axes[i, 1], drawcontours=True)
        axes[i, 1].set_title(f"rho at phi={phi_deg:.1f}°", fontsize=11, fontweight='bold')
        axes[i, 1].set_xlabel("R [m]")
        axes[i, 1].set_ylabel("Z [m]")
        axes[i, 1].set_aspect('equal')

    plt.tight_layout()

    # Save to PNG
    print(f"\n4. Saving plot to: {OUTPUT_PNG}")
    plt.savefig(OUTPUT_PNG, dpi=300, bbox_inches='tight')
    print(f"   PNG saved successfully!")

    # Show interactive plot
    print(f"\n5. Displaying interactive plot...")
    print(f"   (Close the plot window to exit)")
    plt.show()

    # Clean up
    a5.input_free()

    print("\n" + "=" * 60)
    print("Plotting complete!")
    print("=" * 60)
    print(f"\nOutputs:")
    print(f"  - PNG plot: {OUTPUT_PNG}")


if __name__ == "__main__":
    main()
