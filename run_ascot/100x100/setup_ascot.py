#!/usr/bin/env python3
"""
Setup ASCOT5 HDF5 file with VMEC magnetic field.

This script:
1. Creates an ASCOT5 HDF5 file
2. Loads VMEC magnetic field data from wout.nc
3. Stores the field in the HDF5 file for later use

Uses the same VMEC file as SIMPLE orbit tracer for consistency.
"""
import numpy as np
from pathlib import Path
from a5py import Ascot

# Configuration
VMEC_FILE = "wout.nc"
OUTPUT_H5 = "ascot.h5"

# Grid parameters for VMEC field import
VMEC_PARAMS = {
    "phimin": 0,
    "phimax": 361,
    "nphi": 181,
    "ntheta": 120,
    "nr": 100,
    "nz": 100,
    "psipad": 0.0,
    "extrapolate": True
}


def main():
    print("=" * 60)
    print("ASCOT5 HDF5 Setup with VMEC Field")
    print("=" * 60)

    # Check if VMEC file exists
    vmec_path = Path(VMEC_FILE)
    if not vmec_path.exists():
        raise FileNotFoundError(f"VMEC file not found: {VMEC_FILE}")
    print(f"\nVMEC file: {vmec_path.resolve()}")

    # Check if HDF5 file already exists
    h5_path = Path(OUTPUT_H5)
    if h5_path.exists():
        print(f"\n1. HDF5 file already exists: {OUTPUT_H5}")
        print(f"   Delete it to regenerate: rm {OUTPUT_H5}")
        print(f"   Opening existing file...")
        a5 = Ascot(OUTPUT_H5)
    else:
        print(f"\n1. Creating new HDF5 file: {OUTPUT_H5}")
        a5 = Ascot(OUTPUT_H5, create=True)

        print(f"\n2. Loading VMEC magnetic field from: {VMEC_FILE}")
        print(f"   Grid parameters:")
        print(f"   - Toroidal: phi=[{VMEC_PARAMS['phimin']}, {VMEC_PARAMS['phimax']}] deg, nphi={VMEC_PARAMS['nphi']}")
        print(f"   - Poloidal: ntheta={VMEC_PARAMS['ntheta']}")
        print(f"   - Spatial: nr={VMEC_PARAMS['nr']}, nz={VMEC_PARAMS['nz']}")
        print(f"   - Extrapolate outside LCFS: {VMEC_PARAMS['extrapolate']}")

        a5.data.create_input(
            "vmec_sts",
            ncfile=VMEC_FILE,
            **VMEC_PARAMS,
            activate=True,
            desc="VMEC stellarator field"
        )

        print("\n3. Field loaded successfully!")
        print(f"   Input name: {a5.data.bfield.active.get_name()}")
        print(f"   QID: {a5.data.bfield.active.get_qid()}")

    # Read field data to get extents
    print(f"\n4. Reading field data...")
    field_data = a5.data.bfield.active.read()
    rmin = field_data["b_rmin"].item() if hasattr(field_data["b_rmin"], 'item') else field_data["b_rmin"]
    rmax = field_data["b_rmax"].item() if hasattr(field_data["b_rmax"], 'item') else field_data["b_rmax"]
    zmin = field_data["b_zmin"].item() if hasattr(field_data["b_zmin"], 'item') else field_data["b_zmin"]
    zmax = field_data["b_zmax"].item() if hasattr(field_data["b_zmax"], 'item') else field_data["b_zmax"]

    print(f"   Field extent: R=[{rmin:.2f}, {rmax:.2f}] m, Z=[{zmin:.2f}, {zmax:.2f}] m")

    # Get field period
    import unyt
    nfp = field_data.get("nfp", 5)
    if isinstance(nfp, np.ndarray):
        nfp = int(nfp.item())
    else:
        nfp = int(nfp)

    field_period = 360.0 / nfp
    print(f"   Field period: {field_period:.1f}° (nfp={nfp})")

    print("\n" + "=" * 60)
    print("Setup complete!")
    print("=" * 60)
    print(f"\nHDF5 file ready: {OUTPUT_H5}")
    print(f"Use 'make plot' to visualize the equilibrium")


if __name__ == "__main__":
    main()
