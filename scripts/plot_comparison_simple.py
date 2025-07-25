#!/usr/bin/env python3
"""
Simplified plot comparison between SIMPLE and firm3d orbit traces.
"""

import os
import glob
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from datetime import datetime


def find_latest_files(run_dir="run"):
    """Find the latest SIMPLE and firm3d orbit files."""
    simple_files = glob.glob(os.path.join(run_dir, "trace_orbit_simple_*.nc"))
    firm3d_files = glob.glob(os.path.join(run_dir, "trace_orbit_firm3d_*.nc"))
    
    if not simple_files:
        raise FileNotFoundError(f"No SIMPLE orbit files found in {run_dir}")
    if not firm3d_files:
        raise FileNotFoundError(f"No firm3d orbit files found in {run_dir}")
    
    # Get the most recent files
    simple_file = max(simple_files, key=os.path.getmtime)
    firm3d_file = max(firm3d_files, key=os.path.getmtime)
    
    return simple_file, firm3d_file


def plot_comparison(simple_file, firm3d_file, output_dir="plot"):
    """Create comparison plots between SIMPLE and firm3d orbits."""
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Load datasets
    ds_simple = xr.open_dataset(simple_file)
    ds_firm3d = xr.open_dataset(firm3d_file)
    
    print(f"Loaded SIMPLE data from: {simple_file}")
    print(f"Loaded firm3d data from: {firm3d_file}")
    
    # Create figure with subplots
    fig = plt.figure(figsize=(15, 10))
    
    # 1. Poloidal projection (R-Z)
    ax1 = plt.subplot(2, 3, 1)
    if 'R' in ds_simple:
        ax1.plot(ds_simple.R / 100, ds_simple.Z / 100, 'b-', label='SIMPLE', alpha=0.7)  # Convert cm to m
    if 'R' in ds_firm3d:
        ax1.plot(ds_firm3d.R, ds_firm3d.Z, 'r--', label='firm3d', alpha=0.7)
    ax1.set_xlabel('R [m]')
    ax1.set_ylabel('Z [m]')
    ax1.set_title('Poloidal Projection')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.axis('equal')
    
    # 2. Flux surface evolution
    ax2 = plt.subplot(2, 3, 2)
    ax2.plot(ds_simple.time * 1e6, ds_simple.s, 'b-', label='SIMPLE')  # Convert to microseconds
    ax2.plot(ds_firm3d.time * 1e6, ds_firm3d.s, 'r--', label='firm3d')
    ax2.set_xlabel('Time [μs]')
    ax2.set_ylabel('s (normalized flux)')
    ax2.set_title('Flux Surface Evolution')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Parallel velocity
    ax3 = plt.subplot(2, 3, 3)
    if 'v_parallel' in ds_simple:
        ax3.plot(ds_simple.time * 1e6, ds_simple.v_parallel, 'b-', label='SIMPLE (normalized)')
    if 'v_parallel' in ds_firm3d:
        ax3.plot(ds_firm3d.time * 1e6, ds_firm3d.v_parallel / 1e6, 'r--', label='firm3d [×10⁶ m/s]')
    ax3.set_xlabel('Time [μs]')
    ax3.set_ylabel('Parallel Velocity')
    ax3.set_title('Parallel Velocity')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. Magnetic field strength
    ax4 = plt.subplot(2, 3, 4)
    if 'B' in ds_simple:
        ax4.plot(ds_simple.time * 1e6, ds_simple.B, 'b-', label='SIMPLE')
    if 'B' in ds_firm3d:
        ax4.plot(ds_firm3d.time * 1e6, ds_firm3d.B, 'r--', label='firm3d')
    ax4.set_xlabel('Time [μs]')
    ax4.set_ylabel('B [T]')
    ax4.set_title('Magnetic Field Strength')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # 5. Toroidal angle evolution
    ax5 = plt.subplot(2, 3, 5)
    if 'phi_vmec' in ds_simple:
        ax5.plot(ds_simple.time * 1e6, ds_simple.phi_vmec, 'b-', label='SIMPLE (VMEC)', alpha=0.7)
    elif 'phi_can' in ds_simple:
        ax5.plot(ds_simple.time * 1e6, ds_simple.phi_can, 'b-', label='SIMPLE (canonical)', alpha=0.7)
    if 'phi_booz' in ds_firm3d:
        ax5.plot(ds_firm3d.time * 1e6, ds_firm3d.phi_booz, 'r--', label='firm3d (Boozer)', alpha=0.7)
    elif 'phi' in ds_firm3d:
        ax5.plot(ds_firm3d.time * 1e6, ds_firm3d.phi, 'r--', label='firm3d', alpha=0.7)
    ax5.set_xlabel('Time [μs]')
    ax5.set_ylabel('Toroidal Angle [rad]')
    ax5.set_title('Toroidal Angle Evolution')
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    
    # 6. Conservation check (if available)
    ax6 = plt.subplot(2, 3, 6)
    if 'mu' in ds_simple:
        ax6.plot(ds_simple.time * 1e6, ds_simple.mu / ds_simple.mu[0], 'b-', label='SIMPLE μ/μ₀')
    if 'v_total' in ds_simple:
        ax6.plot(ds_simple.time * 1e6, ds_simple.v_total / ds_simple.v_total[0], 'b--', label='SIMPLE E/E₀', alpha=0.7)
    ax6.set_xlabel('Time [μs]')
    ax6.set_ylabel('Conservation')
    ax6.set_title('Conservation Properties')
    ax6.legend()
    ax6.grid(True, alpha=0.3)
    ax6.set_ylim(0.98, 1.02)
    
    plt.suptitle('SIMPLE vs firm3d Orbit Comparison', fontsize=16)
    plt.tight_layout()
    
    # Save figure
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = os.path.join(output_dir, f"orbit_comparison_{timestamp}.png")
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nSaved comparison plot to: {output_file}")
    
    # Print comparison statistics
    print("\nComparison Statistics:")
    print("=" * 50)
    
    # Time ranges
    print(f"\nTime ranges:")
    print(f"  SIMPLE: {ds_simple.time.values[0]:.2e} - {ds_simple.time.values[-1]:.2e} s")
    print(f"  firm3d: {ds_firm3d.time.values[0]:.2e} - {ds_firm3d.time.values[-1]:.2e} s")
    
    # Orbit characteristics
    print(f"\nOrbit characteristics:")
    print(f"  SIMPLE: s range = [{ds_simple.s.min().values:.3f}, {ds_simple.s.max().values:.3f}]")
    print(f"  firm3d: s range = [{ds_firm3d.s.min().values:.3f}, {ds_firm3d.s.max().values:.3f}]")
    
    # Field strength ranges
    if 'B' in ds_simple and 'B' in ds_firm3d:
        print(f"\nMagnetic field:")
        print(f"  SIMPLE: B range = [{ds_simple.B.min().values:.3f}, {ds_simple.B.max().values:.3f}] T")
        print(f"  firm3d: B range = [{ds_firm3d.B.min().values:.3f}, {ds_firm3d.B.max().values:.3f}] T")
    
    # Close datasets
    ds_simple.close()
    ds_firm3d.close()
    
    return output_file


if __name__ == "__main__":
    try:
        # Find latest orbit files
        simple_file, firm3d_file = find_latest_files()
        
        # Create comparison plots
        plot_comparison(simple_file, firm3d_file)
        
    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("Please run both 'make run' and 'make run-firm3d' before plotting comparison.")
        exit(1)
    except Exception as e:
        print(f"Error creating comparison plots: {e}")
        exit(1)