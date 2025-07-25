#!/usr/bin/env python3
"""
Plot comparison between SIMPLE and firm3d orbit traces.
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
    fig = plt.figure(figsize=(16, 12))
    
    # 1. Poloidal projection (R-Z)
    ax1 = plt.subplot(2, 3, 1)
    ax1.plot(ds_simple.R / 100, ds_simple.Z / 100, 'b-', label='SIMPLE', alpha=0.7)  # Convert cm to m
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
    # SIMPLE uses normalized units, firm3d uses m/s
    v_th = 1.0  # Thermal velocity normalization for SIMPLE (would need actual value)
    ax3.plot(ds_simple.time * 1e6, ds_simple.v_parallel, 'b-', label='SIMPLE (normalized)')
    ax3.plot(ds_firm3d.time * 1e6, ds_firm3d.v_parallel / 1e6, 'r--', label='firm3d [×10⁶ m/s]')
    ax3.set_xlabel('Time [μs]')
    ax3.set_ylabel('Parallel Velocity')
    ax3.set_title('Parallel Velocity Evolution')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. Magnetic moment conservation
    ax4 = plt.subplot(2, 3, 4)
    ax4.plot(ds_simple.time * 1e6, ds_simple.mu / ds_simple.mu[0], 'b-', label='SIMPLE')
    if 'mu' in ds_firm3d:
        ax4.plot(ds_firm3d.time * 1e6, ds_firm3d.mu / ds_firm3d.mu[0], 'r--', label='firm3d')
    ax4.set_xlabel('Time [μs]')
    ax4.set_ylabel('μ/μ₀')
    ax4.set_title('Magnetic Moment Conservation')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    ax4.set_ylim(0.99, 1.01)
    
    # 5. Magnetic field strength
    ax5 = plt.subplot(2, 3, 5)
    ax5.plot(ds_simple.time * 1e6, ds_simple.B, 'b-', label='SIMPLE')
    ax5.plot(ds_firm3d.time * 1e6, ds_firm3d.B, 'r--', label='firm3d')
    ax5.set_xlabel('Time [μs]')
    ax5.set_ylabel('B [T]')
    ax5.set_title('Magnetic Field Strength')
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    
    # 6. Energy conservation (for SIMPLE only, since firm3d doesn't have total velocity in same units)
    ax6 = plt.subplot(2, 3, 6)
    if 'v_total' in ds_simple:
        ax6.plot(ds_simple.time * 1e6, ds_simple.v_total / ds_simple.v_total[0], 'b-', label='SIMPLE')
    ax6.set_xlabel('Time [μs]')
    ax6.set_ylabel('E/E₀')
    ax6.set_title('Energy Conservation')
    ax6.legend()
    ax6.grid(True, alpha=0.3)
    ax6.set_ylim(0.99, 1.01)
    
    plt.suptitle('SIMPLE vs firm3d Orbit Comparison', fontsize=16)
    plt.tight_layout()
    
    # Save figure
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = os.path.join(output_dir, f"orbit_comparison_{timestamp}.png")
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nSaved comparison plot to: {output_file}")
    
    # Create additional detailed plots
    
    # Phase space plot
    fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # SIMPLE phase space (s, v_parallel)
    ax1.plot(ds_simple.s, ds_simple.v_parallel, 'b-', alpha=0.7)
    ax1.set_xlabel('s (normalized flux)')
    ax1.set_ylabel('v_parallel (normalized)')
    ax1.set_title('SIMPLE: Phase Space')
    ax1.grid(True, alpha=0.3)
    
    # firm3d phase space
    ax2.plot(ds_firm3d.s, ds_firm3d.v_parallel / 1e6, 'r-', alpha=0.7)
    ax2.set_xlabel('s (normalized flux)')
    ax2.set_ylabel('v_parallel [×10⁶ m/s]')
    ax2.set_title('firm3d: Phase Space')
    ax2.grid(True, alpha=0.3)
    
    plt.suptitle('Phase Space Comparison', fontsize=14)
    plt.tight_layout()
    
    output_file2 = os.path.join(output_dir, f"phase_space_comparison_{timestamp}.png")
    plt.savefig(output_file2, dpi=150, bbox_inches='tight')
    print(f"Saved phase space plot to: {output_file2}")
    
    # Print comparison statistics
    print("\nComparison Statistics:")
    print("=" * 50)
    
    # Initial conditions
    print("\nInitial Conditions:")
    print(f"  SIMPLE: s={ds_simple.attrs.get('initial_s', 'N/A'):.3f}, "
          f"theta={ds_simple.attrs.get('initial_theta', 'N/A'):.3f}, "
          f"phi={ds_simple.attrs.get('initial_phi', 'N/A'):.3f}")
    print(f"  firm3d: s={ds_firm3d.attrs.get('initial_s', 'N/A'):.3f}, "
          f"theta_booz={ds_firm3d.attrs.get('initial_theta_booz', 'N/A'):.3f}, "
          f"phi_booz={ds_firm3d.attrs.get('initial_phi_booz', 'N/A'):.3f}")
    
    # Conservation properties
    print("\nConservation Properties:")
    if 'v_total' in ds_simple:
        energy_error_simple = abs(1 - (ds_simple.v_total[-1] / ds_simple.v_total[0]).values)
        print(f"  SIMPLE energy error: {energy_error_simple:.2e}")
    
    mu_error_simple = abs(1 - (ds_simple.mu[-1] / ds_simple.mu[0]).values)
    mu_error_firm3d = abs(1 - (ds_firm3d.mu[-1] / ds_firm3d.mu[0]).values)
    print(f"  SIMPLE μ error: {mu_error_simple:.2e}")
    print(f"  firm3d μ error: {mu_error_firm3d:.2e}")
    
    # Orbit characteristics
    print("\nOrbit Characteristics:")
    print(f"  SIMPLE: s range = [{ds_simple.s.min().values:.3f}, {ds_simple.s.max().values:.3f}]")
    print(f"  firm3d: s range = [{ds_firm3d.s.min().values:.3f}, {ds_firm3d.s.max().values:.3f}]")
    
    # Close datasets
    ds_simple.close()
    ds_firm3d.close()
    
    return output_file, output_file2


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