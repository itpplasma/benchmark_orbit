#!/usr/bin/env python3
"""
Trace a single particle orbit in VMEC equilibrium using SIMPLE.
Saves orbit data to netCDF file using xarray.
"""

import numpy as np
import xarray as xr
import os
import sys
from datetime import datetime

# Add SIMPLE to path if needed
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'codes', 'SIMPLE'))

from pysimple import simple_main, simple, params, orbit_symplectic
from pysimple import get_can_sub as coord

def trace_orbit(vmec_file, output_dir="run", trace_time=1e-5):
    """
    Trace a single particle orbit and save to netCDF.
    
    Parameters:
    -----------
    vmec_file : str
        Path to VMEC wout file
    output_dir : str
        Output directory for results
    trace_time : float
        Total tracing time in seconds
    """
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize tracer
    tracy = simple.Tracer()
    
    # Initialize field from VMEC file
    simple_main.init_field(tracy, vmec_file, 3, 3, 3, 1)
    params.params_init()
    
    # Initial conditions: s=0.5, theta=0.5, phi=0.5, vpar=0, 3.5 MeV alpha
    # Convert 3.5 MeV alpha energy to normalized velocity
    # Alpha particle: mass = 4 amu = 6.644e-27 kg
    # Energy = 3.5 MeV = 3.5 * 1.602e-13 J = 5.607e-13 J
    # v = sqrt(2*E/m) = sqrt(2*5.607e-13/6.644e-27) = 1.299e7 m/s
    # For normalization, we need thermal velocity v_th
    # Assuming v/v_th = 1.0 for simplicity (can be adjusted)
    
    z0_vmec = np.array([0.5, 0.5, 0.5, 1.0, 0.0])  # s, th, ph, v/v_th, v_par/v
    z0_can = z0_vmec.copy()
    
    # Convert VMEC to canonical coordinates
    z0_can[1:3] = coord.vmec_to_can(z0_vmec[0], z0_vmec[1], z0_vmec[2])
    
    # Initialize symplectic integrator
    simple.init_sympl(tracy.si, tracy.f, z0_can, params.dtaumin, params.dtaumin, 1e-13, tracy.integmode)
    
    # Calculate number of timesteps
    # Assuming normalized time units, adjust dtau appropriately
    dtau = params.dtaumin
    nt = int(trace_time / dtau) + 1
    
    print(f"Tracing orbit for {trace_time} seconds with {nt} timesteps")
    print(f"Initial B = {tracy.f.bmod} T")
    
    # Prepare arrays for storing data
    time = np.zeros(nt)
    s = np.zeros(nt)
    theta_vmec = np.zeros(nt)
    phi_vmec = np.zeros(nt)
    theta_can = np.zeros(nt)
    phi_can = np.zeros(nt)
    v_total = np.zeros(nt)
    v_par = np.zeros(nt)
    v_perp = np.zeros(nt)
    pphi = np.zeros(nt)
    mu = np.zeros(nt)
    b_field = np.zeros(nt)
    r_cyl = np.zeros(nt)
    z_cyl = np.zeros(nt)
    
    # Store initial values
    time[0] = 0.0
    s[0] = z0_vmec[0]
    theta_vmec[0] = z0_vmec[1]
    phi_vmec[0] = z0_vmec[2]
    theta_can[0] = z0_can[1]
    phi_can[0] = z0_can[2]
    v_total[0] = z0_vmec[3]
    v_par[0] = tracy.f.vpar
    v_perp[0] = np.sqrt(2 * tracy.f.mu * tracy.f.bmod)
    pphi[0] = tracy.si.z[3]
    mu[0] = tracy.f.mu
    b_field[0] = tracy.f.bmod
    r_cyl[0], z_cyl[0] = coord.vmec_to_cyl(s[0], theta_vmec[0], phi_vmec[0])
    
    # Time integration loop
    for kt in range(nt-1):
        # Advance one timestep
        orbit_symplectic.orbit_timestep_sympl_expl_impl_euler(tracy.si, tracy.f)
        
        # Store time
        time[kt+1] = (kt+1) * dtau
        
        # Get integrated variables
        s[kt+1] = tracy.si.z[0]
        theta_can[kt+1] = tracy.si.z[1]
        phi_can[kt+1] = tracy.si.z[2]
        pphi[kt+1] = tracy.si.z[3]
        
        # Convert canonical to VMEC coordinates
        theta_vmec[kt+1], phi_vmec[kt+1] = coord.can_to_vmec(
            s[kt+1], theta_can[kt+1], phi_can[kt+1])
        
        # Store velocities and invariants
        v_par[kt+1] = tracy.f.vpar
        mu[kt+1] = tracy.f.mu
        b_field[kt+1] = tracy.f.bmod
        v_perp[kt+1] = np.sqrt(2 * mu[kt+1] * b_field[kt+1])
        v_total[kt+1] = np.sqrt(v_par[kt+1]**2 + v_perp[kt+1]**2)
        
        # Convert to cylindrical coordinates
        r_cyl[kt+1], z_cyl[kt+1] = coord.vmec_to_cyl(
            s[kt+1], theta_vmec[kt+1], phi_vmec[kt+1])
    
    # Create xarray dataset
    ds = xr.Dataset(
        {
            "s": (["time"], s, {"long_name": "Normalized toroidal flux", "units": ""}),
            "theta_vmec": (["time"], theta_vmec, {"long_name": "VMEC poloidal angle", "units": "rad"}),
            "phi_vmec": (["time"], phi_vmec, {"long_name": "VMEC toroidal angle", "units": "rad"}),
            "theta_can": (["time"], theta_can, {"long_name": "Canonical poloidal angle", "units": "rad"}),
            "phi_can": (["time"], phi_can, {"long_name": "Canonical toroidal angle", "units": "rad"}),
            "v_total": (["time"], v_total, {"long_name": "Total velocity", "units": "v_th"}),
            "v_parallel": (["time"], v_par, {"long_name": "Parallel velocity", "units": "v_th"}),
            "v_perpendicular": (["time"], v_perp, {"long_name": "Perpendicular velocity", "units": "v_th"}),
            "pphi": (["time"], pphi, {"long_name": "Canonical toroidal momentum", "units": ""}),
            "mu": (["time"], mu, {"long_name": "Magnetic moment", "units": ""}),
            "B": (["time"], b_field, {"long_name": "Magnetic field strength", "units": "T"}),
            "R": (["time"], r_cyl, {"long_name": "Major radius", "units": "cm"}),
            "Z": (["time"], z_cyl, {"long_name": "Vertical coordinate", "units": "cm"}),
        },
        coords={
            "time": (["time"], time, {"long_name": "Time", "units": "s"}),
        },
        attrs={
            "description": "Particle orbit trace from SIMPLE",
            "vmec_file": vmec_file,
            "initial_s": float(z0_vmec[0]),
            "initial_theta": float(z0_vmec[1]),
            "initial_phi": float(z0_vmec[2]),
            "initial_v_parallel": float(z0_vmec[4]),
            "particle_type": "alpha",
            "energy_MeV": 3.5,
            "created": datetime.now().isoformat(),
        }
    )
    
    # Save to netCDF
    output_file = os.path.join(output_dir, f"orbit_trace_{datetime.now().strftime('%Y%m%d_%H%M%S')}.nc")
    ds.to_netcdf(output_file)
    print(f"Saved orbit data to: {output_file}")
    
    return ds

if __name__ == "__main__":
    # Check if VMEC file is provided
    if len(sys.argv) > 1:
        vmec_file = sys.argv[1]
    else:
        # Default to wout.nc in current directory
        vmec_file = "wout.nc"
    
    if not os.path.exists(vmec_file):
        print(f"Error: VMEC file '{vmec_file}' not found!")
        print("Usage: python trace_orbit.py [path/to/wout.nc]")
        sys.exit(1)
    
    # Run orbit tracing
    ds = trace_orbit(vmec_file)
    
    # Print summary
    print(f"\nOrbit summary:")
    print(f"  Time range: {ds.time.values[0]:.2e} - {ds.time.values[-1]:.2e} s")
    print(f"  s range: {ds.s.min().values:.3f} - {ds.s.max().values:.3f}")
    print(f"  Energy conservation: {(ds.v_total[-1]/ds.v_total[0]).values:.6f}")
    print(f"  Mu conservation: {(ds.mu[-1]/ds.mu[0]).values:.6f}")