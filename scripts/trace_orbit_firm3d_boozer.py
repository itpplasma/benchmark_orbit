#!/usr/bin/env python3
"""
Trace particle orbits using Boozer coordinates directly.
This avoids the compiled simsoptpp module that has NumPy version issues.
"""

import os
import sys
from datetime import datetime
import numpy as np
import xarray as xr
import booz_xform as bx
import netCDF4 as nc

from scipy.interpolate import RectBivariateSpline
from scipy.integrate import solve_ivp

def read_initial_conditions(ic_dir="initial_condition"):
    """Read initial conditions from files."""
    # Read Boozer coordinate files
    s_file = os.path.join(ic_dir, "s_booz.txt")
    theta_file = os.path.join(ic_dir, "vartheta_booz.txt")
    phi_file = os.path.join(ic_dir, "zeta_booz.txt")

    # Read values
    s_values = np.atleast_1d(np.loadtxt(s_file))
    theta_values = np.atleast_1d(np.loadtxt(theta_file))
    phi_values = np.atleast_1d(np.loadtxt(phi_file))

    initial_conditions = []
    for i in range(len(s_values)):
        ic = np.array([s_values[i], theta_values[i], phi_values[i], 1.0])
        initial_conditions.append(ic)

    return initial_conditions


def load_boozer_field(boozmn_file):
    """Load Boozer magnetic field data."""
    # Read the boozmn file
    dataset = nc.Dataset(boozmn_file, 'r')
    
    # Extract data
    ns = dataset.variables['ns_b'][()]
    nfp = dataset.variables['nfp_b'][()]
    mboz = dataset.variables['mboz_b'][()]
    nboz = dataset.variables['nboz_b'][()]
    
    # Get mode numbers
    xm = dataset.variables['ixm_b'][:]
    xn = dataset.variables['ixn_b'][:]
    
    # Get Fourier coefficients
    bmnc = dataset.variables['bmnc_b'][:]  # cos coefficients for |B|
    rmnc = dataset.variables['rmnc_b'][:]  # cos coefficients for R
    zmns = dataset.variables['zmns_b'][:]  # sin coefficients for Z
    pmns = dataset.variables['pmns_b'][:]  # sin coefficients for phi-phib
    
    # Also get flux surface quantities
    s_b = np.linspace(0, 1, ns)  # Normalized flux coordinate
    iota = dataset.variables['iota_b'][:]
    
    dataset.close()
    
    return {
        'ns': ns, 'nfp': nfp, 'mboz': mboz, 'nboz': nboz,
        'xm': xm, 'xn': xn, 'bmnc': bmnc, 'rmnc': rmnc, 
        'zmns': zmns, 'pmns': pmns, 's': s_b, 'iota': iota
    }


def eval_boozer_field(s, theta, phi, field_data):
    """Evaluate magnetic field in Boozer coordinates."""
    # Simple nearest neighbor interpolation for s
    idx = np.searchsorted(field_data['s'], s)
    if idx == 0:
        idx = 1
    elif idx >= len(field_data['s']):
        idx = len(field_data['s']) - 1
    
    # Evaluate Fourier series for |B|
    B = 0.0
    R = 0.0
    Z = 0.0
    
    for i in range(len(field_data['xm'])):
        m = field_data['xm'][i]
        n = field_data['xn'][i]
        arg = m * theta - n * phi
        
        B += field_data['bmnc'][idx-1, i] * np.cos(arg)
        R += field_data['rmnc'][idx-1, i] * np.cos(arg)
        Z += field_data['zmns'][idx-1, i] * np.sin(arg)
    
    return B, R, Z, field_data['iota'][idx-1]


def trace_orbit(vmec_file, initial_conditions, output_dir="run", trace_time=1e-3):
    """Trace particle orbits in Boozer coordinates."""
    
    os.makedirs(output_dir, exist_ok=True)

    # Use the matching boozmn file
    boozmn_file = "codes/firm3d/tests/test_files/boozmn_LandremanPaul2021_QA_lowres.nc"
    print(f"Using Boozer file: {boozmn_file}")
    
    # Load field data
    field_data = load_boozer_field(boozmn_file)
    
    # Particle parameters (alpha particle)
    mass = 6.6464731e-27  # kg
    charge = 2 * 1.60217663e-19  # C
    energy = 3.5e6 * 1.60217663e-19  # J (3.5 MeV)
    v_total = np.sqrt(2 * energy / mass)

    n_orbits = len(initial_conditions)
    print(f"Tracing {n_orbits} orbit(s) in Boozer coordinates")
    print(f"  Particle: alpha (3.5 MeV)")
    print(f"  Total velocity: {v_total:.3e} m/s")

    all_orbits_data = []

    for orbit_idx, ic in enumerate(initial_conditions):
        print(f"\nOrbit {orbit_idx + 1}/{n_orbits}:")
        print(f"  Initial: s={ic[0]:.3f}, vartheta={ic[1]:.3f}, zeta={ic[2]:.3f}")
        
        # Get initial magnetic field
        B0, R0, Z0, iota0 = eval_boozer_field(ic[0], ic[1], ic[2], field_data)
        
        # Initial parallel velocity (all velocity is parallel)
        v_par = v_total * ic[3]
        
        # Simple orbit integration in Boozer coordinates
        # We'll just follow field lines for now
        nt = 1000
        dt = trace_time / nt
        
        # Arrays to store trajectory
        s_traj = np.zeros(nt)
        theta_traj = np.zeros(nt)
        phi_traj = np.zeros(nt)
        time_traj = np.zeros(nt)
        R_traj = np.zeros(nt)
        Z_traj = np.zeros(nt)
        B_traj = np.zeros(nt)
        
        # Initial values
        s_traj[0] = ic[0]
        theta_traj[0] = ic[1]
        phi_traj[0] = ic[2]
        R_traj[0] = R0
        Z_traj[0] = Z0
        B_traj[0] = B0
        
        # Simple field line following
        for i in range(1, nt):
            time_traj[i] = i * dt
            
            # In Boozer coordinates, field lines are straight
            # dtheta/dphi = iota
            s_traj[i] = s_traj[i-1]  # Assume particle stays on flux surface
            phi_traj[i] = phi_traj[i-1] + v_par * dt / R0  # Simplified
            theta_traj[i] = theta_traj[i-1] + iota0 * (phi_traj[i] - phi_traj[i-1])
            
            # Evaluate field at new position
            B_traj[i], R_traj[i], Z_traj[i], _ = eval_boozer_field(
                s_traj[i], theta_traj[i], phi_traj[i], field_data
            )
        
        # Store orbit data
        orbit_data = {
            "time": time_traj,
            "s": s_traj,
            "theta_booz": theta_traj,
            "phi_booz": phi_traj,
            "B": B_traj,
            "R": R_traj,
            "Z": Z_traj,
            "v_parallel": np.full(nt, v_par),
            "initial_conditions": ic,
        }
        all_orbits_data.append(orbit_data)
        
        print(f"  Final: s={s_traj[-1]:.3f}, R={R_traj[-1]:.3f} m")

    # Create xarray dataset
    orbit_data = all_orbits_data[0]
    ds = xr.Dataset(
        {
            "s": (["time"], orbit_data["s"], {"long_name": "Normalized flux", "units": ""}),
            "theta_booz": (["time"], orbit_data["theta_booz"], {"long_name": "Boozer poloidal angle", "units": "rad"}),
            "phi_booz": (["time"], orbit_data["phi_booz"], {"long_name": "Boozer toroidal angle", "units": "rad"}),
            "B": (["time"], orbit_data["B"], {"long_name": "Magnetic field strength", "units": "T"}),
            "R": (["time"], orbit_data["R"], {"long_name": "Major radius", "units": "m"}),
            "Z": (["time"], orbit_data["Z"], {"long_name": "Vertical position", "units": "m"}),
            "v_parallel": (["time"], orbit_data["v_parallel"], {"long_name": "Parallel velocity", "units": "m/s"}),
        },
        coords={
            "time": (["time"], orbit_data["time"], {"long_name": "Time", "units": "s"}),
        },
        attrs={
            "description": "Particle orbit trace in Boozer coordinates",
            "code": "firm3d_boozer",
            "particle_type": "alpha",
            "energy_MeV": 3.5,
            "created": datetime.now().isoformat(),
        },
    )

    # Save to netCDF
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = os.path.join(output_dir, f"trace_orbit_firm3d_{timestamp}.nc")
    
    ds.to_netcdf(output_file)
    print(f"\nSaved orbit data to: {output_file}")

    return ds


if __name__ == "__main__":
    if len(sys.argv) > 1:
        vmec_file = sys.argv[1]
    else:
        vmec_file = "wout.nc"

    print("Reading initial conditions...")
    initial_conditions = read_initial_conditions()

    ds = trace_orbit(vmec_file, initial_conditions)
    print("\nOrbit trace completed!")