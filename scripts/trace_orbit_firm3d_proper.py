#!/usr/bin/env python3
"""
Proper firm3d orbit tracing using trace_particles_boozer with full guiding center physics.
Based on the firm3d example.
"""

import os
import sys
import time
from datetime import datetime

import numpy as np
import xarray as xr

# Add firm3d to path
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "codes", "firm3d", "src"))

from simsopt.field.boozermagneticfield import (
    BoozerRadialInterpolant,
    InterpolatedBoozerField,
)
from simsopt.field.tracing import (
    MaxToroidalFluxStoppingCriterion,
    trace_particles_boozer,
)
from simsopt.util.constants import (
    ALPHA_PARTICLE_CHARGE,
    ALPHA_PARTICLE_MASS,
    FUSION_ALPHA_PARTICLE_ENERGY,
)


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

    return s_values[0], theta_values[0], phi_values[0]


def trace_orbit(vmec_file, output_dir="run", trace_time=1e-4):
    """Trace particle orbits using firm3d with full guiding center physics."""
    
    os.makedirs(output_dir, exist_ok=True)

    # Read initial conditions
    s_init, theta_init, phi_init = read_initial_conditions()
    print(f"Initial conditions: s={s_init:.3f}, vartheta={theta_init:.3f}, zeta={phi_init:.3f}")

    # Use the reactor scale boozmn file
    boozmn_filename = "booz_xform/boozmn_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"
    print(f"Using Boozer file: {boozmn_filename}")
    
    # Parameters from firm3d example
    order = 3  # Order for radial interpolation
    reltol = 1e-8  # Relative tolerance for the ODE solver
    abstol = 1e-8  # Absolute tolerance for the ODE solver
    tmax = trace_time  # Time for integration
    resolution = 48  # Resolution for field interpolation
    degree = 3  # Degree for 3d interpolation
    
    # Setup radial interpolation
    print("Setting up BoozerRadialInterpolant...")
    bri = BoozerRadialInterpolant(boozmn_filename, order, no_K=True)
    nfp = bri.nfp
    
    # Setup 3d interpolation
    print("Setting up InterpolatedBoozerField...")
    field = InterpolatedBoozerField(
        bri,
        degree,
        ns_interp=resolution,
        ntheta_interp=resolution,
        nzeta_interp=resolution,
    )
    
    # Particle parameters
    Ekin = FUSION_ALPHA_PARTICLE_ENERGY
    mass = ALPHA_PARTICLE_MASS
    charge = ALPHA_PARTICLE_CHARGE
    vpar0 = np.sqrt(2 * Ekin / mass)
    
    # Initialize trapped particle with zero parallel velocity
    vpar_init = [0]  # Zero parallel velocity for trapped particle
    points = np.zeros((1, 3))
    points[0, 0] = s_init  # s coordinate
    points[0, 1] = theta_init  # Boozer poloidal angle
    points[0, 2] = phi_init  # Boozer toroidal angle
    
    print(f"\nTracing alpha particle (3.5 MeV) with v_par = 0")
    print(f"Total velocity: {vpar0:.3e} m/s")
    
    # Trace particle in Boozer coordinates
    traj_booz, res_hits = trace_particles_boozer(
        field,
        points,
        vpar_init,
        tmax=tmax,
        mass=mass,
        charge=charge,
        Ekin=Ekin,
        stopping_criteria=[MaxToroidalFluxStoppingCriterion(1.0)],
        forget_exact_path=False,
        dt_save=1e-7,
        abstol=abstol,
        reltol=reltol,
    )
    
    print(f"Orbit tracing completed. Trajectory points: {len(traj_booz[0])}")
    
    # Extract trajectory data
    traj = traj_booz[0]
    nt = len(traj)
    
    # Extract data from trajectory
    time_arr = np.array([traj[i][5] for i in range(nt)])  # Time
    s_arr = np.array([traj[i][0] for i in range(nt)])  # s coordinate
    theta_arr = np.array([traj[i][1] for i in range(nt)])  # Boozer poloidal angle
    phi_arr = np.array([traj[i][2] for i in range(nt)])  # Boozer toroidal angle
    vpar_arr = np.array([traj[i][3] for i in range(nt)])  # Parallel velocity
    mu_arr = np.array([traj[i][4] for i in range(nt)])  # Magnetic moment
    
    # Get cylindrical coordinates and B field
    R_arr = np.zeros(nt)
    Z_arr = np.zeros(nt)
    B_arr = np.zeros(nt)
    
    for i in range(nt):
        xyz = field.evaluate_xyz(s_arr[i], theta_arr[i], phi_arr[i])
        R_arr[i] = np.sqrt(xyz[0]**2 + xyz[1]**2)
        Z_arr[i] = xyz[2]
        B_arr[i] = field.modB(s_arr[i], theta_arr[i], phi_arr[i])
    
    # Create xarray dataset
    ds = xr.Dataset(
        {
            "s": (["time"], s_arr, {"long_name": "Normalized toroidal flux", "units": ""}),
            "theta_booz": (["time"], theta_arr, {"long_name": "Boozer poloidal angle", "units": "rad"}),
            "phi_booz": (["time"], phi_arr, {"long_name": "Boozer toroidal angle", "units": "rad"}),
            "v_parallel": (["time"], vpar_arr, {"long_name": "Parallel velocity", "units": "m/s"}),
            "mu": (["time"], mu_arr, {"long_name": "Magnetic moment", "units": "J/T"}),
            "B": (["time"], B_arr, {"long_name": "Magnetic field strength", "units": "T"}),
            "R": (["time"], R_arr, {"long_name": "Major radius", "units": "m"}),
            "Z": (["time"], Z_arr, {"long_name": "Vertical position", "units": "m"}),
        },
        coords={
            "time": (["time"], time_arr, {"long_name": "Time", "units": "s"}),
        },
        attrs={
            "description": "Particle orbit trace from firm3d using trace_particles_boozer",
            "code": "firm3d",
            "particle_type": "alpha",
            "energy_MeV": 3.5,
            "initial_s": float(s_init),
            "initial_theta_booz": float(theta_init),
            "initial_phi_booz": float(phi_init),
            "initial_v_parallel": 0.0,
            "created": datetime.now().isoformat(),
        },
    )
    
    # Save to netCDF
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = os.path.join(output_dir, f"trace_orbit_firm3d_{timestamp}.nc")
    
    ds.to_netcdf(output_file)
    print(f"\nSaved orbit data to: {output_file}")
    
    # Print summary
    print(f"\nOrbit summary:")
    print(f"  Time range: {time_arr[0]:.2e} - {time_arr[-1]:.2e} s")
    print(f"  s range: {s_arr.min():.3f} - {s_arr.max():.3f}")
    print(f"  Mu conservation: {mu_arr[-1] / mu_arr[0]:.6f}")
    
    return ds


if __name__ == "__main__":
    if len(sys.argv) > 1:
        vmec_file = sys.argv[1]
    else:
        vmec_file = "wout.nc"

    try:
        ds = trace_orbit(vmec_file)
        print("\nOrbit trace completed successfully!")
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()