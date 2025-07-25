#!/usr/bin/env python3
"""
Trace particle orbits in VMEC equilibrium using firm3d.
Reads initial conditions from files and saves orbit data to netCDF file using xarray.
"""

import os
import sys
from datetime import datetime

import numpy as np
import xarray as xr

# Add firm3d to path if needed
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "codes", "firm3d"))

from simsopt.field.boozermagneticfield import (
    BoozerRadialInterpolant,
    InterpolatedBoozerField,
)
from simsopt.field.tracing import (
    trace_particles_boozer,
    MaxToroidalFluxStoppingCriterion,
)
from simsopt.util.constants import (
    ALPHA_PARTICLE_MASS,
    ALPHA_PARTICLE_CHARGE,
    FUSION_ALPHA_PARTICLE_ENERGY,
)


def read_initial_conditions(ic_dir="initial_condition"):
    """
    Read initial conditions from files in the initial_condition directory.

    Returns:
    --------
    initial_conditions : list of arrays
        Each array contains [s, theta, phi, v_par/v]
    """
    # Read coordinate files
    s_file = os.path.join(ic_dir, "s_booz.txt")
    theta_file = os.path.join(ic_dir, "theta_vmec.txt")
    phi_file = os.path.join(ic_dir, "phi_vmec.txt")

    # Check if files exist
    if not all(os.path.exists(f) for f in [s_file, theta_file, phi_file]):
        print(f"Warning: Initial condition files not found in {ic_dir}")
        print("Using default values: s=0.5, theta=0.5, phi=0.5")
        return [np.array([0.5, 0.5, 0.5, 1.0])]

    # Read values from files
    s_values = np.loadtxt(s_file)
    theta_values = np.loadtxt(theta_file)
    phi_values = np.loadtxt(phi_file)

    # Ensure arrays have the same shape
    s_values = np.atleast_1d(s_values)
    theta_values = np.atleast_1d(theta_values)
    phi_values = np.atleast_1d(phi_values)

    # Check if all have the same number of values
    n_orbits = len(s_values)
    if not (len(theta_values) == n_orbits and len(phi_values) == n_orbits):
        raise ValueError("Initial condition files must have the same number of lines")

    # Create initial conditions array
    initial_conditions = []
    for i in range(n_orbits):
        # s, theta, phi, v_par/v (purely parallel motion like SIMPLE)
        ic = np.array([s_values[i], theta_values[i], phi_values[i], 1.0])
        initial_conditions.append(ic)

    return initial_conditions


def trace_orbit(vmec_file, initial_conditions, output_dir="run", trace_time=1e-3):
    """
    Trace particle orbits using firm3d and save to netCDF.

    Parameters:
    -----------
    vmec_file : str
        Path to VMEC wout file
    initial_conditions : list of arrays
        Initial conditions for each orbit
    output_dir : str
        Output directory for results
    trace_time : float
        Total tracing time in seconds
    """

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Create boozmn file from VMEC if needed
    boozmn_file = os.path.join(output_dir, "boozmn_temp.nc")
    
    # Convert VMEC to Boozer coordinates
    from simsopt.mhd import Vmec, Boozer
    vmec = Vmec(vmec_file)
    booz = Boozer(vmec, mpol=32, ntor=32)
    booz.run()
    booz.save(boozmn_file)
    
    print(f"Created Boozer file: {boozmn_file}")

    # Setup field interpolation (similar to firm3d example)
    order = 3  # Order for radial interpolation
    resolution = 48  # Resolution for field interpolation
    degree = 3  # Degree for 3d interpolation
    
    # Setup radial interpolation
    bri = BoozerRadialInterpolant(boozmn_file, order, no_K=True)
    nfp = bri.nfp
    
    # Setup 3d interpolation
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

    n_orbits = len(initial_conditions)
    print(f"Tracing {n_orbits} orbit(s) using firm3d")
    print("  Particle parameters:")
    print(f"    Type: alpha particle")
    print(f"    Energy: {Ekin / 1.60217663e-19:.1f} eV ({Ekin / 1.60217663e-13:.1f} MeV)")
    print(f"    Mass: {mass:.3e} kg")
    print(f"    Charge: {charge:.3e} C")
    print(f"    Initial v_parallel: {vpar0:.3e} m/s")

    all_orbits_data = []

    for orbit_idx, ic in enumerate(initial_conditions):
        print(f"\nOrbit {orbit_idx + 1}/{n_orbits}:")
        print(f"  Initial conditions: s={ic[0]:.3f}, theta={ic[1]:.3f}, phi={ic[2]:.3f}")

        # Setup initial position (s, theta, phi in Boozer coordinates)
        points = np.zeros((1, 3))
        points[0, 0] = ic[0]  # s
        points[0, 1] = ic[1]  # theta (Boozer)
        points[0, 2] = ic[2]  # phi (Boozer)
        
        # Initial parallel velocity (v_par/v = 1.0 for purely parallel motion)
        vpar_init = [ic[3] * vpar0]

        # Trace particle
        # Use similar parameters to firm3d example
        reltol = 1e-8
        abstol = 1e-8
        dt_save = 1e-7  # Save every 0.1 microsecond
        
        traj_booz, res_hits = trace_particles_boozer(
            field,
            points,
            vpar_init,
            tmax=trace_time,
            mass=mass,
            charge=charge,
            Ekin=Ekin,
            stopping_criteria=[MaxToroidalFluxStoppingCriterion(1.0)],
            forget_exact_path=False,
            dt_save=dt_save,
            abstol=abstol,
            reltol=reltol,
        )

        # Extract trajectory data
        traj = traj_booz[0]  # First (and only) particle
        
        # Convert trajectory to cylindrical coordinates
        R = []
        Z = []
        phi_cyl = []
        
        for i in range(len(traj)):
            s_val = traj[i][0]
            theta_booz = traj[i][1]
            phi_booz = traj[i][2]
            
            # Get cylindrical coordinates
            # Note: firm3d uses different coordinate system than SIMPLE
            # We need to evaluate the field at these Boozer coordinates
            xyz = field.evaluate_xyz(s_val, theta_booz, phi_booz)
            r_cyl = np.sqrt(xyz[0]**2 + xyz[1]**2)
            z_cyl = xyz[2]
            phi_c = np.arctan2(xyz[1], xyz[0])
            
            R.append(r_cyl)
            Z.append(z_cyl)
            phi_cyl.append(phi_c)

        # Store orbit data
        nt = len(traj)
        time = np.arange(nt) * dt_save
        s = np.array([traj[i][0] for i in range(nt)])
        theta_booz = np.array([traj[i][1] for i in range(nt)])
        phi_booz = np.array([traj[i][2] for i in range(nt)])
        v_par = np.array([traj[i][3] for i in range(nt)])
        
        # Calculate perpendicular velocity from energy conservation
        v_total = vpar0  # Total velocity is constant
        v_perp = np.sqrt(v_total**2 - v_par**2)
        
        # Get magnetic field strength at each point
        B = np.zeros(nt)
        for i in range(nt):
            modB_val = field.modB(s[i], theta_booz[i], phi_booz[i])
            B[i] = modB_val
        
        # Calculate magnetic moment (should be conserved)
        mu = 0.5 * mass * v_perp**2 / B

        orbit_data = {
            "time": time,
            "s": s,
            "theta_booz": theta_booz,
            "phi_booz": phi_booz,
            "v_total": np.full(nt, v_total),
            "v_parallel": v_par,
            "v_perpendicular": v_perp,
            "mu": mu,
            "B": B,
            "R": np.array(R),
            "Z": np.array(Z),
            "phi_cyl": np.array(phi_cyl),
            "initial_conditions": ic,
        }
        all_orbits_data.append(orbit_data)

        print(f"  Final s: {s[-1]:.3f}")
        print(f"  Mu conservation: {mu[-1] / mu[0]:.6f}")

    # Create xarray dataset
    if n_orbits == 1:
        # Single orbit
        orbit_data = all_orbits_data[0]
        ds = xr.Dataset(
            {
                "s": (
                    ["time"],
                    orbit_data["s"],
                    {"long_name": "Normalized toroidal flux", "units": ""},
                ),
                "theta_booz": (
                    ["time"],
                    orbit_data["theta_booz"],
                    {"long_name": "Boozer poloidal angle", "units": "rad"},
                ),
                "phi_booz": (
                    ["time"],
                    orbit_data["phi_booz"],
                    {"long_name": "Boozer toroidal angle", "units": "rad"},
                ),
                "v_total": (
                    ["time"],
                    orbit_data["v_total"],
                    {"long_name": "Total velocity", "units": "m/s"},
                ),
                "v_parallel": (
                    ["time"],
                    orbit_data["v_parallel"],
                    {"long_name": "Parallel velocity", "units": "m/s"},
                ),
                "v_perpendicular": (
                    ["time"],
                    orbit_data["v_perpendicular"],
                    {"long_name": "Perpendicular velocity", "units": "m/s"},
                ),
                "mu": (
                    ["time"],
                    orbit_data["mu"],
                    {"long_name": "Magnetic moment", "units": "J/T"},
                ),
                "B": (
                    ["time"],
                    orbit_data["B"],
                    {"long_name": "Magnetic field strength", "units": "T"},
                ),
                "R": (
                    ["time"],
                    orbit_data["R"],
                    {"long_name": "Major radius", "units": "m"},
                ),
                "Z": (
                    ["time"],
                    orbit_data["Z"],
                    {"long_name": "Vertical coordinate", "units": "m"},
                ),
                "phi_cyl": (
                    ["time"],
                    orbit_data["phi_cyl"],
                    {"long_name": "Cylindrical toroidal angle", "units": "rad"},
                ),
            },
            coords={
                "time": (
                    ["time"],
                    orbit_data["time"],
                    {"long_name": "Time", "units": "s"},
                ),
            },
            attrs={
                "description": "Single particle orbit trace from firm3d",
                "vmec_file": vmec_file,
                "initial_s": float(orbit_data["initial_conditions"][0]),
                "initial_theta": float(orbit_data["initial_conditions"][1]),
                "initial_phi": float(orbit_data["initial_conditions"][2]),
                "initial_v_parallel": float(orbit_data["initial_conditions"][3]),
                "particle_type": "alpha",
                "energy_MeV": 3.5,
                "created": datetime.now().isoformat(),
                "code": "firm3d",
            },
        )
    else:
        # Multiple orbits - similar structure but with orbit dimension
        # Implementation for multiple orbits would go here
        pass

    # Save to netCDF
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = os.path.join(output_dir, f"trace_orbit_firm3d_{timestamp}.nc")
    
    ds.to_netcdf(output_file)
    print(f"\nSaved orbit data to: {output_file}")

    # Clean up temporary boozmn file
    if os.path.exists(boozmn_file):
        os.remove(boozmn_file)

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
        print("Usage: python trace_orbit_firm3d.py [path/to/wout.nc]")
        sys.exit(1)

    # Read initial conditions
    print("Reading initial conditions from initial_condition/ directory...")
    initial_conditions = read_initial_conditions()

    # Run orbit tracing
    ds = trace_orbit(vmec_file, initial_conditions)

    # Print summary
    print("\nOrbit summary:")
    print(f"  Time range: {ds.time.values[0]:.2e} - {ds.time.values[-1]:.2e} s")
    print(f"  s range: {ds.s.min().values:.3f} - {ds.s.max().values:.3f}")
    print(f"  Mu conservation: {(ds.mu[-1] / ds.mu[0]).values:.6f}")