#!/usr/bin/env python3
"""
Trace particle orbits in VMEC equilibrium using SIMPLE.
Reads initial conditions from files and saves orbit data to netCDF file using xarray.

Renamed to trace_orbit_simple.py for clarity.
"""

import os
import sys
from datetime import datetime

import numpy as np
import xarray as xr

# Add SIMPLE to path if needed
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "codes", "SIMPLE"))

from pysimple import get_can_sub as coord
from pysimple import orbit_symplectic, params, simple, simple_main


def read_initial_conditions(ic_dir="initial_condition"):
    """
    Read initial conditions from files in the initial_condition directory.

    Returns:
    --------
    initial_conditions : list of arrays
        Each array contains [s, theta, phi, v/v_th, v_par/v]
    """
    # Read coordinate files
    s_file = os.path.join(ic_dir, "s_booz.txt")
    theta_file = os.path.join(ic_dir, "theta_vmec.txt")
    phi_file = os.path.join(ic_dir, "phi_vmec.txt")

    # Check if files exist
    if not all(os.path.exists(f) for f in [s_file, theta_file, phi_file]):
        print(f"Warning: Initial condition files not found in {ic_dir}")
        print("Using default values: s=0.5, theta=0.5, phi=0.5")
        return [np.array([0.5, 0.5, 0.5, 1.0, 0.0])]

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
        # s, theta, phi, v/v_th, v_par/v
        # Using v/v_th = 1.0, v_par/v = 0.0 (trapped particle - purely perpendicular motion)
        ic = np.array([s_values[i], theta_values[i], phi_values[i], 1.0, 0.0])
        initial_conditions.append(ic)

    return initial_conditions


def trace_orbit(vmec_file, initial_conditions, output_dir="run", trace_time=1e-3):
    """
    Trace particle orbits and save to netCDF.

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

    # Initialize tracer
    tracy = simple.Tracer()

    # Initialize field from VMEC file using Boozer coordinates (mode 2)
    simple_main.init_field(tracy, vmec_file, 3, 3, 3, 2)
    params.params_init()

    # Verify we're using the same setup as SIMPLE example
    print(f"  Field integration mode: {tracy.integmode}")
    print(f"  Field period: {tracy.fper:.6f}")

    n_orbits = len(initial_conditions)
    print(f"Tracing {n_orbits} orbit(s) using SIMPLE's default parameters")
    print("  Key parameters:")
    print(f"    RT0 (major radius): {params.rt0:.1f} cm")
    print(f"    v0 (particle velocity): {params.v0:.2e} cm/s")
    print(f"    trace_time: {params.trace_time:.2e} s")
    print(f"    ntimstep: {params.ntimstep}")
    print(f"    npoiper2: {params.npoiper2}")
    print(f"    tau (normalized): {params.tau:.2e}")
    print(f"    dtau: {params.dtau:.2e}")
    print(f"    dtaumin: {params.dtaumin:.3e}")
    print(
        f"    Expected dtaumin_initial: {2 * 3.14159 * params.rt0 / params.npoiper2:.2e}"
    )

    all_orbits_data = []

    for orbit_idx, z0_vmec in enumerate(initial_conditions):
        print(f"\nOrbit {orbit_idx + 1}/{n_orbits}:")
        print(
            f"  Initial conditions: s={z0_vmec[0]:.3f}, theta={z0_vmec[1]:.3f}, phi={z0_vmec[2]:.3f}"
        )

        z0_can = z0_vmec.copy()

        # Convert VMEC to canonical coordinates
        z0_can[1:3] = coord.vmec_to_can(z0_vmec[0], z0_vmec[1], z0_vmec[2])

        # Use SIMPLE's approach: fixed number of timesteps (exactly like example.py)
        nt = 10001  # Same as SIMPLE example

        print(f"  Orbit {orbit_idx + 1}: {nt} timesteps, dtaumin: {params.dtaumin:.2e}")
        print(f"  B = {tracy.f.bmod:.3f} T")

        # Initialize symplectic integrator exactly like SIMPLE example
        simple.init_sympl(
            tracy.si,
            tracy.f,
            z0_can,
            params.dtaumin,
            params.dtaumin,
            1e-13,
            tracy.integmode,
        )

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

        # Store initial values (convert normalized time to physical time)
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
        for kt in range(nt - 1):
            # Advance one timestep
            orbit_symplectic.orbit_timestep_sympl_expl_impl_euler(tracy.si, tracy.f)

            # Store time (normalized time units like SIMPLE example)
            # Each step advances by dtaumin in normalized time
            time[kt + 1] = (kt + 1) * params.dtaumin

            # Get integrated variables
            s[kt + 1] = tracy.si.z[0]
            theta_can[kt + 1] = tracy.si.z[1]
            phi_can[kt + 1] = tracy.si.z[2]
            pphi[kt + 1] = tracy.si.z[3]

            # Convert canonical to VMEC coordinates
            theta_vmec[kt + 1], phi_vmec[kt + 1] = coord.can_to_vmec(
                s[kt + 1], theta_can[kt + 1], phi_can[kt + 1]
            )

            # Store velocities and invariants
            v_par[kt + 1] = tracy.f.vpar
            mu[kt + 1] = tracy.f.mu
            b_field[kt + 1] = tracy.f.bmod
            v_perp[kt + 1] = np.sqrt(2 * mu[kt + 1] * b_field[kt + 1])
            v_total[kt + 1] = np.sqrt(v_par[kt + 1] ** 2 + v_perp[kt + 1] ** 2)

            # Convert to cylindrical coordinates
            r_cyl[kt + 1], z_cyl[kt + 1] = coord.vmec_to_cyl(
                s[kt + 1], theta_vmec[kt + 1], phi_vmec[kt + 1]
            )

        # Store orbit data
        orbit_data = {
            "time": time,
            "s": s,
            "theta_vmec": theta_vmec,
            "phi_vmec": phi_vmec,
            "theta_can": theta_can,
            "phi_can": phi_can,
            "v_total": v_total,
            "v_parallel": v_par,
            "v_perpendicular": v_perp,
            "pphi": pphi,
            "mu": mu,
            "B": b_field,
            "R": r_cyl,
            "Z": z_cyl,
            "initial_conditions": z0_vmec,
        }
        all_orbits_data.append(orbit_data)

        print(
            f"  Final s: {s[-1]:.3f}, Energy conservation: {v_total[-1] / v_total[0]:.6f}"
        )

    # Create combined xarray dataset
    if n_orbits == 1:
        # Single orbit - use original format
        orbit_data = all_orbits_data[0]
        ds = xr.Dataset(
            {
                "s": (
                    ["time"],
                    orbit_data["s"],
                    {"long_name": "Normalized toroidal flux", "units": ""},
                ),
                "theta_vmec": (
                    ["time"],
                    orbit_data["theta_vmec"],
                    {"long_name": "VMEC poloidal angle", "units": "rad"},
                ),
                "phi_vmec": (
                    ["time"],
                    orbit_data["phi_vmec"],
                    {"long_name": "VMEC toroidal angle", "units": "rad"},
                ),
                "theta_can": (
                    ["time"],
                    orbit_data["theta_can"],
                    {"long_name": "Canonical poloidal angle", "units": "rad"},
                ),
                "phi_can": (
                    ["time"],
                    orbit_data["phi_can"],
                    {"long_name": "Canonical toroidal angle", "units": "rad"},
                ),
                "v_total": (
                    ["time"],
                    orbit_data["v_total"],
                    {"long_name": "Total velocity", "units": "v_th"},
                ),
                "v_parallel": (
                    ["time"],
                    orbit_data["v_parallel"],
                    {"long_name": "Parallel velocity", "units": "v_th"},
                ),
                "v_perpendicular": (
                    ["time"],
                    orbit_data["v_perpendicular"],
                    {"long_name": "Perpendicular velocity", "units": "v_th"},
                ),
                "pphi": (
                    ["time"],
                    orbit_data["pphi"],
                    {"long_name": "Canonical toroidal momentum", "units": ""},
                ),
                "mu": (
                    ["time"],
                    orbit_data["mu"],
                    {"long_name": "Magnetic moment", "units": ""},
                ),
                "B": (
                    ["time"],
                    orbit_data["B"],
                    {"long_name": "Magnetic field strength", "units": "T"},
                ),
                "R": (
                    ["time"],
                    orbit_data["R"],
                    {"long_name": "Major radius", "units": "cm"},
                ),
                "Z": (
                    ["time"],
                    orbit_data["Z"],
                    {"long_name": "Vertical coordinate", "units": "cm"},
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
                "description": "Single particle orbit trace from SIMPLE",
                "vmec_file": vmec_file,
                "initial_s": float(orbit_data["initial_conditions"][0]),
                "initial_theta": float(orbit_data["initial_conditions"][1]),
                "initial_phi": float(orbit_data["initial_conditions"][2]),
                "initial_v_parallel": float(orbit_data["initial_conditions"][4]),
                "particle_type": "alpha",
                "energy_MeV": 3.5,
                "created": datetime.now().isoformat(),
            },
        )
    else:
        # Multiple orbits - create 2D arrays
        nt = len(all_orbits_data[0]["time"])

        # Stack all orbit data
        s_all = np.array([orbit["s"] for orbit in all_orbits_data])
        theta_vmec_all = np.array([orbit["theta_vmec"] for orbit in all_orbits_data])
        phi_vmec_all = np.array([orbit["phi_vmec"] for orbit in all_orbits_data])
        theta_can_all = np.array([orbit["theta_can"] for orbit in all_orbits_data])
        phi_can_all = np.array([orbit["phi_can"] for orbit in all_orbits_data])
        v_total_all = np.array([orbit["v_total"] for orbit in all_orbits_data])
        v_par_all = np.array([orbit["v_parallel"] for orbit in all_orbits_data])
        v_perp_all = np.array([orbit["v_perpendicular"] for orbit in all_orbits_data])
        pphi_all = np.array([orbit["pphi"] for orbit in all_orbits_data])
        mu_all = np.array([orbit["mu"] for orbit in all_orbits_data])
        b_field_all = np.array([orbit["B"] for orbit in all_orbits_data])
        r_cyl_all = np.array([orbit["R"] for orbit in all_orbits_data])
        z_cyl_all = np.array([orbit["Z"] for orbit in all_orbits_data])

        # Initial conditions for each orbit
        initial_s = np.array(
            [orbit["initial_conditions"][0] for orbit in all_orbits_data]
        )
        initial_theta = np.array(
            [orbit["initial_conditions"][1] for orbit in all_orbits_data]
        )
        initial_phi = np.array(
            [orbit["initial_conditions"][2] for orbit in all_orbits_data]
        )

        ds = xr.Dataset(
            {
                "s": (
                    ["orbit", "time"],
                    s_all,
                    {"long_name": "Normalized toroidal flux", "units": ""},
                ),
                "theta_vmec": (
                    ["orbit", "time"],
                    theta_vmec_all,
                    {"long_name": "VMEC poloidal angle", "units": "rad"},
                ),
                "phi_vmec": (
                    ["orbit", "time"],
                    phi_vmec_all,
                    {"long_name": "VMEC toroidal angle", "units": "rad"},
                ),
                "theta_can": (
                    ["orbit", "time"],
                    theta_can_all,
                    {"long_name": "Canonical poloidal angle", "units": "rad"},
                ),
                "phi_can": (
                    ["orbit", "time"],
                    phi_can_all,
                    {"long_name": "Canonical toroidal angle", "units": "rad"},
                ),
                "v_total": (
                    ["orbit", "time"],
                    v_total_all,
                    {"long_name": "Total velocity", "units": "v_th"},
                ),
                "v_parallel": (
                    ["orbit", "time"],
                    v_par_all,
                    {"long_name": "Parallel velocity", "units": "v_th"},
                ),
                "v_perpendicular": (
                    ["orbit", "time"],
                    v_perp_all,
                    {"long_name": "Perpendicular velocity", "units": "v_th"},
                ),
                "pphi": (
                    ["orbit", "time"],
                    pphi_all,
                    {"long_name": "Canonical toroidal momentum", "units": ""},
                ),
                "mu": (
                    ["orbit", "time"],
                    mu_all,
                    {"long_name": "Magnetic moment", "units": ""},
                ),
                "B": (
                    ["orbit", "time"],
                    b_field_all,
                    {"long_name": "Magnetic field strength", "units": "T"},
                ),
                "R": (
                    ["orbit", "time"],
                    r_cyl_all,
                    {"long_name": "Major radius", "units": "cm"},
                ),
                "Z": (
                    ["orbit", "time"],
                    z_cyl_all,
                    {"long_name": "Vertical coordinate", "units": "cm"},
                ),
                "initial_s": (
                    ["orbit"],
                    initial_s,
                    {"long_name": "Initial s coordinate", "units": ""},
                ),
                "initial_theta": (
                    ["orbit"],
                    initial_theta,
                    {"long_name": "Initial theta coordinate", "units": "rad"},
                ),
                "initial_phi": (
                    ["orbit"],
                    initial_phi,
                    {"long_name": "Initial phi coordinate", "units": "rad"},
                ),
            },
            coords={
                "time": (
                    ["time"],
                    all_orbits_data[0]["time"],
                    {"long_name": "Time", "units": "s"},
                ),
                "orbit": (
                    ["orbit"],
                    np.arange(n_orbits),
                    {"long_name": "Orbit index", "units": ""},
                ),
            },
            attrs={
                "description": f"Multiple particle orbit traces from SIMPLE ({n_orbits} orbits)",
                "vmec_file": vmec_file,
                "n_orbits": n_orbits,
                "particle_type": "alpha",
                "energy_MeV": 3.5,
                "created": datetime.now().isoformat(),
            },
        )

    # Save to netCDF
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    if n_orbits == 1:
        output_file = os.path.join(output_dir, f"trace_orbit_simple_{timestamp}.nc")
    else:
        output_file = os.path.join(
            output_dir, f"trace_orbit_simple_{n_orbits}orbits_{timestamp}.nc"
        )

    ds.to_netcdf(output_file)
    print(f"\nSaved orbit data to: {output_file}")

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
        print("Usage: python trace_orbit_simple.py [path/to/wout.nc]")
        sys.exit(1)

    # Read initial conditions
    print("Reading initial conditions from initial_condition/ directory...")
    initial_conditions = read_initial_conditions()

    # Run orbit tracing
    ds = trace_orbit(vmec_file, initial_conditions)

    # Print summary
    print("\nOrbit summary:")
    print(f"  Time range: {ds.time.values[0]:.2e} - {ds.time.values[-1]:.2e} s")

    if len(initial_conditions) == 1:
        # Single orbit
        print(f"  s range: {ds.s.min().values:.3f} - {ds.s.max().values:.3f}")
        print(f"  Energy conservation: {(ds.v_total[-1] / ds.v_total[0]).values:.6f}")
        print(f"  Mu conservation: {(ds.mu[-1] / ds.mu[0]).values:.6f}")
    else:
        # Multiple orbits
        print(f"  Number of orbits: {len(initial_conditions)}")
        print(
            f"  s range (all orbits): {ds.s.min().values:.3f} - {ds.s.max().values:.3f}"
        )
        # Energy conservation for each orbit
        energy_conservation = ds.v_total.isel(time=-1) / ds.v_total.isel(time=0)
        mu_conservation = ds.mu.isel(time=-1) / ds.mu.isel(time=0)
        print(
            f"  Energy conservation (min/max): {energy_conservation.min().values:.6f} / {energy_conservation.max().values:.6f}"
        )
        print(
            f"  Mu conservation (min/max): {mu_conservation.min().values:.6f} / {mu_conservation.max().values:.6f}"
        )
