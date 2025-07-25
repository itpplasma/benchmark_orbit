#!/usr/bin/env python3
"""Check velocity values in SIMPLE output."""

import xarray as xr
import numpy as np

# Load the latest SIMPLE output
ds = xr.open_dataset('run/trace_orbit_simple_20250725_122010.nc')

print("Initial Conditions:")
print(f"  s = {ds.attrs.get('initial_s', 'N/A')}")
print(f"  theta = {ds.attrs.get('initial_theta', 'N/A')}")
print(f"  phi = {ds.attrs.get('initial_phi', 'N/A')}")
print(f"  initial_v_parallel = {ds.attrs.get('initial_v_parallel', 'N/A')}")

print("\nVelocity Components at t=0:")
print(f"  v_total[0] = {ds.v_total[0].values}")
print(f"  v_parallel[0] = {ds.v_parallel[0].values}")
print(f"  v_perpendicular[0] = {ds.v_perpendicular[0].values}")

print("\nRatios:")
print(f"  v_par/v_total = {ds.v_parallel[0].values / ds.v_total[0].values}")
print(f"  v_perp/v_total = {ds.v_perpendicular[0].values / ds.v_total[0].values}")
print(f"  v_perp = {ds.v_perpendicular[0].values} (normalized units)")

print("\nFirst 10 timesteps:")
print("Time     v_total    v_par      v_perp     v_par/v_total")
print("-" * 60)
for i in range(min(10, len(ds.time))):
    vt = ds.v_total[i].values
    vpar = ds.v_parallel[i].values
    vperp = ds.v_perpendicular[i].values
    ratio = vpar / vt if vt != 0 else 0
    print(f"{i:4d}    {vt:8.6f}   {vpar:8.6f}   {vperp:8.6f}   {ratio:8.6f}")

# Check if magnetic moment is conserved
print(f"\nMagnetic moment conservation:")
print(f"  mu[0] = {ds.mu[0].values}")
print(f"  mu[-1] = {ds.mu[-1].values}")
print(f"  Ratio = {ds.mu[-1].values / ds.mu[0].values}")

# Check the normalization
print(f"\nNormalization check:")
print(f"  sqrt(v_par^2 + v_perp^2) = {np.sqrt(ds.v_parallel[0].values**2 + ds.v_perpendicular[0].values**2)}")
print(f"  v_total = {ds.v_total[0].values}")
print(f"  v_perp/sqrt(2) = {ds.v_perpendicular[0].values / np.sqrt(2)}")