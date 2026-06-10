# W7-X high mirror: symplectic vs RK loss discrepancy

Externally contributed debug case for itpplasma/SIMPLE#370. SIMPLE's
symplectic integrators (canonical and Boozer coordinates) produce
non-prompt energetic-particle losses on this configuration that the RK
schemes do not show. The discrepancy survives the timestep resolution
scans below.

## Equilibrium

`wout_w7x_high_mirror.nc`: W7-X high-mirror configuration, R0 = 5.49 m,
mirror ratio B01/B00 = 0.096 on axis, beta = 3 percent, ns = 501. The
VMEC input is given in the appendix of Gonzalez-Jerez et al., J. Plasma
Phys. 88, 905880310 (2022), doi:10.1017/S0022377822000393.

## Runs

1000 alpha particles (3.5 MeV, `facE_al = 1.0`) started on s = 0.25,
traced 0.1 s. Spline orders `ns_s = ns_tp = 3`, `multharm = 5` in all
runs. Each directory holds `simple.in`, the shared `start.dat`, and the
outputs (`confined_fraction.dat`, `times_lost.dat`, `fort.6601` Newton
diagnostics, `log`).

| Directory prefix | Integrator | Scan |
|---|---|---|
| 000 | symplectic, canonical (`integmode = 1`, `isw_field_type = 0`) | npoiper2 = 1024, 2048, 5096 |
| 001 | symplectic, Boozer (`integmode = 1`, `isw_field_type = 2`) | npoiper2 = 1024, 2048, 5096 |
| 002 | RK, VMEC coordinates (`integmode = -1`) | relerr = 1e-10, 1e-12 |
| 003 | RK, Boozer coordinates (`integmode = 0`, `isw_field_type = 2`) | relerr = 1e-10, 1e-12, 1e-14 |

The run directories reference the equilibrium as `wout.nc`. To rerun,
link the shared file first:

```
cd 000-canonical_symplectic_npoinper2_1024_npart_1000
ln -s ../wout_w7x_high_mirror.nc wout.nc
simple.x
```

`plot_combined_losses_vs_time.py` reproduces
`confined_fraction_plot.png` from the stored outputs: confined fraction
over time for all integrators. The symplectic curves keep losing
particles on the 1e-3 to 1e-1 s timescale; the RK curves flatten after
the prompt phase.

## Field representation check

`field_representation/` contains a standalone test
(`test_field_representation_error.f90`) comparing SIMPLE's splined
field against the VMEC Fourier representation, with the resulting data
and plots. modB spline interpolation agrees in all coordinate systems;
a consistency test for the field derivatives across coordinate systems
is still open.
