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

## Fixed-start follow-up runs (100-112)

The 000-003 runs do not share initial conditions: `startmode = 1` with
`deterministic = .False.` regenerates `start.dat` per run, so the
confined fractions carry sampling noise of about 1.5 percentage points
(1 sigma at N = 1000). Runs 100-112 all trace the same 1000 particles
(`startmode = 2`, start.dat copied from run 001-1024, converted to
Boozer angles for the chartmap runs) for 0.01 s.

| Run | Field | Integrator | Confined at 0.01 s |
|---|---|---|---|
| 100 | VMEC->Boozer splines, mh5 ns3 | symplectic Euler1, npoiper2 1024 | 0.679 |
| 101 | same | RK Boozer, relerr 1e-12 | 0.729 |
| 102 | mh5 ns5 | symplectic | 0.679 |
| 103 | mh7 ns5 | symplectic | 0.678 |
| 104 | mh3 ns3 | symplectic | 0.685 |
| 105 | canonical-flux field, mh5 ns3 | symplectic | 0.682 |
| 106 | mh5 ns3 | symplectic, npoiper2 256 | 0.681 |
| 107-109 | mh5 ns3 | symplectic, near-axis clamp 0.01 / abs / 1e-4 | 0.679 / 0.693 / 0.689 |
| 110 | booz_xform chartmap (128x97x97) | symplectic | 0.703 |
| 111 | same | RK Boozer, relerr 1e-12 | 0.710 |
| 112 | SIMPLE-export chartmap (501x55x300) | symplectic | 0.679 |

Per particle (identical starts), runs 100 vs 101: 270 particles are
lost in both, with loss times agreeing to 0.13 percent median. The
whole discrepancy is 51 particles lost only in the symplectic run
(41 trapped) against 1 lost only in RK. On the booz_xform chartmap
field (110 vs 111) the same comparison gives 290 shared losses and
7 symplectic-only losses.

Mechanism: the symplectic integrators have no near-axis treatment.
When the Newton solve drives r = s below zero, `orbit_symplectic.f90`
sets r = 0.01 and continues; the RK path instead switches to a
pseudo-Cartesian map inside s = 0.01 (`orbit_timestep_axis`). On this
field the clamp fires 606 times per 1000 particles (`r_negative`
counter), and replacing the constant by reflection r = |r| moves the
confined fraction by 1.4 points (runs 107-109, env knob in an
experiment build). On the booz_xform chartmap field the clamp fires 21
times and symplectic agrees with RK; on the SIMPLE-export chartmap,
which routes SIMPLE's own Boozer transform through the identical
chartmap reader and splines, the pathology returns (run 112,
r_negative 528). The near-axis quality of the internal VMEC-to-Boozer
field data drives the spurious non-prompt losses; the loss level is
insensitive to multharm, spline order, and npoiper2 (runs 102-106).

To reproduce: `boozmn_w7x_high_mirror.nc` from booz_xform (mboz =
nboz = 64) on the shared wout, converted with SIMPLE
`tools/booz_xform_to_boozer_chartmap.py --nrho 128 --ntheta 96
--nzeta 96`; `start_boozer_w7x.dat` and `chartmap_simple_mh5_ns3.nc`
from `export_boozer_chartmap_tool.x`. The .nc files are not committed
(165 MB boozmn, 38 and 266 MB chartmaps). `plot_issue370_comparison.py`
makes `issue370_comparison.png` and `issue370_perparticle.png`.
