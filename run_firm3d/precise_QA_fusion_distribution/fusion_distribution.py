import time

import numpy as np

from firm3d.field.boozermagneticfield import (
    BoozerRadialInterpolant,
    InterpolatedBoozerField,
)
from firm3d.field.tracing import (
    MaxToroidalFluxStoppingCriterion,
    trace_particles_boozer,
)
from firm3d.field.tracing_helpers import (
    initialize_position_profile,
    initialize_velocity_uniform,
)
from firm3d.util.constants import (
    ALPHA_PARTICLE_CHARGE,
    ALPHA_PARTICLE_MASS,
    FUSION_ALPHA_PARTICLE_ENERGY,
)
from firm3d.util.functions import proc0_print, setup_logging
from firm3d.util.mpi import comm_size, comm_world, verbose
from firm3d.field.coordinates import boozer_to_vmec

time1 = time.time()

resolution = 48  # Resolution for field interpolation
nParticles = 5000  # Number of particles to trace
reltol = 1e-8  # Relative tolerance for the ODE solver
abstol = 1e-8  # Absolute tolerance for the ODE solver
order = 3  # Order for radial interpolation
degree = 3  # Degree for 3d interpolation
boozmn_filename = "boozmn_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"
wout_filename = "wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"
tmax = 1e-2  # Time for integration
ns_interp = resolution
ntheta_interp = resolution
nzeta_interp = resolution

# Setup logging to redirect output to file
setup_logging(f"stdout_{nParticles}_{resolution}_{comm_size}.txt")

## Setup radial interpolation
bri = BoozerRadialInterpolant(boozmn_filename, order, no_K=True, comm=comm_world)

## Setup 3d interpolation
field = InterpolatedBoozerField(
    bri,
    degree,
    ns_interp=ns_interp,
    ntheta_interp=ntheta_interp,
    nzeta_interp=nzeta_interp,
)

# Define fusion birth distribution
# Bader, A., et al. "Modeling of energetic particle transport in optimized
# stellarators." Nuclear Fusion 61.11 (2021): 116060.
nD = lambda s: (1 - s**5)  # Normalized density
nT = nD
T = lambda s: 11.5 * (1 - s)  # Temperature in keV


# D-T cross-section
def sigmav(T):
    if T > 0:
        return T ** (-2 / 3) * np.exp(-19.94 * T ** (-1 / 3))
    else:
        return 0


# Reactivity profile
reactivity = lambda s: nD(s) * nT(s) * sigmav(T(s))

points = initialize_position_profile(field, nParticles, reactivity, comm=comm_world, seed=0)

# Translate points to VMEC coordinates for simple 
points_vmec = boozer_to_vmec(wout_filename, field, points)
s_vmec = points_vmec[:, 0]
theta_vmec = points_vmec[:, 1]
zeta_vmec = points_vmec[:, 2]

Ekin = 3.5e6 * 1.6022e-19
mass = 4 * 1.6726e-27
charge = 2 * 4.8032e-10 / 2.9979e9
# Initialize uniformly distributed parallel velocities
vpar0 = np.sqrt(2 * Ekin / mass)
vpar_init = initialize_velocity_uniform(vpar0, nParticles, comm=comm_world, seed=0)

vpar_scaled = vpar_init / vpar0

simple_data = np.column_stack((s_vmec, theta_vmec, zeta_vmec, np.ones_like(s_vmec), vpar_scaled))
np.savetxt("simple_data.txt", simple_data)

## Trace alpha particles in Boozer coordinates until they hit the s = 1 surface
res_tys, res_zeta_hits = trace_particles_boozer(
    field,
    points,
    vpar_init,
    tmax=tmax,
    mass=mass,
    charge=charge,
    comm=comm_world,
    Ekin=Ekin,
    stopping_criteria=[MaxToroidalFluxStoppingCriterion(1.0)],
    forget_exact_path=True,
    abstol=abstol,
    reltol=reltol,
)

time2 = time.time()
proc0_print("Elapsed time for tracing = ", time2 - time1)

## Post-process results to obtain lost particles
if verbose:
    from firm3d.field.trajectory_helpers import compute_loss_fraction

    times, loss_frac = compute_loss_fraction(res_tys, tmin=1e-5, tmax=1e-2)

    np.savetxt("loss_frac.txt", np.column_stack((times, loss_frac)))
