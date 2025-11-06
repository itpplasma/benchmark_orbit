import sys
import numpy as np
import time
from firm3d.field.coordinates import boozer_to_vmec

from firm3d.field.boozermagneticfield import (
    BoozerRadialInterpolant,
    InterpolatedBoozerField,
)
from firm3d.field.tracing import (
    trace_particles_boozer,
    MaxToroidalFluxStoppingCriterion,
)
from firm3d.util.constants import (
    ALPHA_PARTICLE_MASS,
    ALPHA_PARTICLE_CHARGE,
    FUSION_ALPHA_PARTICLE_ENERGY,
)
from firm3d.util.functions import proc0_print

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
except ImportError:
    comm = None

boozmn_filename = "boozmn_LandremanPaul2021_QH_reactorScale_lowres_reference.nc"
wout_filename = "wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc"
order = 3  # Order for radial interpolation
degree = 3  # Degree for 3d interpolation
Ekin = 3.5e6 * 1.6022e-19
mass = 4 * 1.6726e-27
charge = 2 * 4.8032e-10 / 2.9979e9
tmax = 1e-3
resolution = 96
tol = 1e-10

# Initialize single trapped particle 
v0 = np.sqrt(2 * Ekin / mass)  # Initial speed
vpar_init = [0.0 * v0, v0]
s0 = np.loadtxt('s_booz.txt')
theta0 = np.loadtxt('vartheta_booz.txt')
zeta0 = np.loadtxt('zeta_booz.txt')
points = np.zeros((2, 3))
points[:, 0] = s0  # s values
points[:, 1] = theta0  # theta values
points[:, 2] = zeta0  # zeta values

## Setup radial interpolation
bri = BoozerRadialInterpolant(boozmn_filename, order, no_K=False, comm=comm)
nfp = bri.nfp

sys.stdout = open(f"stdout_trajectory_resolution_scan_trapped.txt", "a", buffering=1)

time1 = time.time()

ns_interp = resolution
ntheta_interp = resolution
nzeta_interp = resolution

## Setup 3d interpolation
field = InterpolatedBoozerField(
    bri,
    degree,
    ns_interp=ns_interp,
    ntheta_interp=ntheta_interp,
    nzeta_interp=nzeta_interp,
)

time2 = time.time()
proc0_print('Time to set up field: ', time2 - time1, ' seconds')

reltol = tol
abstol = tol

proc0_print('Running firm3d with tol =', tol, ', tmax =', tmax, ', resolution =', resolution)

time2 = time.time()

## Trace alpha particle in Boozer coordinates until it hits the s = 1 surface
## Set forget_exact_path=False to save the trajectory information.
## Set the dt_save parameter to the time interval for trajectory data
## to be saved.
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

time3 = time.time()

proc0_print('Time for integration: ', time3 - time2, ' seconds')

np.savetxt(f'trajectory_data_tol_{tol}_resolution_{resolution}_tmax_{tmax}_trapped.txt', traj_booz[0])

points_booz = np.zeros((np.shape(traj_booz[0])[0], 3))
points_booz[:, 0] = traj_booz[0][:, 1]
points_booz[:, 1] = traj_booz[0][:, 2]
points_booz[:, 2] = traj_booz[0][:, 3]

points_vmec = boozer_to_vmec(wout_filename, field, points_booz)

np.savetxt(f'trajectory_data_vmec_tol_{tol}_resolution_{resolution}_tmax_{tmax}_trapped.txt', np.column_stack((traj_booz[0][:, 0], traj_booz[0][:,1], points_vmec[:, 1], points_vmec[:, 2], traj_booz[0][:, 4])))

points_booz = np.zeros((np.shape(traj_booz[1])[0], 3))
points_booz[:, 0] = traj_booz[1][:, 1]
points_booz[:, 1] = traj_booz[1][:, 2]
points_booz[:, 2] = traj_booz[1][:, 3]

points_vmec = boozer_to_vmec(wout_filename, field, points_booz)

np.savetxt(f'trajectory_data_vmec_tol_{tol}_resolution_{resolution}_tmax_{tmax}_passing.txt', np.column_stack((traj_booz[1][:, 0], traj_booz[1][:,1], points_vmec[:, 1], points_vmec[:, 2], traj_booz[1][:, 4])))