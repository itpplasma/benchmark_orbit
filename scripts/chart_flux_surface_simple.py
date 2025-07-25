# %%
import matplotlib.pyplot as plt
import numpy as np
from pysimple import get_can_sub as coord
from pysimple import orbit_symplectic, params, simple, simple_main

# %%
tracy = simple.Tracer()

simple_main.init_field(
    tracy,
    "booz_xform/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc",
    3,
    3,
    3,
    1,
)
params.params_init()

# %% Initial conditions
z0_vmec = np.array([0.8, 1.0, 0.2, 1.0, 0.5])  # s, th, ph, v/v_th, v_par/v
z0_can = z0_vmec.copy()  # s, th_c, ph_c, v/v_th, v_par/v

z0_can[1:3] = coord.vmec_to_can(z0_vmec[0], z0_vmec[1], z0_vmec[2])

simple.init_sympl(
    tracy.si, tracy.f, z0_can, params.dtaumin, params.dtaumin, 1e-13, tracy.integmode
)

print(f"B = {tracy.f.bmod}")

call trace_orbit_with_classifiers(tracy, )
