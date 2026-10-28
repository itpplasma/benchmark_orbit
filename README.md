# benchmark_orbit

Benchmark of guiding-center orbit tracers in stellarator equilibria. The same
initial conditions are traced with [SIMPLE](https://github.com/itpplasma/SIMPLE),
[firm3d](https://github.com/ejpaul/firm3d) and ASCOT5, and the resulting orbits
and alpha-particle loss fractions are compared.

Status: research scripts and run data; not a packaged tool.

## Contents

- `booz_xform/`: VMEC equilibria (`wout_*.nc`) of the Landreman-Paul 2021
  precise QA and QH configurations, their Boozer transforms (`boozmn_*.nc`) and
  `run_booz_xform.py` to regenerate them
- `initial_condition/`: initial conditions and their conversion from VMEC to
  Boozer coordinates (`translate_initial_condition.py`, uses simsopt)
- `run_simple/`, `run_firm3d/`, `run_ascot/`: input files, run scripts, plotting
  scripts and stored results for each code
- `compare_orbits/`: Jupyter notebooks comparing single orbits
  (`orbit/compare_traj.ipynb`) and loss fractions
  (`loss_frac/plot_loss_fraction.ipynb`)
- `scripts/`, `src/`: standalone orbit tracing with SIMPLE (Python via
  `pysimple` and Fortran) and firm3d
- `initial_to_simple.py`: writes `start.dat` for SIMPLE from firm3d trajectory
  files in `run_firm3d/trapped/`

## Setup

`make build` installs the Python requirements (`codes/requirements.txt`),
clones SIMPLE and firm3d into `codes/` and installs both with `pip install -e`:

    make build

ASCOT5 is not installed by these scripts.

## Usage

    make run                    # SIMPLE orbit via scripts/trace_orbit_simple.py
    make run-fortran            # Fortran tracer src/trace_orbit_simple.f90
    make help

`VMEC_FILE=path/to/wout.nc` selects another equilibrium. The Fortran target
links against the SIMPLE build in `codes/SIMPLE/build`; the library paths in
the `Makefile` assume Homebrew on macOS and may need adjusting.

The per-code directories have their own Makefiles, for example
`run_ascot/100x100` (set up, run and extract ASCOT5 orbits from SIMPLE's
`start.dat`) and `run_simple/100x100` and `run_firm3d/trapped` (plotting).
