"""
Plot number of lost particles as a function of time.

Reads confined_fraction.dat produced by SIMPLE.

Column layout (0-based):
  0  time         [s]
  1  confpart_pass  (confined passing particles)
  2  confpart_trap  (confined trapped particles)
  3  ntestpart      (total particles)

Usage:
    python plot_losses_vs_time.py [confined_fraction.dat]
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import pathlib


fig, ax = plt.subplots(figsize=(8, 5))
fig.tight_layout()
c = ['b', 'b', 'b', 'r', 'r', 'r', 'c', 'c', 'g', 'g', 'g']
l = ['--', '-.', '-', '--', '-.', '-', '--', '-.', '--', '-.', '-']
#for i, path in enumerate(["000-canonical_symplectic_npoinper2_512", "000-canonical_symplectic_npoinper2_1024", "000-canonical_symplectic_npoinper2_2048", "001-boozer_symplectic_npoinper2_512", "001-boozer_symplectic_npoinper2_1028", "001-boozer_symplectic_npoinper2_2048", "002-w7x_runge_kutta_relerr_1e-8", "002-w7x_runge_kutta_relerr_1e-10", "002-w7x_runge_kutta_relerr_1e-12", "002-w7x_runge_kutta_relerr_1e-14"]):
for i, path in enumerate(["000-canonical_symplectic_npoinper2_1024_npart_1000",
                          "000-canonical_symplectic_npoinper2_2048_npart_1000",
                          "000-canonical_symplectic_npoinper2_5096_npart_1000",
                          "001-boozer_symplectic_npoinper2_1024_npart_1000",
                          "001-boozer_symplectic_npoinper2_2048_npart_1000",
                          "001-boozer_symplectic_npoinper2_5096_npart_1000",
                          "002-runge_kutta_relerr_1e-10_npart_1000",
                          "002-runge_kutta_relerr_1e-12_npart_1000",
                          "003-runge_kutta_relerr_1e-10_boozer_npart_1000",
                          "003-runge_kutta_relerr_1e-12_boozer_npart_1000",
                          "003-runge_kutta_relerr_1e-14_boozer_npart_1000",
                          ]):

    
    fname = "confined_fraction.dat"
    data = np.loadtxt(pathlib.Path(path) / fname)
    
    time         = data[:, 0]
    conf_pass    = data[:, 1]
    conf_trap    = data[:, 2]
    
    lost_total   = 1.0 - conf_pass - conf_trap
    lost_pass    = 1.0 - conf_pass
    lost_trap    = 1.0 - conf_trap
    
    
    ax.plot(time, lost_total, c=c[i], linestyle=l[i],   label=path)
    
    ax.set_xlabel("Time $t$ [s]")
    ax.set_ylabel("Lost fraction")
    ax.set_xscale("log")
    ax.set_xlim(0.1*time[time > 0].min(), time.max())
    ax.set_ylim([0.0, 0.4])
    ax.legend()
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.set_title("Lost particles vs. time")
    
out = fname.replace(".dat", "_plot.png")
fig.savefig(out, dpi=150)
print(f"Saved: {out}")
plt.show()
