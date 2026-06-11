"""
Companion to plot_combined_losses_vs_time.py (the contributed 000-003 runs):
the same lost-fraction-vs-time plot for the fixed-start runs 100-114, which
all trace the identical 1000 particles for 0.01 s.

Produces confined_fraction_plot_fixed_start.png.
"""

import pathlib

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent

RUNS = [
    ("100-sympl_boozer_mh5_ns3", "b", "-", "sympl Boozer, old clamp"),
    ("105-sympl_canonical_mh5_ns3", "b", "--", "sympl canonical, old clamp"),
    ("113-sympl_boozer_mh5_ns3_axisfix", "r", "-",
     "sympl Boozer, axis chart switch (PR 384)"),
    ("114-sympl_canonical_mh5_ns3_axisfix", "r", "--",
     "sympl canonical, axis chart switch (PR 384)"),
    ("110-chartmap_booz_xform_sympl", "c", "-",
     "sympl Boozer, booz_xform chartmap field"),
    ("111-chartmap_booz_xform_rk", "c", "--",
     "RK Boozer, booz_xform chartmap field"),
    ("101-rk_boozer_relerr_1e-12_mh5_ns3", "g", "-", "RK Boozer 1e-12"),
]

fig, ax = plt.subplots(figsize=(8, 5))
for path, c, ls, label in RUNS:
    f = HERE / path / "confined_fraction.dat"
    if not f.exists():
        print(f"skip {path}")
        continue
    data = np.loadtxt(f)
    time = data[:, 0]
    lost_total = 1.0 - data[:, 1] - data[:, 2]
    ax.plot(time, lost_total, c=c, linestyle=ls, label=label)

ax.set_xlabel("Time $t$ [s]")
ax.set_ylabel("Lost fraction")
ax.set_xscale("log")
ax.set_xlim(1e-6, 1e-2)
ax.set_ylim([0.0, 0.4])
ax.legend(fontsize=8)
ax.grid(True, which="both", ls="--", alpha=0.4)
ax.set_title("Lost particles vs. time, fixed start.dat (runs 100-114)")
fig.tight_layout()
out = HERE / "confined_fraction_plot_fixed_start.png"
fig.savefig(out, dpi=150)
print(f"Saved: {out}")
