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

fname = sys.argv[1] if len(sys.argv) > 1 else "confined_fraction.dat"
data = np.loadtxt(fname)

time         = data[:, 0]
conf_pass    = data[:, 1]
conf_trap    = data[:, 2]

lost_total   = 1.0 - conf_pass - conf_trap
lost_pass    = 1.0 - conf_pass
lost_trap    = 1.0 - conf_trap

fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(time, lost_total, color="tab:red",    label="Total lost")

ax.set_xlabel("Time $t$ [s]")
ax.set_ylabel("Lost fraction")
ax.set_xscale("log")
ax.set_xlim(time[time > 0].min(), time.max())
ax.set_ylim(bottom=0)
ax.legend()
ax.grid(True, which="both", ls="--", alpha=0.4)
ax.set_title("Lost particles vs. time")

fig.tight_layout()
out = fname.replace(".dat", "_plot.png")
fig.savefig(out, dpi=150)
print(f"Saved: {out}")
plt.show()
