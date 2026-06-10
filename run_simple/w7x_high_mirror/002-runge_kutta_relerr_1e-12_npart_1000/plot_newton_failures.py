"""
Plot Newton convergence failures per particle vs. loss time.

Reads times_lost.dat produced by SIMPLE (with n_newton_failures column).

Column layout (0-based):
  0  particle index
  1  times_lost   [s]
  2  trap_par     (trapping parameter)
  3  s_start      (initial flux label)
  4  perp_inv     (perpendicular adiabatic invariant)
  5-9 zend(1..5)  (final phase-space coords)
  10 n_newton_failures

Usage:
    python plot_newton_failures.py [times_lost.dat] [trace_time]

    trace_time defaults to the maximum times_lost value found in the file.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

# ── load data ────────────────────────────────────────────────────────────────
fname = sys.argv[1] if len(sys.argv) > 1 else "times_lost.dat"
data = np.loadtxt(fname)

times_lost      = data[:, 1]
n_failures      = data[:, 10].astype(int)

# Particles confined to the end have times_lost == trace_time
trace_time = float(sys.argv[2]) if len(sys.argv) > 2 else times_lost.max()
lost_mask     = times_lost < trace_time
confined_mask = ~lost_mask

# ── summary ──────────────────────────────────────────────────────────────────
n_total    = len(times_lost)
n_lost     = lost_mask.sum()
n_confined = confined_mask.sum()
n_any_fail = (n_failures > 0).sum()

print(f"Particles total:       {n_total}")
print(f"  lost:                {n_lost}")
print(f"  confined:            {n_confined}")
print(f"With ≥1 Newton failure:{n_any_fail}  ({100*n_any_fail/n_total:.1f}%)")
print(f"Max failures (single): {n_failures.max()}")

# ── plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left: scatter of failures vs loss time (lost particles only)
ax = axes[0]
if lost_mask.sum() > 0:
    sc = ax.scatter(
        times_lost[lost_mask],
        n_failures[lost_mask],
        c=n_failures[lost_mask],
        cmap="copper",
        vmin=0, vmax=max(n_failures.max(), 1),
        s=10, alpha=0.7, linewidths=0,
    )
    plt.colorbar(sc, ax=ax, label="Newton failures")
ax.set_xlabel("Loss time $t$ [s]")
ax.set_ylabel("Newton failures (particle lifetime)")
ax.set_xscale("log")
ax.set_title("Newton failures vs. loss time\n(lost particles)")
ax.grid(True, which="both", ls="--", alpha=0.4)

# Right: histogram of failure counts across all particles
ax = axes[1]
bins = np.arange(0, n_failures.max() + 2) - 0.5
ax.hist(n_failures[lost_mask],     bins=bins, alpha=0.7, label=f"Lost ({n_lost})", color="tab:red")
ax.hist(n_failures[confined_mask], bins=bins, alpha=0.7, label=f"Confined ({n_confined})", color="tab:blue")
ax.set_xlabel("Newton failures per particle")
ax.set_ylabel("Count")
ax.set_title("Distribution of Newton failures")
ax.legend()
ax.grid(True, ls="--", alpha=0.4)
if n_failures.max() > 20:
    ax.set_yscale("log")

fig.tight_layout()
out = fname.replace(".dat", "_newton_failures.png")
fig.savefig(out, dpi=150)
print(f"Saved: {out}")
plt.show()
