"""
Visualise field representation errors across phi planes.

Reads field_representation_data.dat produced by test_field_representation_error.x.

Columns:
  0  r (=sqrt(s))
  1  theta   [rad]
  2  phi     [rad]
  3  R_cm    [cm]
  4  Z_cm    [cm]
  5  B_Fourier  [G]
  6  B_VMEC_spl [G]
  7  B_Boozer   [G]
  8  B_Canonical[G]

Usage:
    python plot_field_representation_error.py [field_representation_data.dat]
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.colors import Normalize

fname = sys.argv[1] if len(sys.argv) > 1 else 'field_representation_data.dat'
data = np.loadtxt(fname, comments='#')

r_arr     = data[:, 0]
theta_arr = data[:, 1]
phi_arr   = data[:, 2]
R_cm      = data[:, 3]
Z_cm      = data[:, 4]

labels = ['Fourier ref', 'VMEC spline', 'Boozer', 'Canonical']
bmods  = [data[:, 5], data[:, 6], data[:, 7], data[:, 8]]

# Relative error of each spline representation vs Fourier reference
B_ref  = data[:, 5]
errors = [
    np.abs(data[:, 6] - B_ref) / np.maximum(B_ref, 1e-30),
    np.abs(data[:, 7] - B_ref) / np.maximum(B_ref, 1e-30),
    np.abs(data[:, 8] - B_ref) / np.maximum(B_ref, 1e-30),
]
err_labels = ['VMEC spl - Fourier', 'Boozer - Fourier', 'Canonical - Fourier']

# Identify unique phi planes (group within 1e-8 rad)
phi_unique = np.unique(np.round(phi_arr, 8))
n_phi = len(phi_unique)

# ── Figure 1: |B| in each phi plane for all 4 representations ────────────────
ncols = min(n_phi, 6)
nrows_per_rep = 1
fig1, axes1 = plt.subplots(len(labels), ncols,
                            figsize=(3 * ncols, 3 * len(labels)),
                            squeeze=False)
fig1.suptitle('|B| [G] in each phi plane — all representations', fontsize=11)

bmod_all = np.concatenate(bmods)
vmin_b, vmax_b = bmod_all.min(), bmod_all.max()

for i_phi, phi_val in enumerate(phi_unique[:ncols]):
    mask = np.abs(phi_arr - phi_val) < 1e-7
    Rp, Zp = R_cm[mask], Z_cm[mask]

    for i_rep, (bmod, lab) in enumerate(zip(bmods, labels)):
        ax = axes1[i_rep, i_phi]
        sc = ax.scatter(Rp / 100, Zp / 100, c=bmod[mask],
                        cmap='plasma', vmin=vmin_b, vmax=vmax_b,
                        s=15, linewidths=0)
        ax.set_aspect('equal')
        ax.set_xlabel('R [m]', fontsize=7)
        ax.set_ylabel('Z [m]', fontsize=7)
        ax.tick_params(labelsize=6)
        if i_phi == 0:
            ax.set_title(lab, fontsize=8, loc='left')
        if i_rep == 0:
            ax.set_title(f'φ={phi_val:.3f} rad', fontsize=8)
        if i_phi == ncols - 1:
            plt.colorbar(sc, ax=ax, label='|B| [G]', fraction=0.046)

fig1.tight_layout()
out1 = fname.replace('.dat', '_bmod.png')
fig1.savefig(out1, dpi=150)
print(f'Saved: {out1}')

# ── Figure 2: relative error in each phi plane for 3 spline representations ──
fig2, axes2 = plt.subplots(len(errors), ncols,
                            figsize=(3 * ncols, 3 * len(errors)),
                            squeeze=False)
fig2.suptitle('Relative |B| error vs Fourier reference', fontsize=11)

err_all  = np.concatenate(errors)
err_pos  = err_all[err_all > 0]
vmin_e   = 10 ** np.floor(np.log10(np.percentile(err_pos, 1)))
vmax_e   = 10 ** np.ceil(np.log10(np.percentile(err_pos, 99)))
log_norm = plt.matplotlib.colors.LogNorm(vmin=vmin_e, vmax=vmax_e)

for i_phi, phi_val in enumerate(phi_unique[:ncols]):
    mask = np.abs(phi_arr - phi_val) < 1e-7
    Rp, Zp = R_cm[mask], Z_cm[mask]

    for i_rep, (err, lab) in enumerate(zip(errors, err_labels)):
        ax = axes2[i_rep, i_phi]
        sc = ax.scatter(Rp / 100, Zp / 100, c=np.maximum(err[mask], vmin_e),
                        cmap='copper', norm=log_norm,
                        s=15, linewidths=0)
        ax.set_aspect('equal')
        ax.set_xlabel('R [m]', fontsize=7)
        ax.set_ylabel('Z [m]', fontsize=7)
        ax.tick_params(labelsize=6)
        if i_phi == 0:
            ax.set_title(lab, fontsize=8, loc='left')
        if i_rep == 0:
            ax.set_title(f'φ={phi_val:.3f} rad', fontsize=8)
        if i_phi == ncols - 1:
            plt.colorbar(sc, ax=ax, label='rel. err.', fraction=0.046)

fig2.tight_layout()
out2 = fname.replace('.dat', '_error.png')
fig2.savefig(out2, dpi=150)
print(f'Saved: {out2}')
plt.show()
