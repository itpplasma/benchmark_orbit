#!/usr/bin/env python3
"""Comparison plots for itpplasma/SIMPLE#370 follow-up runs.

Panel A: confined fraction vs time, integrator and field-resolution sweep
         (fixed start.dat, trace_time 1e-2 s).
Panel B: near-axis clamp variants and chartmap (booz_xform / SIMPLE export)
         same-field cross-checks.
Scatter: per-particle loss times, symplectic vs RK on identical starts.
"""

import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent


def conf(run):
    d = np.loadtxt(HERE / run / "confined_fraction.dat")
    return d[:, 0], d[:, 1] + d[:, 2]


def maybe_plot(ax, run, label, **kw):
    path = HERE / run / "confined_fraction.dat"
    if not path.exists():
        print(f"skip {run}")
        return
    t, c = conf(run)
    mask = t > 0
    ax.semilogx(t[mask], c[mask], label=label, **kw)


def panel_a(ax):
    maybe_plot(ax, "100-sympl_boozer_mh5_ns3", "sympl Boozer mh5 ns3 (baseline)",
               color="tab:red", lw=2)
    maybe_plot(ax, "102-sympl_boozer_mh5_ns5", "sympl Boozer mh5 ns5",
               color="tab:orange")
    maybe_plot(ax, "103-sympl_boozer_mh7_ns5", "sympl Boozer mh7 ns5",
               color="tab:brown")
    maybe_plot(ax, "104-sympl_boozer_mh3_ns3", "sympl Boozer mh3 ns3",
               color="tab:pink")
    maybe_plot(ax, "106-sympl_boozer_mh5_ns3_npoiper2_256",
               "sympl Boozer npoiper2=256", color="tab:red", ls=":")
    maybe_plot(ax, "105-sympl_canonical_mh5_ns3", "sympl canonical mh5 ns3",
               color="tab:purple")
    maybe_plot(ax, "101-rk_boozer_relerr_1e-12_mh5_ns3", "RK Boozer 1e-12 mh5 ns3",
               color="tab:blue", lw=2)
    ax.set_title("Field resolution and integrator (shared start.dat)")


def panel_b(ax):
    maybe_plot(ax, "100-sympl_boozer_mh5_ns3", "sympl, VMEC->Boozer splines",
               color="tab:red", lw=2)
    maybe_plot(ax, "107-sympl_boozer_mh5_ns3_clampdefault",
               "sympl, clamp r=0.01 (exp binary)", color="tab:red", ls="--")
    maybe_plot(ax, "108-sympl_boozer_mh5_ns3_clampabs",
               "sympl, clamp r=|r|", color="tab:green")
    maybe_plot(ax, "109-sympl_boozer_mh5_ns3_clamptiny",
               "sympl, clamp r=1e-4", color="tab:olive")
    maybe_plot(ax, "110-chartmap_booz_xform_sympl",
               "sympl, booz_xform chartmap", color="tab:cyan", lw=2)
    maybe_plot(ax, "111-chartmap_booz_xform_rk",
               "RK, booz_xform chartmap", color="tab:gray")
    maybe_plot(ax, "112-chartmap_simple_sympl",
               "sympl, SIMPLE-export chartmap", color="black", ls="--")
    maybe_plot(ax, "101-rk_boozer_relerr_1e-12_mh5_ns3", "RK Boozer 1e-12",
               color="tab:blue", lw=2)
    ax.set_title("Near-axis clamp and same-field (chartmap) cross-checks")


def scatter_lost(ax, run_a, run_b, label_a, label_b):
    pa = HERE / run_a / "times_lost.dat"
    pb = HERE / run_b / "times_lost.dat"
    if not (pa.exists() and pb.exists()):
        return
    a = np.loadtxt(pa)[:, 1]
    b = np.loadtxt(pb)[:, 1]
    a = np.where(a < 0, 1e-2, a)
    b = np.where(b < 0, 1e-2, b)
    ax.loglog(a, b, ".", ms=3, alpha=0.5)
    ax.plot([1e-6, 1e-2], [1e-6, 1e-2], "k-", lw=0.5)
    ax.set_xlabel(f"loss time {label_a} [s]")
    ax.set_ylabel(f"loss time {label_b} [s]")
    n_sympl_only = np.sum((a < 1e-2) & (b >= 1e-2))
    n_rk_only = np.sum((a >= 1e-2) & (b < 1e-2))
    ax.set_title(f"lost only in {label_a}: {n_sympl_only}, "
                 f"only in {label_b}: {n_rk_only}")


def main():
    fig, axs = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    panel_a(axs[0])
    panel_b(axs[1])
    for ax in axs:
        ax.set_xlabel("time [s]")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, loc="lower left")
    axs[0].set_ylabel("confined fraction")
    fig.tight_layout()
    fig.savefig(HERE / "issue370_comparison.png", dpi=160)
    print("wrote issue370_comparison.png")

    fig2, ax2 = plt.subplots(figsize=(5.5, 5))
    scatter_lost(ax2, "100-sympl_boozer_mh5_ns3",
                 "101-rk_boozer_relerr_1e-12_mh5_ns3",
                 "symplectic", "RK")
    fig2.tight_layout()
    fig2.savefig(HERE / "issue370_perparticle.png", dpi=160)
    print("wrote issue370_perparticle.png")


if __name__ == "__main__":
    main()
