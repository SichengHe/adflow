"""
Plot CL history comparison from saved BDF2 and TSTHETA results.

Usage:
    python plot_comparison.py
"""

import os
import pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Problem parameters (for alpha reference curve)
k = 0.0808
M = 0.6
gamma = 1.4
R = 287.085
T = 280.0
c = 1.0
alpha_m = 2.77
alpha_0 = 2.34
omega = 2 * M * np.sqrt(gamma * R * T) * k / c

f = 10.0
period = 1.0 / f
nStepPerPeriod = 8
dt = period / nStepPerPeriod
t_final = period

baseDir = os.path.dirname(os.path.abspath(__file__))
outputDir = os.path.join(baseDir, "output")

bdf2_pkl = os.path.join(outputDir, "bdf2_results.pkl")
theta_pkl = os.path.join(outputDir, "tstheta_results.pkl")

if not os.path.exists(bdf2_pkl):
    print(f"ERROR: {bdf2_pkl} not found. Run run_bdf2.py first.")
    exit(1)
if not os.path.exists(theta_pkl):
    print(f"ERROR: {theta_pkl} not found. Run run_tstheta.py first.")
    exit(1)

with open(bdf2_pkl, "rb") as fh:
    bdf2 = pickle.load(fh)
with open(theta_pkl, "rb") as fh:
    theta = pickle.load(fh)

# Reference alpha(t)
t_ref = np.linspace(0, t_final, 500)
alpha_ref = alpha_m - alpha_0 * np.sin(omega * t_ref)

fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

# CL
ax = axes[0]
ax.plot(bdf2["time"], bdf2["cl"], "b-o", ms=5, lw=1.5, label="ADflow BDF2 (2nd order)")
ax.plot(theta["time"], theta["cl"], "r-s", ms=5, lw=1.5,
        label=r"PETSc TSTHETA ($\theta$=1, 1st order)")
ax.set_ylabel(r"$C_L$", fontsize=13)
ax.legend(fontsize=11, loc="best")
ax.set_title(
    f"Pitching NACA 0012  |  M={M}, Re=4.8M, "
    rf"$\alpha$={alpha_m}$\pm${alpha_0}$^\circ$, k={k}  |  "
    f"dt={dt:.4e} ({nStepPerPeriod} steps/period)",
    fontsize=11,
)
ax.grid(True, alpha=0.3)

# CD
ax = axes[1]
ax.plot(bdf2["time"], bdf2["cd"], "b-o", ms=5, lw=1.5, label="ADflow BDF2")
ax.plot(theta["time"], theta["cd"], "r-s", ms=5, lw=1.5, label="PETSc TSTHETA")
ax.set_ylabel(r"$C_D$", fontsize=13)
ax.legend(fontsize=11, loc="best")
ax.grid(True, alpha=0.3)

# Alpha
ax = axes[2]
ax.plot(t_ref, alpha_ref, "k-", lw=1.0)
ax.set_ylabel(r"$\alpha$ [deg]", fontsize=13)
ax.set_xlabel("Time [s]", fontsize=13)
ax.grid(True, alpha=0.3)

for a in axes:
    a.axvline(period, color="gray", ls="--", lw=0.8, alpha=0.5)

plt.tight_layout()
fig_path = os.path.join(outputDir, "cl_comparison.png")
plt.savefig(fig_path, dpi=150)
print(f"Figure saved to {fig_path}")

# Summary
print(f"\nBDF2    wall={bdf2['wall_time']:.1f}s")
print(f"TSTHETA wall={theta['wall_time']:.1f}s  reason={theta.get('reason','?')}")
for label, d in [("BDF2", bdf2), ("TSTHETA", theta)]:
    print(f"  {label:8s} final CL={d['cl'][-1]:.10f}  CD={d['cd'][-1]:.10f}")
