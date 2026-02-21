"""
Compare CL history: ADflow native BDF2 vs PETSc TSTHETA.

Uses the same test case as tests/reg_tests/test_time_accurate_naca0012.py:
  Pitching NACA 0012, RANS, M=0.6, Re=4.8M, k=0.0808
  8 steps/period, 1 period, dt=0.0125 s

Usage:
    mpirun -np 2 python run_pitching_comparison.py
"""

import os
import sys
import time
import copy
import pickle

import numpy as np
from mpi4py import MPI

# Force line-buffered stdout so output appears promptly even in file redirects
sys.stdout.reconfigure(line_buffering=True)

comm = MPI.COMM_WORLD
rank = comm.rank

# ---- Problem parameters (same as reg test) ----
k = 0.0808
M = 0.6
gamma = 1.4
R = 287.085
T = 280.0
c = 1.0
alpha_m = 2.77
alpha_0 = 2.34

omega = 2 * M * np.sqrt(gamma * R * T) * k / c
deltaAlpha = -alpha_0 * np.pi / 180.0

f = 10.0           # [Hz]
period = 1.0 / f   # [s]
nStepPerPeriod = 8
nPeriods = 1
nfineSteps = nStepPerPeriod * nPeriods
dt = period / nStepPerPeriod  # 0.0125 s
t_final = period * nPeriods

# Paths
baseDir = os.path.dirname(os.path.abspath(__file__))
repoDir = os.path.join(baseDir, "../../..")
gridFile = os.path.join(repoDir, "input_files/naca0012_rans-L2.cgns")
outputDir = os.path.join(baseDir, "output")

if not os.path.exists(gridFile):
    if rank == 0:
        print(f"ERROR: {gridFile} not found")
    sys.exit(1)
os.makedirs(outputDir, exist_ok=True)

if rank == 0:
    print(f"Pitching NACA 0012: BDF2 vs TSTHETA")
    print(f"  dt={dt:.4e}, nSteps={nfineSteps}, t_final={t_final:.4f}")
    print(f"  M={M}, Re=4.8M, alpha={alpha_m}+/-{alpha_0} deg, k={k}")
    print()


def create_ap():
    from baseclasses import AeroProblem
    return AeroProblem(
        name="0012pitching",
        alpha=alpha_m, mach=M, machRef=M,
        reynolds=4800000.0, reynoldsLength=c, T=T, R=R,
        areaRef=1.0, chordRef=c,
        evalFuncs=["cl", "cd", "cmz"],
        xRef=0.25, xRot=0.25,
        degreePol=0, coefPol=[0.0],
        degreeFourier=1, omegaFourier=omega,
        cosCoefFourier=[0.0, 0.0], sinCoefFourier=[deltaAlpha],
    )


def get_options():
    """Options matching test_time_accurate_naca0012.py."""
    return {
        "gridfile": gridFile,
        "outputdirectory": outputDir,
        "writevolumesolution": False,
        "writesurfacesolution": False,
        "vis4": 0.025,
        "vis2": 0.5,
        "restrictionrelaxation": 0.5,
        "smoother": "DADI",
        "equationtype": "RANS",
        "equationmode": "unsteady",
        "timeIntegrationscheme": "BDF",
        "ntimestepsfine": nfineSteps,
        "deltat": dt,
        "nsubiterturb": 10,
        "nsubiter": 5,
        "useale": False,
        "usegridmotion": True,
        "cfl": 2.5,
        "cflcoarse": 1.2,
        "ncycles": 2000,
        "mgcycle": "3w",
        "mgstartlevel": 1,
        "monitorvariables": ["cpu", "resrho", "cl", "cd", "cmz"],
        "usenksolver": False,
        "useanksolver": False,
        "l2convergence": 1e-6,
        "l2convergencecoarse": 1e-4,
        "qmode": True,
        "alphafollowing": False,
        "blockSplitting": True,
        "useblockettes": False,
        "printAllOptions": False,
        "printIterations": False,
    }


# ==============================================================
# BDF2
# ==============================================================
def run_bdf2():
    from adflow import ADFLOW

    if rank == 0:
        print("=" * 70)
        print("BDF2: Native ADflow unsteady solver")
        print("=" * 70)

    options = get_options()
    options["coupledsolution"] = True  # return after init

    ap = create_ap()
    solver = ADFLOW(options=options, debug=False)
    solver(ap)  # setup + solverunsteadyinit, returns immediately

    # Initial state (freestream — CL ≈ 0)
    times = [0.0]
    cl_hist = [0.0]
    cd_hist = [0.0]
    cmz_hist = [0.0]
    if rank == 0:
        print(f"  BDF2   0/{nfineSteps} | t=0.000e+00 | CL=0 (freestream)", flush=True)

    t_start = time.time()
    for i in range(nfineSteps):
        curTime, _ = solver.advanceTimeStepCounter()
        solver.adflow.preprocessingapi.shiftcoorandvolumes()
        solver.adflow.solvers.updateunsteadygeometry()
        solver.solveTimeStep()

        funcs = {}
        solver.evalFunctions(ap, funcs, evalFuncs=["cl", "cd", "cmz"])
        times.append(curTime)
        cl_hist.append(funcs["0012pitching_cl"])
        cd_hist.append(funcs["0012pitching_cd"])
        cmz_hist.append(funcs["0012pitching_cmz"])

        if rank == 0:
            print(f"  BDF2 {i+1:3d}/{nfineSteps} | t={curTime:.4e} | "
                  f"CL={cl_hist[-1]:.8f} | CD={cd_hist[-1]:.8f}", flush=True)

    wall = time.time() - t_start
    if rank == 0:
        print(f"\nBDF2 done in {wall:.1f} s\n", flush=True)

    return {"time": np.array(times), "cl": np.array(cl_hist),
            "cd": np.array(cd_hist), "cmz": np.array(cmz_hist),
            "wall_time": wall}


# ==============================================================
# TSTHETA
# ==============================================================
def run_theta():
    from adflow import ADFLOW, ADflowTS

    if rank == 0:
        print("=" * 70)
        print("TSTHETA: PETSc implicit theta (theta=1.0, backward Euler)")
        print("=" * 70)

    options = get_options()

    ap = create_ap()
    solver = ADFLOW(options=options, debug=False)
    solver.adflow.solvers.solverunsteadyinit()

    ts = ADflowTS(
        solver, ap, dt=dt, t_final=t_final,
        theta=1.0, grid_motion=True,
        snes_rtol=1e-8, snes_max_it=200,
        ksp_rtol=1e-4, ksp_max_it=200,
        jac_type="ad",
    )
    ts.setup()
    ts.ts.setMaxSNESFailures(-1)

    # Add SNES/KSP monitors for live progress
    snes = ts.ts.getSNES()
    ksp = snes.getKSP()

    def snes_mon(snes, its, rnorm):
        if rank == 0:
            print(f"    SNES {its:3d}  ||F||={rnorm:.6e}", flush=True)

    def ksp_mon(ksp, its, rnorm):
        if rank == 0 and (its <= 3 or its % 10 == 0):
            print(f"      KSP {its:4d}  ||r||={rnorm:.6e}", flush=True)

    snes.setMonitor(snes_mon)
    ksp.setMonitor(ksp_mon)

    t_start = time.time()
    reason = ts.solve()
    wall = time.time() - t_start

    if rank == 0:
        print(f"\nTSTHETA done in {wall:.1f} s  (reason={reason})\n")

    hist = ts.get_history()
    hist["wall_time"] = wall
    hist["reason"] = int(reason)
    return hist


# ==============================================================
# Plot
# ==============================================================
def plot_comparison(bdf2, theta):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

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
    plt.close()


# ==============================================================
if __name__ == "__main__":
    bdf2 = run_bdf2()
    theta = run_theta()

    if rank == 0:
        pkl = os.path.join(outputDir, "comparison_results.pkl")
        with open(pkl, "wb") as fh:
            pickle.dump({"bdf2": bdf2, "theta": theta,
                         "params": {"dt": dt, "t_final": t_final,
                                    "nSteps": nfineSteps}}, fh)
        print(f"Results saved to {pkl}")

        plot_comparison(bdf2, theta)

        print("\n" + "=" * 70)
        print("Summary")
        print("=" * 70)
        print(f"  BDF2    wall={bdf2['wall_time']:.1f}s")
        print(f"  TSTHETA wall={theta['wall_time']:.1f}s  reason={theta.get('reason','?')}")
        for label, d in [("BDF2", bdf2), ("TSTHETA", theta)]:
            print(f"  {label:8s} final CL={d['cl'][-1]:.10f}  CD={d['cd'][-1]:.10f}")
        print("=" * 70)
