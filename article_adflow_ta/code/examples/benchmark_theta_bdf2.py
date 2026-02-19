"""
Benchmark: native ADflow BDF2 on pitching NACA 0012.

Runs the unsteady pitching airfoil case with manual time stepping to
record CL, CD, CMz at every physical time step, then generates a CL
history plot.

Usage:
    mpirun -np 2 python benchmark_theta_bdf2.py

Requires:
    - ADflow built and importable
    - Input mesh: input_files/naca0012_rans-L2.cgns
      (download via input_files/get-input-files.sh)
"""

import os
import sys
import time
import pickle

import numpy as np
from mpi4py import MPI
from baseclasses import AeroProblem
from adflow import ADFLOW

sys.stdout.reconfigure(line_buffering=True)

comm = MPI.COMM_WORLD
rank = comm.rank

# --- Problem parameters (same as reg_tests/test_time_accurate_naca0012) ---
freq = 10.0  # [Hz] Forcing frequency
period = 1.0 / freq  # [sec]
n_steps_per_period = 8
n_periods = 3
n_steps = n_steps_per_period * n_periods
dt = period / n_steps_per_period  # [s]
t_final = period * n_periods

# Airfoil parameters
k = 0.0808  # reduced frequency
mach = 0.6
gamma = 1.4
R_gas = 287.085
T_inf = 280.0
chord = 1.0
alpha_mean = 2.77  # degrees
alpha_amp = 2.34  # degrees

omega = 2 * mach * np.sqrt(gamma * R_gas * T_inf) * k / chord
delta_alpha = -alpha_amp * np.pi / 180.0

# Paths
base_dir = os.path.dirname(os.path.abspath(__file__))
repo_dir = os.path.join(base_dir, "../../..")
grid_file = os.path.join(repo_dir, "input_files/naca0012_rans-L2.cgns")
output_dir = os.path.join(base_dir, "output")

if not os.path.exists(grid_file):
    if rank == 0:
        print(f"ERROR: Mesh file not found: {grid_file}")
        print("Run input_files/get-input-files.sh to download test meshes.")
    sys.exit(1)

os.makedirs(output_dir, exist_ok=True)


def create_aeroproblem():
    """Create the NACA 0012 pitching AeroProblem."""
    return AeroProblem(
        name="0012pitching",
        alpha=alpha_mean,
        mach=mach,
        machRef=mach,
        reynolds=4800000.0,
        reynoldsLength=chord,
        T=T_inf,
        R=R_gas,
        areaRef=1.0,
        chordRef=chord,
        evalFuncs=["cl", "cd", "cmz"],
        xRef=0.25,
        xRot=0.25,
        degreePol=0,
        coefPol=[0.0],
        degreeFourier=1,
        omegaFourier=omega,
        cosCoefFourier=[0.0, 0.0],
        sinCoefFourier=[delta_alpha],
    )


def run_bdf2_with_history():
    """Run native ADflow BDF2 with manual stepping, recording CL at each step."""
    if rank == 0:
        print("=" * 70)
        print("Running native ADflow BDF2 — manual stepping")
        print(f"  dt = {dt:.4e}, n_steps = {n_steps}, t_final = {t_final:.4f}")
        print(f"  n_periods = {n_periods}, steps/period = {n_steps_per_period}")
        print("=" * 70)

    options = {
        "gridfile": grid_file,
        "outputdirectory": output_dir,
        "writevolumesolution": False,
        "writesurfacesolution": False,
        "vis4": 0.025,
        "vis2": 0.5,
        "restrictionrelaxation": 0.5,
        "smoother": "DADI",
        "equationtype": "RANS",
        "equationmode": "unsteady",
        "timeIntegrationscheme": "BDF",
        "ntimestepsfine": n_steps,
        "deltat": dt,
        "timeaccuracy": 2,
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

    ap = create_aeroproblem()
    solver = ADFLOW(options=options, debug=False)
    solver.setAeroProblem(ap)
    solver.adflow.solvers.solverunsteadyinit()

    # Record history
    time_hist = [0.0]
    cl_hist = []
    cd_hist = []
    cmz_hist = []

    # Evaluate initial CL
    funcs = {}
    solver.evalFunctions(ap, funcs, evalFuncs=["cl", "cd", "cmz"])
    cl_hist.append(funcs.get(f"{ap.name}_cl", 0.0))
    cd_hist.append(funcs.get(f"{ap.name}_cd", 0.0))
    cmz_hist.append(funcs.get(f"{ap.name}_cmz", 0.0))
    if rank == 0:
        print(f"  Step   0 | t = 0.000000e+00 | "
              f"CL = {cl_hist[0]:.6e} | CD = {cd_hist[0]:.6e} | "
              f"CMz = {cmz_hist[0]:.6e}")

    t_wall_start = time.time()

    for step in range(1, n_steps + 1):
        # Advance time counter
        curTime, _ = solver.advanceTimeStepCounter()

        # Update mesh for prescribed motion
        solver.adflow.preprocessingapi.shiftcoorandvolumes()
        solver.adflow.solvers.updateunsteadygeometry()

        # Converge the implicit system at this time step
        solver.solveTimeStep()

        # Evaluate force coefficients
        funcs = {}
        solver.evalFunctions(ap, funcs, evalFuncs=["cl", "cd", "cmz"])
        cl = funcs.get(f"{ap.name}_cl", 0.0)
        cd = funcs.get(f"{ap.name}_cd", 0.0)
        cmz = funcs.get(f"{ap.name}_cmz", 0.0)

        time_hist.append(curTime)
        cl_hist.append(cl)
        cd_hist.append(cd)
        cmz_hist.append(cmz)

        if rank == 0:
            print(f"  Step {step:3d} | t = {curTime:.6e} | "
                  f"CL = {cl:.6e} | CD = {cd:.6e} | CMz = {cmz:.6e}")

    t_wall = time.time() - t_wall_start

    if rank == 0:
        print(f"\nBDF2 completed in {t_wall:.2f} s")

    history = {
        "time": np.array(time_hist),
        "cl": np.array(cl_hist),
        "cd": np.array(cd_hist),
        "cmz": np.array(cmz_hist),
        "wall_time": t_wall,
        "dt": dt,
        "n_steps": n_steps,
        "n_periods": n_periods,
    }
    return history


def plot_history(history, save_path):
    """Generate CL, CD, CMz history plots."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available — skipping plot")
        return

    t = history["time"]
    cl = history["cl"]
    cd = history["cd"]
    cmz = history["cmz"]

    # Compute alpha(t) for reference
    alpha_t = alpha_mean + alpha_amp * np.sin(omega * t) * (180.0 / np.pi)
    # (Note: the Fourier representation gives alpha in a specific form;
    #  for plotting we just show the sinusoidal variation.)
    alpha_t = alpha_mean + alpha_amp * np.sin(omega * t)

    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

    # CL
    ax = axes[0]
    ax.plot(t, cl, "bo-", markersize=5, label="BDF2")
    ax.set_ylabel("$C_L$")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title(
        f"Pitching NACA 0012 — BDF2, dt={dt:.4e}, {n_steps_per_period} steps/period"
    )

    # CD
    ax = axes[1]
    ax.plot(t, cd, "ro-", markersize=5, label="BDF2")
    ax.set_ylabel("$C_D$")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # CMz
    ax = axes[2]
    ax.plot(t, cmz, "go-", markersize=5, label="BDF2")
    ax.set_ylabel("$C_{Mz}$")
    ax.set_xlabel("Time [s]")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    print(f"Plot saved to {save_path}")
    plt.close()


if __name__ == "__main__":
    history = run_bdf2_with_history()

    if rank == 0:
        # Save raw data
        pkl_path = os.path.join(output_dir, "bdf2_history.pkl")
        with open(pkl_path, "wb") as fh:
            pickle.dump(history, fh)
        print(f"Data saved to {pkl_path}")

        # Plot
        plot_path = os.path.join(output_dir, "cl_history_bdf2.png")
        plot_history(history, plot_path)

        # Print summary table
        print("\n" + "=" * 70)
        print("Time step summary:")
        print(f"{'Step':>5s} {'Time':>12s} {'CL':>14s} {'CD':>14s} {'CMz':>14s}")
        print("-" * 70)
        for i in range(len(history["time"])):
            print(f"{i:5d} {history['time'][i]:12.6e} "
                  f"{history['cl'][i]:14.8e} {history['cd'][i]:14.8e} "
                  f"{history['cmz'][i]:14.8e}")
        print("=" * 70)
