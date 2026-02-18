"""
Benchmark: PETSc TSTHETA (Crank--Nicolson) vs ADflow native BDF2.

Runs the NACA 0012 pitching airfoil test case with both time integration
methods and compares the force coefficient time histories.

Usage:
    mpirun -np 2 python benchmark_theta_bdf2.py

Requires:
    - ADflow built and importable
    - petsc4py installed
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
from adflow import ADFLOW, ADflowTS

comm = MPI.COMM_WORLD
rank = comm.rank

# --- Problem parameters ---
freq = 10.0  # [Hz] Forcing frequency
period = 1.0 / freq  # [sec]
n_steps_per_period = 8
n_periods = 1
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
repo_dir = os.path.join(base_dir, "../..")
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


def get_common_options():
    """Return solver options common to both methods."""
    return {
        "gridfile": grid_file,
        "outputdirectory": output_dir,
        "writevolumesolution": False,
        "writesurfacesolution": False,
        "vis4": 0.025,
        "vis2": 0.5,
        "restrictionrelaxation": 0.5,
        "smoother": "DADI",
        "equationtype": "RANS",
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
    }


def run_bdf2():
    """Run the native ADflow BDF2 unsteady solver."""
    if rank == 0:
        print("\n" + "=" * 70)
        print("Running native ADflow BDF2")
        print("=" * 70)

    options = get_common_options()
    options.update({
        "equationmode": "unsteady",
        "timeIntegrationscheme": "BDF",
        "ntimestepsfine": n_steps,
        "deltat": dt,
        "timeaccuracy": 2,
    })

    ap = create_aeroproblem()
    solver = ADFLOW(options=options, debug=False)

    t_start = time.time()
    solver(ap)
    t_elapsed = time.time() - t_start

    # Evaluate final functions
    funcs = {}
    solver.evalFunctions(ap, funcs)

    if rank == 0:
        print(f"\nBDF2 completed in {t_elapsed:.2f} s")
        for key, val in funcs.items():
            print(f"  {key} = {val:.10e}")

    return funcs, t_elapsed


def run_theta():
    """Run PETSc TSTHETA (Crank--Nicolson) wrapping ADflow."""
    if rank == 0:
        print("\n" + "=" * 70)
        print("Running PETSc TSTHETA (Crank--Nicolson, theta=0.5)")
        print("=" * 70)

    # Set up ADflow in steady mode for spatial residual evaluation.
    # Grid motion is handled by ADflowTS internally.
    options = get_common_options()
    options.update({
        "equationmode": "unsteady",
        "timeIntegrationscheme": "BDF",
        "ntimestepsfine": n_steps,
        "deltat": dt,
        "timeaccuracy": 2,
    })

    ap = create_aeroproblem()
    solver = ADFLOW(options=options, debug=False)

    # Initialize unsteady arrays
    solver.adflow.solvers.solverunsteadyinit()

    # Create and set up the PETSc TS wrapper
    ts_wrapper = ADflowTS(
        solver, ap, dt=dt, t_final=t_final,
        theta=0.5, grid_motion=True,
    )
    ts_wrapper.setup()

    t_start = time.time()
    reason = ts_wrapper.solve()
    t_elapsed = time.time() - t_start

    # Evaluate final functions
    funcs = {}
    solver.evalFunctions(ap, funcs)

    if rank == 0:
        print(f"\nTSTHETA completed in {t_elapsed:.2f} s")
        print(f"  Converged reason: {reason}")
        for key, val in funcs.items():
            print(f"  {key} = {val:.10e}")

    history = ts_wrapper.get_history()
    return funcs, t_elapsed, history


def compare_results(funcs_bdf2, funcs_theta):
    """Compare force coefficients between the two methods."""
    if rank != 0:
        return

    print("\n" + "=" * 70)
    print("Comparison: BDF2 vs TSTHETA (Crank--Nicolson)")
    print("=" * 70)
    print(f"{'Function':<25s} {'BDF2':>15s} {'TSTHETA':>15s} {'Rel. Diff':>12s}")
    print("-" * 70)

    for key in sorted(funcs_bdf2.keys()):
        val_bdf2 = funcs_bdf2[key]
        val_theta = funcs_theta.get(key, float("nan"))
        if abs(val_bdf2) > 1e-15:
            rel_diff = abs(val_theta - val_bdf2) / abs(val_bdf2)
        else:
            rel_diff = abs(val_theta - val_bdf2)
        print(f"  {key:<23s} {val_bdf2:>15.8e} {val_theta:>15.8e} {rel_diff:>12.4e}")

    print("=" * 70)


if __name__ == "__main__":
    # Run BDF2
    funcs_bdf2, time_bdf2 = run_bdf2()

    # Run TSTHETA
    funcs_theta, time_theta, history_theta = run_theta()

    # Compare
    compare_results(funcs_bdf2, funcs_theta)

    # Save results
    if rank == 0:
        results = {
            "bdf2": {"funcs": funcs_bdf2, "time": time_bdf2},
            "theta": {
                "funcs": funcs_theta,
                "time": time_theta,
                "history": history_theta,
            },
            "params": {
                "dt": dt,
                "t_final": t_final,
                "n_steps": n_steps,
                "theta": 0.5,
            },
        }
        pkl_path = os.path.join(output_dir, "benchmark_results.pkl")
        with open(pkl_path, "wb") as fh:
            pickle.dump(results, fh)
        print(f"\nResults saved to {pkl_path}")
