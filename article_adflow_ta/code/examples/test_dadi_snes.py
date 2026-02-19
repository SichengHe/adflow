"""
PETSc TSTHETA with DADI SNES on pitching NACA 0012.

Forward solve via ADflow's native DADI+MG while PETSc TS
manages the trajectory for future adjoint use.

Usage:
    mpirun -np 2 python test_dadi_snes.py
"""

import os, sys, time
import numpy as np
from mpi4py import MPI
from baseclasses import AeroProblem
from adflow import ADFLOW
from adflow.pyADflow_TA import ADflowTS

sys.stdout.reconfigure(line_buffering=True)
comm = MPI.COMM_WORLD
rank = comm.rank

# --- Problem parameters (pitching NACA 0012, same as benchmark_theta_bdf2.py) ---
freq = 10.0
period = 1.0 / freq
n_steps_per_period = 8
n_periods = 3
n_steps = n_steps_per_period * n_periods
dt = period / n_steps_per_period
t_final = period * n_periods

k = 0.0808
mach = 0.6
gamma = 1.4
R_gas = 287.085
T_inf = 280.0
chord = 1.0
alpha_mean = 2.77
alpha_amp = 2.34
omega = 2 * mach * np.sqrt(gamma * R_gas * T_inf) * k / chord
delta_alpha = -alpha_amp * np.pi / 180.0

base_dir = os.path.dirname(os.path.abspath(__file__))
grid_file = os.path.join(base_dir, "../../../input_files/naca0012_rans-L2.cgns")
output_dir = os.path.join(base_dir, "output")

if not os.path.exists(grid_file):
    if rank == 0:
        print(f"ERROR: Mesh not found: {grid_file}")
    sys.exit(1)
os.makedirs(output_dir, exist_ok=True)


def create_aeroproblem():
    return AeroProblem(
        name="0012pitching", alpha=alpha_mean,
        mach=mach, machRef=mach, reynolds=4800000.0,
        reynoldsLength=chord, T=T_inf, R=R_gas,
        areaRef=1.0, chordRef=chord,
        evalFuncs=["cl", "cd", "cmz"],
        xRef=0.25, xRot=0.25,
        degreePol=0, coefPol=[0.0],
        degreeFourier=1, omegaFourier=omega,
        cosCoefFourier=[0.0, 0.0], sinCoefFourier=[delta_alpha],
    )


def main():
    if rank == 0:
        print(f"dt = {dt:.4e}, n_steps = {n_steps}, t_final = {t_final:.4f}")

    options = {
        "gridfile": grid_file,
        "outputdirectory": output_dir,
        "writevolumesolution": False,
        "writesurfacesolution": False,
        "vis4": 0.025, "vis2": 0.5,
        "restrictionrelaxation": 0.5,
        "smoother": "DADI",
        "equationtype": "RANS",
        "equationmode": "unsteady",
        "timeIntegrationscheme": "BDF",
        "ntimestepsfine": n_steps, "deltat": dt,
        "timeaccuracy": 1,
        "nsubiterturb": 10, "nsubiter": 5,
        "useale": False, "usegridmotion": True,
        "cfl": 2.5, "cflcoarse": 1.2, "ncycles": 2000,
        "mgcycle": "3w", "mgstartlevel": 1,
        "monitorvariables": ["cpu", "resrho", "cl", "cd", "cmz"],
        "usenksolver": False, "useanksolver": False,
        "l2convergence": 1e-6, "l2convergencecoarse": 1e-4,
        "qmode": True, "alphafollowing": False,
        "blockSplitting": True, "useblockettes": False,
        "printAllOptions": False, "printIterations": True,
    }

    solver = ADFLOW(options=options, debug=False)
    solver.adflow.solvers.solverunsteadyinit()

    ts = ADflowTS(
        solver, create_aeroproblem(),
        dt=dt, t_final=t_final, theta=1.0,
        grid_motion=True, snes_type="dadi",
        save_trajectory=False,
    )
    ts.setup()

    t0 = time.time()
    reason = ts.solve()
    wall = time.time() - t0

    if rank == 0:
        h = ts.get_history()
        print(f"\nDone in {wall:.1f}s  (reason={int(reason)})")
        print(f"{'Step':>4s} {'Time':>11s} {'CL':>13s} {'CD':>13s} {'CMz':>13s}")
        for i in range(len(h["time"])):
            print(f"{i:4d} {h['time'][i]:11.5e} "
                  f"{h['cl'][i]:13.7e} {h['cd'][i]:13.7e} {h['cmz'][i]:13.7e}")


if __name__ == "__main__":
    main()
