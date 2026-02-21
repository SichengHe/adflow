"""
Run PETSc TSTHETA on pitching NACA 0012 and save CL history.

Same test case as tests/reg_tests/test_time_accurate_naca0012.py.

Usage:
    mpirun -np 2 python run_tstheta.py
"""

import os
import sys
import time
import pickle

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

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

f = 10.0
period = 1.0 / f
nStepPerPeriod = 8
nPeriods = 1
nfineSteps = nStepPerPeriod * nPeriods
dt = period / nStepPerPeriod
t_final = period * nPeriods

baseDir = os.path.dirname(os.path.abspath(__file__))
repoDir = os.path.join(baseDir, "../../..")
gridFile = os.path.join(repoDir, "input_files/naca0012_rans-L2.cgns")
outputDir = os.path.join(baseDir, "output")

if not os.path.exists(gridFile):
    if rank == 0:
        print(f"ERROR: {gridFile} not found")
    sys.exit(1)
os.makedirs(outputDir, exist_ok=True)

from baseclasses import AeroProblem
from adflow import ADFLOW
from adflow.pyADflow_TA import ADflowTS

ap = AeroProblem(
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

options = {
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

if rank == 0:
    print(f"Pitching NACA 0012 — TSTHETA (theta=1.0)")
    print(f"  dt={dt:.4e}, nSteps={nfineSteps}, t_final={t_final:.4f}")
    print("=" * 60)

# --- Initialize via native BDF to get away from freestream ---
# ADflow starts from freestream; PETSc Newton can't converge from such
# a far initial guess.  Run one physical time step using ADflow's native
# DADI+MG solver first, then switch to PETSc TS.
solver = ADFLOW(options=options, debug=False)
solver.setAeroProblem(ap)
solver.adflow.solvers.solverunsteadyinit()
nWarmup = 1
for step in range(nWarmup):
    curTime, _ = solver.advanceTimeStepCounter()
    solver.adflow.preprocessingapi.shiftcoorandvolumes()
    solver.adflow.solvers.updateunsteadygeometry()
    solver.solveTimeStep()
    funcs = {}
    solver.evalFunctions(ap, funcs, evalFuncs=["cl", "cd", "cmz"])
    if rank == 0:
        print(f"  Warmup BDF step {step+1}: t={curTime:.6e}  "
              f"CL={funcs.get(ap.name+'_cl',0):.6e}")

if rank == 0:
    print(f"Warmup done ({nWarmup} BDF steps).  Starting PETSc TS at t={curTime:.6e}.")

# --- Set up PETSc TS from the post-warmup state ---
# Save the post-warmup state before setup() (which calls setAeroProblem
# and would reset back to freestream).
t_start_ts = curTime
t_remaining = t_final - t_start_ts
w_warm = solver.getStates().copy()

ts_wrapper = ADflowTS(
    solver, ap, dt=dt, t_final=t_final,
    theta=1.0, grid_motion=True,
    snes_rtol=1e-4, snes_max_it=2000,
    ksp_type="preonly",  # just apply PC once per SNES iteration (no Krylov)
    jac_type="ad",
    pc_type="asm_ilu",
    pc_shift_factor=1.0,  # Fortran uses CFL-based per-cell diagonal
    save_trajectory=True,  # required for TS adjoint
)
# skip_set_ap=True: solver already has the correct AP and state
# from the warmup step — don't let setup() call setAeroProblem
ts_wrapper.setup(skip_set_ap=True)
ts_wrapper.ts.setTime(t_start_ts)
ts_wrapper.ts.setMaxSteps(int(round(t_remaining / dt)))

# CRITICAL: Set _grid_time so the first IFunction call doesn't spuriously
# call shiftcoorandvolumes (the warmup already advanced the mesh to t_start_ts).
ts_wrapper._grid_time = t_start_ts

# SNES/KSP monitors
snes = ts_wrapper.ts.getSNES()
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
reason = ts_wrapper.solve()
wall = time.time() - t_start

# Example unsteady adjoint: final-time dCL/dw(tf) backward to lambda(t0)
adj = ts_wrapper.solve_adjoint(["cl"])

if rank == 0:
    print(f"\nTSTHETA done in {wall:.1f} s  (reason={reason})")
    print(f"Adjoint complete. ||lambda_cl(t0)||_2 = {np.linalg.norm(adj['cl']):.6e}")
    hist = ts_wrapper.get_history()
    hist["wall_time"] = wall
    hist["reason"] = int(reason)
    hist["lambda0_cl"] = adj["cl"]
    pkl = os.path.join(outputDir, "tstheta_results.pkl")
    with open(pkl, "wb") as fh:
        pickle.dump(hist, fh)
    print(f"Saved to {pkl}")
