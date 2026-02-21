"""
Simplest ADflow + PETSc TSTHETA test: Euler on static mesh.

Time-marches from freestream toward steady state using PETSc TSTHETA
(backward Euler).  No grid motion, no turbulence, no oscillation.

Goal: verify that PETSc SNES can converge the implicit system at each
step when dt is small enough (large temporal shift dominates the Jacobian).

Usage:
    mpirun -np 1 python test_euler_tstheta.py
"""

import os
import sys
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

sys.stdout.reconfigure(line_buffering=True)

comm = MPI.COMM_WORLD
rank = comm.rank

# ---- Paths ----
base_dir = os.path.dirname(os.path.abspath(__file__))
repo_dir = os.path.join(base_dir, "../..")
grid_file = os.path.join(repo_dir, "input_files/naca64A010_euler-L2.cgns")

if not os.path.exists(grid_file):
    if rank == 0:
        print(f"ERROR: {grid_file} not found")
    sys.exit(1)

from baseclasses import AeroProblem
from adflow import ADFLOW
from adflow.pyADflow_TA import ADflowTS

# ---- AeroProblem: simple transonic Euler ----
ap = AeroProblem(
    name="euler_test",
    alpha=1.0,
    mach=0.8,
    T=300.0,
    P=100000.0,
    areaRef=1.0,
    chordRef=1.0,
    evalFuncs=["cl", "cd"],
)

# ---- Solver options: Euler, steady ----
# We initialize in unsteady mode so the time-stepping arrays are allocated,
# but PETSc TS handles the actual time integration.
dt = 1e-5  # small dt so shift = 1/dt = 1e5 >> spectral radius of dR/dw
n_steps = 5
t_final = dt * n_steps

options = {
    "gridfile": grid_file,
    "outputdirectory": os.path.join(base_dir, "output"),
    "writevolumesolution": False,
    "writesurfacesolution": False,
    "equationtype": "Euler",
    "equationmode": "unsteady",
    "timeIntegrationscheme": "BDF",
    "ntimestepsfine": n_steps,
    "deltat": dt,
    "smoother": "DADI",
    "cfl": 1.0,
    "cflcoarse": 0.5,
    "mgcycle": "sg",
    "ncycles": 500,
    "l2convergence": 1e-10,
    "monitorvariables": ["resrho", "cl", "cd"],
    "usenksolver": False,
    "useanksolver": False,
    "printAllOptions": False,
    "printIterations": False,
    "useblockettes": False,
}

if rank == 0:
    print("=" * 70)
    print("Simple Euler + PETSc TSTHETA test")
    print(f"  dt = {dt:.2e}, n_steps = {n_steps}, t_final = {t_final:.2e}")
    print(f"  shift = 1/(theta*dt) = {1.0/dt:.2e}")
    print("=" * 70)

# ---- Initialize ADflow ----
solver = ADFLOW(options=options, debug=False)
solver.setAeroProblem(ap)

n_local = solver.getStateSize()
if rank == 0:
    print(f"  State size (local): {n_local}")

# Initialize unsteady arrays (needed for temporal source terms)
solver.adflow.solvers.solverunsteadyinit()

# ---- Check spatial residual at freestream ----
w0 = solver.getStates().copy()
if rank == 0:
    print(f"  ||w0|| = {np.linalg.norm(w0):.6e}")

# Evaluate spatial residual in steady mode
adflow = solver.adflow
orig_mode = adflow.inputphysics.equationmode
adflow.inputphysics.equationmode = 1  # STEADY
res = np.zeros(n_local)
res = adflow.nksolver.getres(res)
adflow.inputphysics.equationmode = orig_mode
if rank == 0:
    print(f"  ||R(w0)|| at freestream = {np.linalg.norm(res):.6e}")

# ---- Evaluate CL at freestream ----
funcs0 = {}
solver.evalFunctions(ap, funcs0, evalFuncs=["cl", "cd"])
if rank == 0:
    print(f"  CL(w0) = {funcs0.get('euler_test_cl', float('nan')):.6e}")
    print(f"  CD(w0) = {funcs0.get('euler_test_cd', float('nan')):.6e}")

# ---- Set up PETSc TSTHETA ----
ts_wrapper = ADflowTS(
    solver, ap, dt=dt, t_final=t_final,
    theta=1.0,  # backward Euler
    grid_motion=False,
    snes_rtol=1e-6,
    snes_max_it=50,
    ksp_type="gmres",
    ksp_rtol=1e-4,
    ksp_max_it=200,
    jac_type="ad",
    pc_type="asm_ilu",
    pc_shift_factor=1.0,
    save_trajectory=True,
)
ts_wrapper.setup(skip_set_ap=True)

# SNES/KSP monitors
snes = ts_wrapper.ts.getSNES()
ksp = snes.getKSP()

def snes_mon(snes, its, rnorm):
    if rank == 0:
        print(f"    SNES {its:3d}  ||F||={rnorm:.6e}", flush=True)

def ksp_mon(ksp, its, rnorm):
    if rank == 0 and (its <= 3 or its % 20 == 0):
        print(f"      KSP {its:4d}  ||r||={rnorm:.6e}", flush=True)

snes.setMonitor(snes_mon)
ksp.setMonitor(ksp_mon)

# ---- Forward solve ----
if rank == 0:
    print("\nStarting PETSc TSTHETA forward solve...")
reason = ts_wrapper.solve()

if rank == 0:
    print(f"\nForward solve reason: {reason}")
    hist = ts_wrapper.get_history()
    print(f"  Time history: {hist['time']}")
    print(f"  CL history:   {hist['cl']}")
    print(f"  CD history:   {hist['cd']}")

# ---- Check final state ----
w_final = solver.getStates().copy()
delta_w = w_final - w0
if rank == 0:
    print(f"\n  ||w_final - w0|| = {np.linalg.norm(delta_w):.6e}")
    print(f"  ||w_final - w0|| / ||w0|| = {np.linalg.norm(delta_w)/np.linalg.norm(w0):.6e}")

funcs_final = {}
solver.evalFunctions(ap, funcs_final, evalFuncs=["cl", "cd"])
if rank == 0:
    cl_final = funcs_final.get('euler_test_cl', float('nan'))
    print(f"  CL(w_final) = {cl_final:.6e}")
    print(f"  CD(w_final) = {funcs_final.get('euler_test_cd', float('nan')):.6e}")

if rank == 0:
    print("\nDone.")
