"""
Run native ADflow BDF2 on pitching NACA 0012 and save CL history.

Same test case as tests/reg_tests/test_time_accurate_naca0012.py.

Usage:
    mpirun -np 2 python run_bdf2.py
"""

import os
import sys
import time
import pickle

import numpy as np
from mpi4py import MPI

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
    "coupledsolution": True,
}

if rank == 0:
    print(f"Pitching NACA 0012 — BDF2")
    print(f"  dt={dt:.4e}, nSteps={nfineSteps}, t_final={t_final:.4f}")
    print("=" * 60)

solver = ADFLOW(options=options, debug=False)
solver(ap)  # setup + solverunsteadyinit, returns immediately

# Step through time, recording CL at each step
times = [0.0]
cl_hist = [0.0]
cd_hist = [0.0]
cmz_hist = [0.0]

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
    print(f"\nBDF2 done in {wall:.1f} s")
    result = {
        "time": np.array(times), "cl": np.array(cl_hist),
        "cd": np.array(cd_hist), "cmz": np.array(cmz_hist),
        "wall_time": wall,
    }
    pkl = os.path.join(outputDir, "bdf2_results.pkl")
    with open(pkl, "wb") as fh:
        pickle.dump(result, fh)
    print(f"Saved to {pkl}")
