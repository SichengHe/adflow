"""
Minimal test for setupTSPreconditioner crash.
Tests the Fortran PC assembly in isolation.

Usage:
    python debug_pc_crash.py          # np=1
    mpirun -np 2 python debug_pc_crash.py  # np=2
"""
import os, sys
import numpy as np
from mpi4py import MPI
from baseclasses import AeroProblem
from adflow import ADFLOW

sys.stdout.reconfigure(line_buffering=True)
comm = MPI.COMM_WORLD
rank = comm.rank

base_dir = os.path.dirname(os.path.abspath(__file__))
grid_file = os.path.join(base_dir, "../../../input_files/naca0012_rans-L2.cgns")

ap = AeroProblem(
    name="test", alpha=2.77,
    mach=0.6, machRef=0.6, reynolds=4800000.0,
    reynoldsLength=1.0, T=280.0, R=287.085,
    areaRef=1.0, chordRef=1.0,
    evalFuncs=["cl"],
)

options = {
    "gridfile": grid_file,
    "outputdirectory": os.path.join(base_dir, "output"),
    "writevolumesolution": False,
    "writesurfacesolution": False,
    "equationtype": "RANS",
    "equationmode": "unsteady",
    "smoother": "DADI",
    "usenksolver": False, "useanksolver": False,
    "printAllOptions": False, "printIterations": False,
}

solver = ADFLOW(options=options, debug=False)
solver.setAeroProblem(ap)

adflow = solver.adflow
if rank == 0:
    print(f"nCellsLocal(1) = {adflow.adjointvars.ncellslocal[0]}")
    print(f"nw = {adflow.flowvarrefstate.nw}")
    print(f"nTimeIntervalsSpectral = {adflow.inputtimespectral.ntimeintervalsspectral}")
    print(f"equationMode = {adflow.inputphysics.equationmode}")

# Test 1: call setupTSPreconditioner right after init (no forward solve)
if rank == 0:
    print("\n--- Test 1: setupTSPreconditioner right after init ---", flush=True)
try:
    orig_mode = adflow.inputphysics.equationmode
    adflow.inputphysics.equationmode = 1  # STEADY
    adflow.nksolver.setuptspreconditioner(80.0)
    adflow.inputphysics.equationmode = orig_mode
    if rank == 0:
        print("  OK!")
except Exception as e:
    if rank == 0:
        print(f"  FAILED: {e}")

# Test 2: call applyTSPreconditioner
if rank == 0:
    print("\n--- Test 2: applyTSPreconditioner ---", flush=True)
try:
    n = solver.getStateSize()
    in_vec = np.random.randn(n)
    out_vec = np.zeros(n)
    adflow.nksolver.applytspreconditioner(in_vec, out_vec, n)
    if rank == 0:
        print(f"  OK! ||out|| = {np.linalg.norm(out_vec):.6e}")
except Exception as e:
    if rank == 0:
        print(f"  FAILED: {e}")

# Test 3: call applyTSPreconditionerTranspose
if rank == 0:
    print("\n--- Test 3: applyTSPreconditionerTranspose ---", flush=True)
try:
    n = solver.getStateSize()
    np.random.seed(42 + rank)
    in_vec = np.random.randn(n)
    out_fwd = np.zeros(n)
    out_tr = np.zeros(n)
    adflow.nksolver.applytspreconditioner(in_vec, out_fwd, n)
    adflow.nksolver.applytspreconditionertranspose(in_vec, out_tr, n)
    if rank == 0:
        print(f"  Forward: ||out|| = {np.linalg.norm(out_fwd):.6e}")
        print(f"  Transp:  ||out|| = {np.linalg.norm(out_tr):.6e}")
        # They should be different but both nonzero
        if np.linalg.norm(out_tr) < 1e-30:
            print("  WARNING: transpose output is zero!")
        else:
            print("  OK!")
except Exception as e:
    if rank == 0:
        import traceback
        print(f"  FAILED: {e}")
        traceback.print_exc()

if rank == 0:
    print("\nDone.")
