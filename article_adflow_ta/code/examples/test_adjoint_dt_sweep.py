"""
Sweep dt to find where the adjoint KSP transitions from converging to diverging.

This helps identify the relationship between physical shift a=1/dt,
CFL-based shift in the PC, and adjoint convergence.

Usage:
    mpirun -np 2 python test_adjoint_dt_sweep.py
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

base_dir = os.path.dirname(os.path.abspath(__file__))
grid_file = os.path.join(base_dir, "../../../input_files/naca0012_rans-L2.cgns")
output_dir = os.path.join(base_dir, "output")
os.makedirs(output_dir, exist_ok=True)

# Static RANS case (no motion)
ap = AeroProblem(
    name="dtsweep", alpha=2.77,
    mach=0.6, machRef=0.6, reynolds=4800000.0,
    reynoldsLength=1.0, T=280.0, R=287.085,
    areaRef=1.0, chordRef=1.0,
    evalFuncs=["cl", "cd"],
)

options = {
    "gridfile": grid_file,
    "outputdirectory": output_dir,
    "writevolumesolution": False,
    "writesurfacesolution": False,
    "smoother": "DADI",
    "equationtype": "RANS",
    "equationmode": "unsteady",
    "timeIntegrationscheme": "BDF",
    "ntimestepsfine": 2, "deltat": 1e-8,  # Will be overridden
    "timeaccuracy": 1,
    "nsubiterturb": 10, "nsubiter": 5,
    "cfl": 2.5, "cflcoarse": 1.2, "ncycles": 2000,
    "mgcycle": "3w", "mgstartlevel": 1,
    "monitorvariables": ["cpu", "resrho", "cl"],
    "usenksolver": False, "useanksolver": False,
    "l2convergence": 1e-6, "l2convergencecoarse": 1e-4,
    "blockSplitting": True, "useblockettes": False,
    "printAllOptions": False, "printIterations": False,
}

solver = ADFLOW(options=options, debug=False)

# Sweep dt from very small (easy) to physical dt (hard)
# Physical dt = 0.0125 → a = 80
# Spectral radius of dR/dw ~ 3.5e6
dt_values = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.0125]

results = []

for dt in dt_values:
    n_steps = 2
    t_final = n_steps * dt
    a = 1.0 / dt

    if rank == 0:
        print(f"\n{'='*70}")
        print(f"  dt = {dt:.1e}  |  a = 1/dt = {a:.1e}  |  n_steps = {n_steps}")
        print(f"{'='*70}", flush=True)

    ts = ADflowTS(
        solver, ap,
        dt=dt, t_final=t_final, theta=1.0,
        grid_motion=False, snes_type="dadi",
        save_trajectory=True,
        ksp_type="gmres", ksp_rtol=1e-4, ksp_max_it=100,
        pc_type="asm_ilu",
    )
    ts.setup()

    # Sparse KSP monitor
    from petsc4py import PETSc as _PETSc
    snes = ts.ts.getSNES()
    ksp = snes.getKSP()
    ksp_iters = []

    def _ksp_monitor(ksp, its, rnorm, _dt=dt, _iters=ksp_iters):
        _iters.append((its, rnorm))
        if comm.rank == 0 and (its <= 5 or its % 50 == 0):
            print(f"    KSP it {its:4d}  rnorm = {rnorm:.6e}", flush=True)
    ksp.setMonitor(_ksp_monitor)

    # Forward
    reason = ts.solve()

    # Adjoint
    adj_ok = False
    try:
        result = ts.solve_adjoint(["cl"])
        adj_ok = True
        lam_norm = np.linalg.norm(result["cl"])
    except Exception:
        lam_norm = float("nan")

    # Report
    if ksp_iters:
        final_rnorm = ksp_iters[-1][1]
        n_ksp = ksp_iters[-1][0]
    else:
        final_rnorm = float("nan")
        n_ksp = 0

    results.append((dt, a, adj_ok, n_ksp, final_rnorm, lam_norm))
    if rank == 0:
        status = "PASS" if adj_ok else "FAIL"
        print(f"  -> {status}: KSP iters={n_ksp}, final_rnorm={final_rnorm:.2e}, "
              f"||lambda||={lam_norm:.4e}")

# Summary table
if rank == 0:
    print(f"\n{'='*70}")
    print(f"  SUMMARY")
    print(f"{'='*70}")
    print(f"  {'dt':>10s} | {'a=1/dt':>10s} | {'Status':>6s} | {'KSP its':>8s} | {'rnorm':>10s} | {'||lambda||':>12s}")
    print(f"  {'-'*10}-+-{'-'*10}-+-{'-'*6}-+-{'-'*8}-+-{'-'*10}-+-{'-'*12}")
    for dt, a, ok, n_ksp, rnorm, lam_norm in results:
        status = "PASS" if ok else "FAIL"
        print(f"  {dt:10.1e} | {a:10.1e} | {status:>6s} | {n_ksp:8d} | {rnorm:10.2e} | {lam_norm:12.4e}")
    print(f"{'='*70}")
