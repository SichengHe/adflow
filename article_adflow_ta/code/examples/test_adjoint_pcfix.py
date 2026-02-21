"""
Test adjoint at physical dt with different PC strategies.

The core issue: setupTSPreconditioner uses max(shift, CFL_shift) for the
diagonal, but at physical dt, shift=80 while CFL_shift~1e5. This makes
the PC diagonal ~1000x too large, so P^{-1}*J is far from identity.

Strategy A: Brute-force GMRES (many iterations, large restart)
Strategy B: Use pc_shift_factor to scale up shift, reducing the mismatch

Usage:
    mpirun -np 2 python test_adjoint_pcfix.py
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

# Physical dt (same as pitching NACA test)
freq = 10.0
period = 1.0 / freq
dt = period / 8  # = 0.0125 → a = 80

n_steps = 2
t_final = n_steps * dt

base_dir = os.path.dirname(os.path.abspath(__file__))
grid_file = os.path.join(base_dir, "../../../input_files/naca0012_rans-L2.cgns")
output_dir = os.path.join(base_dir, "output")
os.makedirs(output_dir, exist_ok=True)

# Static RANS (no grid motion for simplicity)
ap = AeroProblem(
    name="pcfix", alpha=2.77,
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
    "ntimestepsfine": n_steps, "deltat": dt,
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


def run_test(label, ksp_max_it, ksp_restart, pc_type, pc_shift_factor=1.0,
             ksp_rtol=1e-4):
    if rank == 0:
        print(f"\n{'='*70}")
        print(f"  {label}")
        print(f"  dt={dt:.4f}, a=1/dt={1/dt:.0f}, pc_type={pc_type}, "
              f"factor={pc_shift_factor}, max_it={ksp_max_it}, "
              f"restart={ksp_restart}, rtol={ksp_rtol:.0e}")
        print(f"{'='*70}", flush=True)

    ts = ADflowTS(
        solver, ap,
        dt=dt, t_final=t_final, theta=1.0,
        grid_motion=False, snes_type="dadi",
        save_trajectory=True,
        ksp_type="gmres", ksp_rtol=ksp_rtol, ksp_max_it=ksp_max_it,
        ksp_gmres_restart=ksp_restart,
        pc_type=pc_type, pc_shift_factor=pc_shift_factor,
    )
    ts.setup()

    # KSP monitor: first 10, then every 100th
    from petsc4py import PETSc as _PETSc
    snes = ts.ts.getSNES()
    ksp = snes.getKSP()
    ksp_data = []

    def _ksp_monitor(ksp, its, rnorm, _data=ksp_data):
        _data.append((its, rnorm))
        if comm.rank == 0 and (its <= 10 or its % 100 == 0):
            print(f"    KSP it {its:5d}  rnorm = {rnorm:.6e}", flush=True)
    ksp.setMonitor(_ksp_monitor)

    # Forward
    if rank == 0:
        print("--- Forward solve ---", flush=True)
    reason = ts.solve()

    # Adjoint
    if rank == 0:
        print("--- Adjoint solve ---", flush=True)
    t0 = time.time()
    try:
        result = ts.solve_adjoint(["cl"])
        wall = time.time() - t0
        lam_norm = np.linalg.norm(result["cl"])
        ok = True
    except Exception as e:
        wall = time.time() - t0
        lam_norm = float("nan")
        ok = False
        if rank == 0:
            print(f"  FAILED: {e}")

    if ksp_data:
        total_ksp = ksp_data[-1][0] + 1
        final_rnorm = ksp_data[-1][1]
    else:
        total_ksp = 0
        final_rnorm = float("nan")

    if rank == 0:
        status = "PASS" if ok else "FAIL"
        print(f"  -> {status} | wall={wall:.1f}s | "
              f"total_ksp_its={total_ksp} | final_rnorm={final_rnorm:.2e} | "
              f"||lambda||={lam_norm:.4e}")
    return ok, total_ksp, final_rnorm, lam_norm


results = []

# Test A: Brute-force GMRES (no PC, many iterations)
ok, n, rn, lam = run_test(
    "A: No PC, 2000 GMRES iters, restart=500",
    ksp_max_it=2000, ksp_restart=500, pc_type="none", ksp_rtol=1e-2,
)
results.append(("A: no PC, 2000 its", ok, n, rn, lam))

# Test B: ASM+ILU with many more GMRES iters
ok, n, rn, lam = run_test(
    "B: ASM+ILU PC, 2000 GMRES iters, restart=500",
    ksp_max_it=2000, ksp_restart=500, pc_type="asm_ilu", ksp_rtol=1e-2,
)
results.append(("B: asm_ilu, 2000 its", ok, n, rn, lam))

# Test C: ASM+ILU with large pc_shift_factor to match CFL scale
# CFL_shift ~ 1e5, a=80, so factor ~ 1e5/80 = 1250
ok, n, rn, lam = run_test(
    "C: ASM+ILU PC, shift_factor=1250, 200 GMRES iters",
    ksp_max_it=200, ksp_restart=200, pc_type="asm_ilu",
    pc_shift_factor=1250.0, ksp_rtol=1e-2,
)
results.append(("C: factor=1250, 200 its", ok, n, rn, lam))

# Test D: scale PC (y = x/a) with many iters
ok, n, rn, lam = run_test(
    "D: Scale PC, 2000 GMRES iters, restart=500",
    ksp_max_it=2000, ksp_restart=500, pc_type="scale", ksp_rtol=1e-2,
)
results.append(("D: scale, 2000 its", ok, n, rn, lam))

if rank == 0:
    print(f"\n{'='*70}")
    print(f"  SUMMARY  (dt={dt}, a=1/dt={1/dt:.0f})")
    print(f"{'='*70}")
    for label, ok, n, rn, lam in results:
        s = "PASS" if ok else "FAIL"
        print(f"  {s} | {label:35s} | KSP={n:5d} | rnorm={rn:10.2e} | ||lam||={lam:.4e}")
    print(f"{'='*70}")
