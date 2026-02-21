"""
Well-conditioned adjoint test: verify transpose infrastructure.

Uses a very small dt so that shift a = 1/dt >> spectral_radius(dR/dw).
This makes J = a*I - dR/dw trivially well-conditioned, so even a scale
PC (y = x/a) works perfectly. If the adjoint KSP still diverges here,
the issue is in the transpose dispatch, not PC quality.

Tests:
1. scale PC (y = x/a) — should converge in ~1 iteration
2. asm_ilu PC — should also converge easily
3. no PC — should still converge (well-conditioned)

Usage:
    mpirun -np 2 python test_adjoint_wellcond.py
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

# Very small dt → large shift a = 1e8 >> spectral radius ~3.5M
dt = 1e-8
n_steps = 2
t_final = n_steps * dt

base_dir = os.path.dirname(os.path.abspath(__file__))
grid_file = os.path.join(base_dir, "../../../input_files/naca0012_rans-L2.cgns")
output_dir = os.path.join(base_dir, "output")
os.makedirs(output_dir, exist_ok=True)

# Static (no motion) RANS case
ap = AeroProblem(
    name="wellcond", alpha=2.77,
    mach=0.6, machRef=0.6, reynolds=4800000.0,
    reynoldsLength=1.0, T=280.0, R=287.085,
    areaRef=1.0, chordRef=1.0,
    evalFuncs=["cl", "cd", "cmz"],
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
    "monitorvariables": ["cpu", "resrho", "cl", "cd"],
    "usenksolver": False, "useanksolver": False,
    "l2convergence": 1e-6, "l2convergencecoarse": 1e-4,
    "blockSplitting": True, "useblockettes": False,
    "printAllOptions": False, "printIterations": True,
}

solver = ADFLOW(options=options, debug=False)


def run_adjoint_test(pc_type, label):
    """Run forward + adjoint with given PC type."""
    if rank == 0:
        print(f"\n{'='*70}")
        print(f"  TEST: {label}  (pc_type={pc_type}, dt={dt:.1e}, a=1/dt={1/dt:.1e})")
        print(f"{'='*70}", flush=True)

    ts = ADflowTS(
        solver, ap,
        dt=dt, t_final=t_final, theta=1.0,
        grid_motion=False, snes_type="dadi",
        save_trajectory=True,
        ksp_type="gmres", ksp_rtol=1e-4, ksp_max_it=50,
        pc_type=pc_type,
    )
    ts.setup()

    # KSP monitor — print ALL iterations for this small test
    from petsc4py import PETSc as _PETSc
    snes = ts.ts.getSNES()
    ksp = snes.getKSP()

    def _ksp_monitor(ksp, its, rnorm):
        if comm.rank == 0:
            print(f"    KSP it {its:4d}  rnorm = {rnorm:.6e}", flush=True)
    ksp.setMonitor(_ksp_monitor)

    # Forward solve
    if rank == 0:
        print("\n--- Forward solve ---", flush=True)
    t0 = time.time()
    reason = ts.solve()
    wall_fwd = time.time() - t0
    if rank == 0:
        print(f"Forward done in {wall_fwd:.1f}s (reason={int(reason)})")
        h = ts.get_history()
        for i in range(len(h["time"])):
            print(f"  step {i:2d} | t={h['time'][i]:.5e} | CL={h['cl'][i]:.7e}")

    # Adjoint solve
    if rank == 0:
        print("\n--- Adjoint solve (CL at final time) ---", flush=True)
    t0 = time.time()
    try:
        result = ts.solve_adjoint(["cl"])
        wall_adj = time.time() - t0
        if rank == 0:
            print(f"Adjoint done in {wall_adj:.1f}s")
            for obj, lam in result.items():
                lam_norm = np.linalg.norm(lam)
                print(f"  {obj}: ||lambda(t0)||_2 = {lam_norm:.6e}")
            return True
    except Exception as e:
        if rank == 0:
            import traceback
            print(f"Adjoint FAILED: {e}")
            traceback.print_exc()
        return False


# Test 1: scale PC (trivial: y = x/a)
ok1 = run_adjoint_test("scale", "Scale PC (y = x/a)")

# Test 2: no PC
ok2 = run_adjoint_test("none", "No PC (unpreconditioned GMRES)")

# Test 3: ASM+ILU PC
ok3 = run_adjoint_test("asm_ilu", "ASM+ILU Fortran PC")

if rank == 0:
    print(f"\n{'='*70}")
    print(f"  SUMMARY")
    print(f"{'='*70}")
    print(f"  Scale PC:   {'PASS' if ok1 else 'FAIL'}")
    print(f"  No PC:      {'PASS' if ok2 else 'FAIL'}")
    print(f"  ASM+ILU PC: {'PASS' if ok3 else 'FAIL'}")
    print(f"{'='*70}")
