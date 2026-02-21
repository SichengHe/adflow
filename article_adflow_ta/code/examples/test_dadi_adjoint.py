"""
Test TSAdjointSolve with DADI-based forward trajectory.

Runs a short forward solve (3 time steps) with save_trajectory=True,
then runs the discrete adjoint for CL at the final time.

Usage:
    mpirun -np 2 python test_dadi_adjoint.py
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

# Short test: 3 steps only
freq = 10.0
period = 1.0 / freq
dt = period / 8  # = 0.0125
n_steps = 3
t_final = n_steps * dt

k = 0.0808; mach = 0.6; gamma = 1.4; R_gas = 287.085; T_inf = 280.0
chord = 1.0; alpha_mean = 2.77; alpha_amp = 2.34
omega = 2 * mach * np.sqrt(gamma * R_gas * T_inf) * k / chord
delta_alpha = -alpha_amp * np.pi / 180.0

base_dir = os.path.dirname(os.path.abspath(__file__))
grid_file = os.path.join(base_dir, "../../../input_files/naca0012_rans-L2.cgns")
output_dir = os.path.join(base_dir, "output")
os.makedirs(output_dir, exist_ok=True)

ap = AeroProblem(
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

ts = ADflowTS(
    solver, ap,
    dt=dt, t_final=t_final, theta=1.0,
    grid_motion=True, snes_type="dadi",
    save_trajectory=True,  # Required for adjoint
    ksp_type="gmres", ksp_rtol=1e-4, ksp_max_it=200,
    pc_type="asm_ilu",
)
ts.setup()

# Add KSP monitor (sparse: first 20, then every 500th)
from petsc4py import PETSc as _PETSc
snes = ts.ts.getSNES()
ksp = snes.getKSP()
def _ksp_monitor(ksp, its, rnorm):
    if comm.rank == 0 and (its <= 20 or its % 500 == 0):
        print(f"    KSP it {its:4d}  rnorm = {rnorm:.6e}", flush=True)
ksp.setMonitor(_ksp_monitor)

# --- Forward solve ---
if rank == 0:
    print("\n=== Forward solve ===", flush=True)
t0 = time.time()
reason = ts.solve()
wall_fwd = time.time() - t0
if rank == 0:
    print(f"Forward done in {wall_fwd:.1f}s (reason={int(reason)})")
    h = ts.get_history()
    for i in range(len(h["time"])):
        print(f"  step {i:2d} | t={h['time'][i]:.5e} | CL={h['cl'][i]:.7e}")

# --- MatShell sanity check ---
if rank == 0:
    print("\n=== MatShell transpose sanity check ===", flush=True)
# Set state to final state and evaluate Jacobian at a=80
u_final = ts.U_vec.getArray(readonly=True).copy()
solver.setStates(u_final)
ts.jac_ctx.shift = 1.0 / dt  # a = 80

# Random vectors x, y
np.random.seed(42 + rank)
n_local = ts.U_vec.getLocalSize()
x_test = ts.U_vec.duplicate()
y_test = ts.U_vec.duplicate()
Jx = ts.U_vec.duplicate()
JTy = ts.U_vec.duplicate()
x_test.getArray()[:] = np.random.randn(n_local)
y_test.getArray()[:] = np.random.randn(n_local)

# Compute Jx and J^T y
ts.J_shell.mult(x_test, Jx)
ts.J_shell.multTranspose(y_test, JTy)

# Check <Jx, y> vs <x, J^T y>
dot1 = Jx.dot(y_test)
dot2 = x_test.dot(JTy)
if rank == 0:
    print(f"  <J*x, y>   = {dot1:.10e}")
    print(f"  <x, J^T*y> = {dot2:.10e}")
    rel_err = abs(dot1 - dot2) / max(abs(dot1), abs(dot2), 1e-30)
    print(f"  rel error   = {rel_err:.4e}")
    if rel_err < 1e-8:
        print("  PASS: transpose is consistent")
    else:
        print(f"  WARNING: transpose inconsistency (rel_err = {rel_err:.4e})")

# Also check matvec magnitude
Jx_norm = Jx.norm()
JTy_norm = JTy.norm()
x_norm = x_test.norm()
y_norm = y_test.norm()
if rank == 0:
    print(f"  ||x|| = {x_norm:.4e}, ||J*x|| = {Jx_norm:.4e}, ratio = {Jx_norm/x_norm:.4e}")
    print(f"  ||y|| = {y_norm:.4e}, ||J^T*y|| = {JTy_norm:.4e}, ratio = {JTy_norm/y_norm:.4e}")

# --- Adjoint solve ---
if rank == 0:
    print("\n=== Adjoint solve (objective: CL at final time) ===", flush=True)
t0 = time.time()
try:
    result = ts.solve_adjoint(["cl"])
    wall_adj = time.time() - t0
    if rank == 0:
        print(f"Adjoint done in {wall_adj:.1f}s")
        for obj, lam in result.items():
            print(f"  {obj}: ||lambda(t0)||_2 = {np.linalg.norm(lam):.6e}")
except Exception as e:
    if rank == 0:
        import traceback
        print(f"Adjoint FAILED: {e}")
        traceback.print_exc()
