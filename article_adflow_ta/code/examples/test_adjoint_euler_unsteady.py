"""
Verify TA adjoint dCL/dalpha for Euler with truly UNSTEADY dynamics.

Unlike test_adjoint_euler_simple.py (which starts from steady state and
CL barely changes), this test starts from a converged steady state at
alpha_ss and runs BDF1 at a DIFFERENT alpha_run.  Because w_steady(alpha_ss)
!= w_steady(alpha_run), the flow evolves and CL changes at every step.

The initial condition w0 = w_steady(alpha_ss) does NOT depend on alpha_run,
so dw0/dalpha_run = 0.  The total derivative is:

    dCL/dalpha_run = pCL/palpha + dt * SUM delta_n^T * pR_n/palpha

This is computed by solve_adjoint_manual with dv_sens=["alpha"].
The FD reference runs the full BDF1 trajectory at alpha_run +/- eps.

Usage:
    mpirun -np 2 python test_adjoint_euler_unsteady.py
"""
import os, sys, time
import numpy as np
from mpi4py import MPI
from baseclasses import AeroProblem
from adflow import ADFLOW
from adflow.pyADflow_TA import ADflowTS

_STEADY = 1
_UNSTEADY = 2

sys.stdout.reconfigure(line_buffering=True)
comm = MPI.COMM_WORLD
rank = comm.rank

base_dir = os.path.dirname(os.path.abspath(__file__))
grid_file = os.path.join(base_dir, "../../../input_files/mdo_tutorial_euler.cgns")
output_dir = os.path.join(base_dir, "output")
os.makedirs(output_dir, exist_ok=True)

if not os.path.exists(grid_file):
    if rank == 0:
        print(f"ERROR: Mesh file not found: {grid_file}")
    sys.exit(1)

alpha_ss = 1.8     # steady-state alpha (for IC)
alpha_run = 3.0    # alpha for the unsteady run (different from alpha_ss)
eps_alpha = 1e-2   # FD step (degrees)

ap = AeroProblem(
    name="euler_unst", alpha=alpha_ss,
    mach=0.80,
    P=20000.0, T=220.0, R=287.87,
    areaRef=45.5, chordRef=3.25,
    xRef=0.0, yRef=0.0, zRef=0.0,
    evalFuncs=["cl"],
)
ap.addDV("alpha", name="alpha")

options = {
    "gridfile": grid_file,
    "outputdirectory": output_dir,
    "writevolumesolution": False,
    "writesurfacesolution": False,
    "smoother": "DADI",
    "equationtype": "Euler",
    "equationmode": "unsteady",
    "timeIntegrationscheme": "BDF",
    "ntimestepsfine": 10,
    "deltat": 1e-4,
    "timeaccuracy": 1,
    "nsubiter": 200,
    "nsubiterturb": 10,
    "cfl": 2.5,
    "cflcoarse": 1.2,
    "ncycles": 5000,
    "mgcycle": "3w",
    "mgstartlevel": 1,
    "monitorvariables": ["cpu", "resrho", "cl"],
    "usenksolver": True,
    "useanksolver": True,
    "l2convergence": 1e-10,
    "l2convergencecoarse": 1e-4,
    "blockSplitting": True,
    "useblockettes": False,
    "printAllOptions": False,
    "printIterations": False,
}

solver = ADFLOW(options=options, debug=False)

# Monkey-patch to preserve equationmode
_orig_setAPData = solver._setAeroProblemData


def _patched_setAPData(ap_obj, firstCall=False):
    save_mode = solver.adflow.inputphysics.equationmode
    _orig_setAPData(ap_obj, firstCall=firstCall)
    solver.adflow.inputphysics.equationmode = save_mode


solver._setAeroProblemData = _patched_setAPData


def reset_unsteady(solver):
    solver.adflow.monitor.timeunsteady = 0.0
    solver.adflow.monitor.timestepunsteady = 0
    solver.adflow.iteration.noldsolvavail = 1


# ==================================================================
# Phase A: Converge steady at alpha_ss to get initial condition
# ==================================================================
solver.options["equationmode"] = "steady"
solver.adflow.inputphysics.equationmode = _STEADY

if rank == 0:
    print("=" * 70)
    print(f"Phase A: Converge steady Euler at alpha = {alpha_ss} deg")
    print("=" * 70)

solver(ap)
w_steady = solver.getStates().copy()

funcs = {}
solver.evalFunctions(ap, funcs, evalFuncs=["cl"])
cl_ss = funcs[f"{ap.name}_cl"]
if rank == 0:
    print(f"  CL_steady(alpha={alpha_ss}) = {cl_ss:.10e}")

# --- Diagnostic: verify resBar → alpha sensitivity ---
if rank == 0:
    print("\n  --- Diagnostic: resBar→alpha sensitivity ---")

# Get steady adjoint vector
funcsSens = {}
solver.evalFunctionsSens(ap, funcsSens, evalFuncs=["cl"])
dCL_dalpha_ref = float(funcsSens[f"{ap.name}_cl"]["alpha"])

# Explicit partial
solver.setAeroProblem(ap)
solver.setStates(w_steady)
funcsBar = solver._getFuncsBar("cl")
explicit_sens = solver.computeJacobianVectorProductBwd(
    funcsBar=funcsBar, xDvDeriv=True
)
pCL_palpha = float(explicit_sens.get("alpha", 0.0))
implicit_ref = dCL_dalpha_ref - pCL_palpha

# Test: resBar=ones, xDvDeriv=True
solver.setStates(w_steady)
ones_vec = np.ones(solver.getStateSize())
sens_res_only = solver.computeJacobianVectorProductBwd(
    resBar=ones_vec, xDvDeriv=True
)
alpha_res_only = float(sens_res_only.get("alpha", 0.0)) if isinstance(sens_res_only, dict) else None

# Test: both resBar and funcsBar together (like evalFunctionsSens)
solver.setStates(w_steady)
sens_both = solver.computeJacobianVectorProductBwd(
    resBar=ones_vec, funcsBar=funcsBar, xDvDeriv=True
)
alpha_both = float(sens_both.get("alpha", 0.0)) if isinstance(sens_both, dict) else None

# Test: forward mode dR/dalpha
solver.setStates(w_steady)
res_dot = solver.computeJacobianVectorProductFwd(
    xDvDot={"alpha": 1.0}, residualDeriv=True
)
dR_dalpha_norm = np.sqrt(comm.allreduce(np.dot(res_dot, res_dot)))

if rank == 0:
    print(f"  dCL/dalpha (steady adj ref) = {dCL_dalpha_ref:.10e}")
    print(f"  pCL/palpha (explicit)       = {pCL_palpha:.10e}")
    print(f"  implicit part (ref)         = {implicit_ref:.10e}")
    print(f"  resBar=ones → alpha sens    = {alpha_res_only}")
    print(f"  both resBar+funcsBar → alpha= {alpha_both}")
    print(f"  ||dR/dalpha|| (fwd mode)    = {dR_dalpha_norm:.6e}")
    print(f"  resBar type = {type(sens_res_only)}")


# ==================================================================
# Phase B: Helper — run forward BDF1 at given alpha from w0
# ==================================================================
def run_forward(alpha_val, w0, dt, n_steps):
    """Run BDF1 forward from w0 at alpha_val, return CL_final and trajectory."""
    solver.options["equationmode"] = "unsteady"
    solver.adflow.inputphysics.equationmode = _UNSTEADY
    reset_unsteady(solver)
    ap.alpha = alpha_val
    solver.setAeroProblem(ap)
    solver.setStates(w0)

    t_final = n_steps * dt
    ts = ADflowTS(
        solver, ap,
        dt=dt, t_final=t_final, theta=1.0,
        grid_motion=False, snes_type="dadi",
        save_trajectory=True,
    )
    ts.setup(skip_set_ap=True)
    ts.U_vec.setArray(w0)
    ts.solve()

    h = ts.get_history()
    cl_final = h["cl"][-1]
    return cl_final, ts


# ==================================================================
# Phase C: Run forward at alpha_run from w_steady(alpha_ss)
# ==================================================================
dt = 1e-4
n_steps = 10

if rank == 0:
    print("\n" + "=" * 70)
    print(f"Phase C: Forward BDF1 at alpha_run={alpha_run}, dt={dt:.1e}, N={n_steps}")
    print(f"  Initial condition: w_steady(alpha={alpha_ss})")
    print("=" * 70)

cl_base, ts_base = run_forward(alpha_run, w_steady, dt, n_steps)

if rank == 0:
    print(f"\n  CL_final (base) = {cl_base:.10e}")
    print(f"  CL_steady(alpha_ss={alpha_ss}) = {cl_ss:.10e}")
    print(f"  CL changed by {cl_base - cl_ss:.6e}")


# ==================================================================
# Phase D: FD — run forward at alpha_run +/- eps from same w0
# ==================================================================
if rank == 0:
    print("\n" + "=" * 70)
    print(f"Phase D: FD at alpha_run +/- eps (eps={eps_alpha})")
    print("=" * 70)

# alpha_run + eps
if rank == 0:
    print(f"\n  --- Forward at alpha = {alpha_run + eps_alpha} ---")
cl_plus, _ = run_forward(alpha_run + eps_alpha, w_steady, dt, n_steps)

# alpha_run - eps
if rank == 0:
    print(f"\n  --- Forward at alpha = {alpha_run - eps_alpha} ---")
cl_minus, _ = run_forward(alpha_run - eps_alpha, w_steady, dt, n_steps)

dCL_dalpha_fd = (cl_plus - cl_minus) / (2.0 * eps_alpha)

if rank == 0:
    print(f"\n  CL(+eps) = {cl_plus:.10e}")
    print(f"  CL(-eps) = {cl_minus:.10e}")
    print(f"  dCL/dalpha (FD) = {dCL_dalpha_fd:.10e}")


# ==================================================================
# Phase E: Adjoint with dv_sens
# ==================================================================
if rank == 0:
    print("\n" + "=" * 70)
    print("Phase E: TS adjoint with dv_sens (intermediate param accumulation)")
    print("=" * 70)

# Restore state to base trajectory (need ts_base's trajectory)
# Re-run forward to get clean trajectory state
solver.options["equationmode"] = "unsteady"
solver.adflow.inputphysics.equationmode = _UNSTEADY
reset_unsteady(solver)
ap.alpha = alpha_run
solver.setAeroProblem(ap)
solver.setStates(w_steady)

ts = ADflowTS(
    solver, ap,
    dt=dt, t_final=n_steps * dt, theta=1.0,
    grid_motion=False, snes_type="dadi",
    save_trajectory=True,
)
ts.setup(skip_set_ap=True)
ts.U_vec.setArray(w_steady)
ts.solve()

h = ts.get_history()
cl_check = h["cl"][-1]
if rank == 0:
    print(f"  CL_final check = {cl_check:.10e} (should match {cl_base:.10e})")

# Adjoint with parameter sensitivity
result = ts.solve_adjoint_manual(
    ["cl"], gmres_rtol=1e-10, gmres_max_it=30,
    gmres_restart=30, reassemble=True,
    adjoint_pc="ts_ilu", pc_shift_factor=1.0,
    dv_sens=["alpha"],
)

lam_0 = result["cl"]
dJdp = result["cl_dJdp"]

# Since w0 = w_steady(alpha_ss) does NOT depend on alpha_run:
# dw0/dalpha_run = 0
# Total: dCL/dalpha_run = dJdp["alpha"] + lambda_0^T * 0 = dJdp["alpha"]
adj_estimate = dJdp["alpha"]

lam_norm = np.sqrt(comm.allreduce(np.dot(lam_0, lam_0)))

if rank == 0:
    print(f"\n  ||lambda_0||   = {lam_norm:.6e}")
    print(f"  dJdp['alpha']  = {adj_estimate:.10e}")
    print(f"  dCL/dalpha FD  = {dCL_dalpha_fd:.10e}")
    if abs(dCL_dalpha_fd) > 1e-30:
        rel_err = abs(adj_estimate - dCL_dalpha_fd) / abs(dCL_dalpha_fd)
        print(f"  Rel error      = {rel_err:.4e}")


# ==================================================================
# Phase F: Also test the simplified formula (no intermediate terms)
#          to show it's WRONG for unsteady case
# ==================================================================
if rank == 0:
    print("\n" + "=" * 70)
    print("Phase F: Compare simplified vs full formula")
    print("=" * 70)

# Compute pCL/palpha at final state (explicit partial only)
solver.setStates(ts._state_trajectory[-1][1])
ap.alpha = alpha_run
solver.setAeroProblem(ap)
funcsBar = solver._getFuncsBar("cl")
adflow = solver.adflow
orig_mode = adflow.inputphysics.equationmode
adflow.inputphysics.equationmode = _STEADY
explicit_sens = solver.computeJacobianVectorProductBwd(
    funcsBar=funcsBar, xDvDeriv=True
)
adflow.inputphysics.equationmode = orig_mode
pCL_palpha = float(explicit_sens.get("alpha", 0.0))

# Simplified formula (WRONG for unsteady): pCL/palpha + lam_0^T * dw0/dalpha
# Since dw0/dalpha_run = 0: simplified = pCL/palpha + 0 = pCL/palpha
simplified = pCL_palpha

# Full formula: dJdp["alpha"] = pCL/palpha + dt*sum(delta_n^T * pR_n/palpha)
full = adj_estimate

# Intermediate accumulation = full - pCL/palpha
intermediate = full - pCL_palpha

if rank == 0:
    print(f"  pCL/palpha (terminal partial):  {pCL_palpha:.10e}")
    print(f"  Intermediate accumulation:      {intermediate:.10e}")
    print(f"  Full formula (adj):             {full:.10e}")
    print(f"  Simplified (terminal only):     {simplified:.10e}")
    print(f"  FD reference:                   {dCL_dalpha_fd:.10e}")
    print()
    if abs(dCL_dalpha_fd) > 1e-30:
        err_full = abs(full - dCL_dalpha_fd) / abs(dCL_dalpha_fd)
        err_simp = abs(simplified - dCL_dalpha_fd) / abs(dCL_dalpha_fd)
        print(f"  Rel error (full formula):       {err_full:.4e}")
        print(f"  Rel error (simplified, WRONG):  {err_simp:.4e}")
        print(f"  Intermediate term matters:      {err_simp > 10 * err_full}")
    print("=" * 70)
