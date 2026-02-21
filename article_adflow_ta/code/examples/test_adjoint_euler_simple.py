"""
Verify TA adjoint dCL/dalpha for Euler via chain rule decomposition.

For BDF1 from converged steady state with no time-varying forcing:
  CL_final ≈ CL_steady, so dCL_final/dalpha ≈ dCL_steady/dalpha.

Chain rule:
  dCL_final/dalpha = pCL/palpha + lambda_0^T * dw0/dalpha

where:
  - pCL/palpha: explicit partial (from AD reverse, at fixed state)
  - lambda_0: TS adjoint at t=0 (solve_adjoint_manual)
  - dw0/dalpha: FD of converged steady states at alpha +/- eps

Reference: ADflow steady adjoint (evalFunctionsSens).

Usage:
    mpirun -np 2 python test_adjoint_euler_simple.py
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

alpha_base = 1.8  # degrees

ap = AeroProblem(
    name="euler_dcl", alpha=alpha_base,
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
    "ntimestepsfine": 3,
    "deltat": 1e-6,
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
# Phase A: All steady-state computations
# ==================================================================
solver.options["equationmode"] = "steady"
solver.adflow.inputphysics.equationmode = _STEADY

# --- A1: Converge steady at baseline alpha ---
if rank == 0:
    print("=" * 70)
    print(f"A1: Converge steady Euler at alpha = {alpha_base} deg")
    print("=" * 70)

solver(ap)
w_steady = solver.getStates().copy()

funcs = {}
solver.evalFunctions(ap, funcs, evalFuncs=["cl"])
cl_steady = funcs[f"{ap.name}_cl"]
if rank == 0:
    print(f"  CL_steady = {cl_steady:.10e}")

# --- A2: Steady adjoint dCL/dalpha (reference) ---
if rank == 0:
    print("\n" + "=" * 70)
    print("A2: Steady adjoint dCL/dalpha (reference)")
    print("=" * 70)

funcsSens = {}
solver.evalFunctionsSens(ap, funcsSens, evalFuncs=["cl"])
dCL_dalpha_ref = funcsSens[f"{ap.name}_cl"]["alpha"]
if hasattr(dCL_dalpha_ref, "__len__"):
    dCL_dalpha_ref = float(dCL_dalpha_ref)
if rank == 0:
    print(f"  dCL/dalpha (steady adjoint) = {dCL_dalpha_ref:.10e}")

# --- A3: Explicit partial pCL/palpha at fixed state ---
if rank == 0:
    print("\n" + "=" * 70)
    print("A3: Explicit partial pCL/palpha (at fixed w)")
    print("=" * 70)

solver.setAeroProblem(ap)
solver.setStates(w_steady)
funcsBar = solver._getFuncsBar("cl")
explicit_sens = solver.computeJacobianVectorProductBwd(
    funcsBar=funcsBar, xDvDeriv=True
)
pCL_palpha = float(explicit_sens.get("alpha", 0.0))
implicit_part = dCL_dalpha_ref - pCL_palpha
if rank == 0:
    print(f"  pCL/palpha (explicit)       = {pCL_palpha:.10e}")
    print(f"  -psi^T pR/palpha (implicit) = {implicit_part:.10e}")
    print(f"  Sum (= steady adjoint)      = {pCL_palpha + implicit_part:.10e}")

# --- A4: FD of steady states at alpha +/- eps ---
eps_alpha = 1e-2  # degrees

if rank == 0:
    print("\n" + "=" * 70)
    print(f"A4: FD of steady states (eps = {eps_alpha} deg)")
    print("=" * 70)

# alpha + eps
if rank == 0:
    print(f"\n  --- Solving at alpha = {alpha_base + eps_alpha} ---")
solver.setStates(w_steady)
ap.alpha = alpha_base + eps_alpha
solver(ap)
w_plus = solver.getStates().copy()
funcs_p = {}
solver.evalFunctions(ap, funcs_p, evalFuncs=["cl"])
cl_plus = funcs_p[f"{ap.name}_cl"]

# alpha - eps
if rank == 0:
    print(f"\n  --- Solving at alpha = {alpha_base - eps_alpha} ---")
solver.setStates(w_steady)
ap.alpha = alpha_base - eps_alpha
solver(ap)
w_minus = solver.getStates().copy()
funcs_m = {}
solver.evalFunctions(ap, funcs_m, evalFuncs=["cl"])
cl_minus = funcs_m[f"{ap.name}_cl"]

# Restore baseline
ap.alpha = alpha_base
solver.setStates(w_steady)
solver.setAeroProblem(ap)

# FD results
dCL_dalpha_fd = (cl_plus - cl_minus) / (2.0 * eps_alpha)
dw0_dalpha = (w_plus - w_minus) / (2.0 * eps_alpha)
dw0_norm = np.sqrt(comm.allreduce(np.dot(dw0_dalpha, dw0_dalpha)))

if rank == 0:
    print(f"\n  CL(+eps) = {cl_plus:.10e}")
    print(f"  CL(-eps) = {cl_minus:.10e}")
    print(f"  dCL/dalpha (FD steady) = {dCL_dalpha_fd:.10e}")
    print(f"  dCL/dalpha (steady adj) = {dCL_dalpha_ref:.10e}")
    fd_err = abs(dCL_dalpha_fd - dCL_dalpha_ref) / abs(dCL_dalpha_ref)
    print(f"  FD vs steady adj rel error = {fd_err:.4e}")
    print(f"  ||dw0/dalpha|| = {dw0_norm:.6e}")


# ==================================================================
# Phase B: TS forward + adjoint (dt sweep)
# ==================================================================
solver.options["equationmode"] = "unsteady"
solver.adflow.inputphysics.equationmode = _UNSTEADY


def run_ts_test(dt, n_steps, pc_shift_factor=1.0, gmres_max_it=20):
    """Run TS forward + manual adjoint, compare dCL/dalpha via chain rule."""
    t_final = n_steps * dt
    shift = 1.0 / dt

    if rank == 0:
        print("\n" + "#" * 70)
        print(f"# TS TEST: dt={dt:.1e}, n_steps={n_steps}, shift={shift:.1e}")
        print(f"#          pc_shift_factor={pc_shift_factor}")
        print("#" * 70)

    reset_unsteady(solver)
    solver.setAeroProblem(ap)
    solver.setStates(w_steady)

    ts = ADflowTS(
        solver, ap,
        dt=dt, t_final=t_final, theta=1.0,
        grid_motion=False, snes_type="dadi",
        save_trajectory=True,
    )
    ts.setup(skip_set_ap=True)
    ts.U_vec.setArray(w_steady)
    ts.solve()

    h = ts.get_history()
    cl_final = h["cl"][-1]

    # Manual adjoint
    t0 = time.time()
    result = ts.solve_adjoint_manual(
        ["cl"], gmres_rtol=1e-10, gmres_max_it=gmres_max_it,
        gmres_restart=gmres_max_it, reassemble=True,
        adjoint_pc="ts_ilu", pc_shift_factor=pc_shift_factor,
    )
    lam = result["cl"]
    t_adj = time.time() - t0

    # Restore steady state
    solver.setStates(w_steady)

    # Chain rule: dCL/dalpha = pCL/palpha + lambda_0^T * dw0/dalpha
    lam_dot_dw0 = comm.allreduce(np.dot(lam, dw0_dalpha))
    adj_estimate = pCL_palpha + lam_dot_dw0

    lam_norm = np.sqrt(comm.allreduce(np.dot(lam, lam)))

    if rank == 0:
        print(f"\n  CL_final (TS)   = {cl_final:.10e}")
        print(f"  CL_steady       = {cl_steady:.10e}")
        print(f"  ||lambda_0||    = {lam_norm:.6e}")
        print(f"  Adjoint time    = {t_adj:.1f}s")
        print()
        print(f"  --- dCL/dalpha estimates ---")
        print(f"  Steady adjoint (ref): {dCL_dalpha_ref:18.10e}")
        print(f"  FD (steady states):   {dCL_dalpha_fd:18.10e}")
        print(f"  TS adjoint (chain):   {adj_estimate:18.10e}")
        print(f"    pCL/palpha:         {pCL_palpha:18.10e}")
        print(f"    lam^T * dw0/da:     {lam_dot_dw0:18.10e}")
        if abs(dCL_dalpha_ref) > 1e-30:
            rel_err = abs(adj_estimate - dCL_dalpha_ref) / abs(dCL_dalpha_ref)
            print(f"  Rel error (TS vs ref): {rel_err:.4e}")

    return adj_estimate


if rank == 0:
    print("\n" + "=" * 70)
    print("Phase B: TS adjoint dCL/dalpha verification")
    print("=" * 70)

# Test 1: Well-conditioned (dt=1e-6, shift=1e6)
adj1 = run_ts_test(dt=1e-6, n_steps=3, pc_shift_factor=1.0, gmres_max_it=20)

# Test 2: Moderate (dt=1e-4, shift=1e4)
adj2 = run_ts_test(dt=1e-4, n_steps=3, pc_shift_factor=1.0, gmres_max_it=50)

# ==================================================================
# Summary
# ==================================================================
if rank == 0:
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Steady adjoint (ref):  dCL/dalpha = {dCL_dalpha_ref:.10e}")
    print(f"  FD (steady states):    dCL/dalpha = {dCL_dalpha_fd:.10e}")
    print(f"  TS adj (dt=1e-6):      dCL/dalpha = {adj1:.10e}")
    print(f"  TS adj (dt=1e-4):      dCL/dalpha = {adj2:.10e}")
    print()
    print(f"  Decomposition:")
    print(f"    pCL/palpha       = {pCL_palpha:.10e}")
    print(f"    -psi^T pR/palpha = {implicit_part:.10e}")
    print("=" * 70)
