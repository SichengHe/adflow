"""
Verify TS adjoint dCL/dalpha via chain rule decomposition.

Compares three estimates of dCL/dalpha:
  1. Steady adjoint (ADflow evalFunctionsSens) -- reference
  2. TS adjoint chain rule:  dCL/dalpha = pCL/palpha + lambda_0^T * dw0/dalpha
  3. FD of steady states:    dCL/dalpha ~ (CL(a+e) - CL(a-e)) / (2e)

Decomposition for the time-accurate case:
  dCL_final/dalpha = pCL/palpha + lambda_0^T * dw0/dalpha + [IJacobianP terms]

For tiny dt from steady state, IJacobianP ~ 0, so:
  dCL_final/dalpha ~ pCL/palpha + lambda_0^T * dw0/dalpha

This should match the steady adjoint result because CL_final ~ CL_steady.

Usage:
    mpirun -np 2 python test_adjoint_alpha.py
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
grid_file = os.path.join(base_dir, "../../../input_files/naca0012_rans-L2.cgns")
output_dir = os.path.join(base_dir, "output")
os.makedirs(output_dir, exist_ok=True)

alpha_base = 2.77  # baseline angle of attack [degrees]

ap = AeroProblem(
    name="alphatest", alpha=alpha_base,
    mach=0.6, machRef=0.6, reynolds=4800000.0,
    reynoldsLength=1.0, T=280.0, R=287.085,
    areaRef=1.0, chordRef=1.0,
    evalFuncs=["cl"],
)
ap.addDV("alpha", name="alpha")

# Create solver in UNSTEADY mode so all BDF/unsteady arrays are allocated.
# After creation, we switch to STEADY mode for Phase A (convergence + adjoint),
# then back to UNSTEADY for Phase B (TS forward/adjoint).
options = {
    "gridfile": grid_file,
    "outputdirectory": output_dir,
    "writevolumesolution": False,
    "writesurfacesolution": False,
    "smoother": "DADI",
    "equationtype": "RANS",
    "equationmode": "unsteady",
    "timeIntegrationscheme": "BDF",
    "ntimestepsfine": 3,
    "deltat": 1e-8,
    "timeaccuracy": 1,
    "nsubiterturb": 10,
    "nsubiter": 5,
    "cfl": 2.5, "cflcoarse": 1.2,
    "ncycles": 5000,
    "mgcycle": "3w", "mgstartlevel": 1,
    "monitorvariables": ["cpu", "resrho", "cl"],
    "usenksolver": True, "useanksolver": True,
    "l2convergence": 1e-10,
    "blockSplitting": True, "useblockettes": False,
    "printAllOptions": False, "printIterations": True,
}

solver = ADFLOW(options=options, debug=False)

# ---------------------------------------------------------------
# Monkey-patch: prevent _setAeroProblemData from changing equationmode.
# Safety measure to ensure internal setAeroProblem calls never override
# our explicit equationmode control.
# ---------------------------------------------------------------
_orig_setAPData = solver._setAeroProblemData

def _patched_setAPData(ap_obj, firstCall=False):
    save_mode = solver.adflow.inputphysics.equationmode
    _orig_setAPData(ap_obj, firstCall=firstCall)
    solver.adflow.inputphysics.equationmode = save_mode

solver._setAeroProblemData = _patched_setAPData


# ==================================================================
# Phase A: All steady-state computations
# ==================================================================
# Switch BOTH Python option AND Fortran flag to STEADY.
# Python option controls __call__ dispatch; Fortran flag controls solver behavior.
solver.options["equationmode"] = "steady"
solver.adflow.inputphysics.equationmode = _STEADY

# --- A1: Converge steady at baseline alpha ---
if rank == 0:
    print("=" * 70)
    print(f"A1: Converge steady state at alpha = {alpha_base} deg")
    print("=" * 70)

solver(ap)
w_steady = solver.getStates().copy()

funcs = {}
solver.evalFunctions(ap, funcs, evalFuncs=["cl"])
cl_steady = funcs[f"{ap.name}_cl"]
if rank == 0:
    print(f"  CL_steady = {cl_steady:.10e}")

# --- A2: Steady adjoint dCL/dalpha (gold standard) ---
if rank == 0:
    print("\n" + "=" * 70)
    print("A2: Steady adjoint -- dCL/dalpha (reference)")
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
eps_alpha = 1e-2  # [degrees]

if rank == 0:
    print("\n" + "=" * 70)
    print(f"A4: FD of steady states (eps = {eps_alpha} deg)")
    print("=" * 70)

# alpha + eps (warm-start from baseline)
if rank == 0:
    print(f"\n  --- Solving at alpha = {alpha_base + eps_alpha} ---")
solver.setStates(w_steady)
ap.alpha = alpha_base + eps_alpha
solver(ap)
w_plus = solver.getStates().copy()
funcs_p = {}
solver.evalFunctions(ap, funcs_p, evalFuncs=["cl"])
cl_plus = funcs_p[f"{ap.name}_cl"]

# alpha - eps (warm-start from baseline)
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
dw0_norm_sq = comm.allreduce(np.dot(dw0_dalpha, dw0_dalpha))

# Diagnostic: state change norms
dw_plus_sq = comm.allreduce(np.dot(w_plus - w_steady, w_plus - w_steady))
dw_minus_sq = comm.allreduce(np.dot(w_minus - w_steady, w_minus - w_steady))

if rank == 0:
    print(f"\n  CL(alpha + eps) = {cl_plus:.10e}")
    print(f"  CL(alpha - eps) = {cl_minus:.10e}")
    print(f"  dCL/dalpha (FD steady) = {dCL_dalpha_fd:.10e}")
    print(f"  dCL/dalpha (steady adj) = {dCL_dalpha_ref:.10e}")
    fd_err = abs(dCL_dalpha_fd - dCL_dalpha_ref) / abs(dCL_dalpha_ref)
    print(f"  FD vs steady adj rel error = {fd_err:.4e}")
    print(f"  ||dw0/dalpha|| = {np.sqrt(dw0_norm_sq):.6e}")
    print(f"  ||w_plus - w_steady|| = {np.sqrt(dw_plus_sq):.6e}")
    print(f"  ||w_minus - w_steady|| = {np.sqrt(dw_minus_sq):.6e}")


# ==================================================================
# Phase B: TS forward + adjoint tests (switch to UNSTEADY)
# ==================================================================
solver.options["equationmode"] = "unsteady"
solver.adflow.inputphysics.equationmode = _UNSTEADY


def run_ts_test(dt, n_steps):
    """Run TS forward + adjoint and compare dCL/dalpha estimates."""
    t_final = n_steps * dt

    if rank == 0:
        print("\n" + "#" * 70)
        print(f"# TS TEST: dt = {dt:.1e}, n_steps = {n_steps}, "
              f"a = 1/dt = {1/dt:.1e}")
        print("#" * 70)

    solver.setAeroProblem(ap)
    solver.setStates(w_steady)

    ts = ADflowTS(
        solver, ap,
        dt=dt, t_final=t_final, theta=1.0,
        grid_motion=False, snes_type="dadi",
        save_trajectory=True,
        ksp_type="gmres", ksp_rtol=1e-6, ksp_max_it=200,
        pc_type="asm_ilu",
    )
    ts.setup(skip_set_ap=True)
    ts.U_vec.setArray(w_steady)
    ts.solve()

    h = ts.get_history()
    cl_final = h["cl"][-1]

    # Adjoint
    result = ts.solve_adjoint(["cl"])
    lam = result["cl"]

    # Restore steady state in solver
    solver.setStates(w_steady)

    # Chain rule: dCL/dalpha = pCL/palpha + lambda_0^T * dw0/dalpha
    lam_dot_dw0 = comm.allreduce(np.dot(lam, dw0_dalpha))
    adj_estimate = pCL_palpha + lam_dot_dw0

    # Diagnostics
    lam_norm_sq = comm.allreduce(np.dot(lam, lam))
    lam_norm = np.sqrt(lam_norm_sq)

    if rank == 0:
        print(f"\n  CL_final (TS)   = {cl_final:.10e}")
        print(f"  CL_steady       = {cl_steady:.10e}")
        print(f"  ||lambda_0||    = {lam_norm:.6e}")
        print()
        print(f"  --- dCL/dalpha estimates ---")
        print(f"  Steady adjoint (ref): {dCL_dalpha_ref:18.10e}")
        print(f"  FD (steady states):   {dCL_dalpha_fd:18.10e}")
        print(f"  TS adjoint (chain):   {adj_estimate:18.10e}")
        print(f"    pCL/palpha:         {pCL_palpha:18.10e}")
        print(f"    lam^T * dw0/da:     {lam_dot_dw0:18.10e}")
        if abs(dCL_dalpha_ref) > 1e-30:
            rel_err = abs(adj_estimate - dCL_dalpha_ref) / abs(dCL_dalpha_ref)
            print(f"  rel error (TS vs ref): {rel_err:.4e}")

    return adj_estimate


# Test 1: Well-conditioned (dt=1e-8, a=1e8 >> spectral radius)
adj1 = run_ts_test(dt=1e-8, n_steps=3)

# Test 2: Moderate conditioning (dt=1e-5, a=1e5 ~ spectral radius)
adj2 = run_ts_test(dt=1e-5, n_steps=3)

# ==================================================================
# Summary
# ==================================================================
if rank == 0:
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Steady adjoint (ref):  dCL/dalpha = {dCL_dalpha_ref:.10e}")
    print(f"  FD (steady states):    dCL/dalpha = {dCL_dalpha_fd:.10e}")
    print(f"  TS adj (dt=1e-8):      dCL/dalpha = {adj1:.10e}")
    print(f"  TS adj (dt=1e-5):      dCL/dalpha = {adj2:.10e}")
    print()
    print(f"  Decomposition at steady state:")
    print(f"    pCL/palpha       = {pCL_palpha:.10e}")
    print(f"    -psi^T pR/palpha = {implicit_part:.10e}")
    print("=" * 70)
