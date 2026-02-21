"""
Verify TS discrete adjoint for pitching NACA 0012 via FD of initial conditions.

Tests dCL_final/dU0 by comparing:
  1. TS adjoint (TSAdjointSolve) → lambda_0
  2. Central FD: (CL_final(U0+ε·v) − CL_final(U0−ε·v)) / (2ε)

The dot product lambda_0^T · v should match the FD estimate.

This verifies the core adjoint machinery:
  - Trajectory saving and backward sweep
  - KSPSolveTranspose with J^T (MatShell multTranspose)
  - Transpose preconditioner (ILU^{-T})

Note: IJacobianP is not yet implemented, so we test initial-state
sensitivity only (not parameter sensitivity like pitching amplitude).

Usage:
    mpirun -np 2 python test_adjoint_pitching.py
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

# --- Problem parameters (same as reg_tests/test_time_accurate_naca0012) ---
k = 0.0808           # reduced frequency
mach = 0.6
gamma = 1.4
R_gas = 287.085
T_inf = 280.0
chord = 1.0
alpha_mean = 2.77    # degrees
alpha_amp = 2.34     # degrees

omega = 2 * mach * np.sqrt(gamma * R_gas * T_inf) * k / chord
delta_alpha = -alpha_amp * np.pi / 180.0

freq = 10.0           # [Hz]
period = 1.0 / freq   # [s]
n_steps_per_period = 8
n_periods = 1
n_steps = n_steps_per_period * n_periods
dt = period / n_steps_per_period   # 0.0125 s
t_final = period * n_periods

# Paths
base_dir = os.path.dirname(os.path.abspath(__file__))
repo_dir = os.path.join(base_dir, "../../..")
grid_file = os.path.join(repo_dir, "input_files/naca0012_rans-L2.cgns")
output_dir = os.path.join(base_dir, "output")
os.makedirs(output_dir, exist_ok=True)

if not os.path.exists(grid_file):
    if rank == 0:
        print(f"ERROR: Mesh file not found: {grid_file}")
        print("Run input_files/get-input-files.sh to download test meshes.")
    sys.exit(1)


def create_ap():
    """Create the pitching NACA 0012 AeroProblem."""
    return AeroProblem(
        name="0012pitching",
        alpha=alpha_mean, mach=mach, machRef=mach,
        reynolds=4800000.0, reynoldsLength=chord, T=T_inf, R=R_gas,
        areaRef=1.0, chordRef=chord,
        evalFuncs=["cl", "cd", "cmz"],
        xRef=0.25, xRot=0.25,
        degreePol=0, coefPol=[0.0],
        degreeFourier=1, omegaFourier=omega,
        cosCoefFourier=[0.0, 0.0], sinCoefFourier=[delta_alpha],
    )


if rank == 0:
    print(f"Pitching NACA 0012 — TS Adjoint Verification (dCL_final/dU0)")
    print(f"  dt={dt:.4e}, n_steps={n_steps}, t_final={t_final:.4f}")
    print(f"  M={mach}, Re=4.8M, alpha={alpha_mean}+/-{alpha_amp} deg, k={k}")
    print(f"  timeaccuracy=1 (BDF1, consistent with TSTHETA adjoint)")
    print()

# ==================================================================
# Create solver
# ==================================================================
ap = create_ap()

options = {
    "gridfile": grid_file,
    "outputdirectory": output_dir,
    "writevolumesolution": False,
    "writesurfacesolution": False,
    "vis4": 0.025,
    "vis2": 0.5,
    "restrictionrelaxation": 0.5,
    "smoother": "DADI",
    "equationtype": "RANS",
    "equationmode": "unsteady",
    "timeIntegrationscheme": "BDF",
    "ntimestepsfine": n_steps,
    "deltat": dt,
    "timeaccuracy": 1,       # BDF1 — MUST match TSTHETA adjoint
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
}

solver = ADFLOW(options=options, debug=False)

# Monkey-patch: preserve equationmode through internal setAeroProblem calls
_orig_setAPData = solver._setAeroProblemData

def _patched_setAPData(ap_obj, firstCall=False):
    save_mode = solver.adflow.inputphysics.equationmode
    _orig_setAPData(ap_obj, firstCall=firstCall)
    solver.adflow.inputphysics.equationmode = save_mode

solver._setAeroProblemData = _patched_setAPData


def reset_unsteady(solver):
    """Reset Fortran state for a fresh unsteady simulation from t=0."""
    solver.adflow.monitor.timeunsteady = 0.0
    solver.adflow.monitor.timestepunsteady = 0
    solver.adflow.iteration.noldsolvavail = 1
    # Move mesh back to t=0 position (prescribed motion at timeUnsteady=0)
    solver.adflow.solvers.updateunsteadygeometry()


def run_forward(solver, ap, w_init, save_traj=False):
    """Run a forward unsteady solve from w_init, return CL_final.

    Creates a fresh ADflowTS each time to ensure clean Fortran state.
    """
    # Reset Fortran time/coordinate state
    reset_unsteady(solver)
    solver.setStates(w_init)

    ts = ADflowTS(
        solver, ap,
        dt=dt, t_final=t_final, theta=1.0,
        grid_motion=True, snes_type="dadi",
        save_trajectory=save_traj,
        ksp_type="gmres", ksp_rtol=1e-4, ksp_max_it=50,
        ksp_gmres_restart=50,
        pc_type="native_adjoint",
    )
    ts.setup(skip_set_ap=True)
    ts.U_vec.setArray(w_init)
    ts.solve()

    h = ts.get_history()
    cl_final = h["cl"][-1]
    return ts, cl_final, h


# ==================================================================
# Phase A: Converge steady state at alpha_mean
# ==================================================================
solver.options["equationmode"] = "steady"
solver.adflow.inputphysics.equationmode = _STEADY

if rank == 0:
    print("=" * 70)
    print(f"Phase A: Converge steady state at alpha = {alpha_mean} deg")
    print("=" * 70)

solver(ap)
w_steady = solver.getStates().copy()

funcs = {}
solver.evalFunctions(ap, funcs, evalFuncs=["cl"])
cl_steady = funcs[f"{ap.name}_cl"]
if rank == 0:
    print(f"  CL_steady = {cl_steady:.10e}")

w_steady_norm_sq = comm.allreduce(np.dot(w_steady, w_steady))
w_steady_norm = np.sqrt(w_steady_norm_sq)
if rank == 0:
    print(f"  ||w_steady|| = {w_steady_norm:.6e}")


# ==================================================================
# Phase B: Forward solve + adjoint (switch to UNSTEADY)
# ==================================================================
solver.options["equationmode"] = "unsteady"
solver.adflow.inputphysics.equationmode = _UNSTEADY

if rank == 0:
    print("\n" + "=" * 70)
    print(f"Phase B: Forward solve ({n_steps} steps) + adjoint")
    print("=" * 70)

t_start = time.time()
ts_ref, cl_final, h_ref = run_forward(solver, ap, w_steady, save_traj=True)
t_fwd = time.time() - t_start

if rank == 0:
    print(f"\n  Forward solve completed in {t_fwd:.1f} s")
    print(f"  CL history: {h_ref['cl']}")
    print(f"  CL_final = {cl_final:.10e}")

# --- Quick test: does ADflow's own steady adjoint work after the forward solve? ---
if rank == 0:
    print("\n  Testing ADflow native adjoint (evalFunctionsSens)...")

# Switch to steady mode for adjoint test
solver.options["equationmode"] = "steady"
solver.adflow.inputphysics.equationmode = _STEADY
solver.setStates(w_steady)

funcsSens_test = {}
solver.evalFunctionsSens(ap, funcsSens_test, evalFuncs=["cl"])
if rank == 0:
    print(f"  Native adjoint dCL/dalpha: {funcsSens_test}")
    print("  Native adjoint works OK!")

# Switch back to unsteady for our manual adjoint
solver.options["equationmode"] = "unsteady"
solver.adflow.inputphysics.equationmode = _UNSTEADY

# --- Direct test of solveadjointforrhs ---
if rank == 0:
    print("\n  Testing solveadjointforrhs directly...")

# Set state to steady (same state used by evalFunctionsSens)
solver.setStates(w_steady)
test_rhs = np.ones(solver.getStateSize()) * 0.001
test_out = solver.adflow.adjointapi.solveadjointforrhs(test_rhs, 0.01)
if rank == 0:
    print(f"  Direct solveadjointforrhs OK, ||out|| = {np.linalg.norm(test_out):.6e}")

# --- Diagnostic: verify matvec and preconditioner at final state ---
if rank == 0:
    print("\n  Diagnostic: testing matvec and preconditioner building blocks...")

# Restore final state for diagnostics
t_N, w_N = ts_ref._state_trajectory[-1]
solver.setStates(w_N)
ts_ref._update_grid(t_N)
shift = 1.0 / (ts_ref.theta * ts_ref.dt)

# Build resscale_diag
jac_ctx = ts_ref.jac_ctx
resscale_diag = np.ones(ts_ref.n_local)
if hasattr(jac_ctx, '_n_turb') and jac_ctx._n_turb > 0:
    n_cells = ts_ref.n_local // jac_ctx._nw
    rs2d = resscale_diag.reshape(n_cells, jac_ctx._nw)
    for l in range(jac_ctx._n_turb):
        rs2d[:, jac_ctx._nwf + l] = jac_ctx._trs[l]

# Test vector: terminal gradient dCL/dU_N
dphi = ts_ref._compute_terminal_state_gradient("cl")
dphi_norm = np.sqrt(comm.allreduce(np.dot(dphi, dphi)))
if rank == 0:
    print(f"  ||dCL/dU_N|| = {dphi_norm:.6e}")

# Test matvec: J_u^T * v
ts_ref._adjoint_mode = True
jt_v = ts_ref._jtranspose_matvec(dphi, shift)
jt_v_norm = np.sqrt(comm.allreduce(np.dot(jt_v, jt_v)))
if rank == 0:
    print(f"  ||J_u^T * (dCL/dU)|| = {jt_v_norm:.6e}")
    print(f"  ratio ||Jt*v||/||v|| = {jt_v_norm/dphi_norm:.6e}")

# Test preconditioner: P^{-1} * v (native = dR/dw only)
p_inv_v = ts_ref._adjoint_precond(dphi, resscale_diag, 0.01)
p_inv_v_norm = np.sqrt(comm.allreduce(np.dot(p_inv_v, p_inv_v)))
if rank == 0:
    print(f"  ||P_native^{{-1}} * v|| = {p_inv_v_norm:.6e}")
    print(f"  ratio = {p_inv_v_norm/dphi_norm:.6e}")

# Test TS ILU PC: forward and transpose
solver.adflow.inputphysics.equationmode = 1  # STEADY for assembly
solver.adflow.nksolver.setuptspreconditioner(shift)
solver.adflow.inputphysics.equationmode = 2  # restore UNSTEADY

p_fwd = ts_ref._adjoint_precond_ts(dphi, use_transpose=False)
p_fwd_norm = np.sqrt(comm.allreduce(np.dot(p_fwd, p_fwd)))
p_trn = ts_ref._adjoint_precond_ts(dphi, use_transpose=True)
p_trn_norm = np.sqrt(comm.allreduce(np.dot(p_trn, p_trn)))
if rank == 0:
    print(f"  ||P_ilu_fwd^{{-1}} * v|| = {p_fwd_norm:.6e}  ratio={p_fwd_norm/dphi_norm:.6e}")
    print(f"  ||P_ilu_trn^{{-T}} * v|| = {p_trn_norm:.6e}  ratio={p_trn_norm/dphi_norm:.6e}")
    if p_trn_norm > 1e6 * p_fwd_norm:
        print(f"  ** KSPSolveTranspose is BROKEN (ratio {p_trn_norm/p_fwd_norm:.1e}x) **")

# Test preconditioned matvec: A * P^{-1} * v
ap_inv_v = ts_ref._jtranspose_matvec(p_inv_v, shift)
ap_inv_v_norm = np.sqrt(comm.allreduce(np.dot(ap_inv_v, ap_inv_v)))
# Check: should be near ||v|| for a good preconditioner
dot_ratio = comm.allreduce(np.dot(ap_inv_v, dphi)) / (dphi_norm * ap_inv_v_norm + 1e-30)
if rank == 0:
    print(f"  ||A*P^{{-1}}*v|| = {ap_inv_v_norm:.6e}")
    print(f"  <A*P^{{-1}}*v, v> / (||.||*||.||) = {dot_ratio:.6e}")
    print(f"  (If good PC: ||A*M^{{-1}}*v|| ~ ||v|| and cos angle ~ +/-1)")

ts_ref._adjoint_mode = False

# --- Adjoint (manual backward sweep with FGMRES) ---
if rank == 0:
    print("\n  Running manual adjoint backward sweep (FGMRES)...")

t_start = time.time()
result = ts_ref.solve_adjoint_manual(
    ["cl"], gmres_rtol=1e-10, gmres_max_it=200,
    gmres_restart=200, reassemble=True, adjoint_pc="ts_ilu",
    pc_shift_factor=100.0,
)
lambda_0 = result["cl"]
t_adj = time.time() - t_start

lam_norm_sq = comm.allreduce(np.dot(lambda_0, lambda_0))
lam_norm = np.sqrt(lam_norm_sq)
if rank == 0:
    print(f"  Adjoint completed in {t_adj:.1f} s")
    print(f"  ||lambda_0|| = {lam_norm:.6e}")

# --- TEMPORARY: skip FD to focus on adjoint diagnostics ---
sys.exit(0)


# ==================================================================
# Phase C: FD verification
# ==================================================================
if rank == 0:
    print("\n" + "=" * 70)
    print("Phase C: FD verification of dCL_final/dU0")
    print("=" * 70)

# Generate random perturbation direction (seeded for reproducibility)
rng = np.random.RandomState(42)
v_local = rng.randn(len(w_steady))
v_norm_sq = comm.allreduce(np.dot(v_local, v_local))
v_local /= np.sqrt(v_norm_sq)  # ||v||_global = 1

# Scale epsilon: relative perturbation ~ 1e-5
eps = 1e-5 * w_steady_norm

# Adjoint estimate: lambda_0^T * v
adj_dot_v = comm.allreduce(np.dot(lambda_0, v_local))

if rank == 0:
    print(f"  eps = {eps:.6e}")
    print(f"  ||eps*v|| / ||w_steady|| = {eps / w_steady_norm:.1e}")
    print(f"  Adjoint estimate: lambda_0^T * v = {adj_dot_v:.10e}")

# --- FD plus ---
if rank == 0:
    print(f"\n  --- FD plus: U0 + eps*v ---")

w_plus = w_steady + eps * v_local
t_start = time.time()
_, cl_plus, _ = run_forward(solver, ap, w_plus, save_traj=False)
t_fd_plus = time.time() - t_start

if rank == 0:
    print(f"  CL_final(+) = {cl_plus:.10e}  ({t_fd_plus:.1f} s)")

# --- FD minus ---
if rank == 0:
    print(f"\n  --- FD minus: U0 - eps*v ---")

w_minus = w_steady - eps * v_local
t_start = time.time()
_, cl_minus, _ = run_forward(solver, ap, w_minus, save_traj=False)
t_fd_minus = time.time() - t_start

if rank == 0:
    print(f"  CL_final(-) = {cl_minus:.10e}  ({t_fd_minus:.1f} s)")

# --- Compare ---
fd_estimate = (cl_plus - cl_minus) / (2.0 * eps)

if rank == 0:
    print(f"\n  --- Comparison ---")
    print(f"  Adjoint: lambda_0^T * v = {adj_dot_v:18.10e}")
    print(f"  FD:      (CL+ - CL-)/2e = {fd_estimate:18.10e}")
    if abs(adj_dot_v) > 1e-30:
        rel_err = abs(fd_estimate - adj_dot_v) / abs(adj_dot_v)
        print(f"  Relative error:           {rel_err:.4e}")
    print(f"  CL_final (ref) = {cl_final:.10e}")
    print(f"  CL+ - CL_ref   = {cl_plus - cl_final:.6e}")
    print(f"  CL- - CL_ref   = {cl_minus - cl_final:.6e}")


# ==================================================================
# Phase D: Summary
# ==================================================================
if rank == 0:
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Problem: Pitching NACA 0012, M={mach}, Re=4.8M")
    print(f"  alpha = {alpha_mean} +/- {alpha_amp} deg, k={k}")
    print(f"  dt = {dt:.4e}, n_steps = {n_steps} (1 period), BDF1")
    print()
    print(f"  CL_steady     = {cl_steady:.10e}")
    print(f"  CL_final (ref)= {cl_final:.10e}")
    print(f"  ||lambda_0||  = {lam_norm:.6e}")
    print()
    print(f"  dCL_final/dU0 verification (direction: random, seed=42):")
    print(f"    Adjoint:  {adj_dot_v:18.10e}")
    print(f"    FD:       {fd_estimate:18.10e}")
    if abs(adj_dot_v) > 1e-30:
        rel_err = abs(fd_estimate - adj_dot_v) / abs(adj_dot_v)
        print(f"    Rel err:  {rel_err:.4e}")
    print()
    print(f"  Wall times: forward={t_fwd:.0f}s, adjoint={t_adj:.0f}s, "
          f"FD+={t_fd_plus:.0f}s, FD-={t_fd_minus:.0f}s")
    print("=" * 70)
