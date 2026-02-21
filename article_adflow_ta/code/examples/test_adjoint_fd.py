"""
Verify TS adjoint via finite-difference comparison.

Computes:
  adjoint:  dJ/dU0 . v  =  lambda(t0)^T * v
  FD:       [CL(U0 + eps*v) - CL(U0 - eps*v)] / (2*eps)

where J = CL at the final time, U0 is the initial state, and v is a
random perturbation direction.

Tests two dt values:
  - dt=1e-8 (a=1e8 >> spectral radius): well-conditioned, validates FD infrastructure
  - dt=1e-5 (a=1e5 ~ spectral radius): tests adjoint accuracy at moderate conditioning

Usage:
    mpirun -np 2 python test_adjoint_fd.py
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
size = comm.size

base_dir = os.path.dirname(os.path.abspath(__file__))
grid_file = os.path.join(base_dir, "../../../input_files/naca0012_rans-L2.cgns")
output_dir = os.path.join(base_dir, "output")
os.makedirs(output_dir, exist_ok=True)

ap = AeroProblem(
    name="fdtest", alpha=2.77,
    mach=0.6, machRef=0.6, reynolds=4800000.0,
    reynoldsLength=1.0, T=280.0, R=287.085,
    areaRef=1.0, chordRef=1.0,
    evalFuncs=["cl"],
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
    "ntimestepsfine": 3, "deltat": 1e-8,  # Will be overridden per test
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

# ---------------------------------------------------------------
# Helper: run forward solve from a given initial state, return CL
# ---------------------------------------------------------------
def forward_solve(w0, dt, n_steps):
    """Run forward from initial state w0, return final CL."""
    t_final = n_steps * dt
    solver.setAeroProblem(ap)
    solver.setStates(w0)

    ts = ADflowTS(
        solver, ap,
        dt=dt, t_final=t_final, theta=1.0,
        grid_motion=False, snes_type="dadi",
        save_trajectory=False,
        ksp_type="gmres", ksp_rtol=1e-4, ksp_max_it=200,
        pc_type="asm_ilu",
    )
    ts.setup(skip_set_ap=True)
    ts.U_vec.setArray(w0)
    ts.solve()

    h = ts.get_history()
    return h["cl"][-1]


def run_fd_test(dt, n_steps, ksp_rtol=1e-6, eps_list=None):
    """Run adjoint + FD comparison at given dt."""
    t_final = n_steps * dt
    if eps_list is None:
        eps_list = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7]

    if rank == 0:
        print("\n" + "#" * 70)
        print(f"#  FD TEST: dt = {dt:.1e}, a = 1/dt = {1/dt:.1e}, "
              f"n_steps = {n_steps}, ksp_rtol = {ksp_rtol:.0e}")
        print("#" * 70)

    # --- Step 1: Get baseline state ---
    if rank == 0:
        print("\n  Step 1: Baseline forward solve", flush=True)

    solver.setAeroProblem(ap)
    w0_base = solver.getStates().copy()
    cl_base = forward_solve(w0_base, dt, n_steps)
    if rank == 0:
        print(f"  CL(base) = {cl_base:.12e}")

    # --- Step 2: Adjoint solve ---
    if rank == 0:
        print("\n  Step 2: Adjoint solve", flush=True)

    solver.setAeroProblem(ap)
    solver.setStates(w0_base)

    ts_adj = ADflowTS(
        solver, ap,
        dt=dt, t_final=t_final, theta=1.0,
        grid_motion=False, snes_type="dadi",
        save_trajectory=True,
        ksp_type="gmres", ksp_rtol=ksp_rtol, ksp_max_it=200,
        pc_type="asm_ilu",
    )
    ts_adj.setup(skip_set_ap=True)
    ts_adj.U_vec.setArray(w0_base)
    ts_adj.solve()

    # Terminal condition diagnostics (on all ranks)
    solver.setStates(ts_adj.U_vec.getArray(readonly=True).copy())
    dphi_dw = ts_adj._compute_terminal_state_gradient("cl")
    dphi_local_norm = np.linalg.norm(dphi_dw)
    dphi_global_norm_sq = comm.allreduce(np.dot(dphi_dw, dphi_dw))
    dphi_global_norm = np.sqrt(dphi_global_norm_sq)

    # Print per-rank diagnostics
    for r in range(size):
        if rank == r:
            print(f"  [rank {r}] ||dCL/dw|| = {dphi_local_norm:.6e}  "
                  f"dCL/dw[:5] = {dphi_dw[:5]}", flush=True)
        comm.Barrier()
    if rank == 0:
        print(f"  ||dCL/dw|| (global) = {dphi_global_norm:.6e}")

    # Adjoint
    result = ts_adj.solve_adjoint(["cl"])
    lam = result["cl"]

    # Per-rank lambda diagnostics
    lam_local_norm = np.linalg.norm(lam)
    lam_global_norm_sq = comm.allreduce(np.dot(lam, lam))
    lam_global_norm = np.sqrt(lam_global_norm_sq)

    for r in range(size):
        if rank == r:
            print(f"  [rank {r}] ||lambda(t0)|| = {lam_local_norm:.6e}  "
                  f"lambda[:5] = {lam[:5]}", flush=True)
        comm.Barrier()
    if rank == 0:
        print(f"  ||lambda(t0)|| (global) = {lam_global_norm:.6e}")
        # Ratio check: if adjoint is good, lambda should be similar magnitude to dCL/dw
        ratio = lam_global_norm / max(dphi_global_norm, 1e-30)
        print(f"  ||lambda|| / ||dCL/dw|| = {ratio:.6e}  "
              f"(expect ~1 for small dt*n_steps)")

    # --- Step 3: Random perturbation ---
    np.random.seed(12345 + rank)
    v = np.random.randn(len(w0_base))
    v_norm_sq = comm.allreduce(np.dot(v, v))
    v /= np.sqrt(v_norm_sq)

    # Verify normalization
    v_check = comm.allreduce(np.dot(v, v))
    adj_dot = comm.allreduce(np.dot(lam, v))
    dphi_dot = comm.allreduce(np.dot(dphi_dw, v))

    if rank == 0:
        print(f"\n  ||v|| (global) = {np.sqrt(v_check):.10e}  (should be 1)")
        print(f"  lambda(t0)^T * v = {adj_dot:.10e}")
        print(f"  dCL/dw^T * v     = {dphi_dot:.10e}  (terminal cond)")

    # --- Step 4: FD sweep ---
    if rank == 0:
        print(f"\n  {'eps':>12s} | {'FD (central)':>18s} | {'Adjoint':>18s} | "
              f"{'rel error':>12s} | {'CL+':>14s} | {'CL-':>14s}")
        print(f"  {'-'*12}-+-{'-'*18}-+-{'-'*18}-+-{'-'*12}-+-{'-'*14}-+-{'-'*14}")

    for eps in eps_list:
        w_plus = w0_base + eps * v
        w_minus = w0_base - eps * v

        cl_plus = forward_solve(w_plus, dt, n_steps)
        cl_minus = forward_solve(w_minus, dt, n_steps)

        fd_dot = (cl_plus - cl_minus) / (2.0 * eps)
        rel_err = abs(fd_dot - adj_dot) / max(abs(adj_dot), abs(fd_dot), 1e-30)

        if rank == 0:
            print(f"  {eps:12.1e} | {fd_dot:18.10e} | {adj_dot:18.10e} | "
                  f"{rel_err:12.4e} | {cl_plus:14.10e} | {cl_minus:14.10e}")

    if rank == 0:
        print()
    return adj_dot, dphi_dot


# =================================================================
# Run tests
# =================================================================

# Test 1: Well-conditioned case (dt=1e-8, a=1e8 >> spectral radius ~3.5e6)
# This validates the FD infrastructure.
adj1, dphi1 = run_fd_test(dt=1e-8, n_steps=3, ksp_rtol=1e-6,
                          eps_list=[1e-3, 1e-5, 1e-7])

# Test 2: Moderate conditioning (dt=1e-5, a=1e5 < spectral radius)
adj2, dphi2 = run_fd_test(dt=1e-5, n_steps=3, ksp_rtol=1e-6,
                          eps_list=[1e-3, 1e-5, 1e-7])

if rank == 0:
    print("=" * 70)
    print("  SUMMARY")
    print("=" * 70)
    print(f"  dt=1e-8 (well-conditioned): adj_dot = {adj1:.10e}, dphi_dot = {dphi1:.10e}")
    print(f"  dt=1e-5 (moderate):         adj_dot = {adj2:.10e}, dphi_dot = {dphi2:.10e}")
    print("=" * 70)
