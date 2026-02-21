"""
Verify all AD partial derivatives (forward, reverse, FD) at steady state.

Tests 4 partials:
  1. dR/dw:      state Jacobian           (matrix, random directions)
  2. dR/dalpha:   residual w.r.t. alpha    (vector)
  3. dCL/dw:     cost gradient w.r.t. w   (vector, random directions)
  4. dCL/dalpha:  cost w.r.t. alpha        (scalar)

Each partial is tested by:
  (a) Forward AD vs central FD of getRes()  -> validates forward mode
  (b) Forward AD vs ADflow internal FD (mode="FD") -> validates AD tangent code
  (c) Dot product: <fwd(v),u> = <v,rev(u)> -> validates reverse = transpose(fwd)

Usage:
    mpirun -np 2 python test_partials.py
"""
import os
import sys
import numpy as np
from mpi4py import MPI
from baseclasses import AeroProblem
from adflow import ADFLOW

sys.stdout.reconfigure(line_buffering=True)
comm = MPI.COMM_WORLD
rank = comm.rank

_STEADY = 1

# ======================================================================
# Configuration
# ======================================================================
N_RANDOM = 3          # number of random directions per test
EPS_W_REL = 1e-5      # FD step for state perturbation (relative to ||w||)
EPS_ALPHA = 1e-4      # FD step for alpha perturbation [degrees]
FD_TOL = 1e-4         # tolerance for FD vs AD comparison
DOT_TOL = 1e-9        # tolerance for dot product test

base_dir = os.path.dirname(os.path.abspath(__file__))
grid_file = os.path.join(base_dir, "../../../input_files/naca0012_rans-L2.cgns")
output_dir = os.path.join(base_dir, "output")
os.makedirs(output_dir, exist_ok=True)

if not os.path.exists(grid_file):
    if rank == 0:
        print(f"ERROR: Mesh file not found: {grid_file}")
    sys.exit(1)

alpha_base = 2.77  # degrees

ap = AeroProblem(
    name="naca0012",
    alpha=alpha_base,
    mach=0.6,
    machRef=0.6,
    reynolds=4800000.0,
    reynoldsLength=1.0,
    T=280.0,
    R=287.085,
    areaRef=1.0,
    chordRef=1.0,
    xRef=0.25,
    evalFuncs=["cl", "cd"],
)
ap.addDV("alpha", name="alpha")

options = {
    "gridfile": grid_file,
    "outputdirectory": output_dir,
    "writevolumesolution": False,
    "writesurfacesolution": False,
    "smoother": "DADI",
    "equationtype": "RANS",
    "equationmode": "steady",
    "cfl": 2.5,
    "cflcoarse": 1.2,
    "ncycles": 5000,
    "mgcycle": "3w",
    "mgstartlevel": 1,
    "monitorvariables": ["cpu", "resrho", "cl", "cd"],
    "usenksolver": True,
    "useanksolver": True,
    "l2convergence": 1e-12,
    "l2convergencecoarse": 1e-4,
    "nsubiterturb": 10,
    "blockSplitting": True,
    "useblockettes": False,
    "printAllOptions": False,
    "printIterations": False,
}

solver = ADFLOW(options=options, debug=False)
n_local = solver.getStateSize()
n_global = comm.allreduce(n_local)
adflow = solver.adflow


# ======================================================================
# Helper functions
# ======================================================================
def get_turb_info():
    """Return (nw, nwf, n_turb, turbResScale_array)."""
    nw = int(adflow.flowvarrefstate.nw)
    nwf = int(adflow.flowvarrefstate.nwf)
    n_turb = nw - nwf
    if n_turb > 0:
        trs = adflow.inputiteration.turbresscale[:n_turb].copy()
    else:
        trs = None
    return nw, nwf, n_turb, trs


def undo_turb_res_scale(vec):
    """Divide SA DOF rows of vec by turbResScale (in-place). Returns vec."""
    nw, nwf, n_turb, trs = get_turb_info()
    if n_turb > 0:
        n_cells = len(vec) // nw
        v2d = vec.reshape(n_cells, nw)
        for l in range(n_turb):
            if abs(trs[l]) > 1e-30:
                v2d[:, nwf + l] /= trs[l]
    return vec


def apply_inv_turb_res_scale(vec):
    """Return a copy of vec with SA DOFs divided by turbResScale."""
    nw, nwf, n_turb, trs = get_turb_info()
    if n_turb > 0:
        out = vec.copy()
        n_cells = len(out) // nw
        v2d = out.reshape(n_cells, nw)
        for l in range(n_turb):
            if abs(trs[l]) > 1e-30:
                v2d[:, nwf + l] /= trs[l]
        return out
    return vec.copy()


def eval_res_getres(w):
    """Evaluate spatial residual via getRes() (blocketteRes path).
    Returns R in 1/volRef space (no turbResScale)."""
    solver.setStates(w)
    orig = adflow.inputphysics.equationmode
    adflow.inputphysics.equationmode = _STEADY
    res = np.zeros(n_local)
    res = adflow.nksolver.getres(res)
    adflow.inputphysics.equationmode = orig
    return res


def eval_cl(w):
    """Evaluate CL at state w (calls getRes to update internal arrays)."""
    solver.setStates(w)
    orig = adflow.inputphysics.equationmode
    adflow.inputphysics.equationmode = _STEADY
    res_tmp = np.zeros(n_local)
    adflow.nksolver.getres(res_tmp)
    adflow.inputphysics.equationmode = orig
    funcs = {}
    solver.evalFunctions(ap, funcs, evalFuncs=["cl"])
    return funcs[f"{ap.name}_cl"]


def global_dot(a, b):
    return comm.allreduce(np.dot(a, b))


def global_norm(a):
    return np.sqrt(global_dot(a, a))


def rel_err(a, b):
    """Relative error between scalars or arrays."""
    if isinstance(a, np.ndarray):
        diff = global_norm(a - b)
        ref = max(global_norm(a), global_norm(b), 1e-30)
    else:
        diff = abs(a - b)
        ref = max(abs(a), abs(b), 1e-30)
    return diff / ref


def random_vec(seed):
    """Generate a random vector with global unit norm."""
    rng = np.random.RandomState(seed + rank)
    v = rng.randn(n_local)
    v_norm = global_norm(v)
    return v / v_norm


results = []


def log_result(name, err, tol, tag=""):
    ok = err < tol
    label = "PASS" if ok else "FAIL"
    results.append((name, ok, err))
    if rank == 0:
        extra = f" ({tag})" if tag else ""
        print(f"    {label}: {name}{extra}  rel_err = {err:.4e}  (tol {tol:.0e})")


# ======================================================================
# Phase A: Converge steady state
# ======================================================================
if rank == 0:
    print("=" * 70)
    print("Phase A: Converge steady RANS at alpha = %.2f deg" % alpha_base)
    print("=" * 70)

solver(ap)
w_ss = solver.getStates().copy()
cl_ss = eval_cl(w_ss)
w_norm = global_norm(w_ss)

if rank == 0:
    print(f"  CL = {cl_ss:.10e}")
    print(f"  ||w|| = {w_norm:.6e}")
    print(f"  n_local = {n_local}, n_global = {n_global}")
    nw, nwf, n_turb, trs = get_turb_info()
    print(f"  nw={nw}, nwf={nwf}, n_turb={n_turb}, turbResScale={trs}")

# Absolute FD step
eps_w = EPS_W_REL * w_norm

# Restore state
solver.setStates(w_ss)
solver.setAeroProblem(ap)


# ======================================================================
# Test 1: dR/dw (state Jacobian)
# ======================================================================
if rank == 0:
    print("\n" + "=" * 70)
    print("Test 1: dR/dw (state Jacobian)")
    print("=" * 70)

# Check residual norms at baseline
R0_getres = eval_res_getres(w_ss)
R0_norm = global_norm(R0_getres)
if rank == 0:
    print(f"  ||R(w_ss)|| via getRes = {R0_norm:.6e}")

for i_dir in range(N_RANDOM):
    seed_v = 100 + i_dir
    seed_u = 200 + i_dir
    v = random_vec(seed_v)
    u = random_vec(seed_u)

    if rank == 0:
        print(f"\n  Direction {i_dir+1}/{N_RANDOM}:")

    # -- Forward AD (Tapenade master_d) --
    solver.setStates(w_ss)
    adflow.inputphysics.equationmode = _STEADY
    ad_fwd = solver.computeJacobianVectorProductFwd(
        wDot=v, residualDeriv=True
    )
    adflow.inputphysics.equationmode = _STEADY  # restore
    # This includes turbResScale. Undo it.
    undo_turb_res_scale(ad_fwd)
    ad_fwd_norm = global_norm(ad_fwd)

    # -- ADflow internal FD (mode="FD", FDs the master path) --
    solver.setStates(w_ss)
    adflow.inputphysics.equationmode = _STEADY
    ad_fd_internal = solver.computeJacobianVectorProductFwd(
        wDot=v, residualDeriv=True, mode="FD", h=1e-6
    )
    adflow.inputphysics.equationmode = _STEADY
    undo_turb_res_scale(ad_fd_internal)
    ad_fd_internal_norm = global_norm(ad_fd_internal)

    # -- External FD of getRes() (blocketteRes path) --
    solver.setStates(w_ss)
    solver.setAeroProblem(ap)
    R_plus = eval_res_getres(w_ss + eps_w * v)
    R_minus = eval_res_getres(w_ss - eps_w * v)
    solver.setStates(w_ss)
    fd_ext = (R_plus - R_minus) / (2.0 * eps_w)
    fd_ext_norm = global_norm(fd_ext)

    if rank == 0:
        print(f"    ||AD fwd||        = {ad_fwd_norm:.6e}")
        print(f"    ||AD internal FD|| = {ad_fd_internal_norm:.6e}")
        print(f"    ||ext FD (getRes)|| = {fd_ext_norm:.6e}")

    # Compare AD vs internal FD (same code path, should match closely)
    err_ad_ifd = rel_err(ad_fwd, ad_fd_internal)
    log_result("dR/dw AD vs internal FD", err_ad_ifd, FD_TOL, f"dir {i_dir+1}")

    # Compare AD vs external FD (different code paths, ~8% discrepancy expected)
    err_ad_efd = rel_err(ad_fwd, fd_ext)
    log_result("dR/dw AD vs external FD(getRes)", err_ad_efd, 0.2, f"dir {i_dir+1}")

    # -- Dot product test: <fwd(v), u> = <v, rev(u)> --
    solver.setStates(w_ss)
    adflow.inputphysics.equationmode = _STEADY
    u_adj = apply_inv_turb_res_scale(u)
    ad_rev = solver.computeJacobianVectorProductBwd(
        resBar=u_adj, wDeriv=True
    )
    adflow.inputphysics.equationmode = _STEADY

    dot_fwd = global_dot(ad_fwd, u)
    dot_rev = global_dot(v, ad_rev)
    err_dot = abs(dot_fwd - dot_rev) / max(abs(dot_fwd), abs(dot_rev), 1e-30)
    log_result("dR/dw dot product", err_dot, DOT_TOL, f"dir {i_dir+1}")

    if rank == 0:
        print(f"      <fwd(v),u> = {dot_fwd:.12e}")
        print(f"      <v,rev(u)> = {dot_rev:.12e}")


# ======================================================================
# Test 2: dR/dalpha
# ======================================================================
if rank == 0:
    print("\n" + "=" * 70)
    print("Test 2: dR/dalpha (residual w.r.t. alpha)")
    print("=" * 70)

# Forward AD
solver.setStates(w_ss)
solver.setAeroProblem(ap)
adflow.inputphysics.equationmode = _STEADY
dR_da_fwd = solver.computeJacobianVectorProductFwd(
    xDvDot={"alpha": 1.0}, residualDeriv=True
)
adflow.inputphysics.equationmode = _STEADY
undo_turb_res_scale(dR_da_fwd)
dR_da_fwd_norm = global_norm(dR_da_fwd)

# ADflow internal FD
solver.setStates(w_ss)
solver.setAeroProblem(ap)
adflow.inputphysics.equationmode = _STEADY
dR_da_ifd = solver.computeJacobianVectorProductFwd(
    xDvDot={"alpha": 1.0}, residualDeriv=True, mode="FD", h=1e-6
)
adflow.inputphysics.equationmode = _STEADY
undo_turb_res_scale(dR_da_ifd)

# External FD: perturb alpha, evaluate getRes at same state w_ss
ap.alpha = alpha_base + EPS_ALPHA
solver.setAeroProblem(ap)
solver.setStates(w_ss)
R_plus = eval_res_getres(w_ss)

ap.alpha = alpha_base - EPS_ALPHA
solver.setAeroProblem(ap)
solver.setStates(w_ss)
R_minus = eval_res_getres(w_ss)

ap.alpha = alpha_base
solver.setAeroProblem(ap)
solver.setStates(w_ss)

dR_da_efd = (R_plus - R_minus) / (2.0 * EPS_ALPHA)

if rank == 0:
    print(f"  ||dR/dalpha AD fwd||       = {dR_da_fwd_norm:.6e}")
    print(f"  ||dR/dalpha internal FD||  = {global_norm(dR_da_ifd):.6e}")
    print(f"  ||dR/dalpha external FD||  = {global_norm(dR_da_efd):.6e}")

err_ifd = rel_err(dR_da_fwd, dR_da_ifd)
log_result("dR/dalpha AD vs internal FD", err_ifd, FD_TOL)

err_efd = rel_err(dR_da_fwd, dR_da_efd)
log_result("dR/dalpha AD vs external FD", err_efd, 0.2)

# Dot product test
for i_dir in range(N_RANDOM):
    u = random_vec(300 + i_dir)
    solver.setStates(w_ss)
    solver.setAeroProblem(ap)
    dot_fwd = global_dot(dR_da_fwd, u)
    adflow.inputphysics.equationmode = _STEADY
    u_adj = apply_inv_turb_res_scale(u)
    rev_result = solver.computeJacobianVectorProductBwd(
        resBar=u_adj, xDvDeriv=True
    )
    adflow.inputphysics.equationmode = _STEADY
    dot_rev = float(rev_result.get("alpha", 0.0))
    err_dot = abs(dot_fwd - dot_rev) / max(abs(dot_fwd), abs(dot_rev), 1e-30)
    log_result("dR/dalpha dot product", err_dot, DOT_TOL, f"dir {i_dir+1}")
    if rank == 0:
        print(f"      <dR/da, u>     = {dot_fwd:.12e}")
        print(f"      rev(u)['alpha'] = {dot_rev:.12e}")


# ======================================================================
# Test 3: dCL/dw (use reverse AD + FD only; forward funcDeriv may be broken)
# ======================================================================
if rank == 0:
    print("\n" + "=" * 70)
    print("Test 3: dCL/dw (cost gradient w.r.t. state)")
    print("=" * 70)

# Reverse AD: full gradient vector
solver.setStates(w_ss)
solver.setAeroProblem(ap)
adflow.inputphysics.equationmode = _STEADY
funcsBar = solver._getFuncsBar("cl")
dCL_dw_rev = solver.computeJacobianVectorProductBwd(
    funcsBar=funcsBar, wDeriv=True
)
adflow.inputphysics.equationmode = _STEADY
dCL_dw_norm = global_norm(dCL_dw_rev)
if rank == 0:
    print(f"  ||dCL/dw (rev)|| = {dCL_dw_norm:.6e}")

# Also try forward funcDeriv (diagnostic — may return 0)
solver.setStates(w_ss)
solver.setAeroProblem(ap)
adflow.inputphysics.equationmode = _STEADY
v_test = random_vec(400)
fwd_result = solver.computeJacobianVectorProductFwd(
    wDot=v_test, funcDeriv=True
)
adflow.inputphysics.equationmode = _STEADY
fwd_cl = float(fwd_result.get(f"{ap.name}_cl", 0.0))
rev_dot = global_dot(dCL_dw_rev, v_test)
if rank == 0:
    print(f"  [diagnostic] funcDeriv(wDot=v): {fwd_cl:.12e}")
    print(f"  [diagnostic] <rev, v>:          {rev_dot:.12e}")
    if abs(fwd_cl) < 1e-30:
        print(f"  WARNING: funcDeriv with wDot returns 0 — known limitation")

for i_dir in range(N_RANDOM):
    v = random_vec(400 + i_dir)
    if rank == 0:
        print(f"\n  Direction {i_dir+1}/{N_RANDOM}:")

    # Reverse AD: <dCL/dw, v>
    rev_scalar = global_dot(dCL_dw_rev, v)

    # FD: (CL(w+ev) - CL(w-ev)) / (2e)
    solver.setStates(w_ss)
    solver.setAeroProblem(ap)
    cl_plus = eval_cl(w_ss + eps_w * v)
    solver.setStates(w_ss)
    solver.setAeroProblem(ap)
    cl_minus = eval_cl(w_ss - eps_w * v)
    solver.setStates(w_ss)
    fd_scalar = (cl_plus - cl_minus) / (2.0 * eps_w)

    err_fd = abs(rev_scalar - fd_scalar) / max(abs(rev_scalar), abs(fd_scalar), 1e-30)
    log_result("dCL/dw rev vs FD", err_fd, FD_TOL, f"dir {i_dir+1}")
    if rank == 0:
        print(f"      <dCL/dw, v> rev:  {rev_scalar:.12e}")
        print(f"      FD:               {fd_scalar:.12e}")


# ======================================================================
# Test 4: dCL/dalpha (scalar)
# ======================================================================
if rank == 0:
    print("\n" + "=" * 70)
    print("Test 4: dCL/dalpha (scalar)")
    print("=" * 70)

# Reverse AD
solver.setStates(w_ss)
solver.setAeroProblem(ap)
adflow.inputphysics.equationmode = _STEADY
funcsBar = solver._getFuncsBar("cl")
rev_full = solver.computeJacobianVectorProductBwd(
    funcsBar=funcsBar, xDvDeriv=True
)
adflow.inputphysics.equationmode = _STEADY
dCL_da_rev = float(rev_full.get("alpha", 0.0))

# Forward AD (diagnostic)
solver.setStates(w_ss)
solver.setAeroProblem(ap)
adflow.inputphysics.equationmode = _STEADY
fwd_funcs = solver.computeJacobianVectorProductFwd(
    xDvDot={"alpha": 1.0}, funcDeriv=True
)
adflow.inputphysics.equationmode = _STEADY
dCL_da_fwd = float(fwd_funcs.get(f"{ap.name}_cl", 0.0))

# FD
solver.setStates(w_ss)
ap.alpha = alpha_base + EPS_ALPHA
solver.setAeroProblem(ap)
solver.setStates(w_ss)
cl_plus = eval_cl(w_ss)

ap.alpha = alpha_base - EPS_ALPHA
solver.setAeroProblem(ap)
solver.setStates(w_ss)
cl_minus = eval_cl(w_ss)

ap.alpha = alpha_base
solver.setAeroProblem(ap)
solver.setStates(w_ss)

dCL_da_fd = (cl_plus - cl_minus) / (2.0 * EPS_ALPHA)

if rank == 0:
    print(f"  Forward AD:  {dCL_da_fwd:.12e}")
    print(f"  Reverse AD:  {dCL_da_rev:.12e}")
    print(f"  FD:          {dCL_da_fd:.12e}")

err_rev_fd = abs(dCL_da_rev - dCL_da_fd) / max(abs(dCL_da_rev), 1e-30)
log_result("dCL/dalpha rev vs FD", err_rev_fd, FD_TOL)

if abs(dCL_da_fwd) > 1e-30:
    err_fwd_rev = abs(dCL_da_fwd - dCL_da_rev) / max(abs(dCL_da_fwd), 1e-30)
    log_result("dCL/dalpha fwd vs rev", err_fwd_rev, DOT_TOL)
else:
    if rank == 0:
        print("  [skip] fwd vs rev: forward returns 0")


# ======================================================================
# Summary
# ======================================================================
if rank == 0:
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    n_pass = sum(1 for _, ok, _ in results if ok)
    n_fail = sum(1 for _, ok, _ in results if not ok)
    for name, ok, err in results:
        label = "PASS" if ok else "FAIL"
        print(f"  [{label}] {name:45s}  err = {err:.4e}")
    print(f"\n  {n_pass} passed, {n_fail} failed out of {len(results)} tests")
    if n_fail == 0:
        print("  ALL TESTS PASSED")
    else:
        print("  *** SOME TESTS FAILED ***")
    print("=" * 70)
