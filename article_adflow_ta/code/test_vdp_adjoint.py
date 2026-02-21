"""
Van der Pol adjoint verification with PETSc TS.

Tests TSTHETA with theta=1.0 (backward Euler) and theta=0.5 (Crank-Nicolson).
Verifies adjoint sensitivities d[u(tf)]/d[u0] against central finite differences.

This is the foundation test: if adjoint doesn't match FD here, nothing downstream
will work.

Usage:
    python test_vdp_adjoint.py
"""

import sys
import numpy as np

import petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc

# -- Problem: Van der Pol ODE --
#   u1' = u2
#   u2' = mu * ((1 - u1^2) * u2 - u1)
mu = 1.0e3
t_final = 0.5
dt = 0.001

y0 = 2.0
z0 = -2.0 / 3.0 + 10.0 / (81.0 * mu) - 292.0 / (2187.0 * mu**2)


def ifunction(ts, t, U, Udot, F):
    """F(t, U, Udot) = Udot - f(U)"""
    u = U.getArray(readonly=True)
    udot = Udot.getArray(readonly=True)
    f = F.getArray()
    f[0] = udot[0] - u[1]
    f[1] = udot[1] - mu * ((1.0 - u[0] ** 2) * u[1] - u[0])


def ijacobian(ts, t, U, Udot, a, J, P):
    """J = a * I - df/dU"""
    u = U.getArray(readonly=True)
    P.zeroEntries()
    P.setValue(0, 0, a)
    P.setValue(0, 1, -1.0)
    P.setValue(1, 0, mu * (2.0 * u[0] * u[1] + 1.0))
    P.setValue(1, 1, a - mu * (1.0 - u[0] ** 2))
    P.assemble()
    if J != P:
        J.assemble()
    return True


def ijacobianp_dummy(ts, t, U, Udot, a, Jp):
    """Dummy parameter Jacobian (no parameters). Required by TSAdjointSolve."""
    Jp.zeroEntries()
    Jp.assemble()


def forward_solve(y_init, z_init, theta, save_traj=False):
    """Run forward solve. Returns (y_tf, z_tf, ts_or_None)."""
    A = PETSc.Mat().createDense([2, 2])
    A.setUp()
    U = PETSc.Vec().createSeq(2)
    U.getArray()[:] = [y_init, z_init]

    # Dummy parameter Jacobian (1 dummy parameter, 2 state DOFs)
    Jacp = PETSc.Mat().createDense([2, 1])
    Jacp.setUp()

    ts = PETSc.TS().create()
    ts.setType(PETSc.TS.Type.THETA)
    ts.setTheta(theta)
    ts.setIFunction(ifunction, U.duplicate())
    ts.setIJacobian(ijacobian, A, A)
    ts.setIJacobianP(ijacobianp_dummy, Jacp)
    ts.setTime(0.0)
    ts.setMaxTime(t_final)
    ts.setTimeStep(dt)
    ts.setMaxSteps(int(round(t_final / dt)))
    ts.setExactFinalTime(PETSc.TS.ExactFinalTime.MATCHSTEP)

    opts = PETSc.Options()
    opts.setValue("-ts_adapt_type", "none")
    if save_traj:
        ts.setSaveTrajectory()
    ts.setFromOptions()
    ts.setTimeStep(dt)  # re-set after setFromOptions

    ts.solve(U)
    result = U.getArray(readonly=True).copy()

    if save_traj:
        return result[0], result[1], ts, U, A, Jacp
    else:
        A.destroy()
        Jacp.destroy()
        U.destroy()
        ts.destroy()
        return result[0], result[1]


def adjoint_solve(ts, U):
    """Run adjoint for both components of u(tf). Returns 2x2 sensitivity matrix."""
    # lambda0 computes d[u1(tf)]/d[u0]
    lambda0 = U.duplicate()
    lambda0.getArray()[:] = [1.0, 0.0]
    # lambda1 computes d[u2(tf)]/d[u0]
    lambda1 = U.duplicate()
    lambda1.getArray()[:] = [0.0, 1.0]

    # mu vectors for parameter sensitivities (not used here, but API
    # requires same length as lambda list). Size=1 placeholder.
    mu0 = PETSc.Vec().createSeq(1)
    mu0.set(0.0)
    mu1 = PETSc.Vec().createSeq(1)
    mu1.set(0.0)

    ts.setCostGradients([lambda0, lambda1], [mu0, mu1])
    ts.adjointSolve()

    lam0 = lambda0.getArray(readonly=True).copy()
    lam1 = lambda1.getArray(readonly=True).copy()
    lambda0.destroy()
    lambda1.destroy()
    mu0.destroy()
    mu1.destroy()

    # sens[i][j] = d[u_i(tf)] / d[u_j(0)]
    return np.array([[lam0[0], lam0[1]],
                     [lam1[0], lam1[1]]])


def fd_sensitivities(theta, h=1e-7):
    """Central FD for d[u(tf)]/d[u0]. Returns 2x2 matrix."""
    sens = np.zeros((2, 2))
    for j, (dy, dz) in enumerate([(h, 0), (0, h)]):
        yp, zp = forward_solve(y0 + dy, z0 + dz, theta)
        ym, zm = forward_solve(y0 - dy, z0 - dz, theta)
        sens[0, j] = (yp - ym) / (2 * h)
        sens[1, j] = (zp - zm) / (2 * h)
    return sens


def test_theta(theta, label):
    """Test forward + adjoint for a given theta value."""
    print(f"\n{'=' * 70}")
    print(f"TSTHETA  theta = {theta}  ({label})")
    print(f"  dt = {dt}, t_final = {t_final}, steps = {int(t_final/dt)}")
    print(f"{'=' * 70}")

    # Forward + adjoint
    y_tf, z_tf, ts, U, A, Jacp = forward_solve(y0, z0, theta, save_traj=True)
    print(f"  Forward: u(tf) = [{y_tf:.10e}, {z_tf:.10e}]")
    print(f"  Steps completed: {ts.getStepNumber()}")

    adj = adjoint_solve(ts, U)
    print(f"\n  Adjoint sensitivities:")
    print(f"    d[u1(tf)]/d[u1(0)] = {adj[0,0]:+.10e}")
    print(f"    d[u1(tf)]/d[u2(0)] = {adj[0,1]:+.10e}")
    print(f"    d[u2(tf)]/d[u1(0)] = {adj[1,0]:+.10e}")
    print(f"    d[u2(tf)]/d[u2(0)] = {adj[1,1]:+.10e}")

    # FD
    fd = fd_sensitivities(theta)
    print(f"\n  FD sensitivities:")
    print(f"    d[u1(tf)]/d[u1(0)] = {fd[0,0]:+.10e}")
    print(f"    d[u1(tf)]/d[u2(0)] = {fd[0,1]:+.10e}")
    print(f"    d[u2(tf)]/d[u1(0)] = {fd[1,0]:+.10e}")
    print(f"    d[u2(tf)]/d[u2(0)] = {fd[1,1]:+.10e}")

    # Compare
    print(f"\n  {'Sensitivity':<25s} {'Adjoint':>15s} {'FD':>15s} {'Rel Error':>12s}")
    print(f"  {'-' * 67}")
    all_ok = True
    for i in range(2):
        for j in range(2):
            name = f"d[u{i+1}]/d[u{j+1}(0)]"
            a_val = adj[i, j]
            f_val = fd[i, j]
            if abs(f_val) > 1e-30:
                rel_err = abs(a_val - f_val) / abs(f_val)
            else:
                rel_err = abs(a_val - f_val)
            status = "OK" if rel_err < 1e-4 else "FAIL"
            if rel_err >= 1e-4:
                all_ok = False
            print(f"  {name:<25s} {a_val:>15.8e} {f_val:>15.8e} {rel_err:>12.4e}  {status}")

    A.destroy()
    Jacp.destroy()
    U.destroy()
    ts.destroy()

    return all_ok


if __name__ == "__main__":
    ok_cn = test_theta(0.5, "Crank-Nicolson")
    ok_be = test_theta(1.0, "Backward Euler")

    print(f"\n{'=' * 70}")
    print(f"SUMMARY")
    print(f"  Crank-Nicolson (theta=0.5): {'PASS' if ok_cn else 'FAIL'}")
    print(f"  Backward Euler (theta=1.0): {'PASS' if ok_be else 'FAIL'}")
    print(f"{'=' * 70}")
