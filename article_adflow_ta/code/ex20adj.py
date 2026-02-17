"""
Performs adjoint sensitivity analysis for the van der Pol equation using petsc4py,
with finite difference verification.

This program solves the van der Pol ODE:
    u_1' = u_2
    u_2' = mu * ((1 - u_1^2) * u_2 - u_1)

on 0 <= t <= 0.5, with initial conditions:
    u_1(0) = 2
    u_2(0) = -2/3 + 10/(81*mu) - 292/(2187*mu^2)

and computes the sensitivities of the final solution w.r.t. initial conditions
and parameter mu using discrete adjoint via TSAdjoint.

Python translation of PETSc src/ts/tutorials/ex20adj.c
"""

import sys
import petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import numpy as np


# -- Global parameter --
mu = 1.0e3


def ifunction(ts, t, U, Udot, F):
    """Implicit residual: F(t, U, Udot) = Udot - f(U)"""
    u = U.getArray(readonly=True)
    udot = Udot.getArray(readonly=True)
    f = F.getArray()
    f[0] = udot[0] - u[1]
    f[1] = udot[1] - mu * ((1.0 - u[0]**2) * u[1] - u[0])


def ijacobian(ts, t, U, Udot, a, J, P):
    """Jacobian: a * dF/dUdot + dF/dU"""
    u = U.getArray(readonly=True)
    P.zeroEntries()
    P.setValue(0, 0, a)
    P.setValue(0, 1, -1.0)
    P.setValue(1, 0, mu * (2.0 * u[0] * u[1] + 1.0))
    P.setValue(1, 1, a - mu * (1.0 - u[0]**2))
    P.assemble()
    if J != P:
        J.assemble()
    return True  # same nonzero structure


def ijacobianp(ts, t, U, Udot, a, Jp):
    """Jacobian of F w.r.t. parameter mu: dF/dp (note: F = Udot - f, so dF/dp = -df/dp)"""
    u = U.getArray(readonly=True)
    Jp.zeroEntries()
    Jp.setValue(0, 0, 0.0)
    Jp.setValue(1, 0, -((1.0 - u[0]**2) * u[1] - u[0]))
    Jp.assemble()


def forward_solve(y0, z0, mu_val):
    """Run forward solve with given initial conditions and mu. Returns (y_tf, z_tf)."""
    global mu
    mu = mu_val

    A = PETSc.Mat().createDense([2, 2], comm=PETSc.COMM_WORLD)
    A.setUp()
    Jacp = PETSc.Mat().createDense([2, 1], comm=PETSc.COMM_WORLD)
    Jacp.setUp()
    U = PETSc.Vec().createSeq(2, comm=PETSc.COMM_WORLD)

    u0 = U.getArray()
    u0[0] = y0
    u0[1] = z0

    ts = PETSc.TS().create(comm=PETSc.COMM_WORLD)
    ts.setEquationType(PETSc.TS.EquationType.ODE_EXPLICIT)
    ts.setType(PETSc.TS.Type.CN)
    ts.setIFunction(ifunction, U.duplicate())
    ts.setIJacobian(ijacobian, A, A)
    ts.setIJacobianP(ijacobianp, Jacp)
    ts.setTime(0.0)
    ts.setMaxTime(0.5)
    ts.setTimeStep(0.001)
    ts.setExactFinalTime(PETSc.TS.ExactFinalTime.MATCHSTEP)
    ts.setFromOptions()

    ts.solve(U)
    result = U.getArray(readonly=True).copy()

    A.destroy()
    Jacp.destroy()
    U.destroy()
    ts.destroy()
    return result[0], result[1]


# ============================================================
# Baseline initial conditions
# ============================================================
mu_base = 1.0e3
y0_base = 2.0
z0_base = -2.0 / 3.0 + 10.0 / (81.0 * mu_base) - 292.0 / (2187.0 * mu_base * mu_base)

# ============================================================
# Adjoint solve
# ============================================================
mu = mu_base

A = PETSc.Mat().createDense([2, 2], comm=PETSc.COMM_WORLD)
A.setUp()
Jacp = PETSc.Mat().createDense([2, 1], comm=PETSc.COMM_WORLD)
Jacp.setUp()
U = PETSc.Vec().createSeq(2, comm=PETSc.COMM_WORLD)

u0 = U.getArray()
u0[0] = y0_base
u0[1] = z0_base

ts = PETSc.TS().create(comm=PETSc.COMM_WORLD)
ts.setEquationType(PETSc.TS.EquationType.ODE_EXPLICIT)
ts.setType(PETSc.TS.Type.CN)
ts.setIFunction(ifunction, U.duplicate())
ts.setIJacobian(ijacobian, A, A)
ts.setIJacobianP(ijacobianp, Jacp)
ts.setTime(0.0)
ts.setMaxTime(0.5)
ts.setTimeStep(0.001)
ts.setExactFinalTime(PETSc.TS.ExactFinalTime.MATCHSTEP)
ts.setSaveTrajectory()
ts.setFromOptions()

ts.solve(U)
ftime = ts.getSolveTime()
steps = ts.getStepNumber()
y_tf = U.getArray()[0]
z_tf = U.getArray()[1]
PETSc.Sys.Print(f"\nForward solve: ftime = {ftime}, steps = {steps}")
PETSc.Sys.Print(f"Solution at final time: y(tf) = {y_tf}, z(tf) = {z_tf}")

# Adjoint initial conditions
lambda0 = U.duplicate()
lambda1 = U.duplicate()
lam0 = lambda0.getArray()
lam0[0] = 1.0; lam0[1] = 0.0
lam1 = lambda1.getArray()
lam1[0] = 0.0; lam1[1] = 1.0

mup0 = PETSc.Vec().createSeq(1, comm=PETSc.COMM_WORLD)
mup0.set(0.0)
mup1 = PETSc.Vec().createSeq(1, comm=PETSc.COMM_WORLD)
mup1.set(0.0)

ts.setCostGradients([lambda0, lambda1], [mup0, mup1])
ts.adjointSolve()

# Extract adjoint results
adj_dy_dy0 = lambda0.getArray(readonly=True)[0]
adj_dy_dz0 = lambda0.getArray(readonly=True)[1]
adj_dz_dy0 = lambda1.getArray(readonly=True)[0]
adj_dz_dz0 = lambda1.getArray(readonly=True)[1]

# Total sensitivity w.r.t. mu (chain rule: z0 depends on mu)
dz0_dmu = -10.0 / (81.0 * mu_base**2) + 2.0 * 292.0 / (2187.0 * mu_base**3)
adj_dy_dmu = adj_dy_dz0 * dz0_dmu + mup0.getArray(readonly=True)[0]
adj_dz_dmu = adj_dz_dz0 * dz0_dmu + mup1.getArray(readonly=True)[0]

PETSc.Sys.Print("\n===== Adjoint sensitivities =====")
PETSc.Sys.Print(f"  d[y(tf)]/d[y0] = {adj_dy_dy0}")
PETSc.Sys.Print(f"  d[y(tf)]/d[z0] = {adj_dy_dz0}")
PETSc.Sys.Print(f"  d[z(tf)]/d[y0] = {adj_dz_dy0}")
PETSc.Sys.Print(f"  d[z(tf)]/d[z0] = {adj_dz_dz0}")
PETSc.Sys.Print(f"  d[y(tf)]/d[mu] = {adj_dy_dmu}")
PETSc.Sys.Print(f"  d[z(tf)]/d[mu] = {adj_dz_dmu}")

# Cleanup adjoint objects
A.destroy(); Jacp.destroy(); U.destroy()
lambda0.destroy(); lambda1.destroy()
mup0.destroy(); mup1.destroy(); ts.destroy()

# ============================================================
# Finite difference verification
# ============================================================
PETSc.Sys.Print("\n===== Finite difference verification =====")
h = 1.0e-7

# d/dy0
yp, zp = forward_solve(y0_base + h, z0_base, mu_base)
ym, zm = forward_solve(y0_base - h, z0_base, mu_base)
fd_dy_dy0 = (yp - ym) / (2.0 * h)
fd_dz_dy0 = (zp - zm) / (2.0 * h)

# d/dz0
yp, zp = forward_solve(y0_base, z0_base + h, mu_base)
ym, zm = forward_solve(y0_base, z0_base - h, mu_base)
fd_dy_dz0 = (yp - ym) / (2.0 * h)
fd_dz_dz0 = (zp - zm) / (2.0 * h)

# d/dmu (total: perturb mu and also update z0 accordingly)
h_mu = 1.0e-3  # larger step for mu since mu ~ 1e3
mu_p = mu_base + h_mu
z0_p = -2.0 / 3.0 + 10.0 / (81.0 * mu_p) - 292.0 / (2187.0 * mu_p * mu_p)
yp, zp = forward_solve(y0_base, z0_p, mu_p)

mu_m = mu_base - h_mu
z0_m = -2.0 / 3.0 + 10.0 / (81.0 * mu_m) - 292.0 / (2187.0 * mu_m * mu_m)
ym, zm = forward_solve(y0_base, z0_m, mu_m)

fd_dy_dmu = (yp - ym) / (2.0 * h_mu)
fd_dz_dmu = (zp - zm) / (2.0 * h_mu)

PETSc.Sys.Print(f"\n{'Sensitivity':<25s} {'Adjoint':>15s} {'FD':>15s} {'Rel Error':>15s}")
PETSc.Sys.Print(f"{'-'*70}")

for name, adj_val, fd_val in [
    ("d[y(tf)]/d[y0]", adj_dy_dy0, fd_dy_dy0),
    ("d[y(tf)]/d[z0]", adj_dy_dz0, fd_dy_dz0),
    ("d[z(tf)]/d[y0]", adj_dz_dy0, fd_dz_dy0),
    ("d[z(tf)]/d[z0]", adj_dz_dz0, fd_dz_dz0),
    ("d[y(tf)]/d[mu]",  adj_dy_dmu,  fd_dy_dmu),
    ("d[z(tf)]/d[mu]",  adj_dz_dmu,  fd_dz_dmu),
]:
    if abs(fd_val) > 1e-30:
        rel_err = abs(adj_val - fd_val) / abs(fd_val)
    else:
        rel_err = abs(adj_val - fd_val)
    PETSc.Sys.Print(f"{name:<25s} {adj_val:>15.8e} {fd_val:>15.8e} {rel_err:>15.8e}")
