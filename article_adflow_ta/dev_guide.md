# Developer's Guide: Time-Accurate Adjoint for ADflow

## Overview

This guide outlines a staged implementation plan for adding time-accurate (TA) adjoint sensitivity analysis to ADflow.
The approach leverages PETSc's `TS` (time stepper) and `TSAdjoint` infrastructure via `petsc4py`, wrapping ADflow's existing spatial residual and Tapenade-generated AD routines as callbacks.

### Architecture Summary

**Current ADflow unsteady pipeline (Python level):**
```
for each time step:
    advanceTimeStepCounter()
    [update mesh if needed]
    solveTimeStep()         # calls Fortran solverUnsteadyStep()
    writeSolution()
```

**Fortran `solverUnsteadyStep()` (src/solver/solvers.F90:317):**
```
shiftSolution()            # shift wOld arrays for BDF history
setCoefTimeIntegrator()    # set coefTime(0:nOldLevels) based on timeAccuracy & nOldSolAvail
solveState()               # RK/DADI -> ANK -> NK cascade to convergence
```

**`solveState()` inner loop (src/solver/solvers.F90:894):**
```
while not converged:
    if totalR > ANK_switchTol * totalR0:   executeMGCycle()     # RK/DADI
    elif totalR > NK_switchTol * totalR0:  ANKStep()            # approx Newton (1st-order Jacobian, CFL-based)
    else:                                  NKStep()             # exact Newton (2nd-order Jacobian, line search)
```

**New TA pipeline (target):**
```
PETSc TS (TSTHETA, theta=0.5 for Crank-Nicolson)
    IFunction(t, U, Udot, F) = Udot - R(U, t)     # wraps ADflow spatial residual
    IJacobian(t, U, Udot, a, J) = aI - dR/dU       # MatShell wrapping AD products
    [Optional: SNES uses ANK/NK as preconditioner]
Forward solve -> TSSetSaveTrajectory -> TSSolve
Adjoint solve -> TSSetCostGradients -> TSAdjointSolve
```

---

## Stage 1: Forward Unsteady Solver via PETSc TS

### Goal
Replace ADflow's native BDF time loop with PETSc `TSTHETA` (Crank--Nicolson, theta = 0.5), driven from Python via `petsc4py`.
Verify that the PETSc-driven forward solve reproduces ADflow's native unsteady results.

### 1.1 Expose the Spatial Residual

ADflow's residual evaluation is currently embedded inside `solveState()`.
We need a standalone callable that computes `R(w)` (spatial residual only, no temporal terms) given a state vector `w`.

**Key routines:**
- `computeResidualNK()` in `src/NKSolver/NKSolvers.F90:1084` -- calls `blocketteRes()` which evaluates the full spatial residual, boundary conditions, and halo exchanges.
  This is the same residual used by ANK and NK solvers.
- The residual is stored in `dw(i,j,k,l)` in block pointers after the call.

**What to implement:**
1. A new Fortran subroutine `computeSpatialResidual(w_vec, R_vec, n)` that:
   - Takes a flat state vector `w_vec` (same layout as NK solver's `wVec`)
   - Scatters it into ADflow's block data structures (`w` arrays)
   - Calls `computeResidualNK()` (which calls `blocketteRes()`)
   - Gathers the residual `dw` into a flat output vector `R_vec`
   - Does NOT include temporal terms (no `coefTime`, no `wOld` contributions)

2. The existing `setwVec` / `getwVec` patterns in `NKSolvers.F90` already do scatter/gather.
   Reuse or wrap these.

**Important:** ADflow's `dw` stores `R(w)` with sign convention `dw = -R` (the negative of the spatial residual, since `dw/dt = R(w)` but internally `dw` accumulates fluxes with the convention that it is subtracted).
Confirm the sign convention carefully.
The residual in `residuals.F90` adds the temporal term as:
```fortran
dw(i,j,k,l) = coefTime(0)*vol*w + sum(coefTime(m)*vol*wOld(m,...))
             + spatial_fluxes
```
So the spatial part is whatever is in `dw` when `equationMode == steady` or when the temporal source is not added.

### 1.2 IFunction Callback

For PETSc `TSTHETA`, we need an implicit residual:
```
F(t, U, Udot) = Udot - R(U, t)
```

**Python callback:**
```python
def ifunction(ts, t, U, Udot, F):
    # 1. Copy U into ADflow's state arrays
    w_array = U.getArray(readonly=True)
    adflow.set_state_vector(w_array)

    # 2. Update mesh/BCs for current time t if needed
    update_mesh_and_bcs(t)

    # 3. Evaluate spatial residual
    R = adflow.compute_spatial_residual()

    # 4. F = Udot - R
    f = F.getArray()
    udot = Udot.getArray(readonly=True)
    f[:] = udot - R
```

**Mesh motion:** For moving mesh problems (pitching airfoil), the mesh coordinates depend on `t`.
The `surfaceMeshCallback` / `volumeMeshCallback` pattern from `pyADflow.py:1385` can be adapted.

### 1.3 IJacobian Callback (MatShell)

PETSc needs the Jacobian `J = a * dF/dUdot + dF/dU = a*I - dR/dU`.

ADflow already provides matrix-free products via Tapenade AD:
- **Forward:** `computeMatrixFreeProductFwd(wdot, ...) -> dwdot` computes `(dR/dw) * wdot` via `master_d`
- **Reverse:** `computeMatrixFreeProductBwdFast(dwbar) -> wbar` computes `(dR/dw)^T * dwbar` via `master_b`

**MatShell approach (for adjoint we need transpose):**
```python
def dRdWT_matmult(mat, x, y):
    """y = (a*I - dR/dw)^T * x = a*x - (dR/dw)^T * x"""
    # Use computeMatrixFreeProductBwdFast for (dR/dw)^T * x
    wbar = adflow.computeJacobianVectorProductBwdFast(resBar=x.getArray())
    y_arr = y.getArray()
    y_arr[:] = a * x.getArray(readonly=True) - wbar

J_shell = PETSc.Mat().createPython([n, n], context=dRdWT_matmult)
```

**For the forward solve SNES (nonlinear solver within each TS step):**
PETSc `TSTHETA` uses SNES internally.
The Jacobian is `J = a*I - dR/dw`.
For the mat-vec product (GMRES), we can use the forward AD:
```python
def J_matmult(mat, x, y):
    """y = (a*I - dR/dw) * x"""
    dwdot = adflow.computeJacobianVectorProductFwd(wDot=x.getArray(), residualDeriv=True)
    y_arr = y.getArray()
    y_arr[:] = a * x.getArray(readonly=True) - dwdot
```

**Preconditioning:**
ADflow's ANK/NK solvers already build an approximate Jacobian (ILU-based) for preconditioning.
We can reuse this as a PC shell for PETSc's SNES KSP.
The existing `FormJacobianNK()` builds the approximate Jacobian in PETSc Mat format.
The existing `adjointKSP` uses ILU(p) on `dRdWTPreCon` for the steady adjoint.

### 1.4 Time Stepping Configuration

```python
ts = PETSc.TS().create(comm=comm)
ts.setType(PETSc.TS.Type.THETA)
ts.setTheta(0.5)  # Crank-Nicolson
ts.setEquationType(PETSc.TS.EquationType.ODE_EXPLICIT)  # F(Udot, U, t) = 0 form
ts.setIFunction(ifunction, F_template)
ts.setIJacobian(ijacobian, J_shell, J_shell)  # or separate preconditioner matrix
ts.setTime(0.0)
ts.setMaxTime(T_final)
ts.setTimeStep(dt)
ts.setExactFinalTime(PETSc.TS.ExactFinalTime.MATCHSTEP)
ts.setFromOptions()
```

### 1.5 Verification Plan

1. **Steady-state recovery:** Run the TS solver for a steady problem (constant BCs).
   Should converge to the same steady state as ADflow's native solver.

2. **Cylinder vortex shedding:** Compare lift/drag time histories between:
   - ADflow native BDF2 unsteady solver (`equationMode = unsteady`)
   - PETSc TSTHETA (CN) wrapping ADflow's spatial residual

   Both should produce the same Strouhal number, mean Cd, and amplitude of Cl oscillation (to within time discretization differences between BDF2 and CN).

3. **Order of accuracy:** Run with dt, dt/2, dt/4 and verify 2nd-order convergence rate for CN.

### 1.6 Key Files to Modify/Create

| File | Action | Description |
|------|--------|-------------|
| `adflow/pyADflow_TA.py` | Create | New Python module with PETSc TS wrapper class |
| `src/adjoint/adjointAPI.F90` | Modify | Add `computeSpatialResidual` subroutine |
| `src/NKSolver/NKSolvers.F90` | Reference | Reuse `setwVec`/`getwVec` scatter/gather patterns |

---

## Stage 2: ANK-Based Inner Solver for PETSc TS

### Goal
Use ADflow's ANK solver as the nonlinear solver (SNES) within PETSc TS time steps, rather than PETSc's default Newton.

### 2.1 Motivation

ADflow's ANK solver has several advantages for CFD:
- **Globalization:** CFL-based pseudo-time stepping provides superior robustness compared to line search alone
- **Approximate Jacobian:** Uses 1st-order flux approximations for the preconditioner, which is much cheaper than the full 2nd-order Jacobian
- **Proven for CFD:** ADflow's solver cascade (RK/DADI -> ANK -> NK) is well-tuned for compressible flows

### 2.2 Two Integration Strategies

**Strategy A: PETSc SNES wrapping ANK**
- Register a custom SNES type that internally calls ADflow's ANK/NK cascade
- PETSc TS calls this custom SNES at each time step
- Advantage: Clean separation between time integration (PETSc TS) and nonlinear solve (ADflow)
- Challenge: Need to reformulate ANK's pseudo-time stepping to solve the TS nonlinear system `F(U) = 0` where `F` is the theta-method residual

**Strategy B: Use PETSc SNES with ADflow preconditioner**
- Use PETSc's standard SNES Newton with GMRES
- Provide ADflow's approximate Jacobian (from ANK's `FormJacobianNK()`) as the preconditioner
- The mat-vec product uses exact AD (forward mode via `master_d`)
- Advantage: Simpler integration, leverages PETSc's SNES infrastructure
- Challenge: May lose ANK's CFL-based globalization benefits

**Recommendation:** Start with Strategy B (simpler).
If convergence is problematic, implement Strategy A.

### 2.3 Preconditioning with ADflow's Approximate Jacobian

ADflow builds approximate Jacobians in two contexts:
1. **ANK:** Uses `ANKJac` with 1st-order fluxes and ILU preconditioning
2. **Adjoint:** Uses `dRdWTPreCon` (the transposed approximate Jacobian) with ILU

For the TS nonlinear solve, we need `P ≈ aI - dR/dw`:
```python
def setup_preconditioner(a_shift):
    """Build P = aI - dR/dw (approximate)"""
    # Call ADflow's FormJacobianNK to get approximate dR/dw
    adflow.formJacobianNK()
    # Modify diagonal: add a*I
    # Apply ILU factorization
```

### 2.4 Convergence Tolerance for Inner Solves

This is a critical design decision.
In TA simulations, each time step's nonlinear system does NOT need to be converged to machine zero.

**Default ADflow behavior (native BDF):**
- `L2Conv = 1e-6` to `1e-10` (relative convergence for inner iterations)
- In practice, 2--4 orders of magnitude reduction per time step is common

**For PETSc TS with adjoint:**
- PETSc's SNES convergence tolerance directly affects adjoint accuracy
- If the forward nonlinear system is solved to tolerance `eps_fwd`, the adjoint will have O(`eps_fwd`) errors
- For optimization, `eps_fwd = 1e-8` to `1e-10` is typically needed
- **Trade-off:** Tighter tolerance = more inner iterations = more expensive per step, but cleaner adjoint

**Recommendation:** Use PETSc SNES with `rtol = 1e-10`, `atol = 1e-14`.
This ensures the discrete adjoint matches FD to high accuracy.

---

## Stage 3: Adjoint Time Integration

### Goal
Implement TSAdjoint for computing sensitivities of time-dependent objective functions.

### 3.1 Adjoint Terminal Conditions

For a final-time objective `J = g(U(T))`:
```python
# Adjoint seed: lambda(T) = dg/dU(T)
lambda_vec = U.duplicate()
lam = lambda_vec.getArray()
lam[:] = dg_dU  # e.g., [1, 0, 0, ...] for sensitivity of first state component

# Parameter sensitivity seed
mu_vec = PETSc.Vec().createSeq(n_params)
mu_vec.set(0.0)

ts.setCostGradients([lambda_vec], [mu_vec])
```

### 3.2 IJacobianP Callback

Sensitivity with respect to design parameters (mesh coordinates, angle of attack, etc.):
```python
def ijacobianp(ts, t, U, Udot, a, Jp):
    """dF/dp where F = Udot - R(U; p)"""
    # Jp = -dR/dp
    # Use computeMatrixFreeProductBwd with appropriate parameter seeds
    # For mesh coordinates: xvbar
    # For aero DVs: extrabar
```

ADflow already provides `computeMatrixFreeProductBwd` which computes:
```
Given dwbar (residual seeds):
    wbar     = (dR/dw)^T * dwbar       (state sensitivity)
    xvbar    = (dR/dXv)^T * dwbar      (mesh sensitivity)
    extrabar = (dR/dalpha,...)^T * dwbar (aero param sensitivity)
```

### 3.3 Checkpointing (TSTrajectory)

PETSc's `TSSetSaveTrajectory()` enables checkpointing.
For the adjoint backward sweep, the forward solution at each time step is needed.

**Options:**
- `TSTRAJECTORYMEMORY`: Store all time steps in memory (simplest, most memory)
- `TSTRAJECTORYBASIC`: Store to disk
- `TSTRAJECTORYSINGLEFILE`: Efficient single-file storage
- **Revolve checkpointing** (`TSTRAJECTORYVISUALIZATION` or custom): Optimal O(log N) memory with O(N log N) recomputation

**Configuration:**
```python
ts.setSaveTrajectory()
# For revolve checkpointing:
# ts.setFromOptions()  # with -ts_trajectory_type memory -ts_trajectory_max_cps <num_checkpoints>
```

**Memory estimate:** For a 1M cell mesh with 5 state variables, each time step snapshot is ~40 MB.
For 1000 time steps: 40 GB (memory) vs ~1 GB with revolve (50 checkpoints).

### 3.4 Time-Integrated Objectives

For objectives of the form `J = int_0^T f(U(t), t) dt`:

PETSc provides `TSCreateQuadratureTS()` which creates a sub-integrator for the cost function.
This requires an additional RHS function for the integrand.

```python
def cost_integrand(ts, t, U, F_cost):
    """Evaluate f(U, t) -- the integrand of the time-averaged objective"""
    # e.g., for time-averaged drag:
    w = U.getArray(readonly=True)
    adflow.set_state_vector(w)
    funcs = adflow.evalFunctions(...)
    f = F_cost.getArray()
    f[0] = funcs['cd']

qts = ts.createQuadratureTS(forward=True)
qts.setRHSFunction(cost_integrand, cost_vec)
```

The adjoint for the quadrature is handled automatically by PETSc TSAdjoint.

### 3.5 Key Implementation Steps

1. Add `ts.setSaveTrajectory()` to the forward solve
2. Implement `setCostGradients` with appropriate terminal conditions
3. Call `ts.adjointSolve()`
4. Extract `lambda` (state sensitivities) and `mu` (parameter sensitivities)
5. Chain-rule through DVGeo for shape derivatives: `dJ/dXdv = dJ/dXv * dXv/dXdv`

---

## Stage 4: Verification and Validation

### 4.1 Van der Pol Equation (Sanity Check)

Already implemented in `article_adflow_ta/code/ex20adj.py`.
Verifies PETSc TSAdjoint pipeline with finite differences.

### 4.2 ADflow Spatial Residual Verification

**Test: Steady NACA 0012**
1. Run PETSc TS forward to steady state (large T, constant BCs)
2. Compute adjoint for `dCl/dalpha`
3. Compare with ADflow's existing steady adjoint result
4. Should match to ~12+ digits (same code path for spatial residual)

### 4.3 Unsteady FD Verification

For each test case, verify adjoint sensitivities against finite differences:
```
dJ/dp ≈ [J(p + h) - J(p - h)] / (2h)
```

**Test cases (increasing complexity):**

1. **Cylinder vortex shedding (Re=100, laminar)**
   - 2D, no turbulence model, no mesh motion
   - Objective: time-averaged drag coefficient
   - Parameters: inflow velocity (Mach), initial condition perturbation
   - Expected: adjoint matches FD to ~6+ digits (limited by FD truncation)

2. **Pitching NACA 0012 airfoil**
   - 2D, prescribed mesh motion, RANS (SA turbulence)
   - Objective: peak lift coefficient, time-averaged moment
   - Parameters: pitching amplitude, frequency, Mach number
   - Tests: mesh motion (IJacobianP includes dR/dXv * dXv/dp), turbulence adjoint

3. **Transonic buffet (OAT15A airfoil, M=0.73, alpha=3.5 deg)**
   - 2D, no mesh motion, RANS (SA)
   - Objective: time-averaged pressure distribution, buffet frequency
   - Parameters: angle of attack, Mach number
   - Tests: shock-capturing adjoint stability, long time horizons

### 4.4 Convergence Studies

For each test case:
- **Temporal convergence:** Run with dt, dt/2, dt/4; verify 2nd-order convergence of both forward solution and adjoint sensitivities
- **Inner solve tolerance:** Run with `SNES rtol = 1e-6, 1e-8, 1e-10, 1e-12`; verify adjoint accuracy improves with tighter inner tolerance
- **Checkpointing:** Verify adjoint results are identical with and without revolve checkpointing

### 4.5 Partial Derivative Verification (Current Focus)

**Motivation:** Before debugging the adjoint time-stepping and PC (which are
entangled with solver convergence), verify that all AD building blocks are
correct independently. This isolates AD bugs from solver/PC issues.

**Testing methodology for each partial derivative A of function f:**

1. **Forward AD vs FD:** Pick random direction v, compare `A*v` from AD
   with central FD `(f(x+εv) - f(x-εv))/(2ε)`. Validates forward mode.
2. **Reverse AD vs Forward AD (dot product test):** For random u, v,
   compute `y = A*v` (fwd) and `z = A^T*u` (rev), check `<y,u> = <v,z>`.
   Should match to ~12 digits. Validates reverse is exact transpose.

Together these prove both forward and reverse AD are correct.

**Partials to test:**

| # | Partial | Type | Forward API | Reverse API |
|---|---------|------|-------------|-------------|
| 1 | ∂R/∂w | n×n | `computeJacobianVectorProductFwd(wDot=v, residualDeriv=True)` | `computeJacobianVectorProductBwd(resBar=u, wDeriv=True)` |
| 2 | ∂R/∂alpha | n×1 | `computeJacobianVectorProductFwd(xDvDot={"alpha":1}, residualDeriv=True)` | `computeJacobianVectorProductBwd(resBar=u, xDvDeriv=True)["alpha"]` |
| 3 | ∂CL/∂w | 1×n | `computeJacobianVectorProductFwd(wDot=v, funcDeriv=True)["cl"]` | `computeJacobianVectorProductBwd(funcsBar={"cl":1}, wDeriv=True)` |
| 4 | ∂CL/∂alpha | scalar | `computeJacobianVectorProductFwd(xDvDot={"alpha":1}, funcDeriv=True)["cl"]` | `computeJacobianVectorProductBwd(funcsBar={"cl":1}, xDvDeriv=True)["alpha"]` |

Phase 3 (moving mesh) adds:

| 5 | ∂R/∂Xv | n×3N | `Fwd(xVDot=v, residualDeriv=True)` | `Bwd(resBar=u, xVDeriv=True)` |
| 6 | ∂CL/∂Xv | 1×3N | `Fwd(xVDot=v, funcDeriv=True)` | `Bwd(funcsBar={"cl":1}, xVDeriv=True)` |

**turbResScale handling:**  AD routines return `dR_ad = resScale * dR_raw`
where `resScale = (1/volRef) * turbResScale`.  The TA adjoint uses
`R_ta = R_ad / turbResScale`.  Corrections:
- Forward: divide SA DOF rows of AD output by turbResScale
- Reverse: divide SA DOFs of resBar by turbResScale before passing to AD
- For Euler (no turbulence), no correction needed.

**Instance-wise vs full matrix:**  Instance-wise (random vector) testing is
sufficient.  Full matrix requires `n_global` matvecs (~10k-100k), while
3-5 random directions give high statistical confidence.  Only resort to
unit vectors if a random test fails, to isolate the bad DOF.

**Phased plan:**

- **Phase 1:** Fixed mesh, steady state — test all 4 partials at converged
  steady state with 3 random directions each.
  Script: `article_adflow_ta/code/examples/test_partials.py`
- **Phase 2:** Fixed mesh, non-equilibrium — run a few BDF1 steps from
  steady state, repeat tests at intermediate states.
- **Phase 3:** Moving mesh — add ∂R/∂Xv and ∂CL/∂Xv tests at states from
  a pitching simulation.

---

## Stage 5: Optimization Demonstration

### 5.1 Problem Setup

**Pitching airfoil drag minimization:**
- Objective: minimize time-averaged drag
- Subject to: time-averaged lift constraint
- Design variables: airfoil shape (FFD control points via DVGeo)
- Optimizer: SNOPT or IPOPT via pyOptSparse

### 5.2 Integration with MACH Framework

```
DVGeo (shape) -> IDWarp (mesh) -> ADflow+PETSc TS (flow + adjoint) -> pyOptSparse (optimizer)
```

The adjoint provides `dJ/dXs` (surface mesh sensitivity).
Chain rule through IDWarp and DVGeo gives `dJ/dXdv` (design variable sensitivity).

---

## Critical Design Decisions

### Inexact Inner Solves and Adjoint Accuracy

**The problem:** In time-accurate simulations, ADflow's inner solver (solveState) typically converges the residual by 4--8 orders of magnitude per time step, NOT to machine zero.
This means the forward solution `U^{n+1}` satisfies:
```
G(U^{n+1}) = epsilon  (not exactly zero)
```
where `G` is the nonlinear system for the time step and `epsilon` is the residual tolerance.

**Impact on adjoint:**
The discrete adjoint is derived assuming `G(U^{n+1}) = 0` exactly.
If the forward residual is O(epsilon), the adjoint error is also O(epsilon).
This manifests as disagreement between adjoint and finite-difference sensitivities.

**Quantitative analysis:**
- For optimization, we need adjoint accuracy of ~1e-6 to 1e-8 relative
- ADflow's default L2Conv = 1e-6 to 1e-10
- With 1000 time steps, errors can accumulate: total error ~ N_steps * epsilon_per_step
- Recommendation: use `SNES rtol = 1e-10` (10 orders relative convergence per step)

**Mitigation strategies:**
1. **Tight inner tolerance:** Simply converge each time step to high accuracy (safe but expensive)
2. **Adjoint-consistent inexact Newton:** Modify the adjoint equations to account for nonzero forward residuals (complex, research-level)
3. **Practical approach:** Start with tight tolerance; if too expensive, study the sensitivity of adjoint accuracy to inner tolerance empirically

### Choice of TSTHETA vs TSBDF

**TSTHETA (Crank--Nicolson, theta=0.5):**
- 2nd-order accurate, A-stable
- Supported by PETSc TSAdjoint
- One linear system solve per time step
- May exhibit oscillations for stiff problems (no L-stability)

**TSBDF:**
- PETSc `TSBDF` does NOT currently support `TSAdjointSolve` (confirmed in PETSc source)
- Would require implementing custom adjoint for BDF
- BDF2 is L-stable (no oscillations) but same order as CN

**Decision:** Use TSTHETA.
If oscillations are observed, try theta=0.55 (slightly implicit, damps oscillations at cost of formal order).
For problems requiring L-stability, implement adjoint for BDF2 as a future extension, or use small enough dt that CN oscillations are negligible.

### MatShell vs Explicit Jacobian Assembly

**MatShell (matrix-free):**
- Uses AD (Tapenade) for exact Jacobian-vector products
- No memory for storing the full Jacobian (which can be 5N x 5N, huge for 3D)
- Requires a preconditioner (ILU on approximate Jacobian)
- This is what ADflow already does for the steady adjoint

**Explicit Jacobian assembly:**
- ADflow assembles an approximate Jacobian for ANK/NK preconditioning
- Full exact Jacobian is NOT assembled (too expensive for structured multi-block)
- Not a viable option for the TS Jacobian

**Decision:** Use MatShell for both forward (SNES) and adjoint solves.
Reuse ADflow's existing approximate Jacobian (ANK's 1st-order or adjoint's `dRdWTPreCon`) as the preconditioner.

---

## Implementation Roadmap

### Phase 1 (Weeks 1-4): Forward Solver
- [ ] Expose `computeSpatialResidual()` as standalone callable from Python
- [ ] Implement PETSc TS IFunction callback wrapping ADflow spatial residual
- [ ] Implement MatShell IJacobian using `computeMatrixFreeProductFwd`
- [ ] Set up TSTHETA with CN (theta=0.5)
- [ ] Verify steady-state recovery (compare with native ADflow)
- [ ] Verify cylinder vortex shedding (compare time histories with native ADflow BDF2)

### Phase 2 (Weeks 5-8): Adjoint Solver
- [ ] Add `ts.setSaveTrajectory()` for checkpointing
- [ ] Implement adjoint MatShell using `computeMatrixFreeProductBwdFast`
- [ ] Implement IJacobianP callback for parameter sensitivities
- [ ] Implement `setCostGradients` for final-time objectives
- [ ] Verify adjoint vs FD for cylinder (laminar, no mesh motion)
- [ ] Implement `TSCreateQuadratureTS` for time-integrated objectives
- [ ] Verify time-averaged drag adjoint vs FD

### Phase 3 (Weeks 9-12): Advanced Features + V&V
- [ ] Add mesh motion support (pitching airfoil: IFunction and IJacobianP update with mesh)
- [ ] Verify pitching airfoil adjoint vs FD
- [ ] Study inner solve tolerance vs adjoint accuracy
- [ ] Implement revolve checkpointing for memory-efficient long simulations
- [ ] Transonic buffet test case

### Phase 4 (Weeks 13-16): Optimization + Paper
- [ ] Integration with pyOptSparse / DVGeo / IDWarp for shape optimization
- [ ] Pitching airfoil drag minimization demo
- [ ] Performance profiling and scalability study
- [ ] Complete paper with results

---

## Appendix A: Key ADflow Source Files

| File | Contents |
|------|----------|
| `src/solver/solvers.F90` | `solverUnsteadyStep()`, `solveState()` (main solver cascade) |
| `src/solver/solverUtils.F90` | `shiftSolution()` (shift wOld arrays for BDF) |
| `src/solver/residuals.F90` | `residual()` (spatial + temporal residual), `initres()`, `sourceTerms()` |
| `src/utils/utils.F90` | `setCoefTimeIntegrator()` (BDF coefficient selection) |
| `src/NKSolver/NKSolvers.F90` | `ANKStep()`, `NKStep()`, `computeResidualNK()`, `FormJacobianNK()` |
| `src/adjoint/adjointAPI.F90` | `computeMatrixFreeProductFwd()`, `computeMatrixFreeProductBwd()`, `computeMatrixFreeProductBwdFast()`, `solveAdjoint()` |
| `src/modules/ADjointPETSc.F90` | PETSc objects: `dRdWT`, `dRdWTPreCon`, `adjointKSP`, Vec templates |
| `src/modules/inputParam.F90` | `timeIntegrationScheme`, ANK/NK parameters |
| `adflow/pyADflow.py` | `solveTimeStep()` (line 1526), `solveAdjoint()` (line 4052), `computeJacobianVectorProductFwd/Bwd()` (lines 4557, 4814) |

## Appendix B: PETSc TS Adjoint Reference

Minimal working example: `article_adflow_ta/code/ex20adj.py` (van der Pol equation).

**Key PETSc TS adjoint API (petsc4py):**
```python
ts.setSaveTrajectory()                    # enable trajectory storage for adjoint
ts.solve(U)                               # forward solve
ts.setCostGradients([lambda_vecs], [mu_vecs])  # set adjoint seeds
ts.adjointSolve()                         # backward sweep

# For time-integrated objectives:
qts = ts.createQuadratureTS(forward=True)
qts.setRHSFunction(cost_integrand, cost_vec)
```

**Key callbacks:**
- `IFunction(ts, t, U, Udot, F)`: Implicit residual F(t, U, Udot) = 0
- `IJacobian(ts, t, U, Udot, a, J, P)`: Jacobian dF/dUdot * a + dF/dU
- `IJacobianP(ts, t, U, Udot, a, Jp)`: Parameter Jacobian dF/dp

## Appendix C: Sign Conventions

ADflow internal sign convention for residuals:
```
dw = spatial_fluxes + temporal_source
```
where `dw` is used as the RHS of the pseudo-time iteration `w^{n+1} = w^n + dt_pseudo * dw`.

For PETSc TS IFunction, we need `F(Udot, U) = 0`:
```
F = Udot - R(U)
```
where `R(U)` is the spatial RHS (opposite sign to ADflow's internal `dw` when temporal terms are removed).

**Verify sign by checking:** If ADflow computes `dw` (steady residual, no temporal terms), then `R(U) = -dw` typically.
Confirm empirically by comparing `R(U)` against a known solution.
