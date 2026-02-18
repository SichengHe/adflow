"""
PETSc TS wrapper for ADflow time-accurate simulation.

Wraps ADflow's spatial residual as PETSc TS callbacks (IFunction, IJacobian)
to drive time integration via the implicit theta method (Crank--Nicolson).

This replaces ADflow's native BDF time loop with PETSc TSTHETA,
enabling future adjoint sensitivity analysis via TSAdjointSolve.

Sign convention
---------------
ADflow's getRes() in steady mode returns dw/volRef, where dw accumulates
spatial fluxes. At steady state, dw = 0 (fluxes balance). In the
pseudo-time update w += alpha * dt_pseudo * dw/vol, positive dw means
the state should increase. Therefore dw/vol = R(w), the spatial RHS
in dw/dt = R(w).

For PETSc TS IFunction: F(t, U, Udot) = Udot - R(U) = 0.
"""

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc


# ADflow constant for equationMode (from src/modules/constants.F90)
_STEADY = 1


class ADflowJacobian:
    """MatShell context for the Jacobian J = a*I - dR/dw."""

    def __init__(self, solver):
        self.solver = solver
        self.shift = 1.0

    def mult(self, mat, x, y):
        """y = (a*I - dR/dw) * x (forward mode for SNES)."""
        x_arr = x.getArray(readonly=True).copy()
        dRdw_x = self.solver.computeJacobianVectorProductFwd(
            wDot=x_arr, residualDeriv=True
        )
        y_arr = y.getArray()
        y_arr[:] = self.shift * x_arr - dRdw_x

    def multTranspose(self, mat, x, y):
        """y = (a*I - dR/dw)^T * x (reverse mode for adjoint)."""
        x_arr = x.getArray(readonly=True).copy()
        wbar = self.solver.computeJacobianVectorProductBwd(
            resBar=x_arr, wDeriv=True
        )
        y_arr = y.getArray()
        y_arr[:] = self.shift * x_arr - wbar


class ADflowTS:
    """
    Wraps ADflow's spatial residual with PETSc TSTHETA (Crank--Nicolson).

    ADflow is initialized in unsteady mode (so grid motion arrays are
    allocated), but the equationMode is temporarily switched to steady
    during residual evaluation so that getResidual returns only the
    spatial residual (no temporal BDF terms). PETSc TS handles the
    time integration externally.

    Parameters
    ----------
    solver : ADFLOW
        An initialized ADFLOW solver instance (unsteady mode).
    ap : AeroProblem
        The aerodynamic problem definition.
    dt : float
        Physical time step size [s].
    t_final : float
        Final simulation time [s].
    theta : float, optional
        Implicitness parameter. Default 0.5 (Crank--Nicolson).
    grid_motion : bool, optional
        Whether the mesh moves (prescribed motion). Default False.
    """

    def __init__(self, solver, ap, dt, t_final, theta=0.5, grid_motion=False):
        self.solver = solver
        self.ap = ap
        self.dt = dt
        self.t_final = t_final
        self.theta = theta
        self.grid_motion = grid_motion
        self.comm = solver.comm

        # Local state size on this processor
        self.n_local = solver.getStateSize()

        # Storage for force coefficients at each time step
        self.time_history = []
        self.cl_history = []
        self.cd_history = []
        self.cmz_history = []

        # Jacobian context
        self.jac_ctx = ADflowJacobian(solver)

        # Track the last time at which the grid was updated
        self._grid_time = None

        # Track residual evaluation count
        self.n_res_eval = 0

    def setup(self):
        """Create PETSc TS, vectors, and MatShell."""
        comm = self.comm

        # Create PETSc vectors
        self.F_vec = PETSc.Vec().createMPI(self.n_local, comm=comm)
        self.U_vec = PETSc.Vec().createMPI(self.n_local, comm=comm)

        # Set initial state from ADflow
        w0 = self.solver.getStates()
        self.U_vec.setArray(w0)

        # Create MatShell for Jacobian (supports mult and multTranspose)
        self.J_shell = PETSc.Mat().createPython(
            [(self.n_local, None), (self.n_local, None)],
            context=self.jac_ctx,
            comm=comm,
        )
        self.J_shell.setUp()

        # Create TS
        self.ts = PETSc.TS().create(comm=comm)
        self.ts.setType(PETSc.TS.Type.THETA)
        self.ts.setTheta(self.theta)
        self.ts.setIFunction(self._ifunction, self.F_vec)
        self.ts.setIJacobian(self._ijacobian, self.J_shell, self.J_shell)
        self.ts.setTime(0.0)
        self.ts.setMaxTime(self.t_final)
        self.ts.setTimeStep(self.dt)
        self.ts.setMaxSteps(int(round(self.t_final / self.dt)))
        self.ts.setExactFinalTime(PETSc.TS.ExactFinalTime.MATCHSTEP)

        # SNES (nonlinear solver) settings for each time step
        snes = self.ts.getSNES()
        snes.setTolerances(rtol=1e-10, atol=1e-14, max_it=50)

        # KSP (linear solver) settings
        ksp = snes.getKSP()
        ksp.setType("gmres")
        ksp.setTolerances(rtol=1e-6, max_it=100)

        # Set post-step monitor to record force coefficients
        self.ts.setPostStep(self._post_step)

        # Allow command-line options to override
        self.ts.setFromOptions()

    def _update_grid(self, t):
        """
        Update the mesh coordinates and velocities for time t.

        This sets the internal timeUnsteady variable and calls ADflow's
        grid motion routines. Only updates if the time has changed.
        """
        if not self.grid_motion:
            return
        if self._grid_time is not None and abs(t - self._grid_time) < 1e-15:
            return

        # Set the physical time so that grid motion uses the correct time
        self.solver.adflow.monitor.timeunsteady = t

        # Shift coordinate history and update mesh position
        self.solver.adflow.preprocessingapi.shiftcoorandvolumes()
        self.solver.adflow.solvers.updateunsteadygeometry()

        self._grid_time = t

    def _get_spatial_residual(self):
        """
        Evaluate the spatial residual R(w) without temporal terms.

        Temporarily switches equationMode to steady so that ADflow's
        residual evaluation initializes dw to zero (no BDF temporal
        source terms), then restores the original mode.

        Returns
        -------
        res : ndarray
            Spatial residual vector dw/volRef.
        """
        adflow = self.solver.adflow

        # Save the current equation mode
        orig_mode = adflow.inputphysics.equationmode

        # Switch to steady: residual = spatial fluxes only (no temporal terms)
        adflow.inputphysics.equationmode = _STEADY

        # Evaluate residual
        res = self.solver.getResidual(self.ap, releaseAdjointMemory=False)

        # Restore equation mode
        adflow.inputphysics.equationmode = orig_mode

        return res

    def _ifunction(self, ts, t, U, Udot, F):
        """
        IFunction callback: F = Udot - R(w, t).

        For PETSc TS, the implicit residual is F(t, U, Udot) = 0.
        We have F = Udot - R(U), where R is the spatial residual.
        """
        # Update mesh for current time
        self._update_grid(t)

        # Set ADflow state to U
        u_arr = U.getArray(readonly=True).copy()
        self.solver.setStates(u_arr)

        # Evaluate spatial residual R(w) (no temporal terms)
        res = self._get_spatial_residual()

        # F = Udot - R(w)
        f_arr = F.getArray()
        udot_arr = Udot.getArray(readonly=True)
        f_arr[:] = udot_arr - res

        self.n_res_eval += 1

    def _ijacobian(self, ts, t, U, Udot, a, J, P):
        """
        IJacobian callback: J = a*I - dR/dw.

        Stores the shift 'a' and current state for MatShell mult.
        """
        # Update mesh for current time
        self._update_grid(t)

        # Update the shift in the MatShell context
        self.jac_ctx.shift = a

        # Set ADflow state for Jacobian evaluation
        u_arr = U.getArray(readonly=True).copy()
        self.solver.setStates(u_arr)

        # Signal that the matrix structure has not changed
        J.assemble()
        return True

    def _post_step(self, ts):
        """Record force coefficients after each time step."""
        t = ts.getTime()
        U = ts.getSolution()

        # Set state in ADflow
        u_arr = U.getArray(readonly=True).copy()
        self.solver.setStates(u_arr)

        # Evaluate functions
        funcs = {}
        self.solver.evalFunctions(self.ap, funcs, evalFuncs=["cl", "cd", "cmz"])

        ap_name = self.ap.name
        cl = funcs.get(f"{ap_name}_cl", 0.0)
        cd = funcs.get(f"{ap_name}_cd", 0.0)
        cmz = funcs.get(f"{ap_name}_cmz", 0.0)

        self.time_history.append(t)
        self.cl_history.append(cl)
        self.cd_history.append(cd)
        self.cmz_history.append(cmz)

        if self.comm.rank == 0:
            step = ts.getStepNumber()
            print(f"  TS step {step:4d} | t = {t:.6e} | "
                  f"CL = {cl:.6e} | CD = {cd:.6e} | CMz = {cmz:.6e}")

    def solve(self):
        """Run the forward time integration."""
        if self.comm.rank == 0:
            print("=" * 70)
            print(f"ADflowTS: Starting PETSc TSTHETA (theta={self.theta})")
            print(f"  dt = {self.dt}, t_final = {self.t_final}")
            print(f"  grid_motion = {self.grid_motion}")
            n_global = self.comm.allreduce(self.n_local)
            print(f"  n_local = {self.n_local}, n_global = {n_global}")
            print("=" * 70)

        # Record initial state
        self._post_step(self.ts)

        # Solve
        self.ts.solve(self.U_vec)

        # Check convergence
        reason = self.ts.getConvergedReason()
        if self.comm.rank == 0:
            print("=" * 70)
            print(f"ADflowTS: Finished. Converged reason = {reason}")
            print(f"  Total residual evaluations: {self.n_res_eval}")
            print(f"  Final time: {self.ts.getTime():.6e}")
            print(f"  Total steps: {self.ts.getStepNumber()}")
            print("=" * 70)

        # Set final state back into ADflow
        w_final = self.U_vec.getArray(readonly=True).copy()
        self.solver.setStates(w_final)

        return reason

    def get_history(self):
        """Return time history of force coefficients as a dict."""
        return {
            "time": np.array(self.time_history),
            "cl": np.array(self.cl_history),
            "cd": np.array(self.cd_history),
            "cmz": np.array(self.cmz_history),
        }
