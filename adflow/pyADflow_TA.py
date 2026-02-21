"""
PETSc TS wrapper for ADflow time-accurate simulation.

Wraps ADflow's spatial residual as PETSc TS callbacks (IFunction, IJacobian)
to drive time integration via the implicit theta method (Crank--Nicolson).

This replaces ADflow's native BDF time loop with PETSc TSTHETA,
enabling adjoint sensitivity analysis via TSAdjointSolve.

Architecture
------------
- **IFunction**: ``F(t, U, Udot) = Udot - R(U)``, where R is the spatial
  residual computed via the ``master`` code path (AD forward with zero seed)
  followed by ``getResDw`` to read the primal result.
- **IJacobian**: ``J = a*I - dR/dw`` via MatShell (matrix-free AD forward/reverse).
- **Preconditioner**: assembled approximate Jacobian (1st-order PC stencil)
  with ASM + ILU(2), built by Fortran ``setupTSPreconditioner``.

Both IFunction and Jacobian operate in ``1/volRef`` space (turbResScale
removed from SA DOFs).

Usage::

    mpirun -np 2 python run.py \\
        -snes_monitor -ksp_monitor
"""

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc


# ADflow constants for equationMode (from src/modules/constants.F90)
_STEADY = 1
_UNSTEADY = 2
_TIME_SPECTRAL = 3


class ADflowJacobianFD:
    """MatShell context for the Jacobian J = a*I - dR/dw via finite differences.

    Uses central FD of getRes() (blocketteRes path) so the Jacobian is
    exactly consistent with the IFunction residual. This avoids the
    master_d vs blocketteRes code path discrepancy.

    Cost: 2 extra getRes() evaluations per Jacobian-vector product.
    """

    def __init__(self, solver, h_fd=1e-6):
        self.solver = solver
        self.shift = 1.0
        self.h_fd = h_fd
        self._n_local = solver.getStateSize()

    def _get_res_steady(self):
        """Evaluate spatial residual in steady mode via getRes()."""
        adflow = self.solver.adflow
        orig_mode = adflow.inputphysics.equationmode
        adflow.inputphysics.equationmode = _STEADY
        res = np.zeros(self._n_local)
        res = adflow.nksolver.getres(res)
        adflow.inputphysics.equationmode = orig_mode
        return res

    def mult(self, mat, x, y):
        """y = (a*I - dR/dw) * x via central FD of getRes()."""
        x_arr = x.getArray(readonly=True).copy()
        w0 = self.solver.getStates().copy()
        h = self.h_fd

        # FD: dR/dw * v ≈ (R(w+h*v) - R(w-h*v)) / (2h)
        self.solver.setStates(w0 + h * x_arr)
        rp = self._get_res_steady()
        self.solver.setStates(w0 - h * x_arr)
        rm = self._get_res_steady()
        self.solver.setStates(w0)  # restore

        dRdw_x = (rp - rm) / (2.0 * h)
        y_arr = y.getArray()
        y_arr[:] = self.shift * x_arr - dRdw_x


class ADflowJacobianAD:
    """MatShell context for the Jacobian J = a*I - dR/dw via AD.

    Uses Tapenade-generated master_d (forward mode) and master_b (reverse
    mode) for Jacobian-vector products. The equation mode is temporarily
    switched to STEADY so that only the spatial Jacobian is computed.

    The AD routines include resScale_d / resScale_b, which apply
    turbResScale to SA DOFs. This scaling is undone here to match the
    IFunction which uses getRes() (1/volRef only, no turbResScale).

    Note: master_d and blocketteRes use different code paths, leading to
    a ~8% discrepancy in the Jacobian for flow DOFs. For time-accurate
    adjoint, this class should only be used with a future Fortran routine
    that evaluates the residual via the master path.
    """

    def __init__(self, solver):
        self.solver = solver
        self.shift = 1.0

        # Cache turbulence scaling info for SA DOF correction
        adflow = solver.adflow
        self._nw = int(adflow.flowvarrefstate.nw)
        self._nwf = int(adflow.flowvarrefstate.nwf)
        self._n_turb = self._nw - self._nwf
        if self._n_turb > 0:
            self._trs = adflow.inputiteration.turbresscale[:self._n_turb].copy()
            self._inv_trs = 1.0 / self._trs
        else:
            self._trs = None
            self._inv_trs = None

    def _undo_turb_res_scale(self, vec):
        """Divide SA DOF rows of vec by turbResScale (in-place)."""
        if self._n_turb > 0:
            n_cells = len(vec) // self._nw
            v2d = vec.reshape(n_cells, self._nw)
            for l in range(self._n_turb):
                v2d[:, self._nwf + l] *= self._inv_trs[l]

    def _apply_inv_turb_res_scale_to_resbar(self, resbar):
        """Divide SA DOF elements of resBar by turbResScale (returns copy)."""
        if self._n_turb > 0:
            out = resbar.copy()
            n_cells = len(out) // self._nw
            v2d = out.reshape(n_cells, self._nw)
            for l in range(self._n_turb):
                v2d[:, self._nwf + l] *= self._inv_trs[l]
            return out
        return resbar

    def mult(self, mat, x, y):
        """y = (a*I - dR/dw) * x (forward mode via master_d)."""
        adflow = self.solver.adflow
        orig_mode = adflow.inputphysics.equationmode
        adflow.inputphysics.equationmode = _STEADY

        x_arr = x.getArray(readonly=True).copy()
        dRdw_x = self.solver.computeJacobianVectorProductFwd(
            wDot=x_arr, residualDeriv=True
        )

        # Undo turbResScale on SA DOF rows
        self._undo_turb_res_scale(dRdw_x)

        y_arr = y.getArray()
        y_arr[:] = self.shift * x_arr - dRdw_x

        adflow.inputphysics.equationmode = orig_mode

    def multTranspose(self, mat, x, y):
        """y = (a*I - dR/dw)^T * x (reverse mode via master_b)."""
        adflow = self.solver.adflow
        orig_mode = adflow.inputphysics.equationmode
        adflow.inputphysics.equationmode = _STEADY

        x_arr = x.getArray(readonly=True).copy()

        # Divide SA DOFs of resBar by turbResScale before reverse AD
        x_adj = self._apply_inv_turb_res_scale_to_resbar(x_arr)
        wbar = self.solver.computeJacobianVectorProductBwd(
            resBar=x_adj, wDeriv=True
        )

        y_arr = y.getArray()
        y_arr[:] = self.shift * x_arr - wbar

        adflow.inputphysics.equationmode = orig_mode


class ADflowJacobianMaster:
    """MatShell context for J = a*I - dR/dw via AD — NO turbResScale correction.

    WARNING: The current IFunction undoes turbResScale from SA DOFs.
    This class does NOT undo turbResScale, so the Jacobian and IFunction
    are INCONSISTENT.  Use ``ADflowJacobianAD`` instead for normal operation.

    This class is kept for debugging/testing only.
    """

    def __init__(self, solver):
        self.solver = solver
        self.shift = 1.0

    def mult(self, mat, x, y):
        """y = (a*I - dR/dw) * x (forward mode via master_d)."""
        adflow = self.solver.adflow
        orig_mode = adflow.inputphysics.equationmode
        adflow.inputphysics.equationmode = _STEADY

        x_arr = x.getArray(readonly=True).copy()
        dRdw_x = self.solver.computeJacobianVectorProductFwd(
            wDot=x_arr, residualDeriv=True
        )

        y_arr = y.getArray()
        y_arr[:] = self.shift * x_arr - dRdw_x

        adflow.inputphysics.equationmode = orig_mode

    def multTranspose(self, mat, x, y):
        """y = (a*I - dR/dw)^T * x (reverse mode via master_b)."""
        adflow = self.solver.adflow
        orig_mode = adflow.inputphysics.equationmode
        adflow.inputphysics.equationmode = _STEADY

        x_arr = x.getArray(readonly=True).copy()
        wbar = self.solver.computeJacobianVectorProductBwd(
            resBar=x_arr, wDeriv=True
        )

        y_arr = y.getArray()
        y_arr[:] = self.shift * x_arr - wbar

        adflow.inputphysics.equationmode = orig_mode


# Keep backward-compatible name
ADflowJacobian = ADflowJacobianAD


class _TSPreconditioner:
    """PETSc PC 'python' context that wraps the Fortran ASM+ILU preconditioner."""

    def __init__(self, ts_wrapper):
        self.ts_wrapper = ts_wrapper

    def apply(self, pc, x, y):
        """y = P^{-1} * x using Fortran-side ILU factorization."""
        x_arr = x.getArray(readonly=True).copy()
        y_arr = np.zeros(self.ts_wrapper.n_local)
        self.ts_wrapper.solver.adflow.nksolver.applytspreconditioner(x_arr, y_arr)
        y.getArray()[:] = y_arr

    def applyTranspose(self, pc, x, y):
        """y = P^{-T} * x using Fortran-side ILU^{-T} (for adjoint)."""
        x_arr = x.getArray(readonly=True).copy()
        y_arr = np.zeros(self.ts_wrapper.n_local)
        self.ts_wrapper.solver.adflow.nksolver.applytspreconditionertranspose(x_arr, y_arr)
        y.getArray()[:] = y_arr


class _ScalePreconditioner:
    """Trivial PC: y = x / shift.  For J = a*I - dR/dw with large a,
    this approximates J^{-1} ≈ (1/a) * I."""

    def __init__(self, ts_wrapper):
        self.ts_wrapper = ts_wrapper

    def apply(self, pc, x, y):
        a = self.ts_wrapper.jac_ctx.shift
        x_arr = x.getArray(readonly=True)
        y.getArray()[:] = x_arr / a

    def applyTranspose(self, pc, x, y):
        # Scale PC is symmetric: P^{-T} = P^{-1} = (1/a) * I
        a = self.ts_wrapper.jac_ctx.shift
        x_arr = x.getArray(readonly=True)
        y.getArray()[:] = x_arr / a


class _NativeAdjointPC:
    """PC that uses ADflow's native adjoint solver as preconditioner.

    For the TS adjoint, we need to precondition ``J = a*I - dR/dw`` where
    ``a = 1/(theta*dt)`` is the temporal shift.  Since ``a`` is small
    relative to the spectral radius of ``dR/dw``, the dominant part is
    ``-dR/dw``.  We use ``P = -dR/dw`` as preconditioner:

    - ``PCApply(x, y)``:  ``y = P^{-1} x = -(dR/dw)^{-1} x``
      via ``solveDirectForRHS``
    - ``PCApplyTranspose(x, y)``:  ``y = P^{-T} x = -(dR/dw^T)^{-1} x``
      via ``solveAdjointForRHS``

    For the adjoint, PETSc calls ``KSPSolveTranspose`` which uses
    ``PCApplyTranspose``.  The preconditioned adjoint system is::

        P^{-T} J^T = (dR/dw^T)^{-1} (a*I - dR/dw^T)
                    = a*(dR/dw^T)^{-1} - I

    Since ``a << spectral_radius(dR/dw)``, eigenvalues cluster near -1,
    so GMRES converges in O(1) iterations.

    Note: ``PCApplyTranspose`` uses ``solveAdjointForRHS`` which calls
    ``KSPSolve`` (not ``KSPSolveTranspose``), avoiding nested
    ``KSPSolveTranspose`` calls that cause SEGV in PETSc.
    """

    def __init__(self, ts_wrapper, inner_tol=0.01):
        self.ts_wrapper = ts_wrapper
        self.inner_tol = inner_tol
        self._setup_done = False

    def setUp(self, pc):
        """Set up ADflow's native adjoint solver (once)."""
        if not self._setup_done:
            solver = self.ts_wrapper.solver
            # Set up adjoint PETSc vectors, matrices, and KSP
            # Must destroy NK/ANK first (matches pyADflow._setupAdjoint)
            if not solver.adjointSetup:
                solver.adflow.nksolver.destroynksolver()
                solver.adflow.anksolver.destroyanksolver()
                solver.adflow.adjointapi.createpetscvars()
                solver.adflow.adjointapi.setupallresidualmatricesfwd()
                solver.adflow.adjointapi.setuppetscksp()
                solver.adjointSetup = True
            self._setup_done = True

    def apply(self, pc, x, y):
        """y = P^{-1} * x = -(dR/dw)^{-1} * x."""
        x_arr = x.getArray(readonly=True).copy()
        y_arr = self.ts_wrapper.solver.adflow.adjointapi.solvedirectforrhs(
            x_arr, self.inner_tol
        )
        y.getArray()[:] = -y_arr

    def applyTranspose(self, pc, x, y):
        """y = P^{-T} * x = -(dR/dw^T)^{-1} * x."""
        x_arr = x.getArray(readonly=True).copy()
        y_arr = self.ts_wrapper.solver.adflow.adjointapi.solveadjointforrhs(
            x_arr, self.inner_tol
        )
        y.getArray()[:] = -y_arr


class ADflowDADISNES:
    """Custom PETSc SNES that delegates forward solves to ADflow's native
    DADI+MG (plus optional ANK/NK) solver.

    PETSc TSTHETA(theta=1) calls ``SNESSolve`` once per time step.  Instead
    of Newton iterations with IFunction/IJacobian, this SNES calls ADflow's
    ``solverUnsteadyStep()`` which runs DADI + multigrid sub-iterations
    (and optionally ANK/NK) until the BDF1 residual converges.

    The KSP on this SNES is **not** used during the forward solve but **is**
    used by ``TSAdjointStep_Theta`` for the adjoint linear system
    ``K^T lambda = rhs`` via ``KSPSolveTranspose``.

    This avoids the need to differentiate through DADI/MG iterations: the
    adjoint uses the implicit function theorem and only requires the spatial
    Jacobian dR/dw (evaluated by the existing AD routines).
    """

    def __init__(self, ts_wrapper):
        self.ts_wrapper = ts_wrapper
        self.solver = ts_wrapper.solver
        self._step_count = 0

    def setUp(self, snes):
        """Called by PETSc when SNES is set up."""
        pass

    def solve(self, snes, b, x):
        """Advance one BDF1 time step via ADflow's native DADI+MG solver.

        Parameters
        ----------
        snes : PETSc.SNES
            The PETSc SNES context (not used for nonlinear iteration).
        b : PETSc.Vec or None
            RHS vector (typically None for TS problems).
        x : PETSc.Vec
            On entry: current state X_n (initial guess from TS).
            On exit: converged state X_{n+1}.
        """
        adflow = self.solver.adflow
        comm = self.ts_wrapper.comm
        tw = self.ts_wrapper

        # Set ADflow state to the PETSc initial guess (X_n)
        x_arr = x.getArray(readonly=True).copy()
        self.solver.setStates(x_arr)

        # The TS time is at the START of this step; the new time level
        # is t_new = t_current + dt.
        ts = tw.ts
        t_new = ts.getTime() + ts.getTimeStep()

        # Advance Fortran time counters (matches _solverunsteady path)
        adflow.monitor.timestepunsteady += 1
        adflow.monitor.timeunsteady = t_new

        # Shift coordinate/volume history and update geometry.
        # This must happen EVERY step (not just for grid motion) —
        # ADflow's _solverunsteady calls these unconditionally.
        adflow.preprocessingapi.shiftcoorandvolumes()
        adflow.solvers.updateunsteadygeometry()

        # Note: equationMode is already set to _UNSTEADY (=2) by
        # setAeroProblem. Do NOT override it — the Fortran constants are
        # steady=1, unsteady=2, timeSpectral=3.

        # solverUnsteadyStep internally calls:
        #   shiftSolution         — w -> wOld(1), saves X_n
        #   setCoefTimeIntegrator — sets BDF1 coefficients
        #   solveState            — DADI+MG (+ ANK/NK) iterations
        #   nOldSolAvail += 1
        adflow.solvers.solverunsteadystep()

        self._step_count += 1

        # Read converged state back into PETSc vector
        w_conv = self.solver.getStates()
        x.getArray()[:] = w_conv

        if comm.rank == 0:
            print(f"  [DADI-SNES] step {self._step_count} | t = {t_new:.6e}")

        # Signal convergence to PETSc TS
        snes.setConvergedReason(PETSc.SNES.ConvergedReason.CONVERGED_FNORM_RELATIVE)


class ADflowTS:
    """
    Wraps ADflow's spatial residual with PETSc TSTHETA for time-accurate
    simulation and discrete adjoint via ``TSAdjointSolve``.

    Two forward solver modes are available:

    - ``snes_type="newton"`` (default): PETSc Newton--Krylov SNES with
      IFunction/IJacobian.  Supports arbitrary theta but struggles to
      converge at physical time-step sizes due to the small temporal shift.
    - ``snes_type="dadi"``: ADflow's proven DADI + multigrid (+ optional
      ANK/NK) solver drives the forward step.  Forces ``theta=1`` (BDF1).
      PETSc TS manages the trajectory; the adjoint uses the same
      IJacobian (``a*I - dR/dw``) via ``KSPSolveTranspose``.
      No differentiation through DADI/MG is needed — the adjoint relies
      on the implicit function theorem and only requires ``dR/dw``.

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
        Forced to 1.0 when ``snes_type="dadi"``.
    grid_motion : bool, optional
        Whether the mesh moves (prescribed motion). Default False.
    snes_type : str, optional
        ``"newton"`` (PETSc Newton SNES) or ``"dadi"`` (ADflow native
        DADI+MG).  Default ``"newton"``.
    snes_rtol : float, optional
        Relative tolerance for SNES (Newton mode only). Default 1e-8.
    snes_max_it : int, optional
        Maximum SNES iterations per time step (Newton mode). Default 200.
    ksp_type : str, optional
        KSP solver type (Newton mode). Default "preonly".
    ksp_rtol : float, optional
        Relative tolerance for KSP. In DADI mode this applies to the
        adjoint linear solve. Default 1e-4.
    ksp_max_it : int, optional
        Maximum KSP iterations. Default 200.
    ksp_gmres_restart : int, optional
        GMRES restart parameter. Default 200.
    pc_shift_factor : float, optional
        Factor by which the PC shift is multiplied relative to the true
        shift ``a = 1/(theta*dt)``.  Default 1.0.
    """

    def __init__(
        self, solver, ap, dt, t_final, theta=0.5, grid_motion=False,
        snes_rtol=1e-8, snes_max_it=200,
        ksp_type="preonly", ksp_rtol=1e-4, ksp_max_it=200,
        ksp_gmres_restart=200,
        jac_type="ad", fd_step=1e-6, pc_type="asm_ilu",
        pc_shift_factor=1.0, save_trajectory=False,
        snes_type="newton",
    ):
        self.solver = solver
        self.ap = ap
        self.dt = dt
        self.t_final = t_final
        self.theta = theta
        self.grid_motion = grid_motion
        self.comm = solver.comm

        # Solver tolerances
        self.snes_rtol = snes_rtol
        self.snes_max_it = snes_max_it
        self.ksp_type = ksp_type  # "preonly" (PC only) or "gmres"
        self.ksp_rtol = ksp_rtol
        self.ksp_max_it = ksp_max_it
        self.ksp_gmres_restart = ksp_gmres_restart
        self.pc_type = pc_type  # "asm_ilu" (Fortran PC) or "none"

        # The physical shift a = 1/(theta*dt) can be much smaller than the
        # spatial Jacobian spectral radius, making (a*I - dR/dw) poorly
        # conditioned for ILU.  Multiply the PC shift by this factor so
        # that the assembled PC is well-conditioned.  The MatShell operator
        # keeps the exact shift; GMRES resolves the spectral difference.
        self.pc_shift_factor = pc_shift_factor
        self.save_trajectory = save_trajectory
        self.snes_type = snes_type

        # DADI mode requires theta=1.0 (backward Euler = BDF1)
        if snes_type == "dadi" and self.theta != 1.0:
            if solver.comm.rank == 0:
                print(f"  [ADflowTS] snes_type='dadi' requires theta=1.0. "
                      f"Overriding theta={self.theta} -> 1.0")
            self.theta = 1.0

        # Local state size on this processor
        self.n_local = solver.getStateSize()

        # Storage for force coefficients at each time step
        self.time_history = []
        self.cl_history = []
        self.cd_history = []
        self.cmz_history = []

        # State trajectory for manual backward sweep (saved when save_trajectory=True)
        # List of (t_n, w_n) pairs: t_n = time, w_n = local state numpy array
        self._state_trajectory = []

        # Jacobian context.
        # "ad" (default): uses AD from master_d with turbResScale correction.
        #   Consistent with IFunction which also undoes turbResScale.
        # "fd": uses FD of getRes() for debugging/consistency checks.
        # "master": uses AD from master_d WITHOUT turbResScale correction.
        #   Only valid if IFunction also keeps turbResScale (NOT current default).
        if jac_type == "ad":
            self.jac_ctx = ADflowJacobianAD(solver)
        elif jac_type == "fd":
            self.jac_ctx = ADflowJacobianFD(solver, h_fd=fd_step)
        else:  # "master" (default)
            self.jac_ctx = ADflowJacobianMaster(solver)

        # Track the last time at which the grid was updated
        self._grid_time = None

        # Track residual evaluation count
        self.n_res_eval = 0

        # Adjoint outputs (filled by solve_adjoint)
        self.adjoint_init = {}

        # Flag to indicate adjoint backward sweep is active.
        # During adjoint, _update_grid must NOT call shiftCoorAndVolumes
        # (which shifts the coordinate history stack for forward stepping).
        self._adjoint_mode = False

    def setup(self, skip_set_ap=False):
        """Create PETSc TS, vectors, and MatShell.

        Parameters
        ----------
        skip_set_ap : bool
            If True, skip the setAeroProblem call. Use this when the
            solver already has the correct AeroProblem and state (e.g.
            after a warmup BDF step).
        """
        comm = self.comm

        # Ensure the AeroProblem is set so BCs and state are initialized
        if not skip_set_ap:
            self.solver.setAeroProblem(self.ap)

        # Create PETSc vectors (local size, global auto-determined)
        self.F_vec = PETSc.Vec().createMPI((self.n_local, None), comm=comm)
        self.U_vec = PETSc.Vec().createMPI((self.n_local, None), comm=comm)

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
        if self.save_trajectory:
            self.ts.setSaveTrajectory()

        # Disable time step adaptation — use fixed dt
        PETSc.Options().setValue("-ts_adapt_type", "none")

        # -----------------------------------------------------------
        # SNES / KSP / PC setup (depends on snes_type)
        # -----------------------------------------------------------
        if self.snes_type == "newton":
            # Newton-based SNES with backtracking line search
            PETSc.Options().setValue("-snes_linesearch_type", "bt")
            PETSc.Options().setValue("-snes_linesearch_monitor", "")
            snes = self.ts.getSNES()
            snes.setTolerances(
                rtol=self.snes_rtol, atol=1e-14, max_it=self.snes_max_it,
            )

            # KSP (linear solver) settings
            ksp = snes.getKSP()
            ksp.setType(self.ksp_type)
            ksp.setPCSide(PETSc.PC.Side.RIGHT)
            if self.ksp_type == "gmres":
                ksp.setTolerances(rtol=self.ksp_rtol, max_it=self.ksp_max_it)
                ksp.setGMRESRestart(self.ksp_gmres_restart)
            else:
                ksp.setTolerances(max_it=1)

            # Preconditioner setup
            pc = ksp.getPC()
            if self.pc_type == "none":
                pc.setType("none")
            elif self.pc_type == "scale":
                pc.setType("python")
                pc.setPythonContext(_ScalePreconditioner(self))
            elif self.pc_type == "native_adjoint":
                pc.setType("python")
                pc.setPythonContext(_NativeAdjointPC(self))
            else:
                pc.setType("python")
                pc.setPythonContext(_TSPreconditioner(self))
            self._pc_assembled = False

        elif self.snes_type == "dadi":
            # -----------------------------------------------------------
            # DADI mode: set Python SNES type BEFORE setFromOptions so
            # PETSc's internal TSSetUp/SNESSetUp sees the correct type
            # from the start — avoids SEGV from corrupted ops table when
            # the type is changed after the Newton SNES was fully set up.
            # -----------------------------------------------------------
            snes = self.ts.getSNES()
            snes.setType(PETSc.SNES.Type.PYTHON)
            snes.setPythonContext(ADflowDADISNES(self))

            # KSP for adjoint linear system  K^T * lambda = rhs
            # (used by TSAdjointStep_Theta via KSPSolveTranspose)
            ksp = snes.getKSP()
            ksp.setType(self.ksp_type)
            ksp.setPCSide(PETSc.PC.Side.RIGHT)
            ksp.setTolerances(rtol=self.ksp_rtol, max_it=self.ksp_max_it)
            if self.ksp_type == "gmres":
                ksp.setGMRESRestart(self.ksp_gmres_restart)

            # PC for adjoint
            pc = ksp.getPC()
            if self.pc_type == "none":
                pc.setType("none")
            elif self.pc_type == "scale":
                pc.setType("python")
                pc.setPythonContext(_ScalePreconditioner(self))
            elif self.pc_type == "native_adjoint":
                pc.setType("python")
                pc.setPythonContext(_NativeAdjointPC(self))
            else:
                pc.setType("python")
                pc.setPythonContext(_TSPreconditioner(self))
            self._pc_assembled = False

        # Set post-step monitor to record force coefficients
        self.ts.setPostStep(self._post_step)

        # Allow SNES failures so TS can retry
        PETSc.Options().setValue("-ts_max_snes_failures", "-1")

        # Allow command-line options to override all of the above.
        # For DADI mode, the Python SNES type is already set; as long
        # as -snes_type is not passed on the command line, it persists.
        self.ts.setFromOptions()

        # Re-set dt after setFromOptions to ensure it's not overridden
        self.ts.setTimeStep(self.dt)

    def _normalize_obj_name(self, obj):
        """Accept both 'cl' and '<ap.name>_cl' and return bare objective name."""
        key = obj.lower()
        prefix = f"{self.ap.name.lower()}_"
        if key.startswith(prefix):
            return key[len(prefix):]
        return key

    def _compute_terminal_state_gradient(self, obj):
        """Return d(obj)/dw at the current solver state."""
        obj_name = self._normalize_obj_name(obj)
        funcs_bar = self.solver._getFuncsBar(obj_name)
        return self.solver.computeJacobianVectorProductBwd(
            funcsBar=funcs_bar, wDeriv=True
        )

    def solve_adjoint(self, objectives):
        """Run discrete TS adjoint for final-time objectives.

        Parameters
        ----------
        objectives : iterable[str]
            Objective names, e.g. ``["cl"]`` or ``["0012pitching_cl", "cd"]``.

        Returns
        -------
        dict
            Mapping ``objective -> lambda(t0)`` (local vectors).
        """
        if not self.save_trajectory:
            raise RuntimeError(
                "TS trajectory is disabled. Recreate ADflowTS with "
                "save_trajectory=True before running the forward solve."
            )
        if self.ts.getStepNumber() <= 0:
            raise RuntimeError("Forward solve must complete before solve_adjoint().")
        if self.ksp_type == "preonly":
            if self.comm.rank == 0:
                print(
                    "  [Adjoint] WARNING: forward used ksp_type='preonly'. "
                    "For robust adjoint accuracy, prefer ksp_type='gmres'."
                )

        # Ensure ADflow state is synchronized to TS final state.
        u_final = self.U_vec.getArray(readonly=True).copy()
        self.solver.setStates(u_final)

        lambda_vecs = []
        clean_names = []
        for obj in objectives:
            obj_name = self._normalize_obj_name(obj)
            dphi_dw = self._compute_terminal_state_gradient(obj_name)
            lam = self.U_vec.duplicate()
            lam_arr = lam.getArray()
            lam_arr[:] = dphi_dw
            lambda_vecs.append(lam)
            clean_names.append(obj_name)

        # No explicit parameter Jacobian — pass None for mu vectors.
        self.ts.setCostGradients(lambda_vecs, None)

        # Enable adjoint mode so _update_grid skips shiftCoorAndVolumes
        self._adjoint_mode = True
        self._grid_time = None  # force grid update on first adjoint step
        try:
            self.ts.adjointSolve()
        finally:
            self._adjoint_mode = False

        self.adjoint_init = {}
        for i, obj_name in enumerate(clean_names):
            self.adjoint_init[obj_name] = lambda_vecs[i].getArray(readonly=True).copy()

        if self.comm.rank == 0:
            print("=" * 70)
            print("ADflowTS: Adjoint sweep complete")
            for obj_name in clean_names:
                nrm = np.linalg.norm(self.adjoint_init[obj_name])
                print(f"  objective={obj_name:>12s} | ||lambda(t0)||_2 = {nrm:.6e}")
            print("=" * 70)

        return self.adjoint_init

    def _jtranspose_matvec(self, v, shift):
        """Compute y = J_u^T * v = shift*v - (dR_u/dw)^T * v.

        This is the unscaled adjoint Jacobian-vector product. The AD routine
        ``computeJacobianVectorProductBwd`` returns ``(S * dR_u/dw)^T * resBar``
        where S = diag(resScale). To get ``(dR_u/dw)^T * v``, we pass
        ``resBar = S^{-1} * v`` and read back wDeriv directly.

        Parameters
        ----------
        v : ndarray
            Input vector (local part).
        shift : float
            Temporal shift ``a = 1/(theta*dt)``.

        Returns
        -------
        y : ndarray
            ``J_u^T * v = shift*v - (dR_u/dw)^T * v``
        """
        import time as _time
        t0 = _time.time()
        if self.comm.rank == 0:
            print(f"      [Jt-matvec] enter, shift={shift:.4e}, ||v||={np.linalg.norm(v):.4e}", flush=True)

        adflow = self.solver.adflow
        orig_mode = adflow.inputphysics.equationmode
        adflow.inputphysics.equationmode = _STEADY

        # Apply inverse turbResScale to SA DOFs of resBar
        jac_ctx = self.jac_ctx
        if hasattr(jac_ctx, '_apply_inv_turb_res_scale_to_resbar'):
            v_adj = jac_ctx._apply_inv_turb_res_scale_to_resbar(v)
        else:
            v_adj = v.copy()

        wbar = self.solver.computeJacobianVectorProductBwd(
            resBar=v_adj, wDeriv=True
        )

        adflow.inputphysics.equationmode = orig_mode
        y = shift * v - wbar
        if self.comm.rank == 0:
            print(f"      [Jt-matvec] done in {_time.time()-t0:.2f}s, ||y||={np.linalg.norm(y):.4e}", flush=True)
        return y

    def _adjoint_precond(self, v, resscale_diag, inner_tol):
        """Apply preconditioner P^{-1} * v for the adjoint linear system.

        Uses ADflow's native adjoint solver: P = -(S * dR_u/dw)^T.
        ``solveAdjointForRHS`` inverts (S * dR_u/dw)^T, and we post-multiply
        by S to get (dR_u/dw)^T inverse, then negate.

        Parameters
        ----------
        v : ndarray
            Input vector (local part).
        resscale_diag : ndarray
            Diagonal of resScale (turbResScale for SA DOFs, 1 for flow).
        inner_tol : float
            Relative tolerance for the inner DADI+MG solve.

        Returns
        -------
        y : ndarray
            ``-(dR_u/dw)^{-T} * v``
        """
        import time as _time
        t0 = _time.time()
        if self.comm.rank == 0:
            print(f"      [PC-apply] enter, inner_tol={inner_tol}, ||v||={np.linalg.norm(v):.4e}", flush=True)
        z = self.solver.adflow.adjointapi.solveadjointforrhs(v, inner_tol)
        y = -(resscale_diag * z)
        if self.comm.rank == 0:
            print(f"      [PC-apply] done in {_time.time()-t0:.2f}s, ||z||={np.linalg.norm(z):.4e}, ||y||={np.linalg.norm(y):.4e}", flush=True)
        return y

    def _adjoint_precond_ts(self, v, use_transpose=True):
        """Apply TS preconditioner for the adjoint linear system.

        Uses the ILU factorization of ``shift*I - dR/dw`` assembled by
        ``setupTSPreconditioner(shift)``.  This includes the temporal shift,
        making it a much better approximation to the true Jacobian
        ``J = shift*I - dR/dw`` than the native adjoint PC (which only
        approximates ``dR/dw``).

        Must call ``setupTSPreconditioner(shift)`` before first use.

        Parameters
        ----------
        v : ndarray
            Input vector (local part).
        use_transpose : bool
            If True (default), use ``KSPSolveTranspose`` (exact transpose).
            If False, use ``KSPSolve`` (forward ILU as approximate transpose,
            for debugging when the transpose solve is broken).

        Returns
        -------
        y : ndarray
            Approximate ``(shift*I - dR/dw)^{-T} * v`` via ILU.
        """
        import time as _time
        t0 = _time.time()
        tag = "TS-PC-T" if use_transpose else "TS-PC-fwd"
        if self.comm.rank == 0:
            print(f"      [{tag}] enter, ||v||={np.linalg.norm(v):.4e}", flush=True)
        x_arr = v.copy()
        y_arr = np.zeros(self.n_local)
        if use_transpose:
            self.solver.adflow.nksolver.applytspreconditionertranspose(x_arr, y_arr)
        else:
            self.solver.adflow.nksolver.applytspreconditioner(x_arr, y_arr)
        if self.comm.rank == 0:
            print(f"      [{tag}] done in {_time.time()-t0:.3f}s, ||y||={np.linalg.norm(y_arr):.4e}", flush=True)
        return y_arr

    @staticmethod
    def _fgmres_solve(matvec, precond, b, n_local, comm,
                      rtol=1e-10, max_it=50, restart=30, verbose=False):
        """Flexible right-preconditioned GMRES (FGMRES) for solving A*x = b.

        Unlike standard GMRES, FGMRES stores the preconditioned vectors
        Z[j] = M^{-1}*V[j] explicitly. This handles variable/inexact
        preconditioners correctly (e.g. inner iterative solves that give
        slightly different results each time).

        Parameters
        ----------
        matvec : callable(v) -> ndarray
            Matrix-vector product ``A * v``.
        precond : callable(v) -> ndarray
            Preconditioner application ``M^{-1} * v``.
        b : ndarray
            Right-hand side (local part).
        n_local : int
            Local vector size.
        comm : MPI communicator
            For global dot products.
        rtol : float
            Relative tolerance for residual reduction.
        max_it : int
            Maximum total iterations (across restarts).
        restart : int
            FGMRES restart (number of Arnoldi vectors per cycle).
        verbose : bool
            Print convergence info on rank 0.

        Returns
        -------
        x : ndarray
            Approximate solution (local part).
        converged : bool
            Whether the relative tolerance was achieved.
        n_iter : int
            Total number of iterations.
        """
        def global_dot(a, b_vec):
            return comm.allreduce(np.dot(a, b_vec))

        def global_norm(a):
            return np.sqrt(global_dot(a, a))

        x = np.zeros(n_local)
        r = b.copy()
        b_norm = global_norm(b)
        if b_norm < 1e-30:
            return x, True, 0

        total_iter = 0
        converged = False

        while total_iter < max_it:
            r_norm = global_norm(r)
            if r_norm / b_norm < rtol:
                converged = True
                break

            m = min(restart, max_it - total_iter)

            # Arnoldi basis V[0..m], preconditioned vectors Z[0..m-1]
            V = [None] * (m + 1)
            Z = [None] * m  # Z[j] = M^{-1} * V[j] (stored for FGMRES)
            H = np.zeros((m + 1, m))
            V[0] = r / r_norm

            # Givens rotation arrays
            cs = np.zeros(m)
            sn = np.zeros(m)
            g = np.zeros(m + 1)
            g[0] = r_norm

            for j in range(m):
                # FGMRES: store preconditioned vector
                Z[j] = precond(V[j])
                w = matvec(Z[j])

                # Compute debug norms on ALL ranks (allreduce is collective),
                # but only print on rank 0.
                if verbose and total_iter < 3:
                    nV = global_norm(V[j])
                    nZ = global_norm(Z[j])
                    nW = global_norm(w)
                    if comm.rank == 0:
                        print(f"    [debug] j={j}: ||V[j]||={nV:.4e}, "
                              f"||Z[j]||={nZ:.4e}, "
                              f"||A*Z[j]||={nW:.4e}", flush=True)

                # Modified Gram-Schmidt
                for i in range(j + 1):
                    H[i, j] = global_dot(w, V[i])
                    w = w - H[i, j] * V[i]

                H[j + 1, j] = global_norm(w)
                if H[j + 1, j] > 1e-30:
                    V[j + 1] = w / H[j + 1, j]
                else:
                    V[j + 1] = np.zeros(n_local)

                if verbose and total_iter < 3 and comm.rank == 0:
                    print(f"    [debug] j={j}: H[{j},{j}]={H[j,j]:.4e}, "
                          f"H[{j+1},{j}]={H[j+1,j]:.4e}", flush=True)

                # Apply previous Givens rotations to column j
                for i in range(j):
                    temp = cs[i] * H[i, j] + sn[i] * H[i + 1, j]
                    H[i + 1, j] = -sn[i] * H[i, j] + cs[i] * H[i + 1, j]
                    H[i, j] = temp

                # Compute new Givens rotation
                denom = np.sqrt(H[j, j] ** 2 + H[j + 1, j] ** 2)
                if denom > 1e-30:
                    cs[j] = H[j, j] / denom
                    sn[j] = H[j + 1, j] / denom
                else:
                    cs[j] = 1.0
                    sn[j] = 0.0

                # Apply Givens rotation
                H[j, j] = cs[j] * H[j, j] + sn[j] * H[j + 1, j]
                H[j + 1, j] = 0.0
                g[j + 1] = -sn[j] * g[j]
                g[j] = cs[j] * g[j]

                total_iter += 1
                res_est = abs(g[j + 1])

                if verbose and comm.rank == 0:
                    print(f"    FGMRES iter {total_iter:3d} | "
                          f"|r| = {res_est:.4e} | "
                          f"|r|/|b| = {res_est/b_norm:.4e}")

                if res_est / b_norm < rtol:
                    # Converged — solve triangular system and recover x
                    y = np.zeros(j + 1)
                    for i in range(j, -1, -1):
                        y[i] = g[i]
                        for k in range(i + 1, j + 1):
                            y[i] -= H[i, k] * y[k]
                        y[i] /= H[i, i]

                    # FGMRES: use stored Z vectors (no re-application of precond)
                    for i in range(j + 1):
                        x += y[i] * Z[i]

                    converged = True
                    break

            if converged:
                break

            # Not converged within this restart cycle — solve and update
            j_last = min(m, total_iter) - 1
            if j_last < 0:
                break
            y = np.zeros(j_last + 1)
            for i in range(j_last, -1, -1):
                y[i] = g[i]
                for k in range(i + 1, j_last + 1):
                    y[i] -= H[i, k] * y[k]
                y[i] /= H[i, i]

            # FGMRES: use stored Z vectors
            for i in range(j_last + 1):
                x += y[i] * Z[i]

            # Update residual
            r = b - matvec(x)

        return x, converged, total_iter

    def solve_adjoint_manual(self, objectives, inner_tol=0.01,
                             gmres_rtol=1e-10, gmres_max_it=200,
                             gmres_restart=200, reassemble=False,
                             adjoint_pc="ts_ilu", pc_shift_factor=1.0,
                             dv_sens=None):
        """Run discrete adjoint via manual backward sweep with GMRES.

        Avoids PETSc's internal KSPSolveTranspose (which causes SEGV with
        nested native-adjoint preconditioning) by implementing the BEuler
        adjoint update directly in Python with right-preconditioned GMRES.

        For backward Euler (theta=1), each adjoint step solves::

            J_u^T * delta = lambda / dt

        where ``J_u^T = (1/dt)*I - (dR_u/dw)^T`` (unscaled Jacobian).
        Then ``lambda_new = delta``.

        The GMRES uses:
        - **Matvec**: ``J_u^T * v = shift*v - (dR_u/dw)^T * v``
          via AD reverse mode (same as MatShell multTranspose).
        - **Preconditioner** (depends on ``adjoint_pc``):

          - ``"ts_ilu"`` (default): ILU factorization of
            ``pc_shift*I - dR/dw`` via ``setupTSPreconditioner``.
            Uses forward ILU solve (not transpose, which is broken
            in PETSc 3.18 with ASM+ILU).  ``pc_shift`` may be larger
            than the physical shift (controlled by ``pc_shift_factor``)
            to ensure ILU stability on the indefinite Jacobian.
          - ``"native"``: ADflow's DADI+MG adjoint solver for ``dR/dw^T``.
            Does NOT include the shift, so the preconditioned operator has
            eigenvalues near ``1 - shift/eig_i``, which is poor when
            ``eig_i ~ shift``.  Slow per-application (~0.2s).

        The total derivative for a final-time objective is::

            dJ/dalpha = pJ/palpha
                      + dt * SUM_{n=1}^{N} delta_n^T * pR_n/palpha
                      + lambda_0^T * dw_0/dalpha

        where the first two terms are computed when ``dv_sens`` is provided,
        and ``lambda_0^T * dw_0/dalpha`` (the initial-condition contribution)
        must be computed externally.

        Parameters
        ----------
        objectives : list[str]
            Objective names, e.g. ``["cl"]``.
        inner_tol : float
            Relative tolerance for the inner DADI+MG preconditioner
            (only used when ``adjoint_pc="native"``). Default 0.01.
        gmres_rtol : float
            Relative tolerance for the outer GMRES (default 1e-10).
        gmres_max_it : int
            Maximum GMRES iterations per adjoint step (default 200).
        gmres_restart : int
            GMRES restart parameter (default 200).
        reassemble : bool
            If True, reassemble dRdwT at each time step. Default False
            (uses the matrix from the steady state for all steps).
        adjoint_pc : str
            Preconditioner type: ``"ts_ilu"`` (default) or ``"native"``.
        pc_shift_factor : float
            Multiplier for the shift used in PC assembly (default 1.0).
            When the physical shift is too small for ILU stability
            (matrix indefinite), increase this (e.g. 10-100) so the
            PC matrix ``pc_shift*I - dR/dw`` is diagonally dominant.
        dv_sens : list[str] or None
            If provided, accumulate parameter sensitivities during the
            backward sweep.  Returns ``result[obj + "_dJdp"]`` as a dict
            mapping DV names to ``pJ/palpha + dt * SUM delta_n^T * pR/palpha``
            (excludes the initial-condition term ``lambda_0^T * dw0/dalpha``).

        Returns
        -------
        dict
            Mapping ``objective -> lambda(t0)`` (local numpy arrays).
            If ``dv_sens`` is not None, also contains
            ``objective + "_dJdp" -> {dv_name: float}``.
        """
        if not self.save_trajectory or len(self._state_trajectory) < 2:
            raise RuntimeError(
                "State trajectory not available. Run forward solve with "
                "save_trajectory=True first."
            )

        N = len(self._state_trajectory) - 1  # number of time steps
        shift = 1.0 / (self.theta * self.dt)
        comm = self.comm

        solver = self.solver

        if adjoint_pc == "native":
            # --- Set up native adjoint solver (once) ---
            solver._setupAdjoint()
            # --- Build resScale diagonal for PC scaling correction ---
            jac_ctx = self.jac_ctx
            resscale_diag = np.ones(self.n_local)
            if hasattr(jac_ctx, '_n_turb') and jac_ctx._n_turb > 0:
                n_cells = self.n_local // jac_ctx._nw
                rs2d = resscale_diag.reshape(n_cells, jac_ctx._nw)
                for l in range(jac_ctx._n_turb):
                    rs2d[:, jac_ctx._nwf + l] = jac_ctx._trs[l]

        pc_shift = shift * pc_shift_factor

        if comm.rank == 0:
            print("=" * 70)
            print(f"ADflowTS: Manual backward sweep (GMRES, rtol={gmres_rtol})")
            print(f"  {N} adjoint steps, shift = {shift:.4e}")
            print(f"  adjoint_pc = {adjoint_pc}, GMRES max_it = {gmres_max_it}")
            if adjoint_pc == "ts_ilu":
                print(f"  pc_shift_factor = {pc_shift_factor}, pc_shift = {pc_shift:.4e}")
            if adjoint_pc == "native":
                print(f"  inner PC tol = {inner_tol}")
            print(f"  reassemble = {reassemble}")
            print("=" * 70)

        # Enable adjoint mode so _update_grid skips shiftCoorAndVolumes
        self._adjoint_mode = True
        self._grid_time = None
        try:
            # --- Terminal cost gradients ---
            t_N, w_N = self._state_trajectory[N]
            solver.setStates(w_N)
            self._update_grid(t_N)

            results = {}
            for obj in objectives:
                obj_name = self._normalize_obj_name(obj)
                dphi_dw = self._compute_terminal_state_gradient(obj_name)

                # --- Terminal partial pJ/palpha (if dv_sens requested) ---
                if dv_sens is not None:
                    adflow = solver.adflow
                    orig_mode = adflow.inputphysics.equationmode
                    adflow.inputphysics.equationmode = _STEADY
                    funcsBar = solver._getFuncsBar(obj_name)
                    term_sens = solver.computeJacobianVectorProductBwd(
                        funcsBar=funcsBar, xDvDeriv=True
                    )
                    adflow.inputphysics.equationmode = orig_mode
                    dJdp = {key: float(term_sens.get(key, 0.0)) for key in dv_sens}
                    if comm.rank == 0:
                        print(f"  Terminal partials (pJ/p):")
                        for key in dv_sens:
                            print(f"    {key}: {dJdp[key]:.10e}")

                # --- Backward sweep ---
                lam = dphi_dw.copy()  # lambda_N = dCL/dU_N

                for n in range(N, 0, -1):
                    # Restore forward state at step n
                    t_n, w_n = self._state_trajectory[n]
                    if comm.rank == 0:
                        print(f"\n  --- Adjoint step {N-n+1}/{N}: restoring state at t={t_n:.6e} ---", flush=True)
                    solver.setStates(w_n)
                    self._update_grid(t_n)

                    import time as _time

                    if reassemble or n == N:
                        t_asm = _time.time()
                        if adjoint_pc == "native":
                            if comm.rank == 0:
                                print(f"    Reassembling dR/dw matrices (native)...", flush=True)
                            solver.adflow.adjointapi.setupallresidualmatricesfwd()
                        if comm.rank == 0:
                            print(f"    Assembling TS PC (pc_shift={pc_shift:.4e})...", flush=True)
                        adflow = solver.adflow
                        orig_mode = adflow.inputphysics.equationmode
                        adflow.inputphysics.equationmode = _STEADY
                        adflow.nksolver.setuptspreconditioner(pc_shift)
                        adflow.inputphysics.equationmode = orig_mode
                        if comm.rank == 0:
                            print(f"    Assembly done in {_time.time()-t_asm:.2f}s", flush=True)

                        # Diagnostic: compare forward vs transpose ILU on first step
                        if n == N and adjoint_pc == "ts_ilu":
                            test_v = shift * lam
                            tv_norm = np.sqrt(comm.allreduce(np.dot(test_v, test_v)))
                            y_fwd = self._adjoint_precond_ts(test_v, use_transpose=False)
                            y_trn = self._adjoint_precond_ts(test_v, use_transpose=True)
                            yfn = np.sqrt(comm.allreduce(np.dot(y_fwd, y_fwd)))
                            ytn = np.sqrt(comm.allreduce(np.dot(y_trn, y_trn)))
                            if comm.rank == 0:
                                print(f"    [DIAG] ||test_v||={tv_norm:.4e}")
                                print(f"    [DIAG] ILU forward:   ||y||={yfn:.4e}  ratio={yfn/tv_norm:.4e}")
                                print(f"    [DIAG] ILU transpose: ||y||={ytn:.4e}  ratio={ytn/tv_norm:.4e}")
                                print(f"    [DIAG] If transpose >> forward, KSPSolveTranspose is broken", flush=True)

                    # RHS = shift * lambda (= lambda / dt)
                    rhs = shift * lam

                    # Solve J_u^T * delta = rhs via right-preconditioned GMRES
                    def matvec(v):
                        return self._jtranspose_matvec(v, shift)

                    if adjoint_pc == "ts_ilu":
                        def precond(v):
                            # use_transpose=False: forward ILU as approximate
                            # transpose PC (KSPSolveTranspose+ASM+ILU is broken
                            # in PETSc 3.18, producing O(1e23) norms)
                            return self._adjoint_precond_ts(v, use_transpose=False)
                    else:
                        def precond(v):
                            return self._adjoint_precond(v, resscale_diag, inner_tol)

                    import time as _time
                    t_step_start = _time.time()
                    if comm.rank == 0:
                        rhs_norm = np.linalg.norm(rhs)
                        print(f"    [step {N-n+1}/{N}] FGMRES start, ||rhs||={rhs_norm:.4e}, "
                              f"||lam||={np.linalg.norm(lam):.4e}", flush=True)

                    delta, conv, n_iter = self._fgmres_solve(
                        matvec, precond, rhs,
                        self.n_local, comm,
                        rtol=gmres_rtol, max_it=gmres_max_it,
                        restart=gmres_restart, verbose=True,
                    )

                    # For BEuler: lambda_{n-1} = delta
                    lam = delta

                    # Accumulate intermediate parameter sensitivity:
                    # dt * delta_n^T * pR_n/palpha
                    if dv_sens is not None:
                        # Re-set state (GMRES matvec may have modified internals)
                        solver.setStates(w_n)
                        adflow = solver.adflow
                        orig_mode = adflow.inputphysics.equationmode
                        adflow.inputphysics.equationmode = _STEADY
                        jac_ctx = self.jac_ctx
                        if hasattr(jac_ctx, '_apply_inv_turb_res_scale_to_resbar'):
                            d_adj = jac_ctx._apply_inv_turb_res_scale_to_resbar(delta)
                        else:
                            d_adj = delta.copy()
                        step_sens = solver.computeJacobianVectorProductBwd(
                            resBar=d_adj, xDvDeriv=True
                        )
                        adflow.inputphysics.equationmode = orig_mode
                        for key in dv_sens:
                            val = float(step_sens.get(key, 0.0))
                            dJdp[key] += self.dt * val
                        if comm.rank == 0 and n == N:
                            d_norm = np.linalg.norm(d_adj)
                            print(f"    [dv_sens debug] step {N-n+1}: ||d_adj||={d_norm:.4e}")
                            print(f"    [dv_sens debug] step_sens = {step_sens}")
                            print(f"    [dv_sens debug] step_sens type = {type(step_sens)}")
                            # Test with ones vector (same as diagnostic)
                            solver.setStates(w_n)
                            adflow.inputphysics.equationmode = _STEADY
                            ones_test = np.ones(self.n_local)
                            test_ones = solver.computeJacobianVectorProductBwd(
                                resBar=ones_test, xDvDeriv=True
                            )
                            adflow.inputphysics.equationmode = orig_mode
                            print(f"    [dv_sens debug] resBar=ones → alpha = {test_ones.get('alpha', 0.0)}")
                            # Also try forward mode
                            solver.setStates(w_n)
                            adflow.inputphysics.equationmode = _STEADY
                            res_fwd = solver.computeJacobianVectorProductFwd(
                                xDvDot={"alpha": 1.0}, residualDeriv=True
                            )
                            adflow.inputphysics.equationmode = orig_mode
                            print(f"    [dv_sens debug] ||dR/dalpha|| fwd = {np.linalg.norm(res_fwd):.6e}")
                            # Check dot product
                            dot_val = np.dot(d_adj, res_fwd)
                            dot_val_g = comm.allreduce(dot_val)
                            print(f"    [dv_sens debug] d_adj . dR/dalpha = {dot_val_g:.10e}")

                    lam_norm_sq = comm.allreduce(np.dot(lam, lam))
                    if comm.rank == 0:
                        conv_str = "OK" if conv else "FAIL"
                        t_step = _time.time() - t_step_start
                        print(f"  Adjoint step {N-n+1:3d}/{N} | t = {t_n:.6e} "
                              f"| ||lambda|| = {np.sqrt(lam_norm_sq):.6e} "
                              f"| GMRES: {n_iter} it ({conv_str}) | {t_step:.1f}s",
                              flush=True)

                results[obj_name] = lam.copy()
                if dv_sens is not None:
                    results[obj_name + "_dJdp"] = dJdp.copy()
        finally:
            self._adjoint_mode = False

        norms = {}
        for obj_name in results:
            if obj_name.endswith("_dJdp"):
                continue
            norms[obj_name] = np.sqrt(
                comm.allreduce(np.dot(results[obj_name], results[obj_name]))
            )
        if comm.rank == 0:
            print("=" * 70)
            print("ADflowTS: Manual adjoint sweep complete")
            for obj_name in norms:
                print(f"  objective={obj_name:>12s} | ||lambda(t0)||_2 = "
                      f"{norms[obj_name]:.6e}")
            if dv_sens is not None:
                for obj in objectives:
                    on = self._normalize_obj_name(obj)
                    djdp = results.get(on + "_dJdp", {})
                    for key in dv_sens:
                        print(f"  {on}_dJdp[{key}] = {djdp.get(key, 0.0):.10e}"
                              f"  (excl IC term)")
            print("=" * 70)

        self.adjoint_init = results
        return results

    def _update_grid(self, t):
        """
        Update the mesh coordinates and velocities for time t.

        This sets the internal timeUnsteady variable and calls ADflow's
        grid motion routines. Only updates if the time has changed.

        During the adjoint backward sweep (``_adjoint_mode=True``),
        ``shiftCoorAndVolumes`` is skipped because it shifts the
        coordinate history stack (xOld ← x), which is only meaningful
        during forward time integration. The adjoint only needs the mesh
        positioned at the correct time via ``updateUnsteadyGeometry``.
        """
        if not self.grid_motion:
            return
        if self._grid_time is not None and abs(t - self._grid_time) < 1e-15:
            return

        # Set the physical time so that grid motion uses the correct time
        self.solver.adflow.monitor.timeunsteady = t

        if not self._adjoint_mode:
            # Forward: shift coordinate history, then update mesh
            self.solver.adflow.preprocessingapi.shiftcoorandvolumes()
        # Position the mesh at time t (prescribed motion)
        self.solver.adflow.solvers.updateunsteadygeometry()

        self._grid_time = t

    def _get_spatial_residual(self):
        """
        Evaluate the spatial residual R(w) without temporal terms.

        Uses the ``master`` code path (via ``computeJacobianVectorProductFwd``
        with zero wDot) so the residual is computed through exactly the same
        Fortran routines as the Jacobian. After ``master`` runs (including
        ``resScale``), the primal residual dw is read from block arrays via
        ``getResDw``, then turbResScale is undone from SA DOFs.

        This ensures:
        1. Exact consistency between IFunction and IJacobian (both master path)
        2. No artificial 10000x scaling of SA DOFs (turbResScale removed)

        Returns
        -------
        res : ndarray
            Spatial residual R(w), normalized by 1/volRef only.
        """
        adflow = self.solver.adflow

        # Save and switch equation mode
        orig_mode = adflow.inputphysics.equationmode
        adflow.inputphysics.equationmode = _STEADY

        # Run master (primal) by calling the forward AD with zero tangent.
        # This executes master() which writes dw to block arrays with
        # resScale applied (1/volRef + turbResScale for SA DOFs).
        # Side effect: updates stored volumes/metrics for consistency.
        zero_wdot = np.zeros(self.n_local)
        self.solver.computeJacobianVectorProductFwd(
            wDot=zero_wdot, residualDeriv=True
        )

        # Read the primal residual from block arrays (no recomputation,
        # no additional normalization — dw has resScale from master)
        res = np.zeros(self.n_local)
        res = adflow.nksolver.getresdw(res)

        adflow.inputphysics.equationmode = orig_mode

        # Undo turbResScale from SA DOFs to stay in 1/volRef space
        nw = int(adflow.flowvarrefstate.nw)
        nwf = int(adflow.flowvarrefstate.nwf)
        n_turb = nw - nwf
        if n_turb > 0:
            trs = adflow.inputiteration.turbresscale[:n_turb]
            n_cells = self.n_local // nw
            res_2d = res.reshape(n_cells, nw)
            for l in range(n_turb):
                if abs(trs[l]) > 1e-30:
                    res_2d[:, nwf + l] /= trs[l]

        return res

    def _clip_state(self, u_arr):
        """Clip state variables to physical bounds (in-place).

        ADflow stores [rho, vx, vy, vz, rhoE, nuTilde] per cell.
        Negative density or energy causes NaN in flux computation.
        """
        nw = int(self.solver.adflow.flowvarrefstate.nw)
        n_cells = len(u_arr) // nw
        w2d = u_arr.reshape(n_cells, nw)
        # iRho = 0: density > 0
        np.clip(w2d[:, 0], 1e-10, None, out=w2d[:, 0])
        # iRhoE = 4: total energy > 0
        if nw >= 5:
            np.clip(w2d[:, 4], 1e-10, None, out=w2d[:, 4])
        # SA nuTilde >= 0 (if present)
        nwf = int(self.solver.adflow.flowvarrefstate.nwf)
        for l in range(nwf, nw):
            np.clip(w2d[:, l], 0.0, None, out=w2d[:, l])

    def _ifunction(self, ts, t, U, Udot, F):
        """
        IFunction callback: F = Udot - R(w, t).

        For PETSc TS, the implicit residual is F(t, U, Udot) = 0.
        We have F = Udot - R(U), where R is the spatial residual.
        """
        # Check for NaN/Inf in the input state
        u_arr = U.getArray(readonly=True).copy()
        if np.any(~np.isfinite(u_arr)):
            if self.comm.rank == 0:
                n_nan = np.sum(np.isnan(u_arr))
                n_inf = np.sum(np.isinf(u_arr))
                print(f"  [IFunction] WARNING: U has {n_nan} NaN, {n_inf} Inf")
            F.getArray()[:] = 1e30
            self.n_res_eval += 1
            return

        # Clip unphysical states (negative density/energy) that arise
        # from large Newton updates before they cause NaN in fluxes.
        self._clip_state(u_arr)

        # Update mesh for current time
        self._update_grid(t)

        # Set ADflow state to U
        self.solver.setStates(u_arr)

        # Evaluate spatial residual R(w) (no temporal terms)
        res = self._get_spatial_residual()

        # Check for NaN in residual (can happen near unphysical states)
        if np.any(~np.isfinite(res)):
            if self.comm.rank == 0:
                print(f"  [IFunction] WARNING: R(w) has NaN/Inf")
            F.getArray()[:] = 1e30
            self.n_res_eval += 1
            return

        # F = Udot - R(w)
        f_arr = F.getArray()
        udot_arr = Udot.getArray(readonly=True)
        f_arr[:] = udot_arr - res

        self.n_res_eval += 1

    def _ijacobian(self, ts, t, U, Udot, a, J, P):
        """
        IJacobian callback: J = a*I - dR/dw.

        Stores the shift 'a' and current state for MatShell mult.
        Also assembles the approximate Jacobian for the ILU preconditioner.
        """
        # Update mesh for current time
        self._update_grid(t)

        # Update the shift in the MatShell context
        self.jac_ctx.shift = a
        if self.comm.rank == 0:
            print(f"  [IJacobian] shift a = {a:.6e}  (1/dt = {1.0/self.dt:.6e})")

        # Set ADflow state for Jacobian evaluation
        u_arr = U.getArray(readonly=True).copy()
        self.solver.setStates(u_arr)

        # Assemble the approximate Jacobian P = a_pc*I - dR/dw for the PC.
        # Skip assembly when shift ≈ 0 (PETSc adjoint update step uses
        # shift=0 for MatMultTransposeAdd only, no KSP solve needed).
        # Assembling with shift=0 gives P = -dR/dw which can have zero
        # pivots, causing NaN in the ILU factorization.
        if self.pc_type not in ("none", "scale") and a > 1e-10:
            a_pc = a * self.pc_shift_factor
            adflow = self.solver.adflow
            orig_mode = adflow.inputphysics.equationmode
            adflow.inputphysics.equationmode = _STEADY
            adflow.nksolver.setuptspreconditioner(a_pc)
            adflow.inputphysics.equationmode = orig_mode
            self._pc_assembled = True
            if self.comm.rank == 0:
                print(f"  [PC] shift_pc = {a_pc:.6e}  (factor = {self.pc_shift_factor})")

        # Signal that the matrix structure has not changed
        J.assemble()
        return True

    def _post_step(self, ts):
        """Record force coefficients after each time step."""
        t = ts.getTime()

        # Set state in ADflow from the TS solution
        U = ts.getSolution()
        u_arr = U.getArray(readonly=True).copy()
        self.solver.setStates(u_arr)

        # Compute residual to update all internal Fortran arrays
        # (face pressures, surface forces, etc.) before evalFunctions.
        # Without this, evalFunctions reads stale/uninitialized data.
        adflow = self.solver.adflow
        orig_mode = adflow.inputphysics.equationmode
        adflow.inputphysics.equationmode = _STEADY
        res_tmp = np.zeros(self.n_local)
        adflow.nksolver.getres(res_tmp)
        adflow.inputphysics.equationmode = orig_mode

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

        # Save state for manual backward sweep
        if self.save_trajectory:
            self._state_trajectory.append((t, u_arr.copy()))

        if self.comm.rank == 0:
            step = ts.getStepNumber()
            # Diagnostic: state vector stats
            w_min, w_max = u_arr.min(), u_arr.max()
            w_mean = u_arr.mean()
            actual_dt = ts.getTimeStep()
            print(f"  TS step {step:4d} | t = {t:.6e} | dt = {actual_dt:.6e} | "
                  f"CL = {cl:.6e} | CD = {cd:.6e} | CMz = {cmz:.6e}")
            print(f"           w: min={w_min:.6e}  max={w_max:.6e}  mean={w_mean:.6e}")

    def solve(self):
        """Run the forward time integration."""
        adflow = self.solver.adflow
        n_global = self.comm.allreduce(self.n_local)
        if self.comm.rank == 0:
            print("=" * 70)
            print(f"ADflowTS: Starting PETSc TSTHETA (theta={self.theta}, snes={self.snes_type})")
            print(f"  dt = {self.dt}, t_final = {self.t_final}")
            print(f"  grid_motion = {self.grid_motion}")
            print(f"  n_local = {self.n_local}, n_global = {n_global}")
            if self.snes_type == "dadi":
                print(f"  Forward: ADflow native DADI+MG (+ ANK/NK as configured)")
                print(f"  Adjoint KSP: gmres, rtol={self.ksp_rtol}, max_it={self.ksp_max_it}")
            else:
                print(f"  SNES: rtol={self.snes_rtol}, max_it={self.snes_max_it}")
                print(f"  KSP:  type={self.ksp_type}, rtol={self.ksp_rtol}, max_it={self.ksp_max_it}")
            print("=" * 70)

        # Allocate Fortran convergence/time arrays (required by solveState).
        # Without this, solveState writes to unallocated arrays → SEGV.
        # Mirrors solver.__call__ lines 1237-1244 in pyADflow.py.
        if self.snes_type == "dadi":
            adflow.iteration.itertot = 0
            desired = (adflow.inputiteration.nsgstartup
                       + adflow.inputiteration.ncycles)
            adflow.utils.allocconvarrays(desired)
            n_steps = int(round(self.t_final / self.dt))
            adflow.utils.alloctimearrays(n_steps)
            adflow.solvers.solverunsteadyinit()

        # Record initial state (can't use _post_step here because
        # ts.getSolution() is not valid before ts.solve())
        w0 = self.U_vec.getArray(readonly=True).copy()
        self.solver.setStates(w0)
        adflow = self.solver.adflow
        orig_mode = adflow.inputphysics.equationmode
        adflow.inputphysics.equationmode = _STEADY
        res_tmp = np.zeros(self.n_local)
        adflow.nksolver.getres(res_tmp)
        adflow.inputphysics.equationmode = orig_mode
        funcs = {}
        self.solver.evalFunctions(self.ap, funcs, evalFuncs=["cl", "cd", "cmz"])
        ap_name = self.ap.name
        t0 = self.ts.getTime()
        self.time_history.append(t0)
        self.cl_history.append(funcs.get(f"{ap_name}_cl", 0.0))
        self.cd_history.append(funcs.get(f"{ap_name}_cd", 0.0))
        self.cmz_history.append(funcs.get(f"{ap_name}_cmz", 0.0))
        # Save initial state for trajectory
        if self.save_trajectory:
            self._state_trajectory.append((t0, w0.copy()))
        if self.comm.rank == 0:
            print(f"  TS step    0 | t = {t0:.6e} | "
                  f"CL = {self.cl_history[0]:.6e} | "
                  f"CD = {self.cd_history[0]:.6e} | "
                  f"CMz = {self.cmz_history[0]:.6e}")

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
