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


class _ScalePreconditioner:
    """Trivial PC: y = x / shift.  For J = a*I - dR/dw with large a,
    this approximates J^{-1} ≈ (1/a) * I."""

    def __init__(self, ts_wrapper):
        self.ts_wrapper = ts_wrapper

    def apply(self, pc, x, y):
        a = self.ts_wrapper.jac_ctx.shift
        x_arr = x.getArray(readonly=True)
        y.getArray()[:] = x_arr / a


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
            ksp.setType("gmres")
            ksp.setPCSide(PETSc.PC.Side.RIGHT)
            ksp.setTolerances(rtol=self.ksp_rtol, max_it=self.ksp_max_it)
            ksp.setGMRESRestart(self.ksp_gmres_restart)

            # PC for adjoint
            pc = ksp.getPC()
            if self.pc_type == "none":
                pc.setType("none")
            elif self.pc_type == "scale":
                pc.setType("python")
                pc.setPythonContext(_ScalePreconditioner(self))
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

        # No explicit parameter Jacobian is provided yet, so pass empty mu list.
        self.ts.setCostGradients(lambda_vecs, [])
        self.ts.adjointSolve()

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
        # Use a larger shift (a * pc_shift_factor) so that the assembled
        # matrix is diagonally dominant and ILU factorization is stable.
        # The MatShell operator keeps the exact shift a; GMRES handles the
        # spectral mismatch between PC and operator.
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
