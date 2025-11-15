"""
Compare Torus Time Spectral vs Time Accurate with frequency ratio 1/sqrt(2)

This test demonstrates that torus time spectral can match time-accurate solutions
even with incommensurate frequency ratios by using:
1. Spectral interpolation to evaluate torus solution at arbitrary time points
2. Phase adjustment to align the quasi-periodic solutions
"""

import numpy as np
from mpi4py import MPI
from adflow import ADFLOW
from idwarp import USMesh
from baseclasses import AeroProblem
from scipy.optimize import minimize


def parse_cl_from_log(log_file, n_steps):
    """
    Parse CL values from ADflow monitor output in log file.

    Looks for lines like:
    1    <timestep>  <time>  <iter> ... <CL> <CD> ...

    Parameters
    ----------
    log_file : str
        Path to log file
    n_steps : int
        Expected number of time steps

    Returns
    -------
    cl_values : np.ndarray
        CL values for each time step
    """
    import re

    cl_values = []
    timestep_data = {}

    with open(log_file, 'r') as f:
        for line in f:
            # Pattern: Grid=1, timestep, time, iterations, then monitor variables
            # Example: "      1       1  7.85398E-02     50     50     DADI   2.50E+00  1.00   ----   1.2345E-14   3.3118E-03   4.9708E-04   2.2651E-12"
            # Columns: Grid, TimeStep, Time, Iter, TotalIter, Type, CFL, Step, LinRes, ResRho, CL, CD, TotalRes
            match = re.search(r'^\s*1\s+(\d+)\s+([0-9\.E\+\-]+)\s+\d+\s+\d+\s+\w+\s+[0-9\.E\+\-]+\s+[0-9\.]+\s+\S+\s+([0-9\.E\+\-]+)\s+([0-9\.E\+\-]+)\s+([0-9\.E\+\-]+)', line)
            if match:
                timestep = int(match.group(1))
                cl_val = float(match.group(4))  # CL is 4th captured group
                # Store the LAST CL value for each timestep (final converged value)
                timestep_data[timestep] = cl_val

    # Extract in order
    for ts in range(1, n_steps + 1):
        if ts in timestep_data:
            cl_values.append(timestep_data[ts])
        else:
            print(f"Warning: Could not find CL for timestep {ts}")
            cl_values.append(0.0)

    return np.array(cl_values)


def spectral_interp(U, theta1, theta2):
    """
    Interpolate torus solution U(theta1, theta2) at arbitrary phase points.

    Uses 2D FFT to compute Fourier coefficients, then evaluates at given phases.

    Parameters
    ----------
    U : np.ndarray
        Solution on collocation grid, shape (n1, n2)
    theta1 : np.ndarray
        Phase values for first frequency
    theta2 : np.ndarray
        Phase values for second frequency

    Returns
    -------
    np.ndarray
        Interpolated values at (theta1, theta2) points
    """
    n1, n2 = U.shape
    # Compute 2D Fourier coefficients
    C = np.fft.fft2(U) / (n1 * n2)
    # Frequency indices
    k1 = np.fft.fftfreq(n1, d=1.0/n1)
    k2 = np.fft.fftfreq(n2, d=1.0/n2)
    # Ensure arrays are at least 1D
    theta1, theta2 = np.atleast_1d(theta1), np.atleast_1d(theta2)
    # Exponential basis functions
    e1 = np.exp(1j * np.outer(k1, theta1))  # (n1, m)
    e2 = np.exp(1j * np.outer(k2, theta2))  # (n2, m)
    # Evaluate Fourier series
    vals = np.einsum("ij,im,jm->m", C, e1, e2)
    return np.real(vals)


def phase_alignment(U_torus, theta1_time, theta2_time, values_time):
    """
    Find optimal phase shifts to align torus solution with time-accurate data.

    Minimizes mean-square error between torus solution (with phase shifts)
    and time-accurate solution.

    Parameters
    ----------
    U_torus : np.ndarray
        Torus solution on collocation grid, shape (n1, n2)
    theta1_time : np.ndarray
        Phase values for first frequency from time integration
    theta2_time : np.ndarray
        Phase values for second frequency from time integration
    values_time : np.ndarray
        Time-accurate solution values

    Returns
    -------
    tuple
        (phi1, phi2) optimal phase shifts
    """
    def objective(phi):
        phi1, phi2 = phi
        # Evaluate torus solution with phase shifts
        values_pred = spectral_interp(
            U_torus,
            np.mod(theta1_time + phi1, 2 * np.pi),
            np.mod(theta2_time + phi2, 2 * np.pi)
        )
        # Mean square error
        return np.mean((values_pred - values_time) ** 2)

    # Minimize with Powell method
    res = minimize(objective, x0=np.zeros(2), method="Powell")
    if not res.success:
        print(f"Warning: Phase alignment may not have converged: {res.message}")

    return res.x


def run_torus_time_spectral():
    """Run torus time spectral case with omega1/omega2 = 1/sqrt(2)"""

    print("="*70)
    print("RUNNING TORUS TIME SPECTRAL (3x3, omega1=100, omega2=100*sqrt(2))")
    print("="*70)

    gridFile = '../input_files/naca64A010_euler-L2.cgns'

    n1 = 3
    n2 = 3
    ntot = n1 * n2

    omega1 = 10.0  # Match TA frequency
    omega2 = 10.0 * np.sqrt(2.0)  # Incommensurate ratio 1/sqrt(2)
    A1 = 1.0
    A2 = 1.0

    print(f"\nGrid: {n1} × {n2} = {ntot} instances")
    print(f"Frequencies: ω₁={omega1:.4f}, ω₂={omega2:.4f}")
    print(f"Ratio: ω₁/ω₂ = {omega1/omega2:.6f} = 1/√2")
    print(f"Motion: α(θ₁,θ₂) = {A1}°·sin(θ₁) + {A2}°·sin(θ₂)")

    options = {
        'gridfile': gridFile,
        'outputDirectory': './output_compare',
        'equationType': 'Euler',
        'equationMode': 'time spectral',
        'useTorusTimeSpectral': True,
        'nTimeIntervalsSpectral1': n1,
        'nTimeIntervalsSpectral2': n2,
        'omegaFourier1': omega1,
        'omegaFourier2': omega2,
        'useexternaldynamicmesh': True,
        'usetsinterpolatedgridvelocity': True,
        'mgcycle': 'sg',
        'nCycles': 10000,
        'l2convergence': 1e-10,
        'l2convergencecoarse': 1e-2,
        'monitorvariables': ['resrho', 'cl', 'cd'],
        'usenksolver': True,
        'nkswitchtol': 1e-4,
        'useanksolver': True,
        'ankswitchtol': 1e-2,
        'CFL': 1.5,
        'nSubiter': 3,
        'blocksplitting': True,
        'writeVolumeSolution': False,
        'writeSurfaceSolution': False,
    }

    meshOptions = {"gridFile": gridFile}

    print("\nInitializing ADflow (torus mode)...")
    CFDSolver = ADFLOW(options=options, debug=False)

    print("Setting up mesh...")
    mesh = USMesh(options=meshOptions)
    CFDSolver.setMesh(mesh)

    # Generate alpha grid
    theta1 = np.linspace(0.0, 2.0 * np.pi, n1, endpoint=False)
    theta2 = np.linspace(0.0, 2.0 * np.pi, n2, endpoint=False)

    alpha_grid = np.zeros((n1, n2))

    print("\nPhase and alpha distribution:")
    for i1 in range(n1):
        for i2 in range(n2):
            alpha_deg = A1 * np.sin(theta1[i1]) + A2 * np.sin(theta2[i2])
            alpha_grid[i1, i2] = alpha_deg * np.pi / 180.0
            idx = i2 * n1 + i1
            print(f"  Instance {idx+1} ({i1},{i2}): θ₁={theta1[i1]:.4f}, θ₂={theta2[i2]:.4f}, α={alpha_deg:7.3f}°")

    # Deform meshes
    print("\nDeforming meshes...")
    MDGroup = CFDSolver.allWallsGroup
    cfdPts0 = CFDSolver.getSurfaceCoordinates(MDGroup, includeZipper=False)

    xRot = 0.25
    cfdPts_list = []

    for i2 in range(n2):
        for i1 in range(n1):
            pitch_rad = alpha_grid[i1, i2]
            cc = np.cos(pitch_rad)
            ss = np.sin(pitch_rad)

            cfdPoints_deformed = np.zeros_like(cfdPts0)
            for j in range(len(cfdPts0)):
                x0, y0, z0 = cfdPts0[j]
                cfdPoints_deformed[j, 0] = cc * (x0 - xRot) + ss * y0 + xRot
                cfdPoints_deformed[j, 1] = -ss * (x0 - xRot) + cc * y0
                cfdPoints_deformed[j, 2] = z0

            cfdPts_list.append(cfdPoints_deformed)

    # Set warped meshes
    for sps in range(ntot):
        CFDSolver.mesh.setSurfaceCoordinates(cfdPts_list[sps])
        CFDSolver.mesh.warpMesh()
        m = CFDSolver.mesh.getSolverGrid()
        CFDSolver.adflow.warping.setgridforoneinstance(m, sps=sps + 1)

    CFDSolver._updateGeomInfo = True
    CFDSolver.updateGeometryInfo()

    # Solve
    print("\nSolving torus system...")
    ap = AeroProblem(
        name='torus_ts',
        mach=0.5,
        altitude=10000,
        areaRef=1.0,
        chordRef=1.0,
        evalFuncs=['cl', 'cd']
    )

    CFDSolver(ap)

    # Extract CL values
    forces_torus = CFDSolver.getInstanceForces()
    cl_torus = forces_torus['cl']

    print("\nTorus CL values extracted:")
    for i in range(ntot):
        i1 = i % n1
        i2 = i // n1
        alpha_deg = alpha_grid[i1, i2] * 180 / np.pi
        print(f"  Instance {i+1} ({i1},{i2}): θ₁={theta1[i1]:.4f}, θ₂={theta2[i2]:.4f}, α={alpha_deg:7.3f}°, CL={cl_torus[i]:.10f}")

    # Reshape CL to grid for spectral interpolation
    cl_grid = cl_torus.reshape((n2, n1)).T  # Shape (n1, n2)

    return cl_grid, alpha_grid, theta1, theta2, omega1, omega2, n1, n2


def run_time_accurate():
    """Run time-accurate case with dual-frequency forcing"""

    print("\n" + "="*70)
    print("RUNNING TIME ACCURATE (dual-frequency motion)")
    print("="*70)

    gridFile = '../input_files/naca64A010_euler-L2.cgns'

    # Time step setup (inspired by reg test)
    # Run 5x longer with 5x finer time step for better resolution
    n_steps_per_period1 = 40  # Steps per period (8 * 5 for 5x finer resolution)
    omega1 = 10.0  # Reduced frequency for larger time step
    omega2 = 10.0 * np.sqrt(2.0)  # Maintain 1/sqrt(2) ratio
    T1 = 2 * np.pi / omega1
    dt = T1 / n_steps_per_period1  # 5x smaller time step

    # Run for multiple periods to reach quasi-periodic state
    n_periods = 10  # 5x longer duration (2 * 5 = 10 periods)
    t_end = n_periods * T1
    n_steps = int(t_end / dt)

    A1 = 1.0
    A2 = 1.0
    xRot = 0.25

    print(f"\nFrequencies: ω₁={omega1:.4f}, ω₂={omega2:.4f}")
    print(f"Ratio: ω₁/ω₂ = {omega1/omega2:.6f} = 1/√2")
    print(f"Time step: dt={dt:.6f}, Total time: {t_end:.2f}")
    print(f"Number of periods: {n_periods}")
    print(f"Total time steps: {n_steps}")
    print(f"Motion: α(t) = {A1}°·sin(ω₁t) + {A2}°·sin(ω₂t)")

    # Storage for force history (will be populated by callback)
    force_history = {'time': [], 'cl': [], 'cd': []}

    options = {
        'gridfile': gridFile,
        'outputDirectory': './output_compare',
        'equationType': 'Euler',
        'equationMode': 'unsteady',
        'timeIntegrationScheme': 'BDF',
        'useTorusTimeSpectral': False,  # Disable torus mode for time accurate
        'deltaT': dt,
        'nTimeStepsFine': n_steps,  # CRITICAL: Number of time steps to run
        'useALE': False,  # Use external mesh motion like reg test
        'useGridMotion': True,  # Enable grid motion (callback will provide deformed mesh)
        'mgcycle': '3w',  # Multigrid like reg test
        'mgstartlevel': 1,
        'nCycles': 100,  # Max iterations per time step (prevent wasteful iterations)
        'l2convergence': 1e-6,  # Relative convergence tolerance
        'l2convergencecoarse': 1e-4,  # Coarse grid tolerance
        'monitorvariables': ['resrho', 'cl', 'cd'],
        'usenksolver': False,  # Disable NK for TA like reg test
        'useanksolver': False,  # Disable ANK for TA like reg test
        'CFL': 2.5,  # Higher CFL like reg test
        'CFLCoarse': 1.2,
        'nSubiter': 5,  # More subiters like reg test
        'nSubiterTurb': 10,
        'blocksplitting': True,
        'writeVolumeSolution': False,
        'writeSurfaceSolution': False,
        'nTimeStepsCoarse': 0,
    }

    meshOptions = {"gridFile": gridFile}

    print("\nInitializing ADflow (time accurate)...")
    CFDSolver = ADFLOW(options=options, debug=False)

    print("Setting up mesh...")
    mesh = USMesh(options=meshOptions)
    CFDSolver.setMesh(mesh)

    # Define mesh deformation callback that also extracts forces
    def surfaceMeshCallback(refGrid, t, ts):
        """
        Deform mesh for dual-frequency pitching motion and store forces.

        Parameters
        ----------
        refGrid : np.ndarray
            Reference (undeformed) surface mesh coordinates
        t : float
            Current time
        ts : int
            Current time step number (1-indexed)

        Returns
        -------
        np.ndarray
            Deformed surface mesh coordinates
        """
        # Compute pitch angle from dual-frequency motion
        alpha_deg = A1 * np.sin(omega1 * t) + A2 * np.sin(omega2 * t)
        pitch_rad = alpha_deg * np.pi / 180.0

        # Apply rotation about xRot
        newGrid = np.copy(refGrid)
        x = refGrid[:, 0]
        y = refGrid[:, 1]

        c = np.cos(pitch_rad)
        s = np.sin(pitch_rad)

        newGrid[:, 0] = c * (x - xRot) + s * y + xRot
        newGrid[:, 1] = -s * (x - xRot) + c * y

        # Store time (we'll get forces after the solve)
        if MPI.COMM_WORLD.rank == 0:
            force_history['time'].append(t)
            if ts % 4 == 1 or ts == n_steps:
                print(f"  Step {ts}/{n_steps}, t={t:.4f}, α={alpha_deg:.4f}°")

        return newGrid

    # Create aeroproblem
    ap = AeroProblem(
        name='time_accurate',
        mach=0.5,
        altitude=10000,
        areaRef=1.0,
        chordRef=1.0,
        evalFuncs=['cl', 'cd']
    )

    print(f"\nStarting time integration with callback...")
    print(f"ADflow will run {n_steps} time steps internally\n")

    # Run time-accurate simulation
    # ADflow will call surfaceMeshCallback for each time step automatically
    CFDSolver(ap, surfaceMeshCallback=surfaceMeshCallback)

    # After solve, extract final forces and save times
    # Note: ADflow's monitor output prints CL/CD for each time step
    # We'll save the time array here, and parse CL from log file later
    if MPI.COMM_WORLD.rank == 0:
        print("\nTime accurate solve completed")
        print(f"Time points stored: {len(force_history['time'])}")

        time_history = np.array(force_history['time'])

        # We'll parse CL values from the log file in post-processing
        # For now, save times and create placeholder for CL
        cl_history = np.zeros_like(time_history)

        print(f"  Time range: [{time_history[0]:.4f}, {time_history[-1]:.4f}]")
        print("  Note: CL values will be parsed from log file output")
    else:
        # Non-root ranks: create empty arrays (not used in subprocess mode)
        time_history = np.array([])
        cl_history = np.array([])

    return time_history, cl_history, omega1, omega2


def main():
    import sys
    import subprocess
    import os

    # Only print header on rank 0 when not running in subprocess mode
    if len(sys.argv) <= 1 and MPI.COMM_WORLD.rank == 0:
        print("="*70)
        print("COMPARISON: Torus Time Spectral vs Time Accurate (ω ratio = 1/√2)")
        print("="*70)

    # Check if we should run only one mode (for separate processes)
    if len(sys.argv) > 1:
        if sys.argv[1] == 'torus':
            cl_grid, alpha_grid, theta1, theta2, omega1, omega2, n1, n2 = run_torus_time_spectral()
            # Save results (only rank 0)
            if MPI.COMM_WORLD.rank == 0:
                np.savez('torus_results.npz',
                         cl_grid=cl_grid, alpha_grid=alpha_grid,
                         theta1=theta1, theta2=theta2,
                         omega1=omega1, omega2=omega2, n1=n1, n2=n2)
                print("\nTorus results saved to torus_results.npz")
            sys.exit(0)
        elif sys.argv[1] == 'ta':
            time_ta, cl_ta, omega1, omega2 = run_time_accurate()
            # Save results (only rank 0)
            if MPI.COMM_WORLD.rank == 0:
                np.savez('ta_results.npz',
                         time_ta=time_ta, cl_ta=cl_ta,
                         omega1=omega1, omega2=omega2)
                print("\nTime-accurate results saved to ta_results.npz")
            sys.exit(0)

    # Run both in separate processes to avoid ADflow state conflicts
    # Only run on rank 0 (this is the orchestrator, not an MPI job)
    if MPI.COMM_WORLD.rank != 0:
        return None

    print("\n" + "="*70)
    print("Running solvers in separate processes to avoid state conflicts")
    print("="*70)

    # Determine number of cores to use
    n_cores_torus = 4  # Torus TS is relatively fast
    n_cores_ta = 16    # Time accurate needs more parallelization

    print(f"\n[1/2] Running torus time spectral solver with {n_cores_torus} cores...")
    with open('torus_run.log', 'w') as f:
        result = subprocess.run(['mpirun', '-np', str(n_cores_torus), sys.executable, __file__, 'torus'],
                              stdout=f, stderr=subprocess.STDOUT, cwd=os.getcwd())
    if result.returncode != 0:
        print("Torus run failed! Check torus_run.log for details")
        raise RuntimeError("Torus time spectral run failed")
    print("Torus time spectral completed successfully")

    print(f"\n[2/2] Running time-accurate solver with {n_cores_ta} cores...")
    with open('ta_run.log', 'w') as f:
        result = subprocess.run(['mpirun', '-np', str(n_cores_ta), sys.executable, __file__, 'ta'],
                              stdout=f, stderr=subprocess.STDOUT, cwd=os.getcwd())
    if result.returncode != 0:
        print("Time-accurate run failed! Check ta_run.log for details")
        raise RuntimeError("Time accurate run failed")
    print("Time-accurate completed successfully")

    # Parse CL values from TA log file
    print("\nParsing CL values from time-accurate log file...")
    ta_data_temp = np.load('ta_results.npz')
    n_steps_ta = len(ta_data_temp['time_ta'])
    cl_ta_parsed = parse_cl_from_log('ta_run.log', n_steps_ta)

    # Update ta_results.npz with parsed CL values
    np.savez('ta_results.npz',
             time_ta=ta_data_temp['time_ta'],
             cl_ta=cl_ta_parsed,
             omega1=ta_data_temp['omega1'],
             omega2=ta_data_temp['omega2'])
    print(f"Parsed {len(cl_ta_parsed)} CL values from log")

    # Load results
    print("\n" + "="*70)
    print("Loading and comparing results...")
    print("="*70)

    torus_data = np.load('torus_results.npz')
    ta_data = np.load('ta_results.npz')

    cl_grid = torus_data['cl_grid']
    alpha_grid = torus_data['alpha_grid']
    theta1 = torus_data['theta1']
    theta2 = torus_data['theta2']
    omega1_ts = torus_data['omega1']
    omega2_ts = torus_data['omega2']
    n1 = int(torus_data['n1'])
    n2 = int(torus_data['n2'])

    time_ta = ta_data['time_ta']
    cl_ta = ta_data['cl_ta']
    omega1_ta = ta_data['omega1']
    omega2_ta = ta_data['omega2']

    # Verify frequencies match
    assert np.isclose(omega1_ts, omega1_ta), "Frequency mismatch for omega1"
    assert np.isclose(omega2_ts, omega2_ta), "Frequency mismatch for omega2"
    omega1 = omega1_ts
    omega2 = omega2_ts

    # Compute phase values from time-accurate solution
    # Skip initial transient (first 20% of time)
    transient_idx = int(0.2 * len(time_ta))
    time_steady = time_ta[transient_idx:]
    cl_steady = cl_ta[transient_idx:]

    theta1_ta = (omega1 * time_steady) % (2 * np.pi)
    theta2_ta = (omega2 * time_steady) % (2 * np.pi)

    print("\n" + "="*70)
    print("PHASE ALIGNMENT")
    print("="*70)

    # Find optimal phase shifts
    print("\nComputing optimal phase shifts...")
    phi1, phi2 = phase_alignment(cl_grid, theta1_ta, theta2_ta, cl_steady)
    print(f"Optimal phase shifts: φ₁ = {phi1:.6f} rad, φ₂ = {phi2:.6f} rad")

    # Interpolate torus solution with phase adjustment
    cl_torus_interp = spectral_interp(
        cl_grid,
        np.mod(theta1_ta + phi1, 2 * np.pi),
        np.mod(theta2_ta + phi2, 2 * np.pi)
    )

    # Compute errors
    error = cl_torus_interp - cl_steady
    rel_l2 = np.linalg.norm(error) / np.linalg.norm(cl_steady)
    linf = np.max(np.abs(error))
    mean_error = np.mean(error)
    std_error = np.std(error)

    print("\n" + "="*70)
    print("COMPARISON RESULTS")
    print("="*70)

    print(f"\nTime-accurate solution: {len(cl_steady)} points (after transient removal)")
    print(f"Torus grid: {n1}×{n2} = {n1*n2} instances")

    print(f"\nError metrics:")
    print(f"  Relative L₂ error: {rel_l2:.4e}")
    print(f"  L∞ error:          {linf:.4e}")
    print(f"  Mean error:        {mean_error:.4e}")
    print(f"  Std error:         {std_error:.4e}")

    print(f"\nSample comparison (first 10 points after transient):")
    for i in range(min(10, len(time_steady))):
        t = time_steady[i]
        print(f"  t={t:7.4f}: CL_TA={cl_steady[i]:.8f}, CL_Torus={cl_torus_interp[i]:.8f}, Δ={error[i]:+.2e}")

    # Check if errors are within acceptable tolerance
    tolerance = 1e-3  # Relaxed tolerance for TA vs TS comparison
    if rel_l2 < tolerance:
        print(f"\n✓ SUCCESS: Torus TS matches TA within tolerance ({rel_l2:.2e} < {tolerance})")
    else:
        print(f"\n✗ WARNING: Errors may be larger than expected ({rel_l2:.2e} >= {tolerance})")

    # Plot CL history comparison
    print("\n" + "="*70)
    print("PLOTTING CL HISTORY")
    print("="*70)

    try:
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

        # Plot 1: CL vs time
        ax1.plot(time_steady, cl_steady, 'b-', linewidth=2, label='Time Accurate', alpha=0.7)
        ax1.plot(time_steady, cl_torus_interp, 'r--', linewidth=1.5, label='Torus TS (interpolated)', alpha=0.9)
        ax1.set_xlabel('Time [s]')
        ax1.set_ylabel('CL')
        ax1.set_title(f'CL History: Torus TS vs Time Accurate (ω₁/ω₂ = 1/√2)')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Plot 2: Error history
        ax2.plot(time_steady, error, 'k-', linewidth=1)
        ax2.axhline(y=0, color='r', linestyle='--', alpha=0.5)
        ax2.set_xlabel('Time [s]')
        ax2.set_ylabel('Error (Torus - TA)')
        ax2.set_title(f'Error: Rel L₂ = {rel_l2:.2e}, L∞ = {linf:.2e}')
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig('torus_vs_ta_comparison.pdf', dpi=300, bbox_inches='tight')
        plt.savefig('torus_vs_ta_comparison.png', dpi=150, bbox_inches='tight')
        print("\nPlots saved:")
        print("  - torus_vs_ta_comparison.pdf")
        print("  - torus_vs_ta_comparison.png")

    except ImportError:
        print("\nMatplotlib not available - skipping plots")

    print("\n" + "="*70)

    return {
        'time_ta': time_steady,
        'cl_ta': cl_steady,
        'cl_torus': cl_torus_interp,
        'error': error,
        'rel_l2': rel_l2,
        'linf': linf,
        'phi1': phi1,
        'phi2': phi2,
    }


if __name__ == "__main__":
    results = main()