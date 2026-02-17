#!/usr/bin/env python
"""
Transonic test case with shock waves for resolvent mode visualization.

This case computes resolvent modes for transonic flow over NACA 64A010 airfoil:
- Mach number: 0.8 (transonic, shock waves present)
- Angle of attack: 3.0 degrees
- Altitude: 10000 ft

Expected features:
- Normal shock wave on upper surface
- Strong interaction between shock and mode structure
- Interesting response mode patterns due to shock sensitivity
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Add repo root to path (so local adflow is used)
repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, repo_root)

try:
    from adflow import ADFLOW, ResolventAnalysisMatrixFree
    from baseclasses import AeroProblem
    ADFLOW_AVAILABLE = True
except ImportError:
    print("ERROR: Could not import ADFLOW or ResolventAnalysisMatrixFree")
    ADFLOW_AVAILABLE = False
    sys.exit(1)


def resolvent_sweep_transonic():
    """
    Transonic test case with shock waves.

    Returns
    -------
    resolvent : ResolventAnalysisMatrixFree
        Resolvent object with computed modes
    CFDsolver : ADFLOW
        ADflow solver instance
    ap : AeroProblem
        Aerodynamic problem
    """

    # MPI rank (only rank 0 prints)
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
    except Exception:
        comm = None
        rank = 0
    
    if rank==0:
        print("="*80)
        print("TRANSONIC RESOLVENT MODE TEST (M=0.8, AoA=3°)")
        print("="*80)
        print()

    # Check mesh file
    baseDir = os.path.dirname(os.path.abspath(__file__))
    meshFile = os.path.join(baseDir, "./grids/crm.cgns") 

    if rank == 0:
        print(f"Mesh: {meshFile}")
        print()

    # =========================================================================
    # Setup ADflow - transonic settings
    # =========================================================================

    outputDir = './output_crm'
    os.makedirs(outputDir, exist_ok=True)

    aeroOptions = {
        'gridFile': meshFile,
        'outputDirectory': outputDir,

        # Physics
        'equationType': 'Euler',

        # Solver parameters - tighter convergence for shock resolution
        'CFL': 2.0,
        'L2Convergence': 1e-8,  # Tight convergence for shock
        'nCycles': 5000,
        'useNKSolver': True,

        # Output
        'printIterations': True,
        'printTiming': False,
        'printWarnings': True,
        "volumeVariables": ["lstevecvelx","lstevecvely","lstevecvelx","lstevecrho"],
        'writeVolumeSolution': True,
        'writeSurfaceSolution': True,

        # CRITICAL for resolvent analysis
        'useMatrixFreedrdw': False,
    }

    if rank == 0:
        print("Step 1: Solving transonic CFD (M=0.8, AoA=3°)...")
        print("  - Expecting shock waves on upper surface")
    CFDsolver = ADFLOW(options=aeroOptions, debug=False)

    ap = AeroProblem(
        name='crm_transonic',
        mach=0.8,        # Transonic - shock waves
        alpha=3.0,       # 3 degrees AoA
        altitude=10000,
        areaRef=1.0,
        chordRef=1.0,
    )

    CFDsolver(ap)

    # Print converged forces
    funcs = {}
    CFDsolver.evalFunctions(ap, funcs)
    if rank == 0:
        print()
        print("Converged solution:")
        print(f"  CL = {funcs.get('crm_transonic_cl', 'N/A')}")
        print(f"  CD = {funcs.get('crm_transonic_cd', 'N/A')}")
        print()

    # =========================================================================
    # Compute frequency sweep
    # =========================================================================

    omega_min = 0.005
    omega_max = 0.06
    nPoints = 20
    nModes = 1

    if rank == 0:
        print(f"Step 2: Frequency sweep ({nPoints} points) from ω={omega_min} to ω={omega_max}...")

    resolvent = ResolventAnalysisMatrixFree(CFDsolver, ap, omega=omega_min)

    resolvent.setSVDTol(1e-4)     # SVD convergence tolerance
    resolvent.maxSvdIter = 20     # Max SVD iterations (Arnoldi/SVD iterations)
    resolvent.svdNcv = 10         # Krylov subspace size for SVD (ncv)
    resolvent.setLinearSolveTol(1e-6)     # Absolute GMRES tolerance
    resolvent.setLinearSolveRelTol(1e-4)  # Relative GMRES tolerance
    resolvent.maxLinearIter = 300         # Max GMRES iterations
    resolvent.gmresRestart = 200          # GMRES restart size
    resolvent.fortranILUFill = None#2        # ASM+ILU fill level (None uses ADflow default)
    resolvent.fortranASMOverlap = None#1     # ASM overlap (None uses ADflow default)
    resolvent.resolventPCType = 2       # 1=block‑diag, 2=block‑triangular

    if rank == 0:
        print("  - Solving resolvent with Fortran backend")

    omega_vec = np.linspace(omega_min, omega_max, nPoints)
    sigma_mat = np.zeros((nPoints, nModes)) if rank == 0 else None

    for i, omega in enumerate(omega_vec):
        if rank == 0:
            print(f"[{i+1}/{nPoints}] ω = {omega:.6f}")
        resolvent.setOmega(omega)

        # Solve at this frequency (Fortran backend)
        sigma1 = resolvent.solve(nModes=nModes, method='fortran')

        if rank == 0:
            if resolvent.singularValues is not None and resolvent.singularValues.size >= nModes:
                sigma_mat[i, :] = np.asarray(resolvent.singularValues[:nModes], dtype=float)
            else:
                sigma_mat[i, 0] = float(sigma1)
            print(f"  → σ₁ = {sigma_mat[i, 0]:.6f}")

    if rank == 0:
        idx = int(np.argmax(sigma_mat[:, 0]))
        print(f"✓ Sweep complete. Max σ₁ at ω={omega_vec[idx]:.6f} (σ₁={sigma_mat[idx, 0]:.6f})")
        print()


    # =========================================================================
    # Save sweep data
    # =========================================================================

    if rank == 0:
        print("Saving sweep data to NumPy format...")
    npz_file = os.path.join(outputDir, f'sweep_omega_{omega_min}_to_{omega_max}.npz')
    if rank == 0:
        np.savez_compressed(
            npz_file,
            omega=omega_vec,
            sigma=sigma_mat,
            nModes=nModes,
        )
        print(f"✓ Saved sweep data: {npz_file}")
        print()
        print("Plotting sweep (per mode)...")
        plot_files = []
        for j in range(nModes):
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.plot(omega_vec, sigma_mat[:, j], marker='o', linewidth=1.5)
            ax.set_xlabel(r"$\omega$")
            ax.set_ylabel(rf"$\sigma_{{{j+1}}}(\omega)$")
            ax.set_title(f"Resolvent Sweep: Mode {j+1}")
            ax.grid(True, alpha=0.3)
            fig.tight_layout()
            plot_file = os.path.join(
                outputDir, f'sweep_mode{j+1}_omega_{omega_min}_to_{omega_max}.png'
            )
            fig.savefig(plot_file, dpi=200)
            plt.close(fig)
            plot_files.append(plot_file)
            print(f"✓ Saved plot: {plot_file}")
        print()

    if rank == 0:
        print("="*80)
        print("TRANSONIC RESOLVENT SWEEP COMPLETE")
        print("="*80)
        print()
        print("Output files:")
        print(f"  - Sweep: {npz_file}")
        for plot_file in plot_files:
            print(f"  - Plot:  {plot_file}")
        print()

    return resolvent, CFDsolver, ap


if __name__ == "__main__":
    try:
        from mpi4py import MPI
        _rank = MPI.COMM_WORLD.Get_rank()
    except Exception:
        _rank = 0
    # Run test
    resolvent, CFDsolver, ap = resolvent_sweep_transonic()

    if resolvent is not None and _rank == 0:
        print("Sweep complete.")
