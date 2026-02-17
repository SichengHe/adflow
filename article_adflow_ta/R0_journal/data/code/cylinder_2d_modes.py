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


def test_transonic_modes():
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

    # Check mesh file
    baseDir = os.path.dirname(os.path.abspath(__file__))
    meshFile = os.path.join(baseDir, "./grids/cylinder.cgns")

    if not os.path.exists(meshFile):
        print(f"ERROR: Mesh file not found: {meshFile}")
        print("\nPlease run: cd input_files && ./get-input-files.sh")
        return None, None, None

    if rank == 0:
        print(f"Mesh: {meshFile}")
        print()

    # =========================================================================
    # Setup ADflow - transonic settings
    # =========================================================================

    outputDir = './output_cylinder'
    os.makedirs(outputDir, exist_ok=True)

    aeroOptions = {
        'gridFile': meshFile,
        'outputDirectory': outputDir,

        # Physics
        'equationType': 'Euler',

        # Solver parameters - tighter convergence for shock resolution
        'CFL': 2.0,
        'L2Convergence': 1e-12,  # Tight convergence for shock
        'nCycles': 100000,
        'useNKSolver': False,
        'useANKSolver': False,


        # Output
        'printIterations': True,
        'printTiming': False,
        'printWarnings': True,
        #"volumeVariables": ["lstevecvelx","lstevecvely","lstevecvelx","lstevecrho"],
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
        name='cylinder',
        mach=0.1,       
        alpha=0.0,       
        reynolds=100,
        reynoldsLength=1.0,
        #altitude=10000,
        areaRef=1.0,
        chordRef=1.0,
        T=300.0,
    )

    CFDsolver(ap)

    # Print converged forces
    funcs = {}
    CFDsolver.evalFunctions(ap, funcs)
    if rank == 0:
        print()
        print("Converged solution:")
        print(f"  CL = {funcs.get('cylinder_transonic_cl', 'N/A')}")
        print(f"  CD = {funcs.get('cylinder_transonic_cd', 'N/A')}")
        print()

    # =========================================================================
    # Compute modes
    # =========================================================================

    omega = 0.01  # See the sigma vs omega sweep plot. This is the highest sigma's omega.
    nModes = 1  # Three mode pairs for visualization. RK - Changed to 1 temporarily.

    if rank == 0:
        print(f"Step 2: Computing {nModes} resolvent modes at ω = {omega}...")
        print("  - Modes will show interaction with shock wave")
    resolvent = ResolventAnalysisMatrixFree(CFDsolver, ap, omega=omega)

    # Enable ILU for fast convergence
    if rank == 0:
        print("  - Enabling ILU preconditioner...")
    resolvent.enablePreconditioner(precond_type='ilu', drop_tol=1e-3, fill_factor=10)

    # Solve
    if rank == 0:
        print("  - Solving for modes...")
    sigma_max = resolvent.solve(nModes=nModes, method='fortran')
    if rank == 0:
        print(f"✓ Modes computed (σ₁ = {sigma_max:.6f})")
        print()

    # =========================================================================
    # Display mode information
    # =========================================================================

    if rank == 0:
        print("Mode Information:")
        print("-"*80)
    for i in range(nModes):
        sigma = resolvent.getSingularValue(i)
        v_i = resolvent.getForcingMode(i)
        u_i = resolvent.getResponseMode(i)

        if rank == 0:
            print(f"Mode {i+1}:")
            print(f"  σ_{i+1} = {sigma:.6f}")
            print(f"  Forcing mode shape: {v_i.shape}, dtype={v_i.dtype}")
            print(f"  Response mode shape: {u_i.shape}, dtype={u_i.dtype}")
            print(f"  ||v_{i+1}|| = {np.linalg.norm(v_i):.6f}")
            print(f"  ||u_{i+1}|| = {np.linalg.norm(u_i):.6f}")
            print()

    if rank == 0:
        print("="*80)
        print("Mode computation complete!")
        print("="*80)
        print()

    # =========================================================================
    # Save modes
    # =========================================================================

    if rank == 0:
        print("="*80)
        print("SAVING TRANSONIC MODE DATA")
        print("="*80)
        print()

    

    # Write to PETSc files
    if rank == 0:
        print("Writing modes to PETSc format...")
    petsc_dir = os.path.join(outputDir, 'modes_petsc')
    resolvent.writeModes(petsc_dir, load_into_fortran=True, load_kind='forcing_real')

    # Load a PETSc mode into Fortran storage and write a volume CGNS file
    mode_file = os.path.join(petsc_dir, "mode1_response_real.petsc")
    if rank == 0:
        print(f"Loading PETSc mode into Fortran: {mode_file}")
    try:
        CFDsolver.adflow.resolventvistools.readandseteigenvector(mode_file)
        out_vol = os.path.join(outputDir, "out_vol")
        CFDsolver.writeVolumeSolutionFile(out_vol)
        if rank == 0:
            print(f"✓ Wrote volume solution: {out_vol}.cgns")
        CFDsolver.setOption(
            "surfaceVariables",
            ["Deltarho", "Deltavelx", "Deltavely", "Deltavelz", "Deltarhoe"],
        )
        out_surf = os.path.join(outputDir, "out_surf")
        try:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            comm.Barrier()
            ok = True
            try:
                CFDsolver.writeSurfaceSolutionFile(out_surf)
            except Exception as e:
                ok = False
                print(f"[rank {comm.Get_rank()}] writeSurfaceSolutionFile failed: {e}", flush=True)
            ok_all = comm.allreduce(ok, op=MPI.LAND)
            if not ok_all:
                comm.Abort(1)
            comm.Barrier()
            if rank == 0:
                print(f"✓ Wrote surface solution: {out_surf}.cgns")
        except Exception:
            CFDsolver.writeSurfaceSolutionFile(out_surf)
            if rank == 0:
                print(f"✓ Wrote surface solution: {out_surf}.cgns")
    except Exception as e:
        if rank == 0:
            print(f"✗ PETSc->Fortran->CGNS write failed: {e}")

    return resolvent, CFDsolver, ap


if __name__ == "__main__":
    try:
        from mpi4py import MPI
        _rank = MPI.COMM_WORLD.Get_rank()
    except Exception:
        _rank = 0
    # Run test
    resolvent, CFDsolver, ap = test_transonic_modes()

    if resolvent is not None and _rank == 0:
        print("SUCCESS: Transonic test completed!")
