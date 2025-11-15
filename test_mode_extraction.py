#!/usr/bin/env python
"""
Quick test of resolvent mode extraction functionality.

This script verifies that:
1. Multiple modes can be extracted
2. Modes are properly normalized
3. API works correctly
"""

import sys
import os
import numpy as np

# Add parent directory for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from adflow import ADFLOW
    from baseclasses import AeroProblem
    ADFLOW_AVAILABLE = True
except ImportError:
    print("ERROR: Could not import ADFLOW")
    ADFLOW_AVAILABLE = False
    sys.exit(1)

def test_mode_extraction():
    """
    Quick test of mode extraction on NACA 64A010 airfoil.
    """

    print("="*80)
    print("TESTING RESOLVENT MODE EXTRACTION")
    print("="*80)
    print()

    # Check if mesh file exists
    meshFile = './input_files/naca64A010_euler-L2.cgns'
    if not os.path.exists(meshFile):
        meshFile = './tests/input_files/naca64A010_euler-L2.cgns'

    if not os.path.exists(meshFile):
        print(f"ERROR: Mesh file not found")
        print("Please run from repo root or tests/ directory")
        return False

    print(f"✓ Found mesh: {meshFile}\n")

    # =========================================================================
    # Setup CFD
    # =========================================================================

    aeroOptions = {
        'gridFile': meshFile,
        'outputDirectory': './output_test_modes',

        # Physics
        'equationType': 'Euler',

        # Fast convergence for testing
        'CFL': 2.0,
        'L2Convergence': 1e-6,
        'nCycles': 200,
        'useNKSolver': True,

        # Reduce output
        'printIterations': False,
        'printTiming': False,
        'printWarnings': False,
        'writeVolumeSolution': False,
        'writeSurfaceSolution': False,

        # CRITICAL for resolvent analysis
        'useMatrixFreedrdw': False,
    }

    print("Creating CFD solver...")
    CFDsolver = ADFLOW(options=aeroOptions, debug=False)

    ap = AeroProblem(
        name='test',
        mach=0.5,
        alpha=0.0,
        altitude=10000,
        areaRef=1.0,
        chordRef=1.0,
    )

    print("Solving flow...")
    CFDsolver(ap)
    print("✓ Flow converged\n")

    # =========================================================================
    # Test Matrix-Free Mode Extraction
    # =========================================================================

    print("="*80)
    print("TEST 1: Matrix-Free Mode Extraction")
    print("="*80)
    print()

    from adflow import ResolventAnalysisMatrixFree

    omega = 1.0
    nModes = 3

    print(f"Computing {nModes} modes at ω = {omega}...")
    resolvent = ResolventAnalysisMatrixFree(CFDsolver, ap, omega=omega)

    # Enable ILU for fast convergence
    resolvent.enablePreconditioner(precond_type='ilu', drop_tol=1e-3, fill_factor=10)

    # Solve
    sigma_max = resolvent.solve(nModes=nModes, method='scipy')
    print(f"✓ Modes computed\n")

    # Test extracting modes
    print("Extracting modes...")
    print("-"*80)

    all_pass = True

    for i in range(nModes):
        v_i = resolvent.getForcingMode(i)
        u_i = resolvent.getResponseMode(i)
        sigma_i = resolvent.getSingularValue(i)

        print(f"\nMode {i+1}:")
        print(f"  σ_{i+1} = {sigma_i:.6f}")
        print(f"  ||v_{i+1}|| = {np.linalg.norm(v_i):.6f}")
        print(f"  ||u_{i+1}|| = {np.linalg.norm(u_i):.6f}")

        # Check that modes are normalized
        v_norm = np.linalg.norm(v_i)
        u_norm = np.linalg.norm(u_i)

        if not (0.99 < v_norm < 1.01):
            print(f"  ✗ FAIL: Forcing mode not normalized (||v|| = {v_norm:.6f})")
            all_pass = False
        else:
            print(f"  ✓ Forcing mode normalized")

        if not (0.99 < u_norm < 1.01):
            print(f"  ✗ FAIL: Response mode not normalized (||u|| = {u_norm:.6f})")
            all_pass = False
        else:
            print(f"  ✓ Response mode normalized")

        # Check complex type
        if not np.iscomplexobj(v_i):
            print(f"  ✗ FAIL: Forcing mode should be complex")
            all_pass = False
        else:
            print(f"  ✓ Forcing mode is complex")

    # Test orthogonality
    print("\nTesting orthogonality...")
    print("-"*80)

    for i in range(nModes):
        for j in range(i+1, nModes):
            v_i = resolvent.getForcingMode(i)
            v_j = resolvent.getForcingMode(j)

            inner = np.abs(np.vdot(v_i, v_j))

            if inner > 0.01:
                print(f"  ✗ FAIL: Modes {i+1} and {j+1} not orthogonal: <v{i+1},v{j+1}> = {inner:.6f}")
                all_pass = False
            else:
                print(f"  ✓ Modes {i+1} and {j+1} orthogonal: <v{i+1},v{j+1}> = {inner:.6e}")

    # =========================================================================
    # Test Explicit Mode Extraction
    # =========================================================================

    print("\n" + "="*80)
    print("TEST 2: Explicit Mode Extraction")
    print("="*80)
    print()

    from adflow import ResolventAnalysis

    print(f"Computing {nModes} modes at ω = {omega} (explicit method)...")
    resolvent_explicit = ResolventAnalysis(CFDsolver, ap, omega=omega)
    resolvent_explicit.nModes = nModes
    resolvent_explicit.solveExplicit(formJacobian=True, useRealForm=False, useLU=True)
    print(f"✓ Modes computed\n")

    # Test extracting modes
    print("Extracting modes...")
    print("-"*80)

    for i in range(nModes):
        v_i = resolvent_explicit.getForcingMode(i)
        u_i = resolvent_explicit.getResponseMode(i)
        sigma_i = resolvent_explicit.getSingularValue(i)

        print(f"\nMode {i+1}:")
        print(f"  σ_{i+1} = {sigma_i:.6f}")
        print(f"  ||v_{i+1}|| = {np.linalg.norm(v_i):.6f}")
        print(f"  ||u_{i+1}|| = {np.linalg.norm(u_i):.6f}")

        # Check normalized
        v_norm = np.linalg.norm(v_i)
        if 0.99 < v_norm < 1.01:
            print(f"  ✓ Forcing mode normalized")
        else:
            print(f"  ✗ FAIL: Not normalized (||v|| = {v_norm:.6f})")
            all_pass = False

    # =========================================================================
    # Summary
    # =========================================================================

    print("\n" + "="*80)
    if all_pass:
        print("✓ ALL TESTS PASSED")
    else:
        print("✗ SOME TESTS FAILED")
    print("="*80)
    print()

    return all_pass


if __name__ == "__main__":
    success = test_mode_extraction()
    sys.exit(0 if success else 1)
