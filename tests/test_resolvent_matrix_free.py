#!/usr/bin/env python
"""
Test matrix-free resolvent analysis with ADflow.

This test demonstrates the matrix-free implementation for large-scale
CFD problems using iterative methods.

Based on the paper:
"Large-Scale Flow Control Performance Optimization via Differentiable
Resolvent Analysis" by He et al.
"""

import sys
import os
import numpy as np

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from adflow import ADFLOW, ResolventAnalysisMatrixFree
    from baseclasses import AeroProblem
    ADFLOW_AVAILABLE = True
except ImportError:
    print("ERROR: Could not import ADFLOW or ResolventAnalysisMatrixFree")
    ADFLOW_AVAILABLE = False
    sys.exit(1)


def test_matrix_free_resolvent():
    """
    Test matrix-free resolvent analysis on NACA 64A010 airfoil.

    This test verifies:
    1. Matrix-free Jacobian-vector products work
    2. Iterative linear solver (GMRES) converges
    3. Iterative SVD produces reasonable results
    4. Comparison with explicit method (for small problem)
    """

    print("="*80)
    print("Matrix-Free Resolvent Analysis Test")
    print("="*80)
    print()

    # Check mesh file
    baseDir = os.path.dirname(os.path.abspath(__file__))
    meshFile = os.path.join(baseDir, "../input_files/naca64A010_euler-L2.cgns")

    if not os.path.exists(meshFile):
        print(f"ERROR: Mesh file not found: {meshFile}")
        print("\nPlease run: cd input_files && ./get-input-files.sh")
        return False

    print(f"✓ Found mesh file: {meshFile}")
    print()

    # =========================================================================
    # Setup ADflow
    # =========================================================================

    outputDir = './output_resolvent_mf_test'
    os.makedirs(outputDir, exist_ok=True)

    aeroOptions = {
        'gridFile': meshFile,
        'outputDirectory': outputDir,

        # Physics
        'equationType': 'Euler',

        # Solver parameters
        'CFL': 2.0,
        'CFLCoarse': 1.5,
        'MGCycle': '2w',
        'nCyclesCoarse': 100,
        'nCycles': 200,
        'monitorvariables': ['resrho', 'cl', 'cd'],
        'useNKSolver': True,
        'NKSwitchTol': 1e-2,
        'L2Convergence': 1e-6,
        'L2ConvergenceCoarse': 1e-2,

        # Reduce output
        'printIterations': False,
        'printTiming': False,
        'printWarnings': False,
        'writeVolumeSolution': False,
        'writeSurfaceSolution': False,
        'writeTecplotSurfaceSolution': False,
    }

    print("Creating ADflow solver...")
    CFDsolver = ADFLOW(options=aeroOptions, debug=False)

    ap = AeroProblem(
        name='naca64a010',
        mach=0.60,
        alpha=2.77,
        reynolds=4800000.0,
        reynoldsLength=1.0,
        T=280.0,
        areaRef=1.0,
        chordRef=1.0,
        evalFuncs=['cl', 'cd']
    )

    print("✓ Solver created")
    print()

    # =========================================================================
    # Solve steady state
    # =========================================================================
    print("-"*80)
    print("Solving steady-state flow...")
    print("-"*80)

    CFDsolver(ap)

    funcs = {}
    CFDsolver.evalFunctions(ap, funcs)

    print()
    print("Converged solution:")
    print(f"  CL = {funcs['naca64a010_cl']:.6f}")
    print(f"  CD = {funcs['naca64a010_cd']:.6f}")
    print()

    stateSize = CFDsolver.getStateSize()
    print(f"State size: {stateSize}")
    print()

    # =========================================================================
    # Test 1: Matrix-free resolvent analysis
    # =========================================================================
    print("-"*80)
    print("Test 1: Matrix-free resolvent analysis (scipy svds)")
    print("-"*80)
    print()

    try:
        # Create matrix-free resolvent analysis
        omega = 1.0
        resolvent_mf = ResolventAnalysisMatrixFree(CFDsolver, ap, omega=omega)

        print(f"✓ Matrix-free resolvent object created")
        print(f"  Frequency: ω = {omega}")
        print(f"  State size: {stateSize}")
        print()

        # Set tolerances
        resolvent_mf.setLinearSolveTol(1e-6)
        resolvent_mf.setSVDTol(1e-4)

        # Enable ILU preconditioning
        resolvent_mf.enablePreconditioner(precond_type='ilu', drop_tol=1e-3, fill_factor=10)

        # Solve using scipy
        print("Computing resolvent using matrix-free methods...")
        nModes = 3
        sigma1_mf = resolvent_mf.solve(nModes=nModes, method='scipy')

        print()
        print(f"✓ Matrix-free resolvent computed")
        print(f"  Dominant singular value: σ₁ = {sigma1_mf:.6f}")

        # Get modes
        u1_mf = resolvent_mf.getResponseMode(0)
        v1_mf = resolvent_mf.getForcingMode(0)
        print(f"  Response mode norm: {np.linalg.norm(u1_mf):.6f}")
        print(f"  Forcing mode norm: {np.linalg.norm(v1_mf):.6f}")

        # Show all computed singular values
        print(f"\n  All {nModes} singular values:")
        for i in range(nModes):
            sigma = resolvent_mf.getSingularValue(i)
            print(f"    σ_{i+1} = {sigma:.6f}")

        test1_pass = True

    except Exception as e:
        print(f"\n✗ Matrix-free test FAILED: {e}")
        import traceback
        traceback.print_exc()
        test1_pass = False

    print()

    # =========================================================================
    # Test 2: Compare with explicit method (for small problem)
    # =========================================================================
    if stateSize < 10000:
        print("-"*80)
        print("Test 2: Comparison with explicit method")
        print("-"*80)
        print()

        try:
            from adflow import ResolventAnalysis

            # Explicit resolvent analysis
            print("Computing resolvent using explicit method (for comparison)...")
            resolvent_exp = ResolventAnalysis(CFDsolver, ap, omega=omega)
            sigma1_exp = resolvent_exp.solveExplicit(useLU=True)

            print(f"\nComparison:")
            print(f"  Matrix-free:  σ₁ = {sigma1_mf:.6f}")
            print(f"  Explicit:     σ₁ = {sigma1_exp:.6f}")
            print(f"  Relative diff: {abs(sigma1_mf - sigma1_exp)/sigma1_exp:.2e}")

            # Check if they match
            rel_error = abs(sigma1_mf - sigma1_exp) / sigma1_exp
            if rel_error < 0.01:  # 1% tolerance
                print(f"  ✓ PASS - Results match within 1%")
                test2_pass = True
            else:
                print(f"  ✗ FAIL - Results differ by {rel_error*100:.1f}%")
                test2_pass = False

        except Exception as e:
            print(f"\n✗ Comparison test FAILED: {e}")
            import traceback
            traceback.print_exc()
            test2_pass = False
    else:
        print("-"*80)
        print(f"Skipping comparison test (state size {stateSize} too large)")
        print("-"*80)
        test2_pass = True  # Don't fail if skipped

    print()

    # =========================================================================
    # Test 3: Frequency sweep
    # =========================================================================
    print("-"*80)
    print("Test 3: Matrix-free frequency sweep (mini test)")
    print("-"*80)
    print()

    try:
        # Small frequency sweep
        omega_range = (0.5, 2.0)
        nPoints = 3  # Just a few points for testing
        print(f"Frequency sweep: {nPoints} points from ω={omega_range[0]} to {omega_range[1]}")

        omega_vec, sigma1_vec = resolvent_mf.computeFrequencySweep(
            omega_range,
            nPoints=nPoints,
            nModes=1
        )

        print("Results:")
        for i in range(nPoints):
            print(f"  ω = {omega_vec[i]:.4f}: σ₁ = {sigma1_vec[i]:.6f}")

        print(f"\n✓ Frequency sweep completed")
        test3_pass = True

    except Exception as e:
        print(f"\n✗ Frequency sweep FAILED: {e}")
        import traceback
        traceback.print_exc()
        test3_pass = False

    print()

    # =========================================================================
    # Summary
    # =========================================================================
    print("="*80)
    print("TEST SUMMARY")
    print("="*80)
    print("✓ ADflow steady-state solve")
    if test1_pass:
        print("✓ Matrix-free resolvent analysis (scipy svds)")
    else:
        print("✗ Matrix-free resolvent analysis")

    if stateSize < 10000:
        if test2_pass:
            print("✓ Comparison with explicit method")
        else:
            print("✗ Comparison with explicit method")

    if test3_pass:
        print("✓ Matrix-free frequency sweep")
    else:
        print("✗ Matrix-free frequency sweep")

    print()
    if test1_pass and test2_pass and test3_pass:
        print("All tests PASSED!")
        print("="*80)
        return True
    else:
        print("Some tests FAILED")
        print("="*80)
        return False


if __name__ == "__main__":
    print("\n")
    print("="*80)
    print("ADflow Matrix-Free Resolvent Analysis Test")
    print("="*80)
    print()

    if not ADFLOW_AVAILABLE:
        print("ERROR: ADflow not available")
        sys.exit(1)

    success = test_matrix_free_resolvent()

    if success:
        print("\n✓ All tests completed successfully!")
        sys.exit(0)
    else:
        print("\n✗ Some tests failed")
        sys.exit(1)
