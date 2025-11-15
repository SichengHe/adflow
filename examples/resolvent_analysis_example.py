"""
Example script for resolvent analysis with ADflow.

This example demonstrates how to perform resolvent analysis on a converged
steady-state CFD solution. Resolvent analysis is used to study the frequency
response characteristics of the flow and can identify amplification mechanisms.

Based on the paper:
"Large-Scale Flow Control Performance Optimization via Differentiable
Resolvent Analysis" by He et al.

The resolvent operator is: R(ω) = (jω*I - J)^{-1}
where J is the Jacobian at the steady-state solution.

Author: ADflow development team
"""

import numpy as np
from adflow import ADFLOW, ResolventAnalysis
from baseclasses import AeroProblem

# ==============================================================================
# Setup ADflow options
# ==============================================================================
aeroOptions = {
    # Common Parameters
    'gridFile': 'wing.cgns',  # Grid file
    'outputDirectory': './',

    # Physics Parameters
    'equationType': 'RANS',
    'turbulenceModel': 'SA',
    'turbulenceProduction': 'vorticity',

    # Common Parameters
    'CFL': 1.5,
    'CFLCoarse': 1.25,
    'MGCycle': '3w',
    'MGStartLevel': -1,
    'nCyclesCoarse': 250,
    'nCycles': 1000,
    'monitorvariables': ['resrho', 'cl', 'cd'],
    'useNKSolver': True,
    'NKSwitchTol': 1e-4,
    'NKSubSpaceSize': 60,
    'NKLinearSolveTol': 0.3,
    'L2Convergence': 1e-10,
    'L2ConvergenceCoarse': 1e-2,

    # Adjoint Parameters
    'adjointL2Convergence': 1e-10,
}

# ==============================================================================
# Create solver
# ==============================================================================
CFDsolver = ADFLOW(options=aeroOptions)

# ==============================================================================
# Set up aerodynamic problem
# ==============================================================================
ap = AeroProblem(
    name='wing',
    mach=0.8,
    altitude=10000,
    alpha=2.0,
    areaRef=45.5,
    chordRef=3.25,
    evalFuncs=['cl', 'cd']
)

# ==============================================================================
# Solve the steady-state flow
# ==============================================================================
print("\n" + "="*80)
print("Solving steady-state flow...")
print("="*80)

CFDsolver(ap)

# Get converged solution
funcs = {}
CFDsolver.evalFunctions(ap, funcs)
print(f"\nConverged solution:")
print(f"  CL = {funcs[ap.name + '_cl']:.6f}")
print(f"  CD = {funcs[ap.name + '_cd']:.6f}")

# ==============================================================================
# Perform resolvent analysis
# ==============================================================================
print("\n" + "="*80)
print("Performing resolvent analysis...")
print("="*80)

# Create resolvent analysis object
omega0 = 10.0  # Initial excitation frequency
resolvent = ResolventAnalysis(CFDsolver, ap, omega=omega0)

print(f"\nResolvent analysis setup:")
print(f"  State size: {resolvent.stateSize}")
print(f"  Excitation frequency: ω = {omega0}")

# Note: The current implementation requires additional development for
# large-scale problems. The following demonstrates the intended usage:

# Single frequency analysis
try:
    sigma1 = resolvent.solveExplicit()
    print(f"\nDominant singular value σ₁ = {sigma1:.6f}")

    # Get modes
    u1 = resolvent.getResponseMode()
    v1 = resolvent.getForcingMode()
    print(f"  Response mode shape: {u1.shape}")
    print(f"  Forcing mode shape: {v1.shape}")

except NotImplementedError as e:
    print(f"\nNOTE: {e}")
    print("\nThe resolvent analysis module is currently in development.")
    print("Full integration with ADflow's Jacobian infrastructure is required.")

# Frequency sweep analysis (if implemented)
try:
    omega_range = (0.0, 100.0)
    nPoints = 50

    print(f"\nPerforming frequency sweep:")
    print(f"  Range: ω ∈ [{omega_range[0]}, {omega_range[1]}]")
    print(f"  Number of points: {nPoints}")

    omega_vec, sigma1_vec = resolvent.computeFrequencySweep(omega_range, nPoints)

    # Find peak frequency
    idx_max = np.argmax(sigma1_vec)
    omega_peak = omega_vec[idx_max]
    sigma1_peak = sigma1_vec[idx_max]

    print(f"\nPeak response:")
    print(f"  ω_peak = {omega_peak:.4f}")
    print(f"  σ₁_max = {sigma1_peak:.6f}")

    # Save results
    np.savetxt('resolvent_sweep.dat',
               np.column_stack([omega_vec, sigma1_vec]),
               header='omega sigma1')
    print("\nResults saved to 'resolvent_sweep.dat'")

except NotImplementedError as e:
    print(f"\nNOTE: Frequency sweep requires full implementation")

# ==============================================================================
# Matrix-free resolvent analysis (for large problems)
# ==============================================================================
print("\n" + "="*80)
print("Matrix-free resolvent analysis (recommended for large problems)")
print("="*80)

try:
    from adflow import ResolventAnalysisMatrixFree

    resolvent_mf = ResolventAnalysisMatrixFree(CFDsolver, ap, omega=omega0)
    resolvent_mf.setNModes(3)  # Compute top 3 modes

    sigma1 = resolvent_mf.solve()
    print(f"Dominant singular value: σ₁ = {sigma1:.6f}")

except NotImplementedError as e:
    print(f"\nNOTE: {e}")
    print("\nMatrix-free implementation requires:")
    print("  1. Interface with ADflow's linear solver")
    print("  2. Adjoint Jacobian-vector products")
    print("  3. Iterative SVD algorithms")

print("\n" + "="*80)
print("Example complete!")
print("="*80)

# ==============================================================================
# Additional notes and usage guidance
# ==============================================================================
print("\n")
print("IMPLEMENTATION STATUS:")
print("-" * 80)
print("The resolvent analysis framework has been created with:")
print("  ✓ Class structure and API design")
print("  ✓ Integration with ADflow solver infrastructure")
print("  ✓ Documentation based on the resolvent paper")
print("")
print("To complete the implementation, the following is needed:")
print("  • Interface to ADflow's dRdW (Jacobian) assembly")
print("  • Integration with PETSc linear solvers")
print("  • Matrix-free operator implementations")
print("  • Iterative SVD using Arnoldi/Lanczos methods")
print("  • Adjoint Jacobian-vector products for efficiency")
print("")
print("For the algebraic examples in the paper (see article_resolvent_opt/),")
print("the implementation is complete and can be used as reference.")
print("-" * 80)
