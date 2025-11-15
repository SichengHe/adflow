#!/usr/bin/env python
"""
Quick test to verify adjoint preconditioner is being built and used.
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

print("="*80)
print("Testing Adjoint Preconditioner for Matrix-Free Resolvent")
print("="*80)
print()

# Create a simple 2x2 test
print("Creating simple 2x2 resolvent analysis object...")
print()

# Mock objects for testing
class MockADflow:
    def __init__(self):
        self.adjointSetup = False

    def _setupAdjoint(self):
        self.adjointSetup = True

    def getStateSize(self):
        return 2

    def computeJacobianVectorProductFwd(self, wDot, residualDeriv=True):
        J = np.array([[-0.1, 1.0],
                      [-1.0, -0.1]])
        return J @ wDot

    def computeJacobianVectorProductBwd(self, resBar, wDeriv=True):
        J = np.array([[-0.1, 1.0],
                      [-1.0, -0.1]])
        return J.T @ resBar

    def setupResolventJacobian(self, ap):
        pass

    def getJacobianMatrix(self, outputType="dense"):
        # Return explicit Jacobian for ILU
        J = np.array([[-0.1, 1.0],
                      [-1.0, -0.1]])
        return J

class MockAeroProblem:
    pass

# Create mock solver and problem
CFDsolver = MockADflow()
ap = MockAeroProblem()

# Import resolvent class
from adflow import ResolventAnalysisMatrixFree

# Create resolvent object
omega = 1.0
resolvent = ResolventAnalysisMatrixFree(CFDsolver, ap, omega=omega)

print("✓ Matrix-free resolvent object created")
print(f"  Frequency: ω = {omega}")
print(f"  State size: 2")
print()

# Enable ILU preconditioning
print("Enabling ILU preconditioning...")
resolvent.enablePreconditioner(precond_type='ilu', drop_tol=1e-3, fill_factor=10)
print()

# Check if preconditioners exist before solve
print("Checking preconditioners BEFORE first solve:")
print(f"  Forward preconditioner (_preconditioner):     {resolvent._preconditioner is not None}")
print(f"  Adjoint preconditioner (_preconditioner_adj): {resolvent._preconditioner_adj is not None}")
print()

# Do a simple forward solve to trigger preconditioner build
print("Triggering preconditioner build with first linear solve...")
rhs = np.array([1.0 + 0.5j, 0.5 + 1.0j])
x, info = resolvent._solveLinearSystem(rhs, method='gmres')
print(f"✓ Forward solve completed: info = {info}")
print()

# Check if preconditioners were built
print("Checking preconditioners AFTER first solve:")
print(f"  Forward preconditioner (_preconditioner):     {resolvent._preconditioner is not None}")
print(f"  Adjoint preconditioner (_preconditioner_adj): {resolvent._preconditioner_adj is not None}")
print()

if resolvent._preconditioner is not None and resolvent._preconditioner_adj is not None:
    print("✓ SUCCESS: Both forward and adjoint preconditioners were built!")
    print()
    print("This means:")
    print("  - Forward GMRES will use ILU preconditioner")
    print("  - Adjoint GMRES will use separate ILU preconditioner")
    print("  - Both should converge in ~5-10 iterations")
else:
    print("✗ FAILURE: Preconditioners not built correctly!")
    if resolvent._preconditioner is None:
        print("  - Forward preconditioner is None")
    if resolvent._preconditioner_adj is None:
        print("  - Adjoint preconditioner is None")

print()
print("="*80)
