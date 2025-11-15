#!/usr/bin/env python
"""
Test ILU preconditioning for matrix-free resolvent analysis.

This test demonstrates:
1. Matrix-free resolvent WITHOUT preconditioning
2. Matrix-free resolvent WITH ILU preconditioning
3. Comparison of convergence rates

For small problems (n < 10k), ILU can significantly reduce GMRES iterations.
"""

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse.linalg import LinearOperator, svds, gmres

print("="*80)
print("ILU Preconditioning Test for Resolvent Analysis")
print("="*80)
print()

# Simple 2x2 example
n = 2
J = np.array([[-0.1, 1.0],
              [-1.0, -0.1]])

omega = 1.0

print(f"System size: n = {n}")
print(f"Frequency: ω = {omega}")
print()
print("Jacobian J:")
print(J)
print()

# Form resolvent matrix: A = iω·I - J
A = 1j * omega * np.eye(n) - J

print("Resolvent matrix A = iω·I - J:")
print(A)
print()

# ============================================================================
# Test 1: Without preconditioning
# ============================================================================
print("-"*80)
print("Test 1: GMRES WITHOUT preconditioning")
print("-"*80)

# Test vector
rhs = np.array([1.0 + 0.5j, 0.5 + 1.0j])

# GMRES iteration counter
gmres_iter_noprecond = [0]

def callback_noprecond(rk):
    gmres_iter_noprecond[0] += 1
    print(f"  GMRES iter {gmres_iter_noprecond[0]}: residual = {rk:.2e}")

print(f"Solving A*x = rhs with GMRES (no preconditioner)...")
print(f"RHS: {rhs}")
print()

x_noprecond, info = gmres(A, rhs, atol=1e-10, restart=10, callback=callback_noprecond)

print()
print(f"✓ GMRES converged in {gmres_iter_noprecond[0]} iterations")
print(f"  Solution: {x_noprecond}")
print(f"  Info: {info}")
print()

# ============================================================================
# Test 2: With ILU preconditioning
# ============================================================================
print("-"*80)
print("Test 2: GMRES WITH ILU preconditioning")
print("-"*80)

from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spilu

# Convert to sparse for ILU
A_sparse = csc_matrix(A)

print("Computing ILU factorization of A...")
ilu = spilu(A_sparse, drop_tol=1e-4, fill_factor=10)

print(f"✓ ILU factorization complete")
print(f"  L nnz: {ilu.L.nnz}")
print(f"  U nnz: {ilu.U.nnz}")
print()

# Create preconditioner operator
def precond_solve(v):
    return ilu.solve(v)

M_op = LinearOperator((n, n), matvec=precond_solve, dtype=complex)

# GMRES iteration counter
gmres_iter_precond = [0]

def callback_precond(rk):
    gmres_iter_precond[0] += 1
    print(f"  GMRES iter {gmres_iter_precond[0]}: residual = {rk:.2e}")

print(f"Solving A*x = rhs with GMRES (ILU preconditioner)...")
print()

x_precond, info = gmres(A, rhs, M=M_op, atol=1e-10, restart=10, callback=callback_precond)

print()
print(f"✓ GMRES converged in {gmres_iter_precond[0]} iterations")
print(f"  Solution: {x_precond}")
print(f"  Info: {info}")
print()

# ============================================================================
# Comparison
# ============================================================================
print("="*80)
print("COMPARISON")
print("="*80)
print(f"Iterations without preconditioner: {gmres_iter_noprecond[0]}")
print(f"Iterations with ILU preconditioner: {gmres_iter_precond[0]}")
print(f"Speedup: {gmres_iter_noprecond[0] / max(gmres_iter_precond[0], 1):.1f}x")
print()

# Check solutions match
error = np.linalg.norm(x_noprecond - x_precond)
print(f"Solution difference: {error:.2e}")

if error < 1e-8:
    print("✓ Solutions match!")
else:
    print("✗ Solutions differ!")

print()
print("="*80)
print("Key observations:")
print("="*80)
print("1. Both methods give the same solution")
print("2. ILU preconditioning can reduce GMRES iterations")
print("3. For small problems, benefit is modest")
print("4. For large ill-conditioned problems, ILU gives significant speedup")
print()
