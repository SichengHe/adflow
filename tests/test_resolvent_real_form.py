#!/usr/bin/env python
"""
Simple test of resolvent analysis using REAL-VALUED formulation.

This tests the real doubled form where the complex system:
  (iω·I - J) u = f

Is converted to a real system:
  [ -J   -ωI ] [ u_real ]   [ f_real ]
  [  ωI  -J  ] [ u_imag ] = [ f_imag ]

This avoids complex arithmetic at the cost of doubling the system size.
"""

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse.linalg import LinearOperator, svds

print("="*80)
print("Resolvent Analysis: Real-Valued Formulation Test")
print("="*80)
print()

# Simple 2x2 example from the paper
n = 2
J = np.array([[-0.1, 1.0],
              [-1.0, -0.1]])

omega = 1.0

print(f"Testing real-valued formulation")
print(f"System size: n = {n}")
print(f"Frequency: ω = {omega}")
print()
print("Jacobian J:")
print(J)
print()

# ============================================================================
# Method 1: Complex form (for comparison)
# ============================================================================
print("-"*80)
print("Method 1: Complex formulation")
print("-"*80)

A_complex = 1j * omega * np.eye(n) - J
print(f"A = iω·I - J")
print(f"  Shape: {A_complex.shape}")
print(f"  Type: complex")

# LU decomposition
lu_c, piv_c = lu_factor(A_complex)

# Define operators
def matvec_complex(f):
    return lu_solve((lu_c, piv_c), f)

def rmatvec_complex(f):
    return lu_solve((lu_c, piv_c), f, trans=2)

R_op_complex = LinearOperator(
    (n, n),
    matvec=matvec_complex,
    rmatvec=rmatvec_complex,
    dtype=complex
)

# SVD - use full SVD for small matrices
print("Computing SVD (complex)...")
# For small n, use full SVD instead of svds
if n <= 5:
    # Build explicit matrix for small case
    R_explicit = np.zeros((n, n), dtype=complex)
    for i in range(n):
        e = np.zeros(n, dtype=complex)
        e[i] = 1.0
        R_explicit[:, i] = matvec_complex(e)
    U_c, S_c, Vh_c = np.linalg.svd(R_explicit)
    sigma1_complex = S_c[0]
else:
    U_c, S_c, Vh_c = svds(R_op_complex, k=1, which='LM')
    sigma1_complex = S_c[0]

print(f"✓ Dominant singular value: σ₁ = {sigma1_complex:.6f}")
print()

# ============================================================================
# Method 2: Real doubled form
# ============================================================================
print("-"*80)
print("Method 2: Real doubled formulation")
print("-"*80)

# Create real doubled matrix
A_real = np.zeros((2*n, 2*n), dtype=float)
A_real[:n, :n] = -J
A_real[:n, n:] = -omega * np.eye(n)
A_real[n:, :n] = omega * np.eye(n)
A_real[n:, n:] = -J

print(f"A_real = [ -J   -ωI ]")
print(f"         [  ωI  -J  ]")
print(f"  Shape: {A_real.shape}")
print(f"  Type: float")
print()
print("A_real matrix:")
print(A_real)
print()

# LU decomposition
lu_r, piv_r = lu_factor(A_real)

# Define operators
print("Setting up LinearOperator with matvec and rmatvec...")

def matvec_real(f):
    """Apply R*f = A_real^{-1}*f"""
    u = lu_solve((lu_r, piv_r), f)
    return u

def rmatvec_real(f):
    """Apply R^T*f = (A_real^{-1})^T*f by solving A_real^T*u = f"""
    u = lu_solve((lu_r, piv_r), f, trans=1)  # trans=1 = transpose
    return u

R_op_real = LinearOperator(
    (2*n, 2*n),
    matvec=matvec_real,
    rmatvec=rmatvec_real,
    dtype=float
)

print("✓ LinearOperator created")
print(f"  Shape: {R_op_real.shape}")
print(f"  Has rmatvec: {R_op_real.rmatvec is not None}")
print()

# SVD - use full SVD for small matrices
print("Computing SVD (real doubled form)...")
print("This will call matvec and rmatvec multiple times...")
print()

# For small n, use full SVD instead of svds
if n <= 5:
    # Build explicit matrix for small case
    R_explicit_real = np.zeros((2*n, 2*n), dtype=float)
    for i in range(2*n):
        e = np.zeros(2*n, dtype=float)
        e[i] = 1.0
        R_explicit_real[:, i] = matvec_real(e)
    U_r, S_r, Vh_r = np.linalg.svd(R_explicit_real)
    sigma1_real = S_r[0]
else:
    U_r, S_r, Vh_r = svds(R_op_real, k=1, which='LM')
    sigma1_real = S_r[0]

print()
print(f"✓ Dominant singular value: σ₁ = {sigma1_real:.6f}")
print()

# ============================================================================
# Comparison
# ============================================================================
print("="*80)
print("COMPARISON")
print("="*80)
print(f"Complex form:       σ₁ = {sigma1_complex:.6f}")
print(f"Real doubled form:  σ₁ = {sigma1_real:.6f}")
print(f"Relative difference: {abs(sigma1_complex - sigma1_real)/sigma1_complex:.2e}")
print()

if abs(sigma1_complex - sigma1_real)/sigma1_complex < 1e-6:
    print("✓ PASS: Real and complex forms agree!")
else:
    print("✗ FAIL: Results differ!")

print()
print("="*80)
print("Key observations:")
print("="*80)
print("1. Real form doubles system size: 2n × 2n")
print("2. Avoids complex arithmetic")
print("3. Should give identical results to complex form")
print("4. Requires rmatvec for iterative SVD")
print()
