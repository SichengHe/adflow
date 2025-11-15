import numpy as np
from adflow import ADFLOW
from idwarp import USMesh

gridFile = '../input_files/naca64A010_euler-L2.cgns'

print("="*70)
print("VERIFY: dscalar matrix diagonal elements")
print("="*70)

# Test 1: Torus 3x3
print("\n1. Torus 3x3 (omega1=100, omega2=100)")
options_torus = {
    'gridfile': gridFile,
    'equationType': 'Euler',
    'equationMode': 'time spectral',
    'useTorusTimeSpectral': True,
    'nTimeIntervalsSpectral1': 3,
    'nTimeIntervalsSpectral2': 3,
    'omegaFourier1': 100.0,
    'omegaFourier2': 100.0,
    'printIterations': False,
    'printAllOptions': False,
}

solver_torus = ADFLOW(options=options_torus, debug=False)
mesh_torus = USMesh(options={"gridFile": gridFile})
solver_torus.setMesh(mesh_torus)

dscalar_torus = solver_torus.adflow.inputtimespectral.dscalar

print("\nFull 9x9 dscalar matrix:")
print(dscalar_torus[0, :, :])

print("\nDiagonal elements (should ALL be zero):")
diag = np.diag(dscalar_torus[0, :, :])
print(diag)
print(f"Max absolute diagonal value: {np.max(np.abs(diag)):.2e}")

if np.max(np.abs(diag)) < 1e-15:
    print("✅ All diagonal elements are zero (as expected)")
else:
    print("❌ Some diagonal elements are non-zero (BUG!)")

# Test 2: 1D Time Spectral
print("\n2. 1D Time Spectral (3 instances, omega=100)")
from baseclasses import AeroProblem

options_1d = {
    'gridfile': gridFile,
    'equationType': 'Euler',
    'equationMode': 'time spectral',
    'timeIntervals': 3,
    'timePeriod': 2*np.pi/100.0,
    'printIterations': False,
    'printAllOptions': False,
}

solver_1d = ADFLOW(options=options_1d, debug=False)
mesh_1d = USMesh(options={"gridFile": gridFile})
solver_1d.setMesh(mesh_1d)

dscalar_1d = solver_1d.adflow.inputtimespectral.dscalar

print("\nFull 3x3 dscalar matrix:")
print(dscalar_1d[0, :, :])

print("\nDiagonal elements (should ALL be zero):")
diag_1d = np.diag(dscalar_1d[0, :, :])
print(diag_1d)
print(f"Max absolute diagonal value: {np.max(np.abs(diag_1d)):.2e}")

if np.max(np.abs(diag_1d)) < 1e-15:
    print("✅ All diagonal elements are zero (as expected)")
else:
    print("❌ Some diagonal elements are non-zero (BUG!)")

print("\n" + "="*70)
print("EXPLANATION")
print("="*70)
print("""
The spectral differentiation matrix D computes time derivatives:
  ∂u/∂t ≈ D·u

For periodic functions, D is a circulant matrix from FFT differentiation.
The diagonal elements D[i,i] are ALWAYS zero because:
  - Spectral derivative at point i depends on ALL other points
  - Point i itself contributes zero (derivative formula has sin/cos terms)
  - The spatial terms (convection, diffusion) provide the diagonal contribution

So yes, you're correct: the diagonal should be zero, and we rely on
the spatial discretization to fill in the diagonal of the full Jacobian!
""")
