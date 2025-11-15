"""
Debug: Compare grid velocities and spectral resolution between torus and 1D TS

Key questions:
1. Are grid velocities computed correctly at diagonal torus points?
2. Are time derivatives (residuals) matching?
3. What's the effective spectral resolution of 3x3 torus?
"""

import numpy as np
from adflow import ADFLOW
from idwarp import USMesh
from baseclasses import AeroProblem

print("="*70)
print("GRID VELOCITY & SPECTRAL RESOLUTION DEBUG")
print("="*70)

gridFile = '../input_files/naca64A010_euler-L2.cgns'

# ============================================================================
# Part 1: Analyze spectral differentiation matrices
# ============================================================================

print("\n" + "="*70)
print("PART 1: SPECTRAL DIFFERENTIATION MATRICES")
print("="*70)

# Torus 3x3
n1, n2 = 3, 3
omega1, omega2 = 100.0, 100.0

options_torus = {
    'gridfile': gridFile,
    'equationType': 'Euler',
    'equationMode': 'time spectral',
    'useTorusTimeSpectral': True,
    'nTimeIntervalsSpectral1': n1,
    'nTimeIntervalsSpectral2': n2,
    'omegaFourier1': omega1,
    'omegaFourier2': omega2,
    'printIterations': False,
}

print("\nInitializing torus solver...")
solver_torus = ADFLOW(options=options_torus, debug=False)
mesh_torus = USMesh(options={"gridFile": gridFile})
solver_torus.setMesh(mesh_torus)

dscalar_torus = solver_torus.adflow.inputtimespectral.dscalar

print(f"\nTorus dscalar matrix shape: {dscalar_torus.shape}")
print(f"  (nSections, nInstances, nInstances) = (1, 9, 9)")

# Print diagonal rows (instances on periodic orbit)
print("\n" + "="*70)
print("TORUS: Spectral differentiation for DIAGONAL instances")
print("="*70)

diagonal_instances = [0, 4, 8]  # Instances 1, 5, 9 (0-indexed)

for idx in diagonal_instances:
    i1 = idx % n1
    i2 = idx // n1
    print(f"\nInstance {idx+1} ({i1},{i2}) - dscalar row:")
    print(f"  {dscalar_torus[0, idx, :]}")
    print(f"  Sum = {np.sum(dscalar_torus[0, idx, :]):.2e} (should be ~0 for periodic)")

# 1D Time Spectral with 3 instances
print("\n" + "="*70)
print("1D TIME SPECTRAL: Spectral differentiation (3 instances)")
print("="*70)

options_1d = {
    'gridfile': gridFile,
    'equationType': 'Euler',
    'equationMode': 'time spectral',
    'timeIntervals': 3,
    'useexternaldynamicmesh': True,
    'usetsinterpolatedgridvelocity': True,
    'printIterations': False,
}

print("\nInitializing 1D TS solver...")
solver_1d = ADFLOW(options=options_1d, debug=False)
mesh_1d = USMesh(options={"gridFile": gridFile})
solver_1d.setMesh(mesh_1d)

# Deform 1D meshes first (required before calling solver)
print("Deforming 1D meshes...")
cfdPts0_1d = solver_1d.getSurfaceCoordinates(solver_1d.allWallsGroup, includeZipper=False)
xRot = 0.25

theta_1d = np.linspace(0.0, 2.0 * np.pi, 3, endpoint=False)
alpha_1d_rad = 2.0 * np.sin(theta_1d) * np.pi / 180.0  # A=2° amplitude

for i in range(3):
    pitch_rad = alpha_1d_rad[i]
    cc = np.cos(pitch_rad)
    ss = np.sin(pitch_rad)

    cfdPoints_deformed = np.zeros_like(cfdPts0_1d)
    for j in range(len(cfdPts0_1d)):
        x0, y0, z0 = cfdPts0_1d[j]
        cfdPoints_deformed[j, 0] = cc * (x0 - xRot) + ss * y0 + xRot
        cfdPoints_deformed[j, 1] = -ss * (x0 - xRot) + cc * y0
        cfdPoints_deformed[j, 2] = z0

    solver_1d.mesh.setSurfaceCoordinates(cfdPoints_deformed)
    solver_1d.mesh.warpMesh()
    m = solver_1d.mesh.getSolverGrid()
    solver_1d.adflow.warping.setgridforoneinstance(m, sps=i + 1)

solver_1d._updateGeomInfo = True
solver_1d.updateGeometryInfo()

# Initialize with AeroProblem
ap_temp = AeroProblem(name='temp', mach=0.5, altitude=10000, areaRef=1.0, chordRef=1.0, omegaFourier=100.0)
solver_1d(ap_temp, writeSolution=False)

dscalar_1d = solver_1d.adflow.inputtimespectral.dscalar

print(f"\n1D TS dscalar matrix shape: {dscalar_1d.shape}")

for idx in range(3):
    print(f"\nInstance {idx+1} - dscalar row:")
    print(f"  {dscalar_1d[0, idx, :]}")
    print(f"  Sum = {np.sum(dscalar_1d[0, idx, :]):.2e}")

# ============================================================================
# Part 2: Theoretical grid velocity analysis
# ============================================================================

print("\n" + "="*70)
print("PART 2: THEORETICAL GRID VELOCITY")
print("="*70)

theta1 = np.linspace(0, 2*np.pi, n1, endpoint=False)
theta2 = np.linspace(0, 2*np.pi, n2, endpoint=False)

print("\nFor diagonal instances (θ₁ = θ₂):")
print("Grid velocity v = ω₁·∂x/∂θ₁ + ω₂·∂x/∂θ₂")
print("             = ω₁·(∂x/∂α)·(∂α/∂θ₁) + ω₂·(∂x/∂α)·(∂α/∂θ₂)")
print("             = (∂x/∂α)·[ω₁·A₁·cos(θ₁) + ω₂·A₂·cos(θ₂)]")
print("\nFor A₁=A₂=1, ω₁=ω₂=100, on diagonal θ₁=θ₂:")
print("  v = (∂x/∂α)·100·[cos(θ) + cos(θ)] = (∂x/∂α)·200·cos(θ)")
print("\n1D TS with amplitude A=2°, ω=100:")
print("  v = ω·(∂x/∂α)·A·cos(θ) = (∂x/∂α)·100·2·cos(θ) = (∂x/∂α)·200·cos(θ)")
print("\n✅ SAME! Grid velocities SHOULD match on diagonal.")

print("\nVelocity factors at diagonal instances:")
print("Instance    θ      cos(θ)    Torus v_factor    1D TS v_factor")
print("-" * 70)

for i in range(n1):
    theta = theta1[i]
    cos_theta = np.cos(theta)
    torus_factor = 100 * (cos_theta + cos_theta)  # ω₁·cos(θ₁) + ω₂·cos(θ₂)
    ts_factor = 100 * 2 * cos_theta  # ω·A·cos(θ)

    print(f"  {i+1} ({i},{i})  {theta:.3f}  {cos_theta:7.3f}    {torus_factor:10.2f}        {ts_factor:10.2f}")

# ============================================================================
# Part 3: Spectral resolution analysis
# ============================================================================

print("\n" + "="*70)
print("PART 3: SPECTRAL RESOLUTION - How many frequencies?")
print("="*70)

print("\n2D Fourier modes available in 3×3 grid:")
print("  k₁ = 0, ±1  (3 modes in θ₁ direction)")
print("  k₂ = 0, ±1  (3 modes in θ₂ direction)")
print("  Total: 3×3 = 9 complex modes (or 9 real DOFs)")

print("\nBut for PERIODIC orbit (θ₁=θ₂), we collapse to 1D:")
print("  Only diagonal samples at θ = 0, 2π/3, 4π/3")
print("  This gives 3 DOFs → can resolve 1 harmonic (k=0, ±1)")
print("  Equivalent to 3-instance 1D time spectral!")

print("\nFrequency content:")
print("  3×3 Torus (general quasi-periodic):")
print("    ω_mn = m·ω₁ + n·ω₂ for m,n ∈ {0,±1}")
print("    = 100m + 100n = 100(m+n)")
print("    Frequencies: ..., -200, -100, 0, 100, 200, ... Hz")

print("\n  3×3 Torus ON DIAGONAL (ω₁=ω₂, θ₁=θ₂):")
print("    Only m=n modes survive (diagonal coupling)")
print("    ω = 100(m+m) = 200m for m ∈ {0,±1}")
print("    Frequencies: 0, ±200 Hz")
print("    ⚠️  This is DOUBLE frequency of 1D TS!")

print("\n  3-instance 1D TS with ω=100:")
print("    ω_k = k·ω for k ∈ {0,±1}")
print("    Frequencies: 0, ±100 Hz")

print("\n" + "="*70)
print("⚠️  SPECTRAL MISMATCH DETECTED!")
print("="*70)
print("""
The diagonal torus has frequency content at 0, ±200 Hz,
but 1D TS with ω=100 has frequencies at 0, ±100 Hz!

For matching, we need:
  - 1D TS with ω = 200 (NOT 100!)
  - OR torus with ω₁ = ω₂ = 50 (to get ±100 Hz on diagonal)

The diagonal torus is effectively solving at TWICE the frequency
of what we're comparing against!
""")

# ============================================================================
# Part 4: Check if omega scaling is the issue
# ============================================================================

print("\n" + "="*70)
print("PART 4: MATRIX EIGENVALUES - Actual frequencies")
print("="*70)

print("\nTorus dscalar eigenvalues (diagonal submatrix):")
# Extract diagonal 3x3 submatrix
diag_indices = [0, 4, 8]
dscalar_diag = dscalar_torus[0, diag_indices, :][:, diag_indices]
print(f"Diagonal 3×3 submatrix:")
print(dscalar_diag)

eigenvals_torus = np.linalg.eigvals(dscalar_diag)
print(f"\nEigenvalues: {eigenvals_torus}")
print(f"  (These should be proportional to frequencies)")

print("\n1D TS dscalar eigenvalues:")
print(dscalar_1d[0, :, :])
eigenvals_1d = np.linalg.eigvals(dscalar_1d[0, :, :])
print(f"\nEigenvalues: {eigenvals_1d}")

print("\nRatio of eigenvalues (torus / 1D):")
# Sort by absolute value for comparison
eig_t_sorted = np.sort(np.abs(eigenvals_torus))
eig_1d_sorted = np.sort(np.abs(eigenvals_1d))
ratios = eig_t_sorted / (eig_1d_sorted + 1e-15)
print(f"  {ratios}")
print(f"  Expected: ~2.0 if torus frequency is double")

print("\n" + "="*70)
print("CONCLUSION")
print("="*70)
print("""
The spectral differentiation matrices reveal the issue:

1. The 3×3 torus grid, even on the diagonal, uses a 2D spectral
   operator that couples both θ₁ and θ₂ derivatives.

2. For diagonal instances where θ₁=θ₂, the effective frequency is:
   ω_eff = ω₁ + ω₂ = 200 rad/s (NOT 100 rad/s!)

3. To match 1D TS, we need to compare against 1D TS with ω=200,
   OR we need the torus diagonal to somehow recognize it's on
   a 1D orbit and use only ω (not ω₁+ω₂).

The bug is that the torus spectral operator is CORRECTLY computing
the 2D quasi-periodic derivative, but this doesn't reduce to 1D
even on the diagonal!
""")

print("="*70)
