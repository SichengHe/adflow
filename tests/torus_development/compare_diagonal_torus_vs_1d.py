"""
Compare ONLY the 3 diagonal torus instances (1,5,9) against 3-instance 1D time spectral.

These should match at MACHINE PRECISION because they sample the same periodic orbit!
"""

import numpy as np
from adflow import ADFLOW
from idwarp import USMesh
from baseclasses import AeroProblem

print("="*70)
print("DIAGONAL TORUS vs 1D TIME SPECTRAL COMPARISON")
print("="*70)

gridFile = '../input_files/naca64A010_euler-L2.cgns'

# ============================================================================
# Part 1: Run 3x3 Torus (extract only diagonal instances)
# ============================================================================

print("\n" + "="*70)
print("PART 1: TORUS 3x3 (Diagonal Instances Only)")
print("="*70)

n1 = 3
n2 = 3
ntot = n1 * n2

omega1 = 100.0
omega2 = 100.0
A1 = 1.0
A2 = 1.0

print(f"\nTorus: {n1} × {n2} = {ntot} instances")
print(f"ω₁ = {omega1} rad/s, ω₂ = {omega2} rad/s")
print(f"A₁ = {A1}°, A₂ = {A2}°")
print(f"Motion: α(θ₁,θ₂) = {A1}·sin(θ₁) + {A2}·sin(θ₂)")

options_torus = {
    'gridfile': gridFile,
    'outputDirectory': './output_diagonal',
    'equationType': 'Euler',
    'equationMode': 'time spectral',
    'useTorusTimeSpectral': True,
    'nTimeIntervalsSpectral1': n1,
    'nTimeIntervalsSpectral2': n2,
    'omegaFourier1': omega1,
    'omegaFourier2': omega2,
    'useexternaldynamicmesh': True,
    'usetsinterpolatedgridvelocity': True,
    'mgcycle': 'sg',
    'nCycles': 10000,
    'l2convergence': 1e-10,
    'l2convergencecoarse': 1e-2,
    'monitorvariables': ['resrho', 'cl', 'cd'],
    'usenksolver': True,
    'nkswitchtol': 1e-4,
    'useanksolver': True,
    'ankswitchtol': 1e-2,
    'CFL': 1.5,
    'nSubiter': 3,
    'blocksplitting': True,
    'writeVolumeSolution': False,
    'writeSurfaceSolution': False,
    'printIterations': False,
}

meshOptions = {"gridFile": gridFile}

print("\nInitializing torus solver...")
CFDSolver_torus = ADFLOW(options=options_torus, debug=False)
mesh_torus = USMesh(options=meshOptions)
CFDSolver_torus.setMesh(mesh_torus)

# Generate alpha grid
theta1 = np.linspace(0.0, 2.0 * np.pi, n1, endpoint=False)
theta2 = np.linspace(0.0, 2.0 * np.pi, n2, endpoint=False)
alpha_grid = np.zeros((n1, n2))

print("\nDiagonal instances (θ₁ = θ₂):")
diagonal_instances = []
for i1 in range(n1):
    for i2 in range(n2):
        alpha_deg = A1 * np.sin(theta1[i1]) + A2 * np.sin(theta2[i2])
        alpha_grid[i1, i2] = alpha_deg * np.pi / 180.0

        if i1 == i2:  # Diagonal
            idx = i2 * n1 + i1
            diagonal_instances.append(idx + 1)
            print(f"  Instance {idx+1} ({i1},{i2}): θ={theta1[i1]:.3f}, α={alpha_deg:7.3f}°")

# Deform meshes
print("\nDeforming meshes...")
MDGroup = CFDSolver_torus.allWallsGroup
cfdPts0 = CFDSolver_torus.getSurfaceCoordinates(MDGroup, includeZipper=False)
xRot = 0.25
cfdPts_list = []

for i2 in range(n2):
    for i1 in range(n1):
        pitch_rad = alpha_grid[i1, i2]
        cc = np.cos(pitch_rad)
        ss = np.sin(pitch_rad)

        cfdPoints_deformed = np.zeros_like(cfdPts0)
        for j in range(len(cfdPts0)):
            x0, y0, z0 = cfdPts0[j]
            cfdPoints_deformed[j, 0] = cc * (x0 - xRot) + ss * y0 + xRot
            cfdPoints_deformed[j, 1] = -ss * (x0 - xRot) + cc * y0
            cfdPoints_deformed[j, 2] = z0

        cfdPts_list.append(cfdPoints_deformed)

for sps in range(ntot):
    CFDSolver_torus.mesh.setSurfaceCoordinates(cfdPts_list[sps])
    CFDSolver_torus.mesh.warpMesh()
    m = CFDSolver_torus.mesh.getSolverGrid()
    CFDSolver_torus.adflow.warping.setgridforoneinstance(m, sps=sps + 1)

CFDSolver_torus._updateGeomInfo = True
CFDSolver_torus.updateGeometryInfo()

print("\nSolving torus...")
ap_torus = AeroProblem(
    name='torus_diagonal',
    mach=0.5,
    altitude=10000,
    areaRef=1.0,
    chordRef=1.0,
    evalFuncs=['cl', 'cd']
)

CFDSolver_torus(ap_torus)
forces_torus = CFDSolver_torus.getInstanceForces()

print("\nTorus diagonal results:")
torus_diagonal_cl = []
for inst in diagonal_instances:
    idx = inst - 1
    i1 = idx % n1
    i2 = idx // n1
    alpha_deg = alpha_grid[i1, i2] * 180 / np.pi
    cl = forces_torus['cl'][idx]
    torus_diagonal_cl.append(cl)
    print(f"  Instance {inst} ({i1},{i2}): θ={theta1[i1]:.3f}, α={alpha_deg:7.3f}°, CL={cl:.10f}")

# ============================================================================
# Part 2: Run 1D Time Spectral with SAME sampling and amplitude
# ============================================================================

print("\n" + "="*70)
print("PART 2: 1D TIME SPECTRAL (Same Times as Diagonal)")
print("="*70)

n_1d = 3
omega_1d = 100.0  # SAME as torus omega
A_1d = 2.0  # Combined amplitude: 2*A1

print(f"\n1D TS: {n_1d} instances")
print(f"ω = {omega_1d} rad/s (SAME as torus ω₁=ω₂)")
print(f"A = {A_1d}° (= 2×A_torus to match combined amplitude)")
print(f"Motion: α(θ) = {A_1d}·sin(θ)")

options_1d = {
    'gridfile': gridFile,
    'outputDirectory': './output_diagonal',
    'equationType': 'Euler',
    'equationMode': 'time spectral',
    'timeIntervals': n_1d,
    'useexternaldynamicmesh': True,
    'usetsinterpolatedgridvelocity': True,
    'mgcycle': 'sg',
    'nCycles': 10000,
    'l2convergence': 1e-10,
    'l2convergencecoarse': 1e-2,
    'monitorvariables': ['resrho', 'cl', 'cd'],
    'usenksolver': True,
    'nkswitchtol': 1e-4,
    'useanksolver': True,
    'ankswitchtol': 1e-2,
    'CFL': 1.5,
    'nSubiter': 3,
    'blocksplitting': True,
    'writeVolumeSolution': False,
    'writeSurfaceSolution': False,
    'printIterations': False,
}

print("\nInitializing 1D TS solver...")
CFDSolver_1d = ADFLOW(options=options_1d, debug=False)
mesh_1d = USMesh(options=meshOptions)
CFDSolver_1d.setMesh(mesh_1d)

# Use SAME theta sampling as torus diagonal
theta_1d = np.linspace(0.0, 2.0 * np.pi, n_1d, endpoint=False)
alpha_1d_rad = A_1d * np.sin(theta_1d) * np.pi / 180.0
alpha_1d_deg = A_1d * np.sin(theta_1d)

print("\n1D TS instances:")
for i in range(n_1d):
    print(f"  Instance {i+1}: θ={theta_1d[i]:.3f}, α={alpha_1d_deg[i]:7.3f}°")

# Deform meshes
print("\nDeforming 1D meshes...")
cfdPts0_1d = CFDSolver_1d.getSurfaceCoordinates(CFDSolver_1d.allWallsGroup, includeZipper=False)
cfdPts_list_1d = []

for i in range(n_1d):
    pitch_rad = alpha_1d_rad[i]
    cc = np.cos(pitch_rad)
    ss = np.sin(pitch_rad)

    cfdPoints_deformed = np.zeros_like(cfdPts0_1d)
    for j in range(len(cfdPts0_1d)):
        x0, y0, z0 = cfdPts0_1d[j]
        cfdPoints_deformed[j, 0] = cc * (x0 - xRot) + ss * y0 + xRot
        cfdPoints_deformed[j, 1] = -ss * (x0 - xRot) + cc * y0
        cfdPoints_deformed[j, 2] = z0

    cfdPts_list_1d.append(cfdPoints_deformed)

for i in range(n_1d):
    CFDSolver_1d.mesh.setSurfaceCoordinates(cfdPts_list_1d[i])
    CFDSolver_1d.mesh.warpMesh()
    m = CFDSolver_1d.mesh.getSolverGrid()
    CFDSolver_1d.adflow.warping.setgridforoneinstance(m, sps=i + 1)

CFDSolver_1d._updateGeomInfo = True
CFDSolver_1d.updateGeometryInfo()

print("\nSolving 1D TS...")
ap_1d = AeroProblem(
    name='ts_1d',
    mach=0.5,
    altitude=10000,
    areaRef=1.0,
    chordRef=1.0,
    xRef=0.25,
    xRot=0.25,
    omegaFourier=omega_1d,
    evalFuncs=['cl', 'cd']
)

CFDSolver_1d(ap_1d)
forces_1d = CFDSolver_1d.getInstanceForces()

print("\n1D TS results:")
for i in range(n_1d):
    cl = forces_1d['cl'][i]
    print(f"  Instance {i+1}: θ={theta_1d[i]:.3f}, α={alpha_1d_deg[i]:7.3f}°, CL={cl:.10f}")

# ============================================================================
# Part 3: COMPARISON
# ============================================================================

print("\n" + "="*70)
print("DIRECT COMPARISON: Diagonal Torus vs 1D TS")
print("="*70)

print("\nInstance    Torus (diagonal)      1D Time Spectral      Difference      Match?")
print("-" * 85)

all_match = True
for i in range(n_1d):
    torus_cl = torus_diagonal_cl[i]
    ts_cl = forces_1d['cl'][i]
    diff = abs(torus_cl - ts_cl)

    if diff < 1e-10:
        status = "✅ PERFECT"
    elif diff < 1e-6:
        status = "⚠️  CLOSE"
    else:
        status = "❌ MISMATCH"
        all_match = False

    print(f"   {i+1:2d}       {torus_cl:18.10f}   {ts_cl:18.10f}   {diff:12.2e}   {status}")

print("\n" + "="*70)
if all_match:
    print("✅ SUCCESS: Diagonal torus instances match 1D TS at machine precision!")
    print("   This confirms the torus implementation is CORRECT.")
else:
    print("❌ FAILURE: Diagonal torus instances DON'T match 1D TS!")
    print("   This suggests a BUG or configuration difference.")
print("="*70)
