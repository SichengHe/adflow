"""
True Degenerate Torus: omega2 = 0 (or A2 = 0)

This SHOULD reduce to 1D time spectral because:
- alpha(theta1, theta2) = A1 * sin(theta1) + 0
- Motion is purely 1D periodic in theta1
- All instances with same theta1 should have IDENTICAL CL at machine precision
"""

import numpy as np
from adflow import ADFLOW
from idwarp import USMesh
from baseclasses import AeroProblem

print("="*70)
print("TRUE DEGENERATE TORUS: omega2 = 0, A2 = 0")
print("="*70)

gridFile = '../input_files/naca64A010_euler-L2.cgns'

# Torus grid
n1 = 3
n2 = 3
ntot = n1 * n2

# TRUE DEGENERATE: omega2 = 0!
omega1 = 100.0
omega2 = 0.0  # ZERO! This makes it truly 1D
A1 = 1.0
A2 = 0.0  # ZERO!

print(f"\nGrid: {n1} × {n2} = {ntot} instances")
print(f"omega1 = {omega1} rad/s")
print(f"omega2 = {omega2} rad/s  ← ZERO = True 1D!")
print(f"A1 = {A1}°")
print(f"A2 = {A2}°  ← ZERO!")
print(f"\nMotion: α(θ₁,θ₂) = {A1}°·sin(θ₁) + {A2}°·sin(θ₂) = {A1}°·sin(θ₁)")
print("This is TRULY 1D periodic!")

options = {
    'gridfile': gridFile,
    'outputDirectory': './output_true_degenerate',
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
}

meshOptions = {"gridFile": gridFile}

print("\nInitializing ADflow...")
CFDSolver = ADFLOW(options=options, debug=False)

print("Setting up mesh...")
mesh = USMesh(options=meshOptions)
CFDSolver.setMesh(mesh)

# Generate alpha grid
theta1 = np.linspace(0.0, 2.0 * np.pi, n1, endpoint=False)
theta2 = np.linspace(0.0, 2.0 * np.pi, n2, endpoint=False)

alpha_grid = np.zeros((n1, n2))

print("\nAlpha distribution (degrees):")
print(f"{'θ₂→':>10}", end="")
for i2 in range(n2):
    print(f"  θ₂={theta2[i2]:.2f}  ", end="")
print()
print("-" * 70)

for i1 in range(n1):
    print(f"θ₁={theta1[i1]:.2f}:", end="")
    for i2 in range(n2):
        alpha_deg = A1 * np.sin(theta1[i1]) + A2 * np.sin(theta2[i2])
        alpha_grid[i1, i2] = alpha_deg * np.pi / 180.0
        print(f"  {alpha_deg:7.2f}°  ", end="")
    print()

print("\nNote: All instances in same ROW (same θ₁) have IDENTICAL α")
print("      They SHOULD also have IDENTICAL CL at machine precision!")

# Deform meshes
print("\nDeforming meshes...")
MDGroup = CFDSolver.allWallsGroup
cfdPts0 = CFDSolver.getSurfaceCoordinates(MDGroup, includeZipper=False)

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

# Set warped meshes
for sps in range(ntot):
    CFDSolver.mesh.setSurfaceCoordinates(cfdPts_list[sps])
    CFDSolver.mesh.warpMesh()
    m = CFDSolver.mesh.getSolverGrid()
    CFDSolver.adflow.warping.setgridforoneinstance(m, sps=sps + 1)

CFDSolver._updateGeomInfo = True
CFDSolver.updateGeometryInfo()

# Solve
print("\nSolving...")
ap = AeroProblem(
    name='torus_true_degen',
    mach=0.5,
    altitude=10000,
    areaRef=1.0,
    chordRef=1.0,
    evalFuncs=['cl', 'cd']
)

CFDSolver(ap)

# Extract forces using NEW API
forces = CFDSolver.getInstanceForces()
cl_array = forces['cl']
cd_array = forces['cd']
cmz_array = forces['cmz']

print("\n" + "="*70)
print("RESULTS")
print("="*70)

print("\nInstance  (i1,i2)    θ₁       θ₂        α(deg)       CL              CD")
print("-" * 80)

for sps in range(ntot):
    i1 = sps % n1
    i2 = sps // n1
    alpha_deg = alpha_grid[i1, i2] * 180 / np.pi

    print(f"   {sps+1:2d}     ({i1},{i2})   {theta1[i1]:6.3f}  {theta2[i2]:6.3f}   {alpha_deg:7.3f}°   {cl_array[sps]:15.10f}  {cd_array[sps]:15.10e}")

# Group by theta1 (should have identical CL within each group)
print("\n" + "="*70)
print("VERIFICATION: Instances with Same θ₁ Should Match")
print("="*70)

all_match = True
for i1 in range(n1):
    instances = []
    cl_values = []

    for i2 in range(n2):
        sps = i2 * n1 + i1
        instances.append(sps + 1)
        cl_values.append(cl_array[sps])

    cl_diff = max(cl_values) - min(cl_values)
    alpha_deg = alpha_grid[i1, 0] * 180 / np.pi

    print(f"\nθ₁ = {theta1[i1]:.3f}, α = {alpha_deg:7.3f}°")
    print(f"  Instances: {instances}")
    print(f"  CL values: {[f'{cl:.10f}' for cl in cl_values]}")
    print(f"  CL range:  {cl_diff:.2e}", end="")

    if cl_diff < 1e-10:
        print("  ✅ MACHINE PRECISION MATCH!")
    else:
        print(f"  ❌ MISMATCH (expected < 1e-10)")
        all_match = False

print("\n" + "="*70)
if all_match:
    print("✅ SUCCESS: True degenerate torus (ω₂=0) verified!")
    print("   All instances with same θ₁ have IDENTICAL CL at machine precision")
    print("   This confirms the torus implementation is CORRECT!")
else:
    print("❌ FAILURE: Instances with same θ₁ have different CL")
    print("   This indicates a BUG in the torus implementation")
print("="*70)
