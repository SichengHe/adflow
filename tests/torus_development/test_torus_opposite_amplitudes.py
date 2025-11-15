"""
Torus with Opposite Amplitudes: omega1 = omega2, A2 = -A1

This creates: alpha(theta1, theta2) = A * [sin(theta1) - sin(theta2)]
              = 2A * cos((theta1+theta2)/2) * sin((theta1-theta2)/2)

Let's see if instances with matching alpha have matching CL.
"""

import numpy as np
from adflow import ADFLOW
from idwarp import USMesh
from baseclasses import AeroProblem

print("="*70)
print("TORUS WITH OPPOSITE AMPLITUDES: omega1=omega2, A2=-A1")
print("="*70)

gridFile = '../input_files/naca64A010_euler-L2.cgns'

# Torus grid
n1 = 3
n2 = 3
ntot = n1 * n2

# Opposite amplitudes
omega1 = 100.0
omega2 = 100.0  # SAME frequency
A1 = 1.0
A2 = -1.0  # OPPOSITE amplitude!

print(f"\nGrid: {n1} × {n2} = {ntot} instances")
print(f"omega1 = {omega1} rad/s")
print(f"omega2 = {omega2} rad/s")
print(f"A1 = {A1}°")
print(f"A2 = {A2}°  ← NEGATIVE!")
print(f"\nMotion: α(θ₁,θ₂) = {A1}°·sin(θ₁) + ({A2}°)·sin(θ₂)")
print(f"      = {A1}°·[sin(θ₁) - sin(θ₂)]")
print(f"      = {2*A1}°·cos((θ₁+θ₂)/2)·sin((θ₁-θ₂)/2)")

options = {
    'gridfile': gridFile,
    'outputDirectory': './output_opposite_amp',
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
print(f"{'':>10}", end="")
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

print("\nNote: Diagonal symmetry expected")
print("      α(θ₁,θ₂) = -α(θ₂,θ₁)")

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
    name='torus_opposite',
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

# Check grid velocities for instances with same alpha
print("\n" + "="*70)
print("GRID VELOCITY ANALYSIS")
print("="*70)

# Grid velocity: v ∝ omega1 * d/dtheta1 + omega2 * d/dtheta2
#                  = omega * [cos(theta1) - cos(theta2)]

print("\nGrid velocity factor [cos(θ₁) - cos(θ₂)]:")
print(f"{'':>10}", end="")
for i2 in range(n2):
    print(f"  θ₂={theta2[i2]:.2f}  ", end="")
print()
print("-" * 70)

vel_grid = np.zeros((n1, n2))
for i1 in range(n1):
    print(f"θ₁={theta1[i1]:.2f}:", end="")
    for i2 in range(n2):
        vel_factor = np.cos(theta1[i1]) - np.cos(theta2[i2])
        vel_grid[i1, i2] = vel_factor
        print(f"  {vel_factor:7.2f}   ", end="")
    print()

# Group instances by alpha value
print("\n" + "="*70)
print("INSTANCES GROUPED BY ALPHA")
print("="*70)

alpha_groups = {}
for sps in range(ntot):
    i1 = sps % n1
    i2 = sps // n1
    alpha_deg = alpha_grid[i1, i2] * 180 / np.pi
    alpha_key = round(alpha_deg, 3)

    if alpha_key not in alpha_groups:
        alpha_groups[alpha_key] = []
    alpha_groups[alpha_key].append((sps+1, i1, i2))

print("\nChecking if instances with same α have same CL:")
all_match = True

for alpha_key in sorted(alpha_groups.keys()):
    instances = alpha_groups[alpha_key]
    if len(instances) > 1:
        print(f"\nα = {alpha_key:7.3f}°:")

        for inst, i1, i2 in instances:
            sps = inst - 1
            vel_factor = vel_grid[i1, i2]
            print(f"  Instance {inst} ({i1},{i2}): vel_factor={vel_factor:7.3f}, CL={cl_array[sps]:.10f}")

        # Check CL differences
        cl_values = [cl_array[inst-1] for inst, _, _ in instances]
        cl_diff = max(cl_values) - min(cl_values)

        # Check velocity factor differences
        vel_factors = [vel_grid[i1, i2] for _, i1, i2 in instances]
        vel_diff = max(vel_factors) - min(vel_factors)

        print(f"  CL range: {cl_diff:.2e}", end="")
        if cl_diff < 1e-10:
            print("  ✅ MACHINE PRECISION")
        elif cl_diff < 1e-6:
            print("  ⚠️  CLOSE")
        else:
            print("  ❌ MISMATCH")
            all_match = False

        print(f"  Velocity factor range: {vel_diff:.2e}", end="")
        if vel_diff < 1e-10:
            print("  (velocities match)")
        else:
            print("  (velocities differ!)")

print("\n" + "="*70)
if all_match:
    print("✅ Instances with same α AND same velocity have matching CL")
else:
    print("❌ Some instances with same α have different CL")
    print("   (This is expected if grid velocities differ)")
print("="*70)

# Check anti-symmetric pairs: (i1,i2) vs (i2,i1)
print("\n" + "="*70)
print("ANTI-SYMMETRIC PAIRS: (i1,i2) vs (i2,i1)")
print("="*70)
print("\nChecking if α(i1,i2) = -α(i2,i1) and how CL relates:")

for i1 in range(n1):
    for i2 in range(i1+1, n2):  # Only check upper triangle
        sps1 = i2 * n1 + i1  # (i1, i2)
        sps2 = i1 * n1 + i2  # (i2, i1) - swapped!

        alpha1 = alpha_grid[i1, i2] * 180 / np.pi
        alpha2 = alpha_grid[i2, i1] * 180 / np.pi

        cl1 = cl_array[sps1]
        cl2 = cl_array[sps2]

        print(f"\n({i1},{i2}) vs ({i2},{i1}):")
        print(f"  α: {alpha1:7.3f}° vs {alpha2:7.3f}° (sum={alpha1+alpha2:7.3f}°)")
        print(f"  CL: {cl1:15.10f} vs {cl2:15.10f}")
        print(f"  CL sum: {cl1+cl2:15.10f} (should be ~0 if anti-symmetric)")

        if abs(alpha1 + alpha2) < 0.01 and abs(cl1 + cl2) < 1e-6:
            print("  ✅ Anti-symmetric in both α and CL")
        elif abs(alpha1 + alpha2) < 0.01:
            print("  ⚠️  Anti-symmetric in α but NOT in CL")

print("\n" + "="*70)
