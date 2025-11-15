"""
Test the new getInstanceForces() API
"""
import numpy as np
from adflow import ADFLOW
from idwarp import USMesh
from baseclasses import AeroProblem

print("Testing getInstanceForces() API")
print("="*70)

gridFile = '../input_files/naca64A010_euler-L2.cgns'

# Torus configuration
n1 = 3
n2 = 3
ntot = n1 * n2

omega1 = 100.0
omega2 = 100.0
A1 = 1.0
A2 = 1.0

options = {
    'gridfile': gridFile,
    'outputDirectory': './output_api_test',
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
    'nCycles': 100,  # Quick test
    'l2convergence': 1e-6,
    'monitorvariables': ['resrho', 'cl', 'cd'],
    'printIterations': False,
    'writeVolumeSolution': False,
    'writeSurfaceSolution': False,
}

meshOptions = {"gridFile": gridFile}

print("\nInitializing solver...")
CFDSolver = ADFLOW(options=options, debug=False)
mesh = USMesh(options=meshOptions)
CFDSolver.setMesh(mesh)

# Setup alpha grid
theta1 = np.linspace(0.0, 2.0 * np.pi, n1, endpoint=False)
theta2 = np.linspace(0.0, 2.0 * np.pi, n2, endpoint=False)
alpha_grid = np.zeros((n1, n2))

for i1 in range(n1):
    for i2 in range(n2):
        alpha_deg = A1 * np.sin(theta1[i1]) + A2 * np.sin(theta2[i2])
        alpha_grid[i1, i2] = alpha_deg * np.pi / 180.0

# Deform meshes
print("Deforming meshes...")
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

for sps in range(ntot):
    CFDSolver.mesh.setSurfaceCoordinates(cfdPts_list[sps])
    CFDSolver.mesh.warpMesh()
    m = CFDSolver.mesh.getSolverGrid()
    CFDSolver.adflow.warping.setgridforoneinstance(m, sps=sps + 1)

CFDSolver._updateGeomInfo = True
CFDSolver.updateGeometryInfo()

# Solve
print("Solving...")
ap = AeroProblem(
    name='api_test',
    mach=0.5,
    altitude=10000,
    areaRef=1.0,
    chordRef=1.0,
    evalFuncs=['cl', 'cd']
)

CFDSolver(ap)

# Test the NEW API
print("\n" + "="*70)
print("TESTING NEW getInstanceForces() API")
print("="*70)

try:
    forces = CFDSolver.getInstanceForces()

    print("\nSuccessfully called getInstanceForces()!")
    print(f"Returned keys: {list(forces.keys())}")
    print(f"CL array shape: {forces['cl'].shape}")
    print(f"CD array shape: {forces['cd'].shape}")
    print(f"CMz array shape: {forces['cmz'].shape}")

    print("\n" + "="*70)
    print("INSTANCE-SPECIFIC FORCES FROM NEW API")
    print("="*70)
    print("\nInstance  (i1,i2)    α(deg)        CL              CD              CMz")
    print("-" * 80)

    for sps in range(ntot):
        i1 = sps % n1
        i2 = sps // n1
        alpha_deg = alpha_grid[i1, i2] * 180 / np.pi

        cl = forces['cl'][sps]
        cd = forces['cd'][sps]
        cmz = forces['cmz'][sps]

        print(f"   {sps+1:2d}     ({i1},{i2})   {alpha_deg:7.3f}°   {cl:15.10f}  {cd:15.10e}  {cmz:15.10e}")

    # Verify symmetric pairs
    print("\n" + "="*70)
    print("VERIFYING SYMMETRIC PAIRS")
    print("="*70)

    pairs = [(1, 3), (2, 6), (5, 7)]  # 0-indexed
    all_match = True

    for idx1, idx2 in pairs:
        cl1 = forces['cl'][idx1]
        cl2 = forces['cl'][idx2]
        diff = abs(cl1 - cl2)

        i1_1, i2_1 = idx1 % n1, idx1 // n1
        i1_2, i2_2 = idx2 % n1, idx2 // n1

        status = "✅ MATCH" if diff < 1e-10 else "❌ MISMATCH"
        print(f"\nInstance {idx1+1} ({i1_1},{i2_1}) vs Instance {idx2+1} ({i1_2},{i2_2}):")
        print(f"  CL₁ = {cl1:.10f}")
        print(f"  CL₂ = {cl2:.10f}")
        print(f"  Δ   = {diff:.2e} {status}")

        if diff >= 1e-10:
            all_match = False

    print("\n" + "="*70)
    if all_match:
        print("✅ SUCCESS: New API works correctly!")
        print("   All symmetric pairs match at machine precision")
    else:
        print("⚠️  WARNING: Some symmetric pairs don't match")
    print("="*70)

except Exception as e:
    print(f"\n❌ ERROR calling getInstanceForces(): {e}")
    import traceback
    traceback.print_exc()
