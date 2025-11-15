"""
Test Torus Time Spectral with omega1 = omega2 (degenerate case)

When omega1 = omega2, the torus collapses to 1D time spectral.
The CL distribution should show periodicity/pattern.
"""

import numpy as np
from mpi4py import MPI
from adflow import ADFLOW
from idwarp import USMesh
from baseclasses import AeroProblem

print("="*70)
print("TORUS WITH omega1 = omega2 (DEGENERATE TO 1D)")
print("="*70)

gridFile = '../input_files/naca64A010_euler-L2.cgns'

# Torus grid
n1 = 3
n2 = 3
ntot = n1 * n2

# EQUAL frequencies - should collapse to 1D
omega1 = 100.0
omega2 = 100.0  # SAME as omega1
alpha_mean = 0.0
A1 = 1.0
A2 = 0.5

print(f"\nGrid: {n1} × {n2} = {ntot} instances")
print(f"omega1 = {omega1} rad/s")
print(f"omega2 = {omega2} rad/s  ← SAME!")
print(f"Motion: α(θ₁,θ₂) = {alpha_mean}° + {A1}°·sin(θ₁) + {A2}°·sin(θ₂)")

options = {
    'gridfile': gridFile,
    'outputDirectory': './output_torus_equal_omega',
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
for i1 in range(n1):
    for i2 in range(n2):
        alpha_deg = alpha_mean + A1 * np.sin(theta1[i1]) + A2 * np.sin(theta2[i2])
        alpha_grid[i1, i2] = alpha_deg * np.pi / 180.0

print("\nAlpha distribution (degrees):")
print("         θ₂=0      θ₂=2π/3   θ₂=4π/3")
print("-" * 45)
for i1 in range(n1):
    print(f"θ₁={theta1[i1]:.2f}:", end="")
    for i2 in range(n2):
        print(f"  {alpha_grid[i1,i2]*180/np.pi:7.2f}°", end="")
    print()

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
    name='torus_equal',
    mach=0.5,
    altitude=10000,
    areaRef=1.0,
    chordRef=1.0,
    evalFuncs=['cl', 'cd']
)

CFDSolver(ap)

# Extract forces
funcs = {}
CFDSolver.evalFunctions(ap, funcs)

cl_mean = funcs.get('torus_equal_cl', 0.0)
cd_mean = funcs.get('torus_equal_cd', 0.0)

print("\n" + "="*70)
print("RESULTS")
print("="*70)
print(f"\nMean CL = {cl_mean:.8f}")
print(f"Mean CD = {cd_mean:.8f}")

print("\n" + "="*70)
print("EXPECTED BEHAVIOR:")
print("="*70)
print("With omega1 = omega2, the torus is degenerate (1D time spectral).")
print("The motion only depends on (θ₁ + θ₂) combined.")
print()
print("Since we currently only get mean forces (not instance-specific),")
print("we can't verify the periodic pattern yet.")
print()
print("To verify properly, we would need instance-specific CL values,")
print("which should show that CL(i1, i2) depends only on (θ₁[i1] + θ₂[i2]).")
print("="*70)
