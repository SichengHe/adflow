"""
Truly Degenerate Torus Case: omega1 = omega2 AND A1 = A2

This should reduce to 1D time spectral because:
- alpha(theta1, theta2) = A * [sin(theta1) + sin(theta2)]
- With omega1 = omega2 = omega, we get:
  alpha(t) = A * [sin(omega*t + phi1) + sin(omega*t + phi2)]
           = 2*A * sin(omega*t + (phi1+phi2)/2) * cos((phi1-phi2)/2)

For grid points with same (theta1 + theta2), alpha is determined by
the cos((theta1-theta2)/2) factor, but the time-periodic behavior
should show 1D pattern.
"""

import numpy as np
from mpi4py import MPI
from adflow import ADFLOW
from idwarp import USMesh
from baseclasses import AeroProblem
import sys
import io
import re

print("="*70)
print("DEGENERATE TORUS: omega1 = omega2 AND A1 = A2")
print("="*70)

gridFile = '../input_files/naca64A010_euler-L2.cgns'

# Torus grid
n1 = 3
n2 = 3
ntot = n1 * n2

# DEGENERATE CASE: same omega AND same amplitude
omega1 = 100.0
omega2 = 100.0  # SAME
alpha_mean = 0.0
A1 = 1.0
A2 = 1.0  # SAME AS A1!

print(f"\nGrid: {n1} × {n2} = {ntot} instances")
print(f"omega1 = {omega1} rad/s")
print(f"omega2 = {omega2} rad/s  ← SAME!")
print(f"A1 = {A1}°")
print(f"A2 = {A2}°  ← SAME!")
print(f"\nMotion: α(θ₁,θ₂) = {A1}°·[sin(θ₁) + sin(θ₂)]")
print("This should show 1D-like pattern!")

options = {
    'gridfile': gridFile,
    'outputDirectory': './output_torus_degenerate',
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
print("Computing: α = A·[sin(θ₁) + sin(θ₂)]")
print()
print("         θ₂=0      θ₂=2π/3   θ₂=4π/3    |  (θ₁+θ₂)mod2π")
print("-" * 70)

for i1 in range(n1):
    print(f"θ₁={theta1[i1]:.2f}:", end="")
    for i2 in range(n2):
        alpha_deg = A1 * np.sin(theta1[i1]) + A2 * np.sin(theta2[i2])
        alpha_grid[i1, i2] = alpha_deg * np.pi / 180.0
        print(f"  {alpha_deg:7.2f}°", end="")

    # Show combined phase for each row
    print("   |  ", end="")
    for i2 in range(n2):
        combined = (theta1[i1] + theta2[i2]) % (2*np.pi)
        print(f"{combined:.2f} ", end="")
    print()

print("\nNote: With A1=A2, points with same (θ₁+θ₂) have related α values")
print("      α = 2A·sin((θ₁+θ₂)/2)·cos((θ₁-θ₂)/2)")

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

# Solve - capture stdout to extract instance-specific CL values
print("\nSolving...")
ap = AeroProblem(
    name='torus_degen',
    mach=0.5,
    altitude=10000,
    areaRef=1.0,
    chordRef=1.0,
    evalFuncs=['cl', 'cd']
)

# Capture stdout
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

CFDSolver(ap)

# Restore stdout
sys.stdout = old_stdout
output_text = captured_output.getvalue()

# Print the captured output
print(output_text)

# Extract forces
funcs = {}
CFDSolver.evalFunctions(ap, funcs)

cl_mean = funcs.get('torus_degen_cl', 0.0)
cd_mean = funcs.get('torus_degen_cd', 0.0)

print("\n" + "="*70)
print("RESULTS")
print("="*70)
print(f"\nMean CL = {cl_mean:.8f}")
print(f"Mean CD = {cd_mean:.8f}")

# Parse instance-specific CL values from solver output
print("\n" + "="*70)
print("EXTRACTING INSTANCE-SPECIFIC FORCE COEFFICIENTS")
print("="*70)

# Extract the last NK iteration line for each instance
# Format: 1  <sps>  <iter>  <totalIter>  NK  ...  <resrho>  <CL>  <CD>  <residual>
instance_data = {}
pattern = r'^\s*1\s+(\d+)\s+\d+\s+\d+\s+NK\s+[-\w]+\s+[\d\.]+\s+[\d\.]+\s+([\d\.\-E\+]+)\s+([\d\.\-E\+]+)\s+([\d\.\-E\+]+)'

for line in output_text.split('\n'):
    match = re.search(pattern, line)
    if match:
        sps = int(match.group(1))
        cl_val = float(match.group(2))
        cd_val = float(match.group(3))
        instance_data[sps] = {'cl': cl_val, 'cd': cd_val}

print(f"\nExtracted force coefficients for {len(instance_data)} instances\n")

# Display results with alpha comparison
print("Instance  (i1,i2)    α(deg)        CL              CD")
print("-" * 65)

for sps in range(1, ntot + 1):
    idx = sps - 1  # 0-indexed
    i1 = idx % n1
    i2 = idx // n1
    alpha_deg = alpha_grid[i1, i2] * 180 / np.pi

    if sps in instance_data:
        cl = instance_data[sps]['cl']
        cd = instance_data[sps]['cd']
        print(f"   {sps:2d}     ({i1},{i2})   {alpha_deg:7.3f}°   {cl:15.8f}  {cd:15.8e}")
    else:
        print(f"   {sps:2d}     ({i1},{i2})   {alpha_deg:7.3f}°   [NOT FOUND]")

# Group instances by alpha and check for duplicates
print("\n" + "="*70)
print("DEGENERATE CASE VERIFICATION")
print("="*70)

# Round alpha to identify duplicates
alpha_groups = {}
for sps in range(1, ntot + 1):
    idx = sps - 1
    i1 = idx % n1
    i2 = idx // n1
    alpha_deg = alpha_grid[i1, i2] * 180 / np.pi
    alpha_key = round(alpha_deg, 4)

    if alpha_key not in alpha_groups:
        alpha_groups[alpha_key] = []
    alpha_groups[alpha_key].append(sps)

print("\nInstances grouped by angle of attack:")
print()

all_match = True
for alpha_key in sorted(alpha_groups.keys()):
    instances = alpha_groups[alpha_key]
    print(f"α = {alpha_key:7.3f}°  →  Instances: {instances}")

    if len(instances) > 1:
        # Check if CL values match
        cl_values = [instance_data[sps]['cl'] for sps in instances if sps in instance_data]
        if len(cl_values) > 1:
            cl_diff = max(cl_values) - min(cl_values)
            print(f"               CL values: {[f'{cl:.8f}' for cl in cl_values]}")
            print(f"               CL range:  {cl_diff:.2e}", end="")

            if cl_diff < 1e-6:
                print("  ✅ MATCH")
            else:
                print("  ❌ MISMATCH")
                all_match = False
        print()

print("="*70)
if all_match:
    print("✅ SUCCESS: All instances with same α have identical CL!")
    print("   Degenerate torus case VERIFIED - collapses to 1D as expected")
else:
    print("❌ FAILURE: Some instances with same α have different CL")
    print("   This indicates a problem with the torus implementation")
print("="*70)
