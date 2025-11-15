"""
Compare 3×3 degenerate torus vs 6-point 1D time spectral

The 3×3 torus projects to 5 unique angles: 0°, 60°, 120°, 180°, 240°
A 6-point uniform 1D TS has angles: 0°, 60°, 120°, 180°, 240°, 300°

So we can compare the 5 matching points!
"""

import numpy as np
from mpi4py import MPI
from adflow import ADFLOW
from idwarp import USMesh
from baseclasses import AeroProblem


print("="*70)
print("THEORY: 3×3 TORUS vs 6-POINT 1D TS")
print("="*70)
print()

# Verify the matching
n1 = n2 = 3
theta1 = np.linspace(0, 2*np.pi, n1, endpoint=False)
theta2 = np.linspace(0, 2*np.pi, n2, endpoint=False)

# Project all torus points
theta_proj = []
for i in range(n1):
    for j in range(n2):
        theta_proj.append((theta1[i] + theta2[j]) / 2)

theta_unique = np.unique(np.round(np.array(theta_proj) % (2*np.pi), 6))

print("3×3 Torus projected angles:")
print(f"  {theta_unique * 180 / np.pi}")
print()

theta_6pt = np.linspace(0, 2*np.pi, 6, endpoint=False)
print("6-point 1D TS angles:")
print(f"  {theta_6pt * 180 / np.pi}")
print()

# Find matching indices
matching_1d_indices = []
matching_torus_theta = []

for th_torus in theta_unique:
    for idx, th_1d in enumerate(theta_6pt):
        if abs(th_torus - th_1d) < 1e-6:
            matching_1d_indices.append(idx)
            matching_torus_theta.append(th_torus)
            print(f"Match: torus θ={th_torus*180/np.pi:.1f}° ↔ 1D[{idx}] θ={th_1d*180/np.pi:.1f}°")

print()
print(f"Total matches: {len(matching_1d_indices)}/5")
print()

print("="*70)
print("CONCLUSION:")
print("="*70)
print()
print("To compare 3×3 torus with 1D TS:")
print("  1. Run 3×3 torus → get CL at 9 instances")
print("  2. Run 6-point 1D TS → get CL at 6 instances")
print("  3. Compare at the 5 matching angles")
print()
print("Expected: CL should match at matching angles if coupling is similar")
print("(But they likely won't due to different spectral coupling)")
