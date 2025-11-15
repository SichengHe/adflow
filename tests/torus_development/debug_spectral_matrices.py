"""
Debug script to compare spectral differentiation matrices between
degenerate torus and 1D time spectral
"""
import numpy as np
from adflow import ADFLOW
from idwarp import USMesh

gridFile = '../input_files/naca64A010_euler-L2.cgns'

print("="*70)
print("COMPARING SPECTRAL MATRICES: Torus vs 1D TS")
print("="*70)

# Test 1: Degenerate Torus (3x3 with omega1=omega2=100)
print("\n### TORUS MODE (3x3, omega1=omega2=100) ###\n")
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

# Access the dscalar matrix from Fortran
dscalar_torus = solver_torus.adflow.inputtimespectral.dscalar

print(f"Torus dscalar shape: {dscalar_torus.shape}")
print(f"  Expected: (nSections=1, ntot=9, ntot=9)")
print()

# Print matrix for instance 1 (i1=0, i2=0, idx=0)
idx = 0  # Instance 1 (0-indexed in Python, 1-indexed in Fortran)
print(f"Torus dscalar[0, {idx}, :] (derivatives at instance 1):")
print(f"  {dscalar_torus[0, idx, :]}")
print()

# Test 2: Regular 1D Time Spectral (3 instances, omega=100)
print("\n### 1D TIME SPECTRAL MODE (3 instances, omega=200) ###\n")
print("Note: Using omega=200 to match torus omega1+omega2=200")
print()
options_1d = {
    'gridfile': gridFile,
    'equationType': 'Euler',
    'equationMode': 'time spectral',
    'timeIntervals': 3,
    'printIterations': False,
    'printAllOptions': False,
}

# Need to set omegaFourier via AeroProblem - but we just want matrix, so use timePeriod
from baseclasses import AeroProblem
omega_1d = 200.0  # Should match torus omega1+omega2 for degenerate case
timePeriod = 2.0 * np.pi / omega_1d
options_1d['timePeriod'] = timePeriod

solver_1d = ADFLOW(options=options_1d, debug=False)
mesh_1d = USMesh(options={"gridFile": gridFile})
solver_1d.setMesh(mesh_1d)

dscalar_1d = solver_1d.adflow.inputtimespectral.dscalar

print(f"1D TS dscalar shape: {dscalar_1d.shape}")
print(f"  Expected: (nSections=1, nInstances=3, nInstances=3)")
print()

# Print matrix for instance 1
idx = 0
print(f"1D TS dscalar[0, {idx}, :] (derivatives at instance 1):")
print(f"  {dscalar_1d[0, idx, :]}")
print()

# Test 3: Check if degenerate torus reduces to 1D pattern
print("\n### COMPARISON ###\n")
print("For degenerate case (omega1=omega2), torus matrix should show 1D pattern")
print()

# Expected: torus instances along main diagonal (i1, i2) where i1+i2=const
# should have same differentiation behavior
print("Torus instances with same (i1+i2) should have similar dscalar rows:")
print()

# Map torus instances to (i1, i2)
n1, n2 = 3, 3
for sum_idx in range(n1 + n2 - 1):
    instances = []
    for i1 in range(n1):
        for i2 in range(n2):
            if (i1 + i2) == sum_idx:
                idx = i2 * n1 + i1
                instances.append((i1, i2, idx))

    if instances:
        print(f"Sum(i1+i2) = {sum_idx}:")
        for i1, i2, idx in instances:
            print(f"  Instance ({i1},{i2}) [idx={idx}]: dscalar[0,{idx},:] = {dscalar_torus[0, idx, :]}")
        print()

# Check omega scaling
print("\n### OMEGA SCALING CHECK ###\n")
print(f"1D TS omega (from AeroProblem): Not set in options (uses timePeriod)")
print(f"Torus omega1: 100.0")
print(f"Torus omega2: 100.0")
print()
print("Note: 1D TS uses timePeriod, torus uses omegaFourier directly")
print("Need to check if scaling is consistent!")
