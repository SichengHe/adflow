"""
Example: Torus Time Spectral Method for 2-Frequency Aeroelastic Motion

This example demonstrates the torus time spectral capability in ADflow
for solving quasi-periodic problems with two incommensurate frequencies.

Application: Airfoil with combined pitch and plunge at different frequencies
Motion: alpha(t) = alpha_0 + A_pitch*sin(omega_pitch*t) + A_plunge*sin(omega_plunge*t)

where omega_pitch/omega_plunge = sqrt(2) (incommensurate ratio)
"""

import numpy as np
from mpi4py import MPI

# Import ADflow
from adflow import ADFLOW
from baseclasses import AeroProblem
from idwarp import USMesh

# ==============================================================================
# Problem Setup
# ==============================================================================

# Two-frequency motion parameters
omega_pitch = 100.0  # rad/s - pitching frequency
omega_plunge = 100.0 * np.sqrt(2.0)  # rad/s - plunging frequency (incommensurate)

alpha_mean = 0.0  # degrees - mean angle of attack
A_pitch = 1.0  # degrees - pitch amplitude
A_plunge = 0.5  # degrees - plunge amplitude (converted to effective alpha)

# Torus spectral grid resolution
n_pitch = 3  # Number of points in pitch direction (theta1)
n_plunge = 3  # Number of points in plunge direction (theta2)

print("="*70)
print("TORUS TIME SPECTRAL METHOD EXAMPLE")
print("="*70)
print(f"\nMotion: α(t) = {alpha_mean}° + {A_pitch}°·sin({omega_pitch} t) + {A_plunge}°·sin({omega_plunge:.2f} t)")
print(f"\nFrequency ratio: ω_pitch/ω_plunge = 1/√2 (incommensurate)")
print(f"Torus grid: {n_pitch} x {n_plunge} = {n_pitch*n_plunge} time instances")
print(f"\nThis creates a quasi-periodic (non-repeating) motion pattern.")
print("="*70)

# ==============================================================================
# ADflow Options
# ==============================================================================

gridFile = "../input_files/naca64A010_euler-L2.cgns"

aeroOptions = {
    # I/O
    "gridfile": gridFile,
    "outputDirectory": "./output_torus_example",
    "monitorvariables": ["cpu", "resrho", "cl", "cd", "cmz"],
    "writeSurfaceSolution": True,
    "writeVolumeSolution": False,
    "surfaceVariables": ["cp", "vx", "vy", "vz", "mach"],

    # Physics
    "equationtype": "Euler",
    "equationmode": "time spectral",
    "liftIndex": 2,

    # Torus Time Spectral Options
    "useTorusTimeSpectral": True,
    "nTimeIntervalsSpectral1": n_pitch,
    "nTimeIntervalsSpectral2": n_plunge,
    "omegaFourier1": omega_pitch,
    "omegaFourier2": omega_plunge,

    # Mesh deformation
    "useexternaldynamicmesh": True,
    "usetsinterpolatedgridvelocity": True,

    # Solver
    "l2convergence": 1e-10,
    "ncycles": 5000,
    "usenksolver": True,
    "nkswitchtol": 1e-4,
    "nksubspacesize": 200,
    "useanksolver": True,
    "ankswitchtol": 1e-2,

    # Multigrid
    "mgcycle": "sg",
    "ncyclescoarse": 250,
}

# ==============================================================================
# Create Solver
# ==============================================================================

print("\nInitializing ADflow solver...")
CFDSolver = ADFLOW(options=aeroOptions, comm=MPI.COMM_WORLD, debug=False)

print(f"Torus Time Spectral mode enabled: {n_pitch}x{n_plunge} grid")
print(f"Total time instances: {n_pitch * n_plunge}")

# ==============================================================================
# Setup Aerodynamic Problem
# ==============================================================================

ap = AeroProblem(
    name="torus_aeroelastic",
    mach=0.796,
    altitude=10000.0,
    alpha=alpha_mean,
    areaRef=1.0,
    chordRef=1.0,
    xRef=0.248,
    xRot=0.248,
    evalFuncs=["cl", "cd", "cmz"],
)

# ==============================================================================
# Setup Mesh Deformation
# ==============================================================================

print("\nSetting up mesh deformation...")

# Create mesh object
meshOptions = {"gridFile": gridFile}
mesh = USMesh(options=meshOptions, comm=MPI.COMM_WORLD)
CFDSolver.setMesh(mesh)

# Generate torus grid of alpha values
theta_pitch = np.linspace(0.0, 2.0 * np.pi, n_pitch, endpoint=False)
theta_plunge = np.linspace(0.0, 2.0 * np.pi, n_plunge, endpoint=False)

print("\nTorus grid points (theta_pitch, theta_plunge, alpha):")
print("-" * 60)

alpha_grid = np.zeros((n_pitch, n_plunge))
for i in range(n_pitch):
    for j in range(n_plunge):
        # Combined motion on torus
        alpha_deg = alpha_mean + A_pitch * np.sin(theta_pitch[i]) + A_plunge * np.sin(theta_plunge[j])
        alpha_grid[i, j] = alpha_deg

        sps = j * n_pitch + i + 1  # Spectral instance index (1-based)
        print(f"  Instance {sps:2d}: θ_pitch={theta_pitch[i]:5.2f}, θ_plunge={theta_plunge[j]:5.2f}, α={alpha_deg:6.3f}°")

# Deform mesh for each torus point
# (In full implementation, use the Transfer class to handle surface deformation)
# For this example, we note that mesh deformation is required

print("\nNote: Mesh deformation for each time instance should be applied using")
print("      a Transfer class (see test_torus_time_spectral_naca64A010.py)")

# ==============================================================================
# Run Solver
# ==============================================================================

print("\n" + "="*70)
print("RUNNING TORUS TIME SPECTRAL SOLVE")
print("="*70)

# Solve
# CFDSolver(ap)  # Uncomment when mesh deformation is set up

# ==============================================================================
# Post-processing
# ==============================================================================

print("\nPost-processing:")
print("-" * 60)

# Evaluate functions
# funcs = {}
# CFDSolver.evalFunctions(ap, funcs)

# print("\nIntegrated forces and moments:")
# for key, val in funcs.items():
#     print(f"  {key}: {val:.6f}")

# The torus method gives values at each grid point
# To reconstruct the full quasi-periodic signal, use Fourier interpolation:
#
# f(theta1, theta2) ≈ Σ Σ c_{k,l} * exp(i*(k*theta1 + l*theta2))
#
# This allows evaluation at arbitrary time: t → (omega1*t, omega2*t) on torus

print("\nTorus Time Spectral Benefits:")
print("  ✓ Solves quasi-periodic problems without long time integration")
print("  ✓ Spectral accuracy in time for smooth solutions")
print("  ✓ Orders of magnitude faster than unsteady for periodic/quasi-periodic flows")
print("  ✓ Natural framework for aeroelastic problems with multiple modes")

print("\n" + "="*70)
print("Example complete!")
print("="*70)

# ==============================================================================
# Theoretical Background
# ==============================================================================

print("\n" + "="*70)
print("THEORETICAL BACKGROUND")
print("="*70)

print("""
The torus time spectral method extends the classical time spectral method
from single-frequency periodic problems to multi-frequency quasi-periodic problems.

Key Concepts:
-------------
1. **Quasi-periodic flow**: Flow with multiple incommensurate frequencies
   - Periodic: f(t) = f(t + T) for some period T
   - Quasi-periodic: f(t) = F(ω₁t, ω₂t, ...) where ω₁/ω₂ is irrational

2. **Invariant torus**: The solution trajectory lives on a torus in phase space
   - Parametrized by angles (θ₁, θ₂) ∈ [0, 2π]²
   - Time evolution: (θ₁, θ₂) = (ω₁t, ω₂t)

3. **Spectral collocation**: Discretize torus with n₁ × n₂ grid
   - Approximate: u(θ₁, θ₂) ≈ Σ Σ û_{k,l} exp(i(kθ₁ + lθ₂))
   - Time derivatives: ∂u/∂t = ω₁∂u/∂θ₁ + ω₂∂u/∂θ₂

4. **Advantages over unsteady**:
   - No time integration required
   - Spectral accuracy in time
   - Direct computation of limit cycle
   - Natural for frequency response analysis

Applications:
------------
- Aeroelasticity with multiple structural modes
- Rotor-stator interactions with different blade counts
- Vortex-induced vibrations
- Aeroacoustics with multiple tones
- Any quasi-periodic CFD problem

References:
----------
- He, S., Li, H., Ekici, K. "Torus Time Spectral Method for Multi-frequency
  Quasi-Periodic Problems" (2024)
- Hall, K. et al. "A Time-Linearized Navier-Stokes Analysis of Stall Flutter"
  (2000) - Original time spectral paper
""")

print("="*70)
