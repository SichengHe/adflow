#!/usr/bin/env python
"""
Simple test to verify torus time spectral solver can initialize and run.
This test does NOT require idwarp - just tests basic solver functionality.
"""
import numpy as np
from mpi4py import MPI
from adflow import ADFLOW

# Grid dimensions
n1 = 3
n2 = 3

# Two-frequency motion parameters
omega1 = 100.0  # rad/s
omega2 = 100.0 * np.sqrt(2.0)  # rad/s (incommensurate ratio)

print("="*70)
print("Testing Torus Time Spectral Solver Initialization")
print("="*70)
print(f"Grid: {n1} x {n2} = {n1*n2} time instances")
print(f"Frequency 1: {omega1} rad/s")
print(f"Frequency 2: {omega2} rad/s")
print(f"Ratio: {omega2/omega1:.6f} (sqrt(2) = {np.sqrt(2):.6f})")
print("="*70)

options = {
    'gridfile': '../input_files/naca64A010_euler-L2.cgns',
    'equationtype': 'Euler',
    'equationmode': 'time spectral',
    'useTorusTimeSpectral': True,
    'nTimeIntervalsSpectral1': n1,
    'nTimeIntervalsSpectral2': n2,
    'omegaFourier1': omega1,
    'omegaFourier2': omega2,
    'usenksolver': True,
    'l2convergence': 1e-6,
    'ncycles': 10,  # Just a few cycles to test
}

print("\nInitializing solver...")
CFDSolver = ADFLOW(options=options, comm=MPI.COMM_WORLD)

print("\nSolver initialized successfully!")
print(f"Time intervals set to: {CFDSolver.getOption('timeIntervals')}")

# Set simple aeroproblem
ap = CFDSolver.getAeroProblem()
ap.mach = 0.5
ap.alpha = 1.5
ap.T = 293.15
ap.P = 101325.0

print("\nRunning solver for a few iterations...")
try:
    CFDSolver(ap)
    print("\n" + "="*70)
    print("SUCCESS: Torus Time Spectral solver can solve!")
    print("="*70)
except Exception as e:
    print("\n" + "="*70)
    print(f"Solver error: {e}")
    print("="*70)
    raise