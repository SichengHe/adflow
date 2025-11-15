#!/usr/bin/env python
"""
Save Resolvent Modes to NumPy Format

Simpler alternative to CGNS export - saves modes as .npz files
that can be easily loaded and visualized.

This is recommended for:
- Quick analysis and visualization
- Custom post-processing
- Integration with other Python tools
"""

import numpy as np
import sys
import os

# Add parent directory for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from adflow import ADFLOW, ResolventAnalysisMatrixFree
from baseclasses import AeroProblem


def save_modes_numpy(resolvent, output_file='resolvent_modes.npz', nModes=3):
    """
    Save resolvent modes to NumPy .npz file.

    Parameters
    ----------
    resolvent : ResolventAnalysisMatrixFree or ResolventAnalysis
        Resolvent analysis object with computed modes
    output_file : str
        Output filename (.npz)
    nModes : int
        Number of modes to save

    Returns
    -------
    str
        Path to saved file
    """

    print(f"\nSaving {nModes} modes to NumPy format...")
    print(f"Output file: {output_file}")
    print("-" * 80)

    # Extract all modes
    singular_values = np.array([resolvent.getSingularValue(i) for i in range(nModes)])
    forcing_modes = np.column_stack([resolvent.getForcingMode(i) for i in range(nModes)])
    response_modes = np.column_stack([resolvent.getResponseMode(i) for i in range(nModes)])

    # Save to compressed NumPy archive
    np.savez_compressed(
        output_file,
        omega=resolvent.omega,
        nModes=nModes,
        singular_values=singular_values,
        forcing_modes_real=np.real(forcing_modes),
        forcing_modes_imag=np.imag(forcing_modes),
        response_modes_real=np.real(response_modes),
        response_modes_imag=np.imag(response_modes),
        state_size=forcing_modes.shape[0]
    )

    # Print summary
    print(f"\n✓ Saved {nModes} modes")
    print(f"  File: {output_file}")
    print(f"  Size: {os.path.getsize(output_file) / 1024:.1f} KB")
    print()

    print("Contents:")
    print(f"  - omega: {resolvent.omega}")
    print(f"  - nModes: {nModes}")
    print(f"  - singular_values: shape {singular_values.shape}")
    print(f"  - forcing_modes: shape {forcing_modes.shape} (complex)")
    print(f"  - response_modes: shape {response_modes.shape} (complex)")
    print()

    print("To load:")
    print(f"  data = np.load('{output_file}')")
    print("  v1 = data['forcing_modes_real'][:, 0] + 1j*data['forcing_modes_imag'][:, 0]")
    print("  u1 = data['response_modes_real'][:, 0] + 1j*data['response_modes_imag'][:, 0]")
    print("  sigma1 = data['singular_values'][0]")
    print()

    return output_file


def main():
    """
    Complete example: Compute modes and save to NumPy format.
    """

    # =========================================================================
    # 1. Setup CFD Problem
    # =========================================================================

    # Try multiple possible locations for mesh file
    possible_paths = [
        "../tests/input_files/naca64A010_euler-L2.cgns",  # From examples/
        "./tests/input_files/naca64A010_euler-L2.cgns",   # From repo root
        "./input_files/naca64A010_euler-L2.cgns",         # From tests/
    ]

    meshFile = None
    for path in possible_paths:
        if os.path.exists(path):
            meshFile = path
            break

    if meshFile is None:
        print(f"ERROR: Mesh file not found")
        return

    aeroOptions = {
        'gridFile': meshFile,
        'outputDirectory': './output_modes_numpy',

        # Physics
        'equationType': 'Euler',

        # Solver settings
        'CFL': 2.0,
        'L2Convergence': 1e-6,
        'nCycles': 200,
        'useNKSolver': True,

        # Output settings
        'printIterations': False,
        'printTiming': False,
        'printWarnings': False,
        'writeVolumeSolution': False,
        'writeSurfaceSolution': False,

        # CRITICAL for resolvent analysis
        'useMatrixFreedrdw': False,
    }

    print("="*80)
    print("RESOLVENT MODES → NUMPY FILE EXPORT")
    print("="*80)
    print()

    print("Step 1: Solving CFD...")
    CFDsolver = ADFLOW(options=aeroOptions, debug=False)

    ap = AeroProblem(
        name='naca64a010',
        mach=0.5,
        alpha=0.0,
        altitude=10000,
        areaRef=1.0,
        chordRef=1.0,
    )

    CFDsolver(ap)
    print("✓ Flow converged\n")

    # =========================================================================
    # 2. Compute Resolvent Modes
    # =========================================================================

    omega = 1.0
    nModes = 5

    print(f"Step 2: Computing {nModes} resolvent modes at ω = {omega}...")
    resolvent = ResolventAnalysisMatrixFree(CFDsolver, ap, omega=omega)

    # Enable ILU for fast convergence
    resolvent.enablePreconditioner(precond_type='ilu', drop_tol=1e-3, fill_factor=10)

    # Solve
    sigma_max = resolvent.solve(nModes=nModes, method='scipy')
    print(f"✓ Modes computed (σ₁ = {sigma_max:.6f})\n")

    # =========================================================================
    # 3. Save to NumPy
    # =========================================================================

    print("Step 3: Saving modes to NumPy format...")
    os.makedirs('./output_modes_numpy', exist_ok=True)
    output_file = f'./output_modes_numpy/resolvent_modes_omega{omega}.npz'

    save_modes_numpy(resolvent, output_file, nModes=nModes)

    # =========================================================================
    # 4. Demonstrate Loading
    # =========================================================================

    print("="*80)
    print("VERIFICATION: Loading saved modes")
    print("="*80)
    print()

    data = np.load(output_file)

    print("Loaded data:")
    for key in data.files:
        print(f"  {key}: {data[key].shape if hasattr(data[key], 'shape') else 'scalar'}")
    print()

    # Reconstruct first mode
    v1 = data['forcing_modes_real'][:, 0] + 1j*data['forcing_modes_imag'][:, 0]
    u1 = data['response_modes_real'][:, 0] + 1j*data['response_modes_imag'][:, 0]
    sigma1 = data['singular_values'][0]

    print("Reconstructed Mode 1:")
    print(f"  σ₁ = {sigma1:.6f}")
    print(f"  ||v₁|| = {np.linalg.norm(v1):.6f}")
    print(f"  ||u₁|| = {np.linalg.norm(u1):.6f}")
    print()

    # =========================================================================
    # Summary
    # =========================================================================

    print("="*80)
    print("COMPLETE!")
    print("="*80)
    print()
    print(f"Modes saved to: {output_file}")
    print()
    print("Next steps:")
    print("  1. Load modes: data = np.load('" + output_file + "')")
    print("  2. Analyze mode shapes")
    print("  3. Visualize with matplotlib/custom tools")
    print("  4. Use for control design")
    print()


if __name__ == "__main__":
    main()
