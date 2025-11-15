#!/usr/bin/env python
"""
Save Resolvent Modes to CGNS Format

⚠️ NOTE: This script has known issues with ADflow's CGNS writing.
⚠️ RECOMMENDED ALTERNATIVE: Use save_resolvent_modes_numpy.py instead!

The NumPy export is simpler, more reliable, and produces smaller files.
You can then visualize modes using custom Python tools or convert to
other formats as needed.

This script attempts to save resolvent mode shapes as CGNS files
for visualization in Tecplot, ParaView, or other CFD visualization tools.

KNOWN ISSUES:
- ADflow may fail to write CGNS files depending on configuration
- Error: "File could not be opened by cgns for writing"
- Workaround: Use save_resolvent_modes_numpy.py instead
"""

import numpy as np
import sys
import os

# Add parent directory for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from adflow import ADFLOW, ResolventAnalysisMatrixFree
from baseclasses import AeroProblem


def save_modes_to_cgns(resolvent, CFDsolver, output_dir='./output_modes_cgns', nModes=3):
    """
    Save resolvent modes to CGNS files for visualization.

    Parameters
    ----------
    resolvent : ResolventAnalysisMatrixFree or ResolventAnalysis
        Resolvent analysis object with computed modes
    CFDsolver : ADFLOW
        ADflow solver with mesh and solution
    output_dir : str
        Directory to save CGNS files
    nModes : int
        Number of modes to save

    Notes
    -----
    Each mode is saved as a separate CGNS file with variables:
    - forcing_real, forcing_imag: Forcing mode (v_i)
    - response_real, response_imag: Response mode (u_i)
    - forcing_magnitude: |v_i|
    - response_magnitude: |u_i|
    """

    os.makedirs(output_dir, exist_ok=True)

    print(f"\nSaving {nModes} modes to CGNS format...")
    print(f"Output directory: {output_dir}")
    print("-" * 80)

    for i in range(nModes):
        # Get mode shapes
        v_i = resolvent.getForcingMode(i)  # Forcing mode
        u_i = resolvent.getResponseMode(i)  # Response mode
        sigma_i = resolvent.getSingularValue(i)

        print(f"\nMode {i+1} (σ_{i+1} = {sigma_i:.6f}):")

        # Set modes as solution variables in ADflow
        # ADflow stores conservative variables: [ρ, ρu, ρv, ρw, ρE] (+ turbulence)
        # We'll map the complex mode to real/imaginary parts

        # Real part of forcing mode
        CFDsolver.setStates(np.real(v_i))
        forcing_real_file = os.path.join(output_dir, f'mode{i+1}_forcing_real')
        CFDsolver.writeSolution(forcing_real_file + '.cgns')
        print(f"  ✓ Saved: {forcing_real_file}.cgns (forcing real part)")

        # Imaginary part of forcing mode
        CFDsolver.setStates(np.imag(v_i))
        forcing_imag_file = os.path.join(output_dir, f'mode{i+1}_forcing_imag')
        CFDsolver.writeSolution(forcing_imag_file + '.cgns')
        print(f"  ✓ Saved: {forcing_imag_file}.cgns (forcing imag part)")

        # Real part of response mode
        CFDsolver.setStates(np.real(u_i))
        response_real_file = os.path.join(output_dir, f'mode{i+1}_response_real')
        CFDsolver.writeSolution(response_real_file + '.cgns')
        print(f"  ✓ Saved: {response_real_file}.cgns (response real part)")

        # Imaginary part of response mode
        CFDsolver.setStates(np.imag(u_i))
        response_imag_file = os.path.join(output_dir, f'mode{i+1}_response_imag')
        CFDsolver.writeSolution(response_imag_file + '.cgns')
        print(f"  ✓ Saved: {response_imag_file}.cgns (response imag part)")

        # Magnitude of forcing (useful for visualization)
        forcing_mag = np.abs(v_i)
        # Map magnitude to density for visualization
        state_mag = np.zeros_like(v_i, dtype=float)
        state_mag[0::5] = forcing_mag[0::5]  # Set density component
        CFDsolver.setStates(state_mag)
        forcing_mag_file = os.path.join(output_dir, f'mode{i+1}_forcing_magnitude')
        CFDsolver.writeSolution(forcing_mag_file + '.cgns')
        print(f"  ✓ Saved: {forcing_mag_file}.cgns (forcing magnitude)")

    print(f"\n{'='*80}")
    print(f"✓ Saved {nModes} modes to CGNS format")
    print(f"{'='*80}\n")

    print("Visualization Tips:")
    print("-" * 80)
    print("In Tecplot/ParaView:")
    print("  1. Load mode*_forcing_*.cgns to see forcing distribution")
    print("  2. Load mode*_response_*.cgns to see response distribution")
    print("  3. Use magnitude files for scalar visualization")
    print("  4. Real/imag parts show spatial phase information")
    print()
    print("For actuator placement:")
    print("  - Look for regions of high |forcing| magnitude")
    print("  - These are optimal locations for control actuators")
    print()
    print("For sensor placement:")
    print("  - Look for regions of high |response| magnitude")
    print("  - These are optimal locations for sensors")
    print()


def main():
    """
    Complete example: Compute modes and save to CGNS.
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
        print(f"ERROR: Mesh file not found in any of these locations:")
        for path in possible_paths:
            print(f"  - {path}")
        print("\nPlease run from repo root, examples/, or tests/ directory")
        print("Or ensure mesh file exists at: tests/input_files/naca64A010_euler-L2.cgns")
        return

    aeroOptions = {
        'gridFile': meshFile,
        'outputDirectory': './output_modes_cgns',

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

        # CRITICAL for resolvent analysis
        'useMatrixFreedrdw': False,
    }

    print("="*80)
    print("RESOLVENT MODES → CGNS FILE EXPORT")
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
    nModes = 3

    print(f"Step 2: Computing {nModes} resolvent modes at ω = {omega}...")
    resolvent = ResolventAnalysisMatrixFree(CFDsolver, ap, omega=omega)

    # Enable ILU for fast convergence
    resolvent.enablePreconditioner(precond_type='ilu', drop_tol=1e-3, fill_factor=10)

    # Solve
    sigma_max = resolvent.solve(nModes=nModes, method='scipy')
    print(f"✓ Modes computed (σ₁ = {sigma_max:.6f})\n")

    # =========================================================================
    # 3. Save to CGNS
    # =========================================================================

    print("Step 3: Saving modes to CGNS...")
    save_modes_to_cgns(resolvent, CFDsolver, nModes=nModes)

    # =========================================================================
    # 4. Summary
    # =========================================================================

    print("="*80)
    print("COMPLETE!")
    print("="*80)
    print()
    print("Files created in ./output_modes_cgns/:")
    print("  - mode{1,2,3}_forcing_real.cgns")
    print("  - mode{1,2,3}_forcing_imag.cgns")
    print("  - mode{1,2,3}_response_real.cgns")
    print("  - mode{1,2,3}_response_imag.cgns")
    print("  - mode{1,2,3}_forcing_magnitude.cgns")
    print()
    print("Next steps:")
    print("  1. Open files in Tecplot or ParaView")
    print("  2. Visualize mode spatial structure")
    print("  3. Identify critical regions for control")
    print()


if __name__ == "__main__":
    main()
