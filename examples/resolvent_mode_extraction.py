#!/usr/bin/env python
"""
Example: Extracting Resolvent Mode Shapes

This script demonstrates how to extract forcing and response mode shapes
from resolvent analysis for flow control applications.

The resolvent operator R(ω) = (iω·I - J)^{-1} provides the input-output
relationship between forcing f and response u:

    u = R(ω) * f

The SVD of R(ω) gives:
    R(ω) = U * Σ * V^H

Where:
- U columns = Response modes (flow response to forcing)
- V columns = Forcing modes (optimal forcing shapes)
- Σ diagonal = Amplification factors (singular values)

For flow control:
- Forcing mode V_i: Where to apply control
- Response mode U_i: Expected flow response
- Singular value σ_i: Amplification (how effective control is)
"""

import numpy as np
import sys
import os

# Add parent directory for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from adflow import ADFLOW, ResolventAnalysisMatrixFree
from baseclasses import AeroProblem

def extract_and_analyze_modes():
    """
    Example showing how to extract and analyze resolvent modes.
    """

    # =========================================================================
    # 1. Setup and solve CFD (same as before)
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
        return None

    aeroOptions = {
        'gridFile': meshFile,
        'outputDirectory': './output_mode_extraction',

        # Physics
        'equationType': 'Euler',

        # Solver settings
        'CFL': 2.0,
        'CFLCoarse': 1.5,
        'MGCycle': '2w',
        'nCycles': 200,
        'L2Convergence': 1e-6,
        'useNKSolver': True,
        'NKSwitchTol': 1e-2,

        # Reduce output
        'printIterations': False,
        'printTiming': False,
        'printWarnings': False,
        'writeVolumeSolution': False,
        'writeSurfaceSolution': False,

        # CRITICAL: Must use assembled Jacobian
        'useMatrixFreedrdw': False,
    }

    print("Creating CFD solver...")
    CFDsolver = ADFLOW(options=aeroOptions, debug=False)

    ap = AeroProblem(
        name='naca64a010',
        mach=0.5,
        alpha=0.0,
        altitude=10000,
        areaRef=1.0,
        chordRef=1.0,
    )

    print("Solving flow field...")
    CFDsolver(ap)
    print("✓ Flow converged\n")

    # =========================================================================
    # 2. Perform Resolvent Analysis
    # =========================================================================

    omega = 1.0  # Frequency of interest
    nModes = 5   # Number of modes to extract

    print(f"Setting up resolvent analysis at ω = {omega}")
    resolvent = ResolventAnalysisMatrixFree(CFDsolver, ap, omega=omega)

    # Enable ILU preconditioning for fast convergence
    print("Enabling ILU preconditioning...")
    resolvent.enablePreconditioner(precond_type='ilu', drop_tol=1e-3, fill_factor=10)

    # Solve for modes
    print(f"\nComputing {nModes} leading resolvent modes...")
    sigma_max = resolvent.solve(nModes=nModes, method='scipy')

    print(f"\n{'='*80}")
    print("RESOLVENT ANALYSIS COMPLETE")
    print(f"{'='*80}\n")

    # =========================================================================
    # 3. Extract Mode Shapes
    # =========================================================================

    print("Extracting mode shapes...\n")

    # Get singular values (amplification factors)
    print("Singular Values (Amplification Factors):")
    print("-" * 50)
    for i in range(nModes):
        sigma_i = resolvent.getSingularValue(i)
        print(f"  Mode {i+1}: σ_{i+1} = {sigma_i:.6f}")
    print()

    # Extract individual modes
    print("Mode Shape Details:")
    print("-" * 50)
    for i in range(min(3, nModes)):  # Show first 3 modes
        # Get forcing mode (optimal forcing shape)
        v_i = resolvent.getForcingMode(i)

        # Get response mode (flow response)
        u_i = resolvent.getResponseMode(i)

        # Get amplification
        sigma_i = resolvent.getSingularValue(i)

        print(f"\nMode {i+1}:")
        print(f"  Singular value: σ_{i+1} = {sigma_i:.6f}")
        print(f"  Forcing mode:   ||v_{i+1}|| = {np.linalg.norm(v_i):.6f}")
        print(f"  Response mode:  ||u_{i+1}|| = {np.linalg.norm(u_i):.6f}")
        print(f"  Amplification:  ||u_{i+1}|| / ||v_{i+1}|| ≈ {sigma_i:.6f}")

        # Mode properties
        print(f"  Forcing is complex:  {np.iscomplexobj(v_i)}")
        print(f"  Response is complex: {np.iscomplexobj(u_i)}")
        print(f"  Shape: {v_i.shape}")

    # =========================================================================
    # 4. Physical Interpretation
    # =========================================================================

    print(f"\n{'='*80}")
    print("PHYSICAL INTERPRETATION")
    print(f"{'='*80}\n")

    # Dominant mode (most amplified)
    v1 = resolvent.getForcingMode(0)
    u1 = resolvent.getResponseMode(0)
    sigma1 = resolvent.getSingularValue(0)

    print("Dominant Mode (Most Effective for Control):")
    print("-" * 50)
    print(f"  Amplification factor: {sigma1:.2f}x")
    print(f"  Meaning: A forcing with shape v₁ produces a response {sigma1:.2f}x larger")
    print()
    print("For Flow Control:")
    print(f"  - Apply control forcing proportional to v₁")
    print(f"  - Expect flow response proportional to u₁")
    print(f"  - Energy amplification: {sigma1**2:.1f}x")
    print()

    # Energy distribution across modes
    print("Energy Distribution Across Modes:")
    print("-" * 50)
    total_energy = sum(resolvent.getSingularValue(i)**2 for i in range(nModes))
    for i in range(nModes):
        sigma_i = resolvent.getSingularValue(i)
        energy_fraction = (sigma_i**2 / total_energy) * 100
        print(f"  Mode {i+1}: {energy_fraction:5.1f}% of energy")
    print()

    # =========================================================================
    # 5. Example: Using Modes for Control Design
    # =========================================================================

    print(f"{'='*80}")
    print("EXAMPLE: CONTROL DESIGN")
    print(f"{'='*80}\n")

    print("To design a flow control system:")
    print("-" * 50)
    print("1. Choose target modes (e.g., dominant mode for maximum effect)")
    print("2. Extract forcing mode shape: v₁")
    print("3. Extract response mode shape: u₁")
    print("4. Check amplification: σ₁")
    print()
    print("Control Law:")
    print("  - Actuator placement: Align with v₁ (forcing mode)")
    print("  - Expected response: σ₁ * v₁ → u₁")
    print("  - Sensor placement: Align with u₁ (response mode)")
    print()

    print(f"For this case (ω = {omega}):")
    print(f"  - Amplification: {sigma1:.2f}x")
    print(f"  - Dominant forcing requires {v1.shape[0]} DOF")
    print(f"  - Can achieve {sigma1:.2f}x amplification with optimal forcing")
    print()

    # =========================================================================
    # 6. Saving Modes for Visualization
    # =========================================================================

    print(f"{'='*80}")
    print("SAVING MODES FOR VISUALIZATION")
    print(f"{'='*80}\n")

    output_dir = './output_mode_extraction'
    os.makedirs(output_dir, exist_ok=True)

    # Save all modes to file
    modes_file = os.path.join(output_dir, f'resolvent_modes_omega{omega}.npz')

    np.savez(
        modes_file,
        omega=omega,
        singular_values=np.array([resolvent.getSingularValue(i) for i in range(nModes)]),
        forcing_modes=np.array([resolvent.getForcingMode(i) for i in range(nModes)]).T,
        response_modes=np.array([resolvent.getResponseMode(i) for i in range(nModes)]).T,
        state_size=v1.shape[0]
    )

    print(f"✓ Modes saved to: {modes_file}")
    print()
    print("To load modes later:")
    print(f"  data = np.load('{modes_file}')")
    print("  v1 = data['forcing_modes'][:, 0]  # First forcing mode")
    print("  u1 = data['response_modes'][:, 0]  # First response mode")
    print("  sigma1 = data['singular_values'][0]  # First singular value")
    print()

    # =========================================================================
    # Summary
    # =========================================================================

    print(f"{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}\n")
    print(f"✓ Extracted {nModes} resolvent modes at ω = {omega}")
    print(f"✓ Dominant amplification: {sigma1:.2f}x")
    print(f"✓ Mode shapes saved for visualization")
    print(f"✓ Ready for flow control design")
    print()

    return resolvent


if __name__ == "__main__":
    print("="*80)
    print("RESOLVENT MODE EXTRACTION EXAMPLE")
    print("="*80)
    print()
    print("This example shows how to:")
    print("  1. Perform resolvent analysis")
    print("  2. Extract forcing and response mode shapes")
    print("  3. Interpret modes for flow control")
    print("  4. Save modes for visualization")
    print()

    resolvent = extract_and_analyze_modes()

    print("="*80)
    print("EXAMPLE COMPLETE!")
    print("="*80)
