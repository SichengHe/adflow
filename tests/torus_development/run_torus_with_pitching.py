"""
Torus Time Spectral with Two-Frequency Pitching Motion

Motion: alpha(theta1, theta2) = alpha_mean + A1*sin(theta1) + A2*sin(theta2)

Based on test_torus_time_spectral_naca64A010.py
"""

import numpy as np
from mpi4py import MPI
from adflow import ADFLOW
from idwarp import USMesh
from baseclasses import AeroProblem
import pickle


# =============================================================================
# Torus Transfer Class for Mesh Deformation
# =============================================================================

class TorusTransfer:
    """Transfer class for two-frequency pitching motion on torus grid"""

    def __init__(self, alpha_grid, xRot, aeroSolver, n1, n2):
        """
        Parameters
        ----------
        alpha_grid : ndarray (n1, n2)
            Pitch angle in RADIANS at each torus grid point
        xRot : float
            X-coordinate of rotation axis (typically 0.25 for quarter-chord)
        aeroSolver : ADFLOW object
            The CFD solver
        n1, n2 : int
            Torus grid dimensions
        """
        self.alpha_grid = alpha_grid
        self.n1 = n1
        self.n2 = n2
        self.ntot = n1 * n2
        self.aeroSolver = aeroSolver
        self.xRot = xRot

    def getUndeformedSurfaceNodes(self):
        """Get initial undeformed surface nodes"""
        self.MDGroup = self.aeroSolver.allWallsGroup
        self.cfdPts0 = self.aeroSolver.getSurfaceCoordinates(
            self.MDGroup, includeZipper=False
        )

    def setDisplacements(self):
        """Compute deformed surface nodes for all torus instances"""
        xRot = self.xRot
        alpha_grid = self.alpha_grid
        cfdPoints_init = self.cfdPts0

        N_pts = cfdPoints_init.shape[0]
        self.cfdPts = []

        # Loop over torus grid (column-major: iterate i1 first, then i2)
        for i2 in range(self.n2):
            for i1 in range(self.n1):
                cfdPoints_deformed = np.zeros((N_pts, 3))

                # Pitch angle for this instance
                pitch_rad = alpha_grid[i1, i2]

                # Rotation matrix
                cc = np.cos(pitch_rad)
                ss = np.sin(pitch_rad)

                # Rotate each point about (xRot, 0)
                for j in range(N_pts):
                    x0, y0, z0 = cfdPoints_init[j, :]

                    # Rigid rotation:
                    # x' = cos(α)*(x0 - xRot) + sin(α)*y0 + xRot
                    # y' = -sin(α)*(x0 - xRot) + cos(α)*y0
                    # z' = z0 (2D motion)
                    cfdPoints_deformed[j, 0] = cc * (x0 - xRot) + ss * y0 + xRot
                    cfdPoints_deformed[j, 1] = -ss * (x0 - xRot) + cc * y0
                    cfdPoints_deformed[j, 2] = z0

                self.cfdPts.append(cfdPoints_deformed)

    def setVolumeMesh(self):
        """Warp volume mesh for all torus instances"""
        sps = 0
        for i2 in range(self.n2):
            for i1 in range(self.n1):
                # Set surface coordinates for this instance
                self.aeroSolver.mesh.setSurfaceCoordinates(self.cfdPts[sps])

                # Warp volume mesh
                self.aeroSolver.mesh.warpMesh()

                # Get warped volume grid
                m = self.aeroSolver.mesh.getSolverGrid()

                # Set grid for this instance (sps is 1-indexed in Fortran)
                self.aeroSolver.adflow.warping.setgridforoneinstance(m, sps=sps + 1)

                sps += 1

        # Update geometry info
        self.aeroSolver._updateGeomInfo = True
        self.aeroSolver.updateGeometryInfo()


# =============================================================================
# Main Script
# =============================================================================

def main():
    print("="*70)
    print("TORUS TIME SPECTRAL WITH TWO-FREQUENCY PITCHING")
    print("="*70)

    gridFile = '../input_files/naca64A010_euler-L2.cgns'

    # Torus grid dimensions
    n1 = 3
    n2 = 3
    ntot = n1 * n2

    # Two-frequency motion parameters
    omega1 = 100.0  # rad/s
    omega2 = 100.0 * np.sqrt(2.0)  # rad/s (incommensurate)
    alpha_mean = 0.0  # degrees
    A1 = 1.0  # degrees amplitude for frequency 1
    A2 = 0.5  # degrees amplitude for frequency 2

    f1 = omega1 / (2 * np.pi)
    f2 = omega2 / (2 * np.pi)

    print(f"\nMotion parameters:")
    print(f"  α(θ₁,θ₂) = {alpha_mean}° + {A1}°·sin(θ₁) + {A2}°·sin(θ₂)")
    print(f"  ω₁ = {omega1:.1f} rad/s ({f1:.2f} Hz)")
    print(f"  ω₂ = {omega2:.1f} rad/s ({f2:.2f} Hz)")
    print(f"  Ratio ω₂/ω₁ = √2 = {omega2/omega1:.6f} (quasi-periodic)")
    print(f"  Torus grid: {n1} × {n2} = {ntot} instances")

    # ADflow options
    options = {
        'gridfile': gridFile,
        'outputDirectory': './output_torus_pitching',

        # Physics
        'equationType': 'Euler',
        'equationMode': 'time spectral',

        # Torus time spectral options
        'useTorusTimeSpectral': True,
        'nTimeIntervalsSpectral1': n1,
        'nTimeIntervalsSpectral2': n2,
        'omegaFourier1': omega1,
        'omegaFourier2': omega2,
        'timeintervals': ntot,

        # Mesh deformation
        'useexternaldynamicmesh': True,
        'usetsinterpolatedgridvelocity': True,

        # Solver settings
        'mgcycle': 'sg',
        'nCycles': 10000,
        'l2convergence': 1e-10,
        'l2convergencecoarse': 1e-2,
        'monitorvariables': ['resrho', 'cl', 'cd'],

        # NK solver
        'usenksolver': True,
        'nkswitchtol': 1e-4,
        'NKSubSpaceSize': 400,
        'applypcsubspacesize': 400,

        # ANK solver - Start with DADI smoother first to get better initial guess
        'useanksolver': True,
        'ankswitchtol': 1e-2,  # Switch to ANK after some DADI iterations
        'anksubspacesize': 200,
        'ANKCFL0': 1.0,
        'ANKCFLLimit': 10.0,

        # Use basic DADI smoother first
        'CFL': 1.5,
        'nSubiter': 3,  # More subiterations for DADI

        # Other
        'alphafollowing': False,
        'blocksplitting': True,
        'useblockettes': False,

        # Output
        'writeVolumeSolution': True,
        'writeSurfaceSolution': True,
        'numberSolutions': True,
    }

    # Mesh options
    meshOptions = {
        "gridFile": gridFile,
    }

    print("\n" + "="*70)
    print("STEP 1: Initialize ADflow")
    print("="*70)

    CFDSolver = ADFLOW(options=options, debug=False)

    print("\n" + "="*70)
    print("STEP 2: Setup mesh deformation with IDWarp")
    print("="*70)

    mesh = USMesh(options=meshOptions)
    CFDSolver.setMesh(mesh)
    print("✅ IDWarp mesh initialized")

    print("\n" + "="*70)
    print("STEP 3: Generate torus grid of alpha values")
    print("="*70)

    # Torus angle grid
    theta1 = np.linspace(0.0, 2.0 * np.pi, n1, endpoint=False)
    theta2 = np.linspace(0.0, 2.0 * np.pi, n2, endpoint=False)

    # Compute alpha at each torus point
    alpha_grid = np.zeros((n1, n2))

    print("\nAlpha distribution (degrees):")
    print("       ", end="")
    for i2 in range(n2):
        print(f"  θ₂={theta2[i2]:.2f}", end="")
    print()
    print("-" * 50)

    for i1 in range(n1):
        print(f"θ₁={theta1[i1]:.2f}:", end="")
        for i2 in range(n2):
            # Motion: α(θ₁,θ₂) = α_mean + A1·sin(θ₁) + A2·sin(θ₂)
            alpha_deg = alpha_mean + A1 * np.sin(theta1[i1]) + A2 * np.sin(theta2[i2])
            alpha_rad = alpha_deg * np.pi / 180.0
            alpha_grid[i1, i2] = alpha_rad
            print(f"  {alpha_deg:7.2f}°", end="")
        print()

    print(f"\nAlpha range: [{np.min(alpha_grid)*180/np.pi:.2f}°, {np.max(alpha_grid)*180/np.pi:.2f}°]")

    print("\n" + "="*70)
    print("STEP 4: Deform mesh for all torus instances")
    print("="*70)

    xRot = 0.25  # Rotation axis at quarter chord
    TSTransfer = TorusTransfer(alpha_grid, xRot, CFDSolver, n1, n2)

    print("Getting undeformed surface nodes...")
    TSTransfer.getUndeformedSurfaceNodes()
    print(f"✅ Got {TSTransfer.cfdPts0.shape[0]} surface points")

    print("Computing displacements for all instances...")
    TSTransfer.setDisplacements()
    print(f"✅ Computed {len(TSTransfer.cfdPts)} deformed meshes")

    print("Warping volume meshes...")
    TSTransfer.setVolumeMesh()
    print(f"✅ All {ntot} volume meshes warped and set")

    print("\n" + "="*70)
    print("STEP 5: Solve torus time spectral system")
    print("="*70)

    # Create aeroproblem
    ap = AeroProblem(
        name='torus_pitch',
        mach=0.5,
        altitude=10000,
        areaRef=1.0,
        chordRef=1.0,
        xRef=0.25,
        evalFuncs=['cl', 'cd', 'cmz']
    )

    print("\nRunning torus solver...")
    print("(This couples all 9 instances simultaneously)")

    CFDSolver(ap)

    print("\n" + "="*70)
    print("STEP 6: Extract forces from all instances")
    print("="*70)

    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)

    CL_instances = []
    CD_instances = []
    CM_instances = []

    print("\nForce coefficients at each torus instance:")
    print("Instance  (i1,i2)     θ₁      θ₂      α(deg)      CL         CD         CM")
    print("-" * 85)

    for idx in range(ntot):
        i1 = idx % n1
        i2 = idx // n1

        th1 = theta1[i1]
        th2 = theta2[i2]
        alpha_deg = alpha_grid[i1, i2] * 180 / np.pi

        # Try to extract instance-specific forces
        cl_key = f'torus_pitch_cl_{idx}'
        cd_key = f'torus_pitch_cd_{idx}'
        cm_key = f'torus_pitch_cmz_{idx}'

        cl = funcs.get(cl_key, funcs.get('torus_pitch_cl', 0.0))
        cd = funcs.get(cd_key, funcs.get('torus_pitch_cd', 0.0))
        cm = funcs.get(cm_key, funcs.get('torus_pitch_cmz', 0.0))

        CL_instances.append(cl)
        CD_instances.append(cd)
        CM_instances.append(cm)

        print(f"   {idx:2d}    ({i1},{i2})   {th1:6.3f}  {th2:6.3f}  {alpha_deg:7.2f}   {cl:9.6f}  {cd:9.6f}  {cm:9.6f}")

    CL_instances = np.array(CL_instances)
    CD_instances = np.array(CD_instances)
    CM_instances = np.array(CM_instances)

    print("\n" + "="*70)
    print("RESULTS SUMMARY")
    print("="*70)

    print(f"\nCL statistics:")
    print(f"  Mean: {np.mean(CL_instances):.6f}")
    print(f"  Std:  {np.std(CL_instances):.6f}")
    print(f"  Min:  {np.min(CL_instances):.6f}")
    print(f"  Max:  {np.max(CL_instances):.6f}")

    print(f"\nCD statistics:")
    print(f"  Mean: {np.mean(CD_instances):.6f}")
    print(f"  Std:  {np.std(CD_instances):.6f}")
    print(f"  Min:  {np.min(CD_instances):.6f}")
    print(f"  Max:  {np.max(CD_instances):.6f}")

    # Save results
    torus_data = {
        'n1': n1,
        'n2': n2,
        'omega1': omega1,
        'omega2': omega2,
        'f1': f1,
        'f2': f2,
        'alpha_mean': alpha_mean,
        'A1': A1,
        'A2': A2,
        'alpha_grid': alpha_grid,
        'theta1_grid': theta1,
        'theta2_grid': theta2,
        'CL_instances': CL_instances,
        'CD_instances': CD_instances,
        'CM_instances': CM_instances,
    }

    with open('torus_pitching_results.pkl', 'wb') as f:
        pickle.dump(torus_data, f)

    print("\n✅ Results saved to torus_pitching_results.pkl")
    print("\n" + "="*70)
    print("TORUS TIME SPECTRAL COMPLETED SUCCESSFULLY!")
    print("="*70)


if __name__ == "__main__":
    main()
