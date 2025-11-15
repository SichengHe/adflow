"""
Compare 3×3 degenerate torus vs 5-point 1D time spectral

When we project the 3×3 torus onto the diagonal, we get 5 unique time samples.
So the equivalent 1D TS should have 5 instances at those specific times.
"""

import numpy as np
from mpi4py import MPI
from adflow import ADFLOW
from idwarp import USMesh
from baseclasses import AeroProblem


def run_torus_3x3():
    """Run 3×3 degenerate torus"""

    gridFile = '../input_files/naca64A010_euler-L2.cgns'

    n1 = n2 = 3
    omega = 100.0
    A = 1.0

    options = {
        'gridfile': gridFile,
        'equationType': 'Euler',
        'equationMode': 'time spectral',
        'useTorusTimeSpectral': True,
        'nTimeIntervalsSpectral1': n1,
        'nTimeIntervalsSpectral2': n2,
        'omegaFourier1': omega,
        'omegaFourier2': omega,
        'timeIntervals': n1 * n2,
        'useexternaldynamicmesh': True,
        'l2convergence': 1e-12,
        'nCycles': 10000,
        'writeVolumeSolution': False,
        'writeSurfaceSolution': False,
    }

    print("="*70)
    print("RUNNING 3×3 DEGENERATE TORUS")
    print("="*70)

    CFDSolver = ADFLOW(options=options, debug=False)
    mesh = USMesh(options={"gridFile": gridFile})
    CFDSolver.setMesh(mesh)

    # Torus grid
    theta1 = np.linspace(0, 2*np.pi, n1, endpoint=False)
    theta2 = np.linspace(0, 2*np.pi, n2, endpoint=False)

    # Deform meshes for each torus point
    MDGroup = CFDSolver.allWallsGroup
    cfdPts0 = CFDSolver.getSurfaceCoordinates(MDGroup, includeZipper=False)
    xRot = 0.25

    cfdPts_list = []
    alphas_torus = []

    for j in range(n2):
        for i in range(n1):
            alpha_rad = (A * np.sin(theta1[i]) + A * np.sin(theta2[j])) * np.pi / 180
            alphas_torus.append(alpha_rad * 180 / np.pi)

            cc = np.cos(alpha_rad)
            ss = np.sin(alpha_rad)

            cfdPoints_deformed = np.zeros_like(cfdPts0)
            for k in range(len(cfdPts0)):
                x0, y0, z0 = cfdPts0[k]
                cfdPoints_deformed[k, 0] = cc * (x0 - xRot) + ss * y0 + xRot
                cfdPoints_deformed[k, 1] = -ss * (x0 - xRot) + cc * y0
                cfdPoints_deformed[k, 2] = z0

            cfdPts_list.append(cfdPoints_deformed)

    # Set warped meshes
    for sps in range(n1 * n2):
        CFDSolver.mesh.setSurfaceCoordinates(cfdPts_list[sps])
        CFDSolver.mesh.warpMesh()
        m = CFDSolver.mesh.getSolverGrid()
        CFDSolver.adflow.warping.setgridforoneinstance(m, sps=sps + 1)

    CFDSolver._updateGeomInfo = True
    CFDSolver.updateGeometryInfo()

    # Solve
    ap = AeroProblem(name='torus', mach=0.5, altitude=10000, areaRef=1.0, chordRef=1.0, evalFuncs=['cl', 'cd'])
    CFDSolver(ap)

    # Extract forces
    forces = CFDSolver.getInstanceForces()
    cl_torus = forces['cl']

    print("\nTorus results:")
    for idx in range(n1 * n2):
        i = idx % n1
        j = idx // n1
        print(f"  ({i},{j}): α={alphas_torus[idx]:7.3f}°, CL={cl_torus[idx]:.10f}")

    # Compute projected angles
    theta_proj = []
    for i in range(n1):
        for j in range(n2):
            theta_proj.append((theta1[i] + theta2[j]) / 2)

    return cl_torus, alphas_torus, np.array(theta_proj)


def run_1d_5pt():
    """Run 5-point 1D time spectral at projected angles"""

    gridFile = '../input_files/naca64A010_euler-L2.cgns'

    omega = 100.0
    A_1d = 2.0  # Amplitude for 1D to match torus diagonal

    # The 5 unique projected angles from 3×3 torus
    theta_samples = np.array([0.0, 1.047198, 2.094395, 3.141593, 4.188790])
    n_inst = len(theta_samples)

    # But 1D TS requires UNIFORM spacing - this won't work!
    # We need a different approach

    print("\n" + "="*70)
    print("PROBLEM: 1D TS requires uniform time spacing!")
    print("="*70)
    print()
    print(f"The 3×3 torus projects to {n_inst} points at:")
    print(f"  θ = {theta_samples}")
    print()
    print("These are NOT uniformly spaced, so standard 1D TS won't work.")
    print()
    print("OPTIONS:")
    print("  1. Use torus with n2=1 (true 1D)")
    print("  2. Use 6-point 1D TS and compare subset")
    print("  3. Accept that 3×3 torus is fundamentally different")

    return None


if __name__ == "__main__":
    cl_torus, alphas, theta_proj = run_torus_3x3()

    print("\n" + "="*70)
    print("PROJECTED ANGLES FROM 3×3 TORUS")
    print("="*70)

    theta_unique = np.unique(np.round(theta_proj % (2*np.pi), 6))
    print(f"\nUnique projected angles: {theta_unique}")
    print(f"Number of unique: {len(theta_unique)}")
    print()
    print("Standard 1D TS uses UNIFORM spacing: θ = 2πi/n")
    print(f"For n=5: {np.linspace(0, 2*np.pi, 5, endpoint=False)}")
    print()
    print("These don't match! The 3×3 torus is NOT equivalent to uniform 1D TS.")
