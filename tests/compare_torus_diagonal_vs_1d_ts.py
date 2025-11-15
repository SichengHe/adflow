"""
Compare DIAGONAL instances of 3×3 degenerate torus vs 1D Time Spectral

For a degenerate torus with ω₁=ω₂ and A₁=A₂, only the DIAGONAL instances
(i₁=i₂) correspond to true physical time instants.

Diagonal instances:
- (0,0): θ₁=0, θ₂=0 → θ=0, α = 2A·sin(0) = 0
- (1,1): θ₁=2π/3, θ₂=2π/3 → θ=2π/3, α = 2A·sin(2π/3) = √3·A
- (2,2): θ₁=4π/3, θ₂=4π/3 → θ=4π/3, α = 2A·sin(4π/3) = -√3·A

These should match 3-point 1D TS with amplitude 2A at θ = 0, 2π/3, 4π/3.
"""

import numpy as np
from mpi4py import MPI
from adflow import ADFLOW
from idwarp import USMesh
from baseclasses import AeroProblem


def run_degenerate_torus():
    """Run 3×3 degenerate torus and extract DIAGONAL instances"""

    print("="*70)
    print("RUNNING 3×3 DEGENERATE TORUS (extracting diagonal only)")
    print("="*70)

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
        'l2convergence': 1e-10,
        'nCycles': 10000,
        'writeVolumeSolution': False,
        'writeSurfaceSolution': False,
    }

    CFDSolver = ADFLOW(options=options, debug=False)
    mesh = USMesh(options={"gridFile": gridFile})
    CFDSolver.setMesh(mesh)

    # Generate alpha grid
    theta1 = np.linspace(0.0, 2.0 * np.pi, n1, endpoint=False)
    theta2 = np.linspace(0.0, 2.0 * np.pi, n2, endpoint=False)

    print(f"\nDiagonal instances (θ₁ = θ₂):")
    diagonal_indices = []
    diagonal_theta = []
    diagonal_alpha = []

    for i in range(n1):
        theta = theta1[i]  # = theta2[i] on diagonal
        alpha_rad = (A * np.sin(theta1[i]) + A * np.sin(theta2[i]))
        alpha_deg = alpha_rad
        idx = i * n1 + i  # Diagonal index in flattened array

        diagonal_indices.append(idx)
        diagonal_theta.append(theta)
        diagonal_alpha.append(alpha_deg * np.pi / 180.0)

        print(f"  Instance {idx+1} ({i},{i}): θ={theta:.4f} ({theta*180/np.pi:6.1f}°), α={alpha_deg:7.3f}°")

    # Deform meshes for ALL instances
    MDGroup = CFDSolver.allWallsGroup
    cfdPts0 = CFDSolver.getSurfaceCoordinates(MDGroup, includeZipper=False)
    xRot = 0.25

    cfdPts_list = []

    for i2 in range(n2):
        for i1 in range(n1):
            alpha_rad = (A * np.sin(theta1[i1]) + A * np.sin(theta2[i2])) * np.pi / 180
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
    print("\nSolving torus system...")
    ap = AeroProblem(name='torus', mach=0.5, altitude=10000, areaRef=1.0, chordRef=1.0, evalFuncs=['cl', 'cd'])
    CFDSolver(ap)

    # Extract forces
    forces = CFDSolver.getInstanceForces()
    cl_all = forces['cl']

    # Extract ONLY diagonal values
    cl_diag = [cl_all[idx] for idx in diagonal_indices]

    print("\nDiagonal CL values:")
    for i, idx in enumerate(diagonal_indices):
        i1 = i2 = i
        print(f"  Instance {idx+1} ({i1},{i2}): θ={diagonal_theta[i]:.4f} ({diagonal_theta[i]*180/np.pi:6.1f}°), α={diagonal_alpha[i]*180/np.pi:7.3f}°, CL={cl_diag[i]:.10f}")

    return np.array(cl_diag), np.array(diagonal_alpha), np.array(diagonal_theta)


def run_1d_time_spectral():
    """Run 3-point 1D time spectral"""

    print("\n" + "="*70)
    print("RUNNING 3-POINT 1D TIME SPECTRAL")
    print("="*70)

    gridFile = '../input_files/naca64A010_euler-L2.cgns'

    n_instances = 3
    omega = 100.0
    A = 2.0  # 2*A_torus to match diagonal amplitude

    options = {
        'gridfile': gridFile,
        'equationType': 'Euler',
        'equationMode': 'time spectral',
        'timeIntervals': n_instances,
        'useexternaldynamicmesh': True,
        'l2convergence': 1e-10,
        'nCycles': 10000,
        'writeVolumeSolution': False,
        'writeSurfaceSolution': False,
    }

    CFDSolver = ADFLOW(options=options, debug=False)
    mesh = USMesh(options={"gridFile": gridFile})
    CFDSolver.setMesh(mesh)

    # Generate alpha values at SAME phases as diagonal
    theta = np.linspace(0.0, 2.0 * np.pi, n_instances, endpoint=False)
    alpha_deg = A * np.sin(theta)
    alpha_rad = alpha_deg * np.pi / 180.0

    print("\nPhase and alpha distribution:")
    for i in range(n_instances):
        print(f"  Instance {i+1}: θ={theta[i]:.4f} ({theta[i]*180/np.pi:6.1f}°), α={alpha_deg[i]:7.3f}°")

    # Deform meshes
    MDGroup = CFDSolver.allWallsGroup
    cfdPts0 = CFDSolver.getSurfaceCoordinates(MDGroup, includeZipper=False)
    xRot = 0.25

    cfdPts_list = []

    for i in range(n_instances):
        pitch_rad = alpha_rad[i]
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
    for i in range(n_instances):
        CFDSolver.mesh.setSurfaceCoordinates(cfdPts_list[i])
        CFDSolver.mesh.warpMesh()
        m = CFDSolver.mesh.getSolverGrid()
        CFDSolver.adflow.warping.setgridforoneinstance(m, sps=i + 1)

    CFDSolver._updateGeomInfo = True
    CFDSolver.updateGeometryInfo()

    # Solve
    print("\nSolving 1D time spectral system...")
    ap = AeroProblem(
        name='ts_1d',
        mach=0.5,
        altitude=10000,
        areaRef=1.0,
        chordRef=1.0,
        xRef=0.25,
        xRot=0.25,
        omegaFourier=omega,
        evalFuncs=['cl', 'cd']
    )

    CFDSolver(ap)

    # Extract forces
    forces = CFDSolver.getInstanceForces()
    cl_1d = forces['cl']

    print("\n1D TS CL values:")
    for i in range(n_instances):
        print(f"  Instance {i+1}: θ={theta[i]:.4f} ({theta[i]*180/np.pi:6.1f}°), α={alpha_deg[i]:7.3f}°, CL={cl_1d[i]:.10f}")

    return cl_1d, alpha_rad, theta


def main():
    print("="*70)
    print("COMPARISON: Torus DIAGONAL vs 1D Time Spectral")
    print("="*70)

    # Run both solvers
    cl_torus, alpha_torus, theta_torus = run_degenerate_torus()
    cl_1d, alpha_1d, theta_1d = run_1d_time_spectral()

    # Compare
    print("\n" + "="*70)
    print("DIRECT COMPARISON")
    print("="*70)

    print("\nMatching instances:")
    for i in range(len(cl_torus)):
        print(f"\nInstance {i+1}:")
        print(f"  Torus diagonal ({i},{i}): θ={theta_torus[i]:.4f} ({theta_torus[i]*180/np.pi:6.1f}°), α={alpha_torus[i]*180/np.pi:7.3f}°, CL={cl_torus[i]:.10f}")
        print(f"  1D TS:                    θ={theta_1d[i]:.4f} ({theta_1d[i]*180/np.pi:6.1f}°), α={alpha_1d[i]*180/np.pi:7.3f}°, CL={cl_1d[i]:.10f}")

        # Compare
        theta_diff = abs(theta_torus[i] - theta_1d[i])
        alpha_diff = abs(alpha_torus[i] - alpha_1d[i])
        cl_diff = abs(cl_torus[i] - cl_1d[i])

        print(f"  Δθ = {theta_diff:.3e} rad ({theta_diff*180/np.pi:.3e}°)")
        print(f"  Δα = {alpha_diff:.3e} rad ({alpha_diff*180/np.pi:.3e}°)")
        print(f"  ΔCL = {cl_diff:.3e}")

        if theta_diff < 1e-12 and alpha_diff < 1e-12 and cl_diff < 1e-6:
            print("  ✅ PERFECT MATCH")
        elif theta_diff < 1e-6 and alpha_diff < 1e-6 and cl_diff < 1e-3:
            print("  ✅ GOOD MATCH")
        else:
            print("  ❌ MISMATCH")

    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)

    max_cl_diff = np.max(np.abs(cl_torus - cl_1d))
    max_alpha_diff = np.max(np.abs(alpha_torus - alpha_1d))
    max_theta_diff = np.max(np.abs(theta_torus - theta_1d))

    print(f"\nMax differences:")
    print(f"  Δθ_max = {max_theta_diff:.3e} rad ({max_theta_diff*180/np.pi:.3e}°)")
    print(f"  Δα_max = {max_alpha_diff:.3e} rad ({max_alpha_diff*180/np.pi:.3e}°)")
    print(f"  ΔCL_max = {max_cl_diff:.3e}")

    if max_cl_diff < 1e-6:
        print("\n✅ DIAGONAL TORUS MATCHES 1D TS!")
    else:
        print(f"\n❌ MISMATCH: CL differs by {max_cl_diff:.3e}")

    print("="*70)


if __name__ == "__main__":
    main()
