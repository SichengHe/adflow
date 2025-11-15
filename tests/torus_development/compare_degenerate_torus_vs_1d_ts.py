"""
Compare Degenerate Torus (3x3 with A1=A2, omega1=omega2) vs 1D Time Spectral (3 instances)

The degenerate torus with same frequencies and amplitudes should give a subset
of CL values that match the 1D time spectral case.
"""

import numpy as np
from mpi4py import MPI
from adflow import ADFLOW
from idwarp import USMesh
from baseclasses import AeroProblem
import sys
import io
import re


def extract_cl_from_output(output_text, n_instances):
    """Extract CL values from NK iteration lines"""
    instance_data = {}
    # Correct pattern: 1  instance  iter  total  NK/*NK  ----  CFL  ?  resrho  CL  CD  total_residual
    # Note: Use \s+ to match multiple spaces between fields
    pattern = r'^\s*1\s+(\d+)\s+\d+\s+\d+\s+\*?NK\s+----\s+([\d\.]+)\s+([\d\.]+)\s+([-\d\.E\+]+)\s+([-\d\.E\+]+)\s+([-\d\.E\+]+)\s+([-\d\.E\+]+)'

    for line in output_text.split('\n'):
        match = re.search(pattern, line)
        if match:
            sps = int(match.group(1))
            cl_val = float(match.group(5))  # group(5) is CL
            cd_val = float(match.group(6))  # group(6) is CD
            instance_data[sps] = {'cl': cl_val, 'cd': cd_val}

    # Extract in order
    cl_values = []
    for sps in range(1, n_instances + 1):
        if sps in instance_data:
            cl_values.append(instance_data[sps]['cl'])
        else:
            cl_values.append(np.nan)

    return np.array(cl_values)


def run_degenerate_torus():
    """Run degenerate torus case: 3x3 grid with omega1=omega2, A1=A2"""

    print("="*70)
    print("RUNNING DEGENERATE TORUS (3x3, A1=A2=1.0, omega1=omega2=100)")
    print("="*70)

    gridFile = '../input_files/naca64A010_euler-L2.cgns'

    n1 = 3
    n2 = 3
    ntot = n1 * n2

    omega1 = 100.0
    omega2 = 100.0  # SAME
    A1 = 1.0
    A2 = 1.0  # SAME

    print(f"\nGrid: {n1} × {n2} = {ntot} instances")
    print(f"Motion: α(θ₁,θ₂) = {A1}°·[sin(θ₁) + sin(θ₂)]")

    options = {
        'gridfile': gridFile,
        'outputDirectory': './output_compare',
        'equationType': 'Euler',
        'equationMode': 'time spectral',
        'useTorusTimeSpectral': True,
        'nTimeIntervalsSpectral1': n1,
        'nTimeIntervalsSpectral2': n2,
        'omegaFourier1': omega1,
        'omegaFourier2': omega2,
        'useexternaldynamicmesh': True,
        'usetsinterpolatedgridvelocity': True,
        'mgcycle': 'sg',
        'nCycles': 10000,
        'l2convergence': 1e-10,
        'l2convergencecoarse': 1e-2,
        'monitorvariables': ['resrho', 'cl', 'cd'],
        'usenksolver': True,
        'nkswitchtol': 1e-4,
        'useanksolver': True,
        'ankswitchtol': 1e-2,
        'CFL': 1.5,
        'nSubiter': 3,
        'blocksplitting': True,
        'writeVolumeSolution': False,
        'writeSurfaceSolution': False,
    }

    meshOptions = {"gridFile": gridFile}

    print("\nInitializing ADflow (torus mode)...")
    CFDSolver = ADFLOW(options=options, debug=False)

    print("Setting up mesh...")
    mesh = USMesh(options=meshOptions)
    CFDSolver.setMesh(mesh)

    # Generate alpha grid
    theta1 = np.linspace(0.0, 2.0 * np.pi, n1, endpoint=False)
    theta2 = np.linspace(0.0, 2.0 * np.pi, n2, endpoint=False)

    alpha_grid = np.zeros((n1, n2))
    theta_proj = np.zeros((n1, n2))

    print("\nPhase and alpha distribution:")
    for i1 in range(n1):
        for i2 in range(n2):
            alpha_deg = A1 * np.sin(theta1[i1]) + A2 * np.sin(theta2[i2])
            alpha_grid[i1, i2] = alpha_deg * np.pi / 180.0
            theta_proj[i1, i2] = (theta1[i1] + theta2[i2]) / 2
            idx = i2 * n1 + i1
            print(f"  Instance {idx+1} ({i1},{i2}): θ={theta_proj[i1,i2]:.4f} ({theta_proj[i1,i2]*180/np.pi:6.1f}°), α={alpha_deg:7.3f}°")

    # Deform meshes
    print("\nDeforming meshes...")
    MDGroup = CFDSolver.allWallsGroup
    cfdPts0 = CFDSolver.getSurfaceCoordinates(MDGroup, includeZipper=False)

    xRot = 0.25
    cfdPts_list = []

    for i2 in range(n2):
        for i1 in range(n1):
            pitch_rad = alpha_grid[i1, i2]
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
    for sps in range(ntot):
        CFDSolver.mesh.setSurfaceCoordinates(cfdPts_list[sps])
        CFDSolver.mesh.warpMesh()
        m = CFDSolver.mesh.getSolverGrid()
        CFDSolver.adflow.warping.setgridforoneinstance(m, sps=sps + 1)

    CFDSolver._updateGeomInfo = True
    CFDSolver.updateGeometryInfo()

    # Solve
    print("\nSolving torus system...")
    ap = AeroProblem(
        name='torus_degen',
        mach=0.5,
        altitude=10000,
        areaRef=1.0,
        chordRef=1.0,
        evalFuncs=['cl', 'cd']
    )

    CFDSolver(ap)

    # Extract CL values using NEW API
    forces_torus = CFDSolver.getInstanceForces()
    cl_torus = forces_torus['cl']

    print("\nTorus CL values extracted:")
    for i in range(ntot):
        i1 = i % n1
        i2 = i // n1
        alpha_deg = alpha_grid[i1, i2] * 180 / np.pi
        theta = theta_proj[i1, i2]
        print(f"  Instance {i+1} ({i1},{i2}): θ={theta:.4f} ({theta*180/np.pi:6.1f}°), α={alpha_deg:7.3f}°, CL={cl_torus[i]:.10f}")

    return cl_torus, alpha_grid, theta_proj


def run_1d_time_spectral():
    """Run regular 1D time spectral with 6 instances"""

    print("\n" + "="*70)
    print("RUNNING 1D TIME SPECTRAL (6 instances)")
    print("="*70)

    gridFile = '../input_files/naca64A010_euler-L2.cgns'

    n_instances = 6
    omega = 100.0
    A = 2.0  # Amplitude = 2*A_torus to match sin(θ₁)+sin(θ₂) range

    print(f"\nInstances: {n_instances}")
    print(f"omega = {omega} rad/s")
    print(f"Motion: α(θ) = {A}°·sin(θ)")

    options = {
        'gridfile': gridFile,
        'outputDirectory': './output_compare',
        'equationType': 'Euler',
        'equationMode': 'time spectral',
        'timeIntervals': n_instances,
        'useexternaldynamicmesh': True,
        'usetsinterpolatedgridvelocity': True,
        'mgcycle': 'sg',
        'nCycles': 10000,
        'l2convergence': 1e-10,
        'l2convergencecoarse': 1e-2,
        'monitorvariables': ['resrho', 'cl', 'cd'],
        'usenksolver': True,
        'nkswitchtol': 1e-4,
        'useanksolver': True,
        'ankswitchtol': 1e-2,
        'CFL': 1.5,
        'nSubiter': 3,
        'blocksplitting': True,
        'writeVolumeSolution': False,
        'writeSurfaceSolution': False,
    }

    meshOptions = {"gridFile": gridFile}

    print("\nInitializing ADflow (1D time spectral)...")
    CFDSolver = ADFLOW(options=options, debug=False)

    print("Setting up mesh...")
    mesh = USMesh(options=meshOptions)
    CFDSolver.setMesh(mesh)

    # Generate alpha values
    theta = np.linspace(0.0, 2.0 * np.pi, n_instances, endpoint=False)
    alpha_deg = A * np.sin(theta)
    alpha_rad = alpha_deg * np.pi / 180.0

    print("\nPhase and alpha distribution:")
    for i in range(n_instances):
        print(f"  Instance {i+1}: θ={theta[i]:.4f} ({theta[i]*180/np.pi:6.1f}°), α={alpha_deg[i]:7.3f}°")

    # Deform meshes
    print("\nDeforming meshes...")
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
        omegaFourier=omega,  # CRITICAL: Need this for external mesh case!
        evalFuncs=['cl', 'cd']
    )

    CFDSolver(ap)

    # Extract CL values using NEW API
    forces_1d = CFDSolver.getInstanceForces()
    cl_1d = forces_1d['cl']

    print("\n1D Time Spectral CL values extracted:")
    for i in range(n_instances):
        print(f"  Instance {i+1}: θ={theta[i]:.4f} ({theta[i]*180/np.pi:6.1f}°), α={alpha_deg[i]:7.3f}°, CL={cl_1d[i]:.10f}")

    return cl_1d, alpha_deg, theta


def main():
    print("="*70)
    print("COMPARISON: Degenerate Torus vs 1D Time Spectral")
    print("="*70)

    # Run degenerate torus
    cl_torus, alpha_grid, theta_proj_grid = run_degenerate_torus()

    # Run 1D time spectral
    cl_1d, alpha_1d, theta_1d = run_1d_time_spectral()

    # Compare
    print("\n" + "="*70)
    print("COMPARISON RESULTS")
    print("="*70)

    print("\nDegenerate Torus CL values (3x3 = 9 instances):")
    for i in range(len(cl_torus)):
        i1 = i % 3
        i2 = i // 3
        alpha_deg = alpha_grid[i1, i2] * 180 / np.pi
        theta = theta_proj_grid[i1, i2]
        print(f"  Instance {i+1} ({i1},{i2}): θ={theta:.4f} ({theta*180/np.pi:6.1f}°), α={alpha_deg:7.3f}°, CL={cl_torus[i]:.10f}")

    print(f"\n1D Time Spectral CL values ({len(cl_1d)} instances):")
    for i in range(len(cl_1d)):
        print(f"  Instance {i+1}: θ={theta_1d[i]:.4f} ({theta_1d[i]*180/np.pi:6.1f}°), α={alpha_1d[i]:7.3f}°, CL={cl_1d[i]:.10f}")

    # Find matching phase values
    print("\n" + "="*70)
    print("MATCHING PHASE VALUES (θ)")
    print("="*70)

    print("\nLooking for torus instances that match 1D TS phase values...")
    print()

    for j in range(len(cl_1d)):
        theta_1d_val = theta_1d[j]
        print(f"1D TS Instance {j+1}: θ={theta_1d_val:.4f} ({theta_1d_val*180/np.pi:6.1f}°), α={alpha_1d[j]:7.3f}°, CL={cl_1d[j]:.10f}")

        # Find matching torus instances
        matches = []
        for i in range(len(cl_torus)):
            i1 = i % 3
            i2 = i // 3
            theta_torus = theta_proj_grid[i1, i2]
            alpha_torus = alpha_grid[i1, i2] * 180 / np.pi

            # Match by phase (mod 2π)
            theta_diff = min(abs(theta_torus - theta_1d_val),
                           abs(theta_torus - theta_1d_val + 2*np.pi),
                           abs(theta_torus - theta_1d_val - 2*np.pi))

            if theta_diff < 0.01:  # Within 0.01 radians
                matches.append((i+1, i1, i2, theta_torus, alpha_torus, cl_torus[i]))

        if matches:
            for match in matches:
                inst, i1, i2, theta_t, alpha_t, cl_t = match
                diff = abs(cl_t - cl_1d[j])
                status = "✅ MATCH" if diff < 1e-6 else "❌ MISMATCH"
                print(f"  → Torus Instance {inst} ({i1},{i2}): θ={theta_t:.4f} ({theta_t*180/np.pi:6.1f}°), α={alpha_t:7.3f}°, CL={cl_t:.10f}, Δ={diff:.2e} {status}")
        print()

    print("="*70)


if __name__ == "__main__":
    main()
