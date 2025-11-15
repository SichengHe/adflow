"""
Compare Torus Time Spectral (omega1 = omega2) vs Regular Time Spectral

When omega1 = omega2, the torus grid collapses to a 1D time spectral problem.
This test verifies that torus mode produces the same results as regular time spectral.
"""

import numpy as np
from mpi4py import MPI
from adflow import ADFLOW
from idwarp import USMesh
from baseclasses import AeroProblem
import pickle


def run_torus_mode(n_instances, omega, alpha_mean, amplitude):
    """Run torus mode with omega1 = omega2 (degenerate to 1D)"""

    print("="*70)
    print("RUNNING TORUS MODE (omega1 = omega2)")
    print("="*70)

    gridFile = '../input_files/naca64A010_euler-L2.cgns'

    # Torus grid: n x 1 (degenerate)
    n1 = n_instances
    n2 = 1
    ntot = n1 * n2

    print(f"  Grid: {n1} x {n2} = {ntot} instances")
    print(f"  omega1 = omega2 = {omega} rad/s")
    print(f"  Motion: α(θ) = {alpha_mean}° + {amplitude}°·sin(θ)")

    options = {
        'gridfile': gridFile,
        'outputDirectory': './output_torus_vs_ts',
        'equationType': 'Euler',
        'equationMode': 'time spectral',
        'useTorusTimeSpectral': True,
        'nTimeIntervalsSpectral1': n1,
        'nTimeIntervalsSpectral2': n2,
        'omegaFourier1': omega,
        'omegaFourier2': omega,
        'useexternaldynamicmesh': True,
        'usetsinterpolatedgridvelocity': True,
        'mgcycle': 'sg',
        'nCycles': 10000,
        'l2convergence': 1e-10,
        'l2convergencecoarse': 1e-2,
        'monitorvariables': ['resrho', 'cl', 'cd'],
        'usenksolver': True,
        'nkswitchtol': 1e-4,
        'NKSubSpaceSize': 400,
        'applypcsubspacesize': 400,
        'useanksolver': True,
        'ankswitchtol': 1e-2,
        'anksubspacesize': 200,
        'ANKCFL0': 1.0,
        'ANKCFLLimit': 10.0,
        'CFL': 1.5,
        'nSubiter': 3,
        'alphafollowing': False,
        'blocksplitting': True,
        'useblockettes': False,
        'writeVolumeSolution': False,
        'writeSurfaceSolution': False,
        'numberSolutions': False,
    }

    meshOptions = {"gridFile": gridFile}

    print("\nInitializing ADflow (torus mode)...")
    CFDSolver = ADFLOW(options=options, debug=False)

    print("Setting up mesh deformation...")
    mesh = USMesh(options=meshOptions)
    CFDSolver.setMesh(mesh)

    # Generate alpha values
    theta = np.linspace(0.0, 2.0 * np.pi, n1, endpoint=False)
    alpha_deg = alpha_mean + amplitude * np.sin(theta)
    alpha_rad = alpha_deg * np.pi / 180.0

    print("\nAlpha distribution (degrees):")
    for i in range(n1):
        print(f"  Instance {i}: θ={theta[i]:.3f}, α={alpha_deg[i]:.3f}°")

    # Deform meshes
    print("\nDeforming meshes...")
    MDGroup = CFDSolver.allWallsGroup
    cfdPts0 = CFDSolver.getSurfaceCoordinates(MDGroup, includeZipper=False)

    xRot = 0.25
    cfdPts_list = []

    for i1 in range(n1):
        pitch_rad = alpha_rad[i1]
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
    for i1 in range(n1):
        CFDSolver.mesh.setSurfaceCoordinates(cfdPts_list[i1])
        CFDSolver.mesh.warpMesh()
        m = CFDSolver.mesh.getSolverGrid()
        CFDSolver.adflow.warping.setgridforoneinstance(m, sps=i1 + 1)

    CFDSolver._updateGeomInfo = True
    CFDSolver.updateGeometryInfo()

    # Solve
    print("\nSolving torus system...")
    ap = AeroProblem(
        name='torus_test',
        mach=0.5,
        altitude=10000,
        areaRef=1.0,
        chordRef=1.0,
        xRef=0.25,
        evalFuncs=['cl', 'cd', 'cmz']
    )

    CFDSolver(ap)

    # Extract forces
    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)

    cl = funcs.get('torus_test_cl', 0.0)
    cd = funcs.get('torus_test_cd', 0.0)
    cm = funcs.get('torus_test_cmz', 0.0)

    print(f"\nTorus results:")
    print(f"  CL = {cl:.8f}")
    print(f"  CD = {cd:.8f}")
    print(f"  CM = {cm:.8f}")

    return {'cl': cl, 'cd': cd, 'cm': cm, 'alpha': alpha_deg}


def run_time_spectral_mode(n_instances, omega, alpha_mean, amplitude):
    """Run regular time spectral mode"""

    print("\n" + "="*70)
    print("RUNNING REGULAR TIME SPECTRAL MODE")
    print("="*70)

    gridFile = '../input_files/naca64A010_euler-L2.cgns'

    print(f"  Instances: {n_instances}")
    print(f"  omega = {omega} rad/s")
    print(f"  Motion: α(θ) = {alpha_mean}° + {amplitude}°·sin(θ)")

    options = {
        'gridfile': gridFile,
        'outputDirectory': './output_torus_vs_ts',
        'equationType': 'Euler',
        'equationMode': 'time spectral',
        'timeIntervals': n_instances,
        'alphaFollowing': False,
        'useexternaldynamicmesh': True,
        'usetsinterpolatedgridvelocity': True,
        'mgcycle': 'sg',
        'nCycles': 10000,
        'l2convergence': 1e-10,
        'l2convergencecoarse': 1e-2,
        'monitorvariables': ['resrho', 'cl', 'cd'],
        'usenksolver': True,
        'nkswitchtol': 1e-4,
        'NKSubSpaceSize': 400,
        'applypcsubspacesize': 400,
        'useanksolver': True,
        'ankswitchtol': 1e-2,
        'anksubspacesize': 200,
        'ANKCFL0': 1.0,
        'ANKCFLLimit': 10.0,
        'CFL': 1.5,
        'nSubiter': 3,
        'blocksplitting': True,
        'useblockettes': False,
        'writeVolumeSolution': False,
        'writeSurfaceSolution': False,
        'numberSolutions': False,
    }

    meshOptions = {"gridFile": gridFile}

    print("\nInitializing ADflow (time spectral mode)...")
    CFDSolver = ADFLOW(options=options, debug=False)

    print("Setting up mesh deformation...")
    mesh = USMesh(options=meshOptions)
    CFDSolver.setMesh(mesh)

    # Generate alpha values
    theta = np.linspace(0.0, 2.0 * np.pi, n_instances, endpoint=False)
    alpha_deg = alpha_mean + amplitude * np.sin(theta)
    alpha_rad = alpha_deg * np.pi / 180.0

    print("\nAlpha distribution (degrees):")
    for i in range(n_instances):
        print(f"  Instance {i}: θ={theta[i]:.3f}, α={alpha_deg[i]:.3f}°")

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
    print("\nSolving time spectral system...")
    ap = AeroProblem(
        name='ts_test',
        mach=0.5,
        altitude=10000,
        areaRef=1.0,
        chordRef=1.0,
        xRef=0.25,
        evalFuncs=['cl', 'cd', 'cmz']
    )

    CFDSolver(ap)

    # Extract forces
    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)

    cl = funcs.get('ts_test_cl', 0.0)
    cd = funcs.get('ts_test_cd', 0.0)
    cm = funcs.get('ts_test_cmz', 0.0)

    print(f"\nTime spectral results:")
    print(f"  CL = {cl:.8f}")
    print(f"  CD = {cd:.8f}")
    print(f"  CM = {cm:.8f}")

    return {'cl': cl, 'cd': cd, 'cm': cm, 'alpha': alpha_deg}


def main():
    print("="*70)
    print("TORUS vs TIME SPECTRAL COMPARISON TEST")
    print("Testing that omega1 = omega2 gives same results as 1D time spectral")
    print("="*70)

    # Test parameters
    n_instances = 5
    omega = 100.0  # rad/s
    alpha_mean = 0.0  # degrees
    amplitude = 1.0  # degrees

    # Run both modes
    torus_results = run_torus_mode(n_instances, omega, alpha_mean, amplitude)
    ts_results = run_time_spectral_mode(n_instances, omega, alpha_mean, amplitude)

    # Compare
    print("\n" + "="*70)
    print("COMPARISON")
    print("="*70)

    print("\n{:20s} {:>15s} {:>15s} {:>15s}".format("Variable", "Torus", "Time Spectral", "Difference"))
    print("-"*70)

    cl_diff = abs(torus_results['cl'] - ts_results['cl'])
    cd_diff = abs(torus_results['cd'] - ts_results['cd'])
    cm_diff = abs(torus_results['cm'] - ts_results['cm'])

    print("{:20s} {:15.8f} {:15.8f} {:15.2e}".format("CL", torus_results['cl'], ts_results['cl'], cl_diff))
    print("{:20s} {:15.8f} {:15.8f} {:15.2e}".format("CD", torus_results['cd'], ts_results['cd'], cd_diff))
    print("{:20s} {:15.8f} {:15.8f} {:15.2e}".format("CM", torus_results['cm'], ts_results['cm'], cm_diff))

    # Check tolerance
    tol = 1e-8
    all_match = cl_diff < tol and cd_diff < tol and cm_diff < tol

    print("\n" + "="*70)
    if all_match:
        print("✅ SUCCESS: Torus mode matches time spectral mode!")
        print(f"   All differences < {tol}")
    else:
        print("❌ MISMATCH: Torus mode differs from time spectral mode")
        print(f"   Tolerance: {tol}")
        if cl_diff >= tol:
            print(f"   CL difference: {cl_diff:.2e} (FAIL)")
        if cd_diff >= tol:
            print(f"   CD difference: {cd_diff:.2e} (FAIL)")
        if cm_diff >= tol:
            print(f"   CM difference: {cm_diff:.2e} (FAIL)")
    print("="*70)

    # Save results
    results = {
        'torus': torus_results,
        'ts': ts_results,
        'differences': {'cl': cl_diff, 'cd': cd_diff, 'cm': cm_diff},
        'match': all_match,
        'n_instances': n_instances,
        'omega': omega,
        'alpha_mean': alpha_mean,
        'amplitude': amplitude,
    }

    with open('torus_vs_ts_comparison.pkl', 'wb') as f:
        pickle.dump(results, f)

    print("\n✅ Results saved to torus_vs_ts_comparison.pkl")


if __name__ == "__main__":
    main()
