"""
Verification test: Compare torus time spectral vs unsteady time marching

This test runs a simple 2-frequency pitching motion using both:
1. Torus time spectral method (quasi-periodic)
2. Unsteady time marching (long integration)

The unsteady solution should converge to the torus after sufficient time.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI
import os

# Test parameters
omega1 = 100.0  # rad/s
omega2 = 100.0 * np.sqrt(2.0)  # incommensurate
alpha_mean = 0.0  # deg
A1 = 1.0  # deg
A2 = 0.5  # deg

n1 = 3  # torus grid points in theta1
n2 = 3  # torus grid points in theta2

def alpha_motion(t):
    """Two-frequency pitch motion"""
    return alpha_mean + A1 * np.sin(omega1 * t) + A2 * np.sin(omega2 * t)

def run_torus_spectral():
    """Run torus time spectral simulation"""
    from adflow import ADFLOW
    from idwarp import USMesh
    from baseclasses import AeroProblem

    print("\n" + "="*60)
    print("RUNNING TORUS TIME SPECTRAL SIMULATION")
    print("="*60)

    # Grid file
    gridFile = "../input_files/naca64A010_euler-L2.cgns"

    options = {
        "gridfile": gridFile,
        "outputDirectory": "./output_torus",
        "writeVolumeSolution": False,
        "writeSurfaceSolution": False,
        "equationtype": "Euler",
        "equationmode": "time spectral",
        "l2convergence": 1e-10,
        "ncycles": 10000,
        "monitorvariables": ["resrho", "cl", "cd"],
        "usenksolver": True,
        "nkswitchtol": 1e-4,
        "useTorusTimeSpectral": True,
        "nTimeIntervalsSpectral1": n1,
        "nTimeIntervalsSpectral2": n2,
        "omegaFourier1": omega1,
        "omegaFourier2": omega2,
        "useexternaldynamicmesh": True,
        "usetsinterpolatedgridvelocity": True,
    }

    # Create solver
    CFDSolver = ADFLOW(options=options, comm=MPI.COMM_WORLD)

    # Create aeroproblem
    ap = AeroProblem(
        name="torus_pitch",
        mach=0.796,
        reynolds=12.56e6,
        reynoldsLength=1.0,
        T=305.0,
        areaRef=1.0,
        chordRef=1.0,
        xRef=0.248,
        evalFuncs=["cl", "cd", "cmz"],
    )

    # Setup mesh deformation (simplified - assume small angle)
    mesh = USMesh(options={"gridFile": gridFile})
    CFDSolver.setMesh(mesh)

    # Generate torus grid
    theta1 = np.linspace(0.0, 2.0 * np.pi, n1, endpoint=False)
    theta2 = np.linspace(0.0, 2.0 * np.pi, n2, endpoint=False)

    # For each torus point, compute alpha and deform mesh
    # (In practice, use the Transfer class from the test)

    # Solve
    CFDSolver(ap)

    # Extract forces at each torus point
    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)

    print("\nTorus Time Spectral Results:")
    print(f"  CL = {funcs[ap.name + '_cl']:.6f}")
    print(f"  CD = {funcs[ap.name + '_cd']:.6f}")

    return funcs

def run_unsteady_time_marching():
    """Run unsteady time marching simulation"""
    from adflow import ADFLOW
    from baseclasses import AeroProblem

    print("\n" + "="*60)
    print("RUNNING UNSTEADY TIME MARCHING SIMULATION")
    print("="*60)

    gridFile = "../input_files/naca64A010_euler-L2.cgns"

    # Time marching parameters
    # Need to resolve both frequencies
    dt_max = 2.0 * np.pi / omega2 / 20.0  # 20 points per fastest period
    dt = 0.0001  # s
    n_periods_1 = 5  # integrate over 5 periods of slower frequency
    T_period_1 = 2.0 * np.pi / omega1
    t_final = n_periods_1 * T_period_1
    n_steps = int(t_final / dt)

    print(f"Time step: {dt:.6f} s")
    print(f"Total time: {t_final:.6f} s")
    print(f"Number of steps: {n_steps}")
    print(f"Period 1 (omega1): {T_period_1:.6f} s")
    print(f"Period 2 (omega2): {2.0*np.pi/omega2:.6f} s")

    options = {
        "gridfile": gridFile,
        "outputDirectory": "./output_unsteady",
        "writeVolumeSolution": False,
        "writeSurfaceSolution": False,
        "equationtype": "Euler",
        "equationmode": "unsteady",
        "timeintegrationscheme": "BDF",
        "timeaccuracy": 2,
        "nTimeStepsFine": n_steps,
        "deltaT": dt,
        "l2convergence": 1e-10,
        "ncycles": 1000,
        "monitorvariables": ["resrho", "cl", "cd"],
    }

    # Create solver
    CFDSolver = ADFLOW(options=options, comm=MPI.COMM_WORLD)

    # Time history
    cl_history = []
    cd_history = []
    t_history = []
    alpha_history = []

    # Time integration loop
    for step in range(n_steps):
        t = step * dt
        alpha_t = alpha_motion(t)

        # Create aeroproblem for this time step
        ap = AeroProblem(
            name="unsteady_pitch",
            mach=0.796,
            alpha=alpha_t,
            reynolds=12.56e6,
            reynoldsLength=1.0,
            T=305.0,
            areaRef=1.0,
            chordRef=1.0,
            xRef=0.248,
            evalFuncs=["cl", "cd"],
        )

        # Solve this time step
        CFDSolver(ap)

        # Extract forces
        funcs = {}
        CFDSolver.evalFunctions(ap, funcs)

        cl_history.append(funcs[ap.name + '_cl'])
        cd_history.append(funcs[ap.name + '_cd'])
        t_history.append(t)
        alpha_history.append(alpha_t)

        if step % 100 == 0:
            print(f"  Step {step}/{n_steps}, t={t:.4f}s, alpha={alpha_t:.3f}°, CL={funcs[ap.name + '_cl']:.6f}")

    results = {
        't': np.array(t_history),
        'alpha': np.array(alpha_history),
        'cl': np.array(cl_history),
        'cd': np.array(cd_history),
    }

    print("\nUnsteady Time Marching Complete!")
    print(f"  Mean CL = {np.mean(cl_history):.6f}")
    print(f"  Mean CD = {np.mean(cd_history):.6f}")

    return results

def compare_results(torus_funcs, unsteady_results):
    """Compare torus vs unsteady results"""
    print("\n" + "="*60)
    print("COMPARISON: TORUS vs UNSTEADY")
    print("="*60)

    # Plot time history
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))

    # Alpha history
    axes[0].plot(unsteady_results['t'], unsteady_results['alpha'], 'b-', label='Alpha')
    axes[0].set_ylabel('Alpha (deg)')
    axes[0].set_title('Two-Frequency Pitch Motion: α(t) = %.1f° + %.1f°sin(%.1f t) + %.1f°sin(%.1f t)' %
                       (alpha_mean, A1, omega1, A2, omega2))
    axes[0].grid(True)
    axes[0].legend()

    # CL history
    axes[1].plot(unsteady_results['t'], unsteady_results['cl'], 'b-', label='Unsteady CL')
    axes[1].axhline(y=np.mean(unsteady_results['cl']), color='r', linestyle='--', label='Unsteady Mean')
    # Note: Torus gives values at grid points, would need to reconstruct full signal
    axes[1].set_ylabel('CL')
    axes[1].set_title('Lift Coefficient')
    axes[1].grid(True)
    axes[1].legend()

    # CD history
    axes[2].plot(unsteady_results['t'], unsteady_results['cd'], 'b-', label='Unsteady CD')
    axes[2].axhline(y=np.mean(unsteady_results['cd']), color='r', linestyle='--', label='Unsteady Mean')
    axes[2].set_xlabel('Time (s)')
    axes[2].set_ylabel('CD')
    axes[2].set_title('Drag Coefficient')
    axes[2].grid(True)
    axes[2].legend()

    plt.tight_layout()
    plt.savefig('torus_vs_unsteady_comparison.png', dpi=150)
    print("\nPlot saved to: torus_vs_unsteady_comparison.png")

    # Quantitative comparison
    cl_mean_unsteady = np.mean(unsteady_results['cl'])
    cd_mean_unsteady = np.mean(unsteady_results['cd'])

    print(f"\nQuantitative Comparison:")
    print(f"  Unsteady Mean CL: {cl_mean_unsteady:.6f}")
    print(f"  Unsteady Mean CD: {cd_mean_unsteady:.6f}")
    # Would compare to torus reconstruction here

if __name__ == "__main__":
    # For now, just print the setup
    print("\n" + "="*60)
    print("TORUS TIME SPECTRAL vs UNSTEADY VERIFICATION TEST")
    print("="*60)
    print(f"\nMotion: α(t) = {alpha_mean}° + {A1}°·sin({omega1} rad/s · t) + {A2}°·sin({omega2:.2f} rad/s · t)")
    print(f"\nTorus grid: {n1} x {n2} = {n1*n2} points")
    print(f"Frequency ratio: ω1/ω2 = 1/√2 (incommensurate)")
    print("\nThis test requires:")
    print("  1. Compiled ADflow with torus time spectral support")
    print("  2. NACA 64A010 grid file")
    print("  3. Significant computational resources for unsteady simulation")
    print("\nTo run:")
    print("  mpirun -np 4 python test_torus_vs_unsteady.py")
    print("="*60)

    # Uncomment to actually run:
    # torus_funcs = run_torus_spectral()
    # unsteady_results = run_unsteady_time_marching()
    # compare_results(torus_funcs, unsteady_results)
