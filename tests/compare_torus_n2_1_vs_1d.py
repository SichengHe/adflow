"""
Compare Torus with n2=1 (degenerate to 1D) vs classical 1D Time Spectral

With n2=1, the torus should be identical to 1D TS.
"""

import numpy as np
from mpi4py import MPI
from adflow import ADFLOW
from idwarp import USMesh
from baseclasses import AeroProblem


def compare_operators():
    """Compare spectral operators between torus n2=1 and 1D TS"""

    gridFile = '../input_files/naca64A010_euler-L2.cgns'
    n_inst = 3
    omega = 100.0

    print("="*70)
    print("COMPARING TORUS (n1=3, n2=1) vs 1D TS (n=3)")
    print("="*70)

    # Torus with n2=1
    opts_torus = {
        'gridfile': gridFile,
        'equationType': 'Euler',
        'equationMode': 'time spectral',
        'useTorusTimeSpectral': True,
        'nTimeIntervalsSpectral1': n_inst,
        'nTimeIntervalsSpectral2': 1,
        'omegaFourier1': omega,
        'omegaFourier2': omega,
        'timeIntervals': n_inst,
    }

    # 1D TS
    opts_1d = {
        'gridfile': gridFile,
        'equationType': 'Euler',
        'equationMode': 'time spectral',
        'timeIntervals': n_inst,
        'useExternalDynamicMesh': True,
        'omegaFourier1': omega,
    }

    print("\nInitializing solvers...")
    solver_torus = ADFLOW(options=opts_torus, comm=MPI.COMM_WORLD, debug=False)
    dscalar_torus = np.array(solver_torus.adflow.inputtimespectral.dscalar[0, :, :])

    solver_1d = ADFLOW(options=opts_1d, comm=MPI.COMM_WORLD, debug=False)
    dscalar_1d = np.array(solver_1d.adflow.inputtimespectral.dscalar[0, :, :])

    print("\nTorus (n2=1) dscalar:")
    print(dscalar_torus)

    print("\n1D TS dscalar:")
    print(dscalar_1d)

    print("\nDifference:")
    print(dscalar_torus - dscalar_1d)

    print("\nDifference (with sign flip for convention):")
    print(dscalar_torus - (-dscalar_1d))

    # Normalize and compare
    dscalar_torus_norm = dscalar_torus / np.max(np.abs(dscalar_torus))
    dscalar_1d_norm = dscalar_1d / np.max(np.abs(dscalar_1d))

    diff_direct = np.max(np.abs(dscalar_torus_norm - dscalar_1d_norm))
    diff_flipped = np.max(np.abs(dscalar_torus_norm - (-dscalar_1d_norm)))

    print(f"\nMax normalized difference (direct): {diff_direct:.3e}")
    print(f"Max normalized difference (sign flipped): {diff_flipped:.3e}")

    if diff_flipped < 1e-12:
        print("\n✅ Matrices match (with sign flip)!")
    elif diff_direct < 1e-12:
        print("\n✅ Matrices match (directly)!")
    else:
        print(f"\n❌ Matrices DON'T match (diff={min(diff_direct, diff_flipped):.3e})")


if __name__ == "__main__":
    compare_operators()
