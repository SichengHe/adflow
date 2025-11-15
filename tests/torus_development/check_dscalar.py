import numpy as np
from adflow import ADFLOW
from idwarp import USMesh

gridFile = '../input_files/naca64A010_euler-L2.cgns'

options_torus = {
    'gridfile': gridFile,
    'equationType': 'Euler',
    'equationMode': 'time spectral',
    'useTorusTimeSpectral': True,
    'nTimeIntervalsSpectral1': 3,
    'nTimeIntervalsSpectral2': 3,
    'omegaFourier1': 100.0,
    'omegaFourier2': 100.0,
    'printIterations': False,
}

solver_torus = ADFLOW(options=options_torus, debug=False)
mesh_torus = USMesh(options={"gridFile": gridFile})
solver_torus.setMesh(mesh_torus)

dscalar_torus = solver_torus.adflow.inputtimespectral.dscalar

print("Full torus dscalar matrix:")
print(dscalar_torus[0, :, :])

print("\nDiagonal instances: 0, 4, 8")
diag_indices = [0, 4, 8]
print("\nExtract diagonal 3x3 submatrix:")
dscalar_diag = dscalar_torus[0, diag_indices, :][:, diag_indices]
print(dscalar_diag)

print("\nEigenvalues:")
eigs = np.linalg.eigvals(dscalar_diag)
print(eigs)
print(f"Max abs eigenvalue: {np.max(np.abs(eigs)):.6f}")
