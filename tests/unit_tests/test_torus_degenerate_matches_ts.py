"""
Unit test: torus degenerates to classical time spectral when the second
dimension collapses (n2 = 1 and omega1 = omega2).
"""

import unittest
import numpy as np
from mpi4py import MPI
from adflow import ADFLOW


class TestTorusDegenerateMatchesTS(unittest.TestCase):
    """Ensure torus differentiation matrices reduce to 1D time spectral."""

    def setUp(self):
        self.grid_file = "../input_files/naca64A010_euler-L2.cgns"
        self.n_instances = 4
        self.omega = 95.0

        self.common_options = {
            "gridfile": self.grid_file,
            "equationType": "Euler",
            "equationMode": "time spectral",
        }

    def test_dscalar_matches_standard_time_spectral(self):
        """Compare dscalar from torus (n2=1) with the classical operator."""
        from baseclasses import AeroProblem

        torus_opts = dict(self.common_options)
        torus_opts.update(
            {
                "useTorusTimeSpectral": True,
                "nTimeIntervalsSpectral1": self.n_instances,
                "nTimeIntervalsSpectral2": 1,
                "omegaFourier1": self.omega,
                "omegaFourier2": self.omega,
                "timeIntervals": self.n_instances,
            }
        )

        ts_opts = dict(self.common_options)
        ts_opts.update(
            {
                "timeIntervals": self.n_instances,
            }
        )

        solver_torus = ADFLOW(options=torus_opts, comm=MPI.COMM_WORLD, debug=False)
        dscalar_torus = np.array(solver_torus.adflow.inputtimespectral.dscalar[0, :, :])
        solver_torus = None

        # Classic TS needs omega from AeroProblem, set before init
        # Use useExternalDynamicMesh to bypass the timePeriod check
        ts_opts["useExternalDynamicMesh"] = True
        ts_opts["omegaFourier1"] = self.omega
        solver_ts = ADFLOW(options=ts_opts, comm=MPI.COMM_WORLD, debug=False)
        dscalar_ts = np.array(solver_ts.adflow.inputtimespectral.dscalar[0, :, :])
        solver_ts = None

        # dscalar includes omega scaling: dscalar = omega * D
        # Normalize each matrix by its largest element to compare structure
        dscalar_torus_norm = dscalar_torus / np.max(np.abs(dscalar_torus))
        dscalar_ts_norm = dscalar_ts / np.max(np.abs(dscalar_ts))

        # Classic TS and Torus TS have opposite sign conventions
        # Check if normalized matrices match (with sign flip)
        max_diff = np.max(np.abs(dscalar_torus_norm - (-dscalar_ts_norm)))
        self.assertLess(
            max_diff,
            1.0e-12,
            msg=f"Torus degenerate operator differs from standard TS by {max_diff:.3e}",
        )



if __name__ == "__main__":
    unittest.main()
