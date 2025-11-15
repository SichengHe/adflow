# built-ins
import unittest
import os
import copy
import numpy as np

# MACH testing class
from adflow import ADFLOW
from idwarp import USMesh

# import the testing utilities
import reg_test_utils as utils

from reg_default_options import adflowDefOpts

from reg_aeroproblems import ap_naca64A010_time_spectral
import reg_test_classes


baseDir = os.path.dirname(os.path.abspath(__file__))


class TestTorusTimeSpectral(reg_test_classes.RegTest):
    """
    Tests for torus time spectral case for an airfoil with 2-frequency motion.

    Motion: alpha(t) = alpha_mean + A1*sin(omega1*t) + A2*sin(omega2*t)
    where omega1/omega2 = sqrt(2) (incommensurate frequencies)
    """

    N_PROCS = 2
    ref_file = "solve_euler_torus_time_spectral_naca64A010.json"

    def setUp(self):
        # Introduce a transfer class for displacement transfer from struct to aero.
        class TorusTransfer:
            # Torus time spectral transfer class
            # Converting 2-frequency motion to CFD surface nodes on 2D theta grid

            def __init__(self, alpha_grid, xRot, aeroSolver, n1, n2):
                # alpha_grid: [n1, n2] array of alpha values on torus grid

                self.alpha_grid = alpha_grid
                self.n1 = n1
                self.n2 = n2
                self.ntot = n1 * n2

                self.aeroSolver = aeroSolver
                self.xRot = xRot

            def getUndeformedSurfaceNodes(self):
                self.MDGroup = self.aeroSolver.allWallsGroup
                self.cfdPts0 = self.aeroSolver.getSurfaceCoordinates(self.MDGroup, includeZipper=False)

            def setDisplacements(self):
                xRot = self.xRot
                alpha_grid = self.alpha_grid
                cfdPoints_init = self.cfdPts0

                N_pts = cfdPoints_init.shape[0]
                self.cfdPts = []

                # Loop over torus grid (flattened)
                for i2 in range(self.n2):
                    for i1 in range(self.n1):
                        cfdPoints_deformed = np.zeros((N_pts, 3))

                        pitch_loc = alpha_grid[i1, i2]

                        cc = np.cos(pitch_loc)
                        ss = np.sin(pitch_loc)

                        for j in range(N_pts):
                            cfdPoints_deformed[j, 0] = cc * (cfdPoints_init[j, 0] - xRot) + ss * cfdPoints_init[j, 1] + xRot
                            cfdPoints_deformed[j, 1] = -ss * (cfdPoints_init[j, 0] - xRot) + cc * cfdPoints_init[j, 1]
                            cfdPoints_deformed[j, 2] = cfdPoints_init[j, 2]

                        self.cfdPts.append(cfdPoints_deformed)

            def setVolumeMesh(self):
                sps = 0
                for i2 in range(self.n2):
                    for i1 in range(self.n1):
                        self.aeroSolver.mesh.setSurfaceCoordinates(self.cfdPts[sps])
                        self.aeroSolver.mesh.warpMesh()
                        m = self.aeroSolver.mesh.getSolverGrid()
                        self.aeroSolver.adflow.warping.setgridforoneinstance(m, sps=sps + 1)
                        sps += 1

                self.aeroSolver._updateGeomInfo = True
                self.aeroSolver.updateGeometryInfo()

        # Set up the test for the parent RegTest class.
        super().setUp()

        gridFile = os.path.join(baseDir, "../../input_files/naca64A010_euler-L2.cgns")

        # Torus grid dimensions
        n1 = 3
        n2 = 3
        ntot = n1 * n2

        # Two-frequency motion parameters
        omega1 = 100.0  # rad/s
        # omega2 = 100.0 * np.sqrt(2.0)  # rad/s (incommensurate ratio)
        omega2 = 100.0
        alpha_mean = 0.0  # degrees
        A1 = 1.0  # degrees
        A2 = 0.5  # degrees

        options = copy.copy(adflowDefOpts)
        options.update(
            {
                "gridfile": gridFile,
                "outputDirectory": os.path.join(baseDir, "../output_files"),
                "writeVolumeSolution": False,
                "writeSurfaceSolution": False,
                "blocksplitting": True,
                "useblockettes": False,
                "equationtype": "Euler",
                "equationmode": "time spectral",
                "mgcycle": "sg",
                "l2convergence": 1e-15,
                "ncycles": 200000,
                "monitorvariables": ["resrho", "cl"],
                "usenksolver": True,
                "nkswitchtol": 1e-4,
                "NKSubSpaceSize": 400,
                "applypcsubspacesize": 400,
                "useanksolver": True,
                "ankswitchtol": 1e-2,
                "anksubspacesize": 50,
                "alphafollowing": False,
                # Torus time spectral options
                "useTorusTimeSpectral": True,
                "nTimeIntervalsSpectral1": n1,
                "nTimeIntervalsSpectral2": n2,
                "omegaFourier1": omega1,
                "omegaFourier2": omega2,
                "timeintervals": ntot,  # Will be set automatically, but explicit for clarity
                "useexternaldynamicmesh": True,
                "usetsinterpolatedgridvelocity": True,
            }
        )

        # Grid option
        meshOptions = {
            "gridFile": gridFile,
        }

        # Setup aeroproblem (use same as regular time spectral for consistency)
        self.ap = copy.copy(ap_naca64A010_time_spectral)

        # Generate torus grid of alpha values
        # theta1, theta2 in [0, 2*pi]
        theta1 = np.linspace(0.0, 2.0 * np.pi, n1, endpoint=False)
        theta2 = np.linspace(0.0, 2.0 * np.pi, n2, endpoint=False)

        alpha_grid = np.zeros((n1, n2))
        for i1 in range(n1):
            for i2 in range(n2):
                # alpha(theta1, theta2) = alpha_mean + A1*sin(theta1) + A2*sin(theta2)
                alpha_rad = (alpha_mean + A1 * np.sin(theta1[i1]) + A2 * np.sin(theta2[i2])) * np.pi / 180.0
                alpha_grid[i1, i2] = alpha_rad

        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=False)

        # Deform the mesh
        mesh = USMesh(options=meshOptions)
        self.CFDSolver.setMesh(mesh)

        # deformation
        xRot = 0.25  # Hard copied from the reference file.
        TSTransfer = TorusTransfer(alpha_grid, xRot, self.CFDSolver, n1, n2)
        TSTransfer.getUndeformedSurfaceNodes()
        TSTransfer.setDisplacements()
        TSTransfer.setVolumeMesh()

    def test_solve(self):
        # do the solve
        self.CFDSolver(self.ap)

        # check if the solution failed
        self.assert_solution_failure()

        # check its accuracy
        utils.assert_functions_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-8)


if __name__ == "__main__":
    unittest.main()
