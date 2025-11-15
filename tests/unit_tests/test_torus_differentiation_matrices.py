"""
Unit Test 1: Verify Torus Spectral Differentiation Matrices

Tests:
1. Matrix dimensions are correct
2. Derivative of constant is zero (consistency check)
3. Derivative of sine/cosine functions matches analytical derivatives
4. Sparsity pattern is correct (cross pattern, not dense)
5. dvector is properly constructed from dscalar
"""

import numpy as np
import unittest
from mpi4py import MPI


class TestTorusDifferentiationMatrices(unittest.TestCase):
    """Test spectral differentiation matrices for torus time spectral"""

    def setUp(self):
        """Setup test configuration"""
        self.n1 = 3
        self.n2 = 3
        self.ntot = self.n1 * self.n2
        self.omega1 = 100.0  # rad/s
        self.omega2 = 100.0 * np.sqrt(2.0)  # rad/s

        # Create theta grids
        self.theta1 = np.linspace(0, 2*np.pi, self.n1, endpoint=False)
        self.theta2 = np.linspace(0, 2*np.pi, self.n2, endpoint=False)

    def test_matrix_dimensions(self):
        """Test 1: Verify matrix dimensions"""
        print("\n" + "="*70)
        print("TEST 1: Differentiation Matrix Dimensions")
        print("="*70)

        # Import ADflow and create solver
        from adflow import ADFLOW

        options = {
            'gridfile': '../input_files/naca64A010_euler-L2.cgns',
            'equationmode': 'time spectral',
            'useTorusTimeSpectral': True,
            'nTimeIntervalsSpectral1': self.n1,
            'nTimeIntervalsSpectral2': self.n2,
            'omegaFourier1': self.omega1,
            'omegaFourier2': self.omega2,
        }

        try:
            solver = ADFLOW(options=options, comm=MPI.COMM_WORLD)

            # Access Fortran differentiation matrices
            dscalar = solver.adflow.inputtimespectral.dscalar
            dvector = solver.adflow.inputtimespectral.dvector

            # Check dscalar dimensions: (nSections, ntot, ntot)
            self.assertEqual(dscalar.shape[1], self.ntot,
                           f"dscalar rows: expected {self.ntot}, got {dscalar.shape[1]}")
            self.assertEqual(dscalar.shape[2], self.ntot,
                           f"dscalar cols: expected {self.ntot}, got {dscalar.shape[2]}")

            # Check dvector dimensions: (nSections, 3*ntot, 3*ntot)
            self.assertEqual(dvector.shape[1], 3*self.ntot,
                           f"dvector rows: expected {3*self.ntot}, got {dvector.shape[1]}")
            self.assertEqual(dvector.shape[2], 3*self.ntot,
                           f"dvector cols: expected {3*self.ntot}, got {dvector.shape[2]}")

            print(f"✓ dscalar shape: {dscalar.shape}")
            print(f"✓ dvector shape: {dvector.shape}")
            print("✓ Matrix dimensions CORRECT")

        except Exception as e:
            self.fail(f"Failed to create solver or access matrices: {e}")

    def test_constant_derivative_is_zero(self):
        """Test 2: Derivative of constant should be zero"""
        print("\n" + "="*70)
        print("TEST 2: Derivative of Constant = 0")
        print("="*70)

        from adflow import ADFLOW

        options = {
            'gridfile': '../input_files/naca64A010_euler-L2.cgns',
            'equationmode': 'time spectral',
            'useTorusTimeSpectral': True,
            'nTimeIntervalsSpectral1': self.n1,
            'nTimeIntervalsSpectral2': self.n2,
            'omegaFourier1': self.omega1,
            'omegaFourier2': self.omega2,
        }

        solver = ADFLOW(options=options, comm=MPI.COMM_WORLD)
        dscalar = solver.adflow.inputtimespectral.dscalar

        # Constant function: f(theta1, theta2) = 1.0
        f_const = np.ones(self.ntot)

        # Compute derivative using dscalar
        # d/dt = sum over all j: dscalar(i,j) * f(j)
        section = 1  # First section (Fortran 1-based, Python uses 0)
        df_dt = np.zeros(self.ntot)
        for i in range(self.ntot):
            for j in range(self.ntot):
                df_dt[i] += dscalar[section-1, i, j] * f_const[j]

        # Check if derivative is zero (within numerical tolerance)
        max_error = np.max(np.abs(df_dt))
        print(f"Max |d(const)/dt| = {max_error:.2e}")

        self.assertLess(max_error, 1e-12,
                       f"Derivative of constant should be ~0, got max={max_error}")
        print("✓ Derivative of constant is zero (spectral accuracy)")

    def test_sine_cosine_derivatives(self):
        """Test 3: Derivative of sine/cosine matches analytical"""
        print("\n" + "="*70)
        print("TEST 3: Spectral Derivatives of Sine/Cosine")
        print("="*70)

        from adflow import ADFLOW

        options = {
            'gridfile': '../input_files/naca64A010_euler-L2.cgns',
            'equationmode': 'time spectral',
            'useTorusTimeSpectral': True,
            'nTimeIntervalsSpectral1': self.n1,
            'nTimeIntervalsSpectral2': self.n2,
            'omegaFourier1': self.omega1,
            'omegaFourier2': self.omega2,
        }

        solver = ADFLOW(options=options, comm=MPI.COMM_WORLD)
        dscalar = solver.adflow.inputtimespectral.dscalar

        # Test function: f(theta1, theta2) = sin(theta1) + cos(theta2)
        # Analytical: df/dt = omega1*cos(theta1) - omega2*sin(theta2)

        f = np.zeros(self.ntot)
        df_analytical = np.zeros(self.ntot)

        idx = 0
        for i2 in range(self.n2):
            for i1 in range(self.n1):
                f[idx] = np.sin(self.theta1[i1]) + np.cos(self.theta2[i2])
                df_analytical[idx] = (self.omega1 * np.cos(self.theta1[i1]) -
                                     self.omega2 * np.sin(self.theta2[i2]))
                idx += 1

        # Compute derivative using dscalar
        section = 1
        df_numerical = np.zeros(self.ntot)
        for i in range(self.ntot):
            for j in range(self.ntot):
                df_numerical[i] += dscalar[section-1, i, j] * f[j]

        # Compare
        error = np.abs(df_numerical - df_analytical)
        max_error = np.max(error)
        rms_error = np.sqrt(np.mean(error**2))

        print(f"Max error: {max_error:.2e}")
        print(f"RMS error: {rms_error:.2e}")

        # For spectral methods, error should be near machine precision for smooth functions
        self.assertLess(max_error, 1e-10,
                       f"Spectral derivative error too large: {max_error}")
        print("✓ Spectral derivatives accurate to machine precision")

    def test_sparsity_pattern(self):
        """Test 4: Verify sparsity pattern (cross, not dense)"""
        print("\n" + "="*70)
        print("TEST 4: Sparsity Pattern (Cross Pattern)")
        print("="*70)

        from adflow import ADFLOW

        options = {
            'gridfile': '../input_files/naca64A010_euler-L2.cgns',
            'equationmode': 'time spectral',
            'useTorusTimeSpectral': True,
            'nTimeIntervalsSpectral1': self.n1,
            'nTimeIntervalsSpectral2': self.n2,
            'omegaFourier1': self.omega1,
            'omegaFourier2': self.omega2,
        }

        solver = ADFLOW(options=options, comm=MPI.COMM_WORLD)
        dscalar = solver.adflow.inputtimespectral.dscalar

        section = 1
        tolerance = 1e-14

        # Count non-zeros
        nnz = 0
        for i in range(self.ntot):
            for j in range(self.ntot):
                if abs(dscalar[section-1, i, j]) > tolerance:
                    nnz += 1

        # Expected: Each row has n1 + n2 - 1 non-zeros (including diagonal=0)
        # But diagonal is 0, so n1 + n2 - 2 actual non-zeros
        # Total: ntot * (n1 + n2 - 2) for off-diagonals + ntot * 0 for diagonal
        expected_nnz_per_row = self.n1 + self.n2 - 1  # Counting diagonal as "present"
        expected_nnz_total = self.ntot * (self.n1 + self.n2 - 2)  # Actual non-zeros

        print(f"Non-zeros found: {nnz}")
        print(f"Expected (cross pattern): {expected_nnz_total}")
        print(f"If dense (all-to-all): {self.ntot * self.ntot} ({self.ntot**2})")
        print(f"Sparsity: {100*(1 - nnz/(self.ntot**2)):.1f}% sparse")

        # Verify sparsity
        self.assertLess(nnz, 0.6 * self.ntot**2,
                       "Matrix is too dense! Should be sparse cross pattern")
        self.assertGreater(nnz, expected_nnz_total * 0.9,
                          "Matrix is too sparse! Missing couplings")

        # Check specific instance (center of grid)
        i1_center, i2_center = self.n1 // 2, self.n2 // 2
        idx_center = i2_center * self.n1 + i1_center

        print(f"\nCenter instance ({i1_center},{i2_center}) = {idx_center}:")
        print("Couples to:")

        for j in range(self.ntot):
            if abs(dscalar[section-1, idx_center, j]) > tolerance:
                j1 = j % self.n1
                j2 = j // self.n1
                same_theta1 = (j1 == i1_center)
                same_theta2 = (j2 == i2_center)

                if same_theta1 or same_theta2:
                    print(f"  Instance {j} ({j1},{j2}) - " +
                          f"{'same θ₁' if same_theta1 else 'same θ₂'}")
                else:
                    self.fail(f"Unexpected coupling to ({j1},{j2}) - neither same θ₁ nor θ₂!")

        print("✓ Sparsity pattern is correct (cross pattern)")

    def test_dvector_construction(self):
        """Test 5: Verify dvector is block-diagonal with dscalar"""
        print("\n" + "="*70)
        print("TEST 5: dvector Construction from dscalar")
        print("="*70)

        from adflow import ADFLOW

        options = {
            'gridfile': '../input_files/naca64A010_euler-L2.cgns',
            'equationmode': 'time spectral',
            'useTorusTimeSpectral': True,
            'nTimeIntervalsSpectral1': self.n1,
            'nTimeIntervalsSpectral2': self.n2,
            'omegaFourier1': self.omega1,
            'omegaFourier2': self.omega2,
        }

        solver = ADFLOW(options=options, comm=MPI.COMM_WORLD)
        dscalar = solver.adflow.inputtimespectral.dscalar
        dvector = solver.adflow.inputtimespectral.dvector

        section = 1
        max_error = 0.0

        # Check block-diagonal structure
        # dvector(3*i+k, 3*j+k) should equal dscalar(i,j) for k=0,1,2
        for i in range(self.ntot):
            for j in range(self.ntot):
                for k in range(3):
                    expected = dscalar[section-1, i, j]
                    actual = dvector[section-1, 3*i+k, 3*j+k]
                    error = abs(actual - expected)
                    max_error = max(max_error, error)

                    if error > 1e-12:
                        print(f"ERROR at ({i},{j},comp={k}): " +
                              f"expected {expected:.6e}, got {actual:.6e}")

        print(f"Max error in dvector blocks: {max_error:.2e}")
        self.assertLess(max_error, 1e-12,
                       "dvector blocks don't match dscalar!")

        # Check off-diagonal components are zero
        off_diag_max = 0.0
        for i in range(3*self.ntot):
            for j in range(3*self.ntot):
                if i % 3 != j % 3:  # Different velocity components
                    val = abs(dvector[section-1, i, j])
                    off_diag_max = max(off_diag_max, val)

        print(f"Max off-diagonal coupling (should be 0): {off_diag_max:.2e}")
        self.assertLess(off_diag_max, 1e-12,
                       "dvector should be block-diagonal!")

        print("✓ dvector correctly constructed as block-diagonal from dscalar")


if __name__ == "__main__":
    print("\n" + "="*70)
    print("TORUS TIME SPECTRAL: Differentiation Matrix Unit Tests")
    print("="*70)

    # Run tests
    suite = unittest.TestLoader().loadTestsFromTestCase(TestTorusDifferentiationMatrices)
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    # Summary
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)
    print(f"Tests run: {result.testsRun}")
    print(f"Successes: {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")

    if result.wasSuccessful():
        print("\n✓ ALL TESTS PASSED!")
    else:
        print("\n✗ SOME TESTS FAILED")
        exit(1)
