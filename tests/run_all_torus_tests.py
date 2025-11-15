"""
Master Test Runner for Torus Time Spectral Implementation

Runs all unit tests and verification tests in order:
1. Differentiation matrices
2. Jacobian pattern
3. Grid velocity
4. Mesh deformation
5. Unsteady interpolation (verification)
6. Full regression test

Usage:
    python run_all_torus_tests.py [--quick] [--verbose]

Options:
    --quick     Skip expensive verification tests
    --verbose   Print detailed output
"""

import sys
import os
import unittest
import argparse
from datetime import datetime


class TorusTestSuite:
    """Manages the complete test suite for torus time spectral"""

    def __init__(self, verbose=False, quick=False):
        self.verbose = verbose
        self.quick = quick
        self.results = {}

    def run_test_file(self, test_file, description):
        """Run a single test file"""
        print("\n" + "="*80)
        print(f"Running: {description}")
        print(f"File: {test_file}")
        print("="*80)

        # Import and run
        try:
            # Add test directory to path
            test_dir = os.path.dirname(os.path.abspath(test_file))
            sys.path.insert(0, test_dir)

            # Load and run
            loader = unittest.TestLoader()
            suite = loader.discover(test_dir, pattern=os.path.basename(test_file))

            verbosity = 2 if self.verbose else 1
            runner = unittest.TextTestRunner(verbosity=verbosity)
            result = runner.run(suite)

            # Record result
            self.results[description] = {
                'passed': result.wasSuccessful(),
                'tests': result.testsRun,
                'failures': len(result.failures),
                'errors': len(result.errors),
            }

            return result.wasSuccessful()

        except Exception as e:
            print(f"\n✗ ERROR running {test_file}: {e}")
            self.results[description] = {
                'passed': False,
                'tests': 0,
                'failures': 0,
                'errors': 1,
                'error_msg': str(e)
            }
            return False

    def run_all(self):
        """Run all tests in sequence"""
        start_time = datetime.now()

        print("\n" + "="*80)
        print("TORUS TIME SPECTRAL - COMPLETE TEST SUITE")
        print("="*80)
        print(f"Start time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Mode: {'QUICK' if self.quick else 'FULL'}")
        print()

        # Test sequence
        tests = [
            ("unit_tests/test_torus_differentiation_matrices.py",
             "Unit Test 1: Differentiation Matrices"),

            ("unit_tests/test_torus_jacobian_pattern.py",
             "Unit Test 2: Jacobian Sparsity Pattern"),

            ("unit_tests/test_torus_grid_velocity.py",
             "Unit Test 3: Grid Velocity Computation"),

            ("unit_tests/test_torus_mesh_deformation.py",
             "Unit Test 4: Mesh Deformation"),
        ]

        # Add expensive tests if not quick mode
        if not self.quick:
            tests.append(
                ("verification/test_torus_unsteady_interpolation.py",
                 "Verification Test: Unsteady Interpolation")
            )

            tests.append(
                ("reg_tests/test_torus_time_spectral_naca64A010.py",
                 "Regression Test: Full Solver")
            )

        # Run each test
        all_passed = True
        for test_file, description in tests:
            passed = self.run_test_file(test_file, description)
            all_passed = all_passed and passed

            if not passed and not self.quick:
                print("\n⚠ Test failed! Continue anyway? (y/n): ", end='')
                response = input().lower()
                if response != 'y':
                    print("Aborting test suite.")
                    break

        # Summary
        end_time = datetime.now()
        duration = end_time - start_time

        print("\n" + "="*80)
        print("TEST SUITE SUMMARY")
        print("="*80)
        print(f"End time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Duration: {duration}")
        print()

        total_tests = 0
        total_failures = 0
        total_errors = 0
        passed_count = 0

        for description, result in self.results.items():
            status = "✓ PASS" if result['passed'] else "✗ FAIL"
            print(f"{status} - {description}")
            print(f"       Tests: {result['tests']}, " +
                  f"Failures: {result['failures']}, " +
                  f"Errors: {result['errors']}")

            total_tests += result['tests']
            total_failures += result['failures']
            total_errors += result['errors']
            if result['passed']:
                passed_count += 1

        print()
        print("-" * 80)
        print(f"Total: {len(self.results)} test suites, {total_tests} tests")
        print(f"Passed: {passed_count}/{len(self.results)} suites")
        print(f"Failed: {total_failures} tests")
        print(f"Errors: {total_errors} tests")
        print("-" * 80)

        if all_passed:
            print("\n" + "🎉 " * 20)
            print("ALL TESTS PASSED! Torus time spectral implementation verified!")
            print("🎉 " * 20)
            return 0
        else:
            print("\n" + "⚠ " * 20)
            print("SOME TESTS FAILED! Review output above for details.")
            print("⚠ " * 20)
            return 1


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="Run torus time spectral test suite"
    )
    parser.add_argument('--quick', action='store_true',
                       help='Skip expensive verification tests')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose output')

    args = parser.parse_args()

    # Create and run test suite
    suite = TorusTestSuite(verbose=args.verbose, quick=args.quick)
    exit_code = suite.run_all()

    sys.exit(exit_code)


if __name__ == "__main__":
    main()
