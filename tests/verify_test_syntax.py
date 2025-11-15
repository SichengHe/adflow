#!/usr/bin/env python
"""
Verify Test Syntax - Check all torus test files for syntax errors

This script performs basic syntax checking on all test files without
requiring ADflow to be compiled. Useful for CI and pre-compilation checks.

Usage:
    python verify_test_syntax.py
"""

import ast
import sys
import os
from pathlib import Path


def check_python_syntax(filepath):
    """Check if a Python file has valid syntax"""
    try:
        with open(filepath, 'r') as f:
            source = f.read()
        ast.parse(source)
        return True, None
    except SyntaxError as e:
        return False, f"Line {e.lineno}: {e.msg}"
    except Exception as e:
        return False, str(e)


def main():
    """Check all torus test files"""
    print("=" * 80)
    print("TORUS TIME SPECTRAL - TEST SYNTAX VERIFICATION")
    print("=" * 80)
    print()

    # Get test directory
    test_dir = Path(__file__).parent

    # List of test files to check
    test_files = [
        "run_all_torus_tests.py",
        "unit_tests/test_torus_differentiation_matrices.py",
        "unit_tests/test_torus_jacobian_pattern.py",
        "unit_tests/test_torus_grid_velocity.py",
        "unit_tests/test_torus_mesh_deformation.py",
        "verification/test_torus_unsteady_interpolation.py",
        "reg_tests/test_torus_time_spectral_naca64A010.py",
    ]

    all_passed = True
    results = []

    for test_file in test_files:
        filepath = test_dir / test_file

        # Check if file exists
        if not filepath.exists():
            all_passed = False
            results.append((test_file, False, "File not found"))
            continue

        # Check syntax
        passed, error = check_python_syntax(filepath)
        results.append((test_file, passed, error))

        if not passed:
            all_passed = False

    # Print results
    print("SYNTAX CHECK RESULTS:")
    print("-" * 80)

    for test_file, passed, error in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status} - {test_file}")

        if not passed and error:
            print(f"       Error: {error}")

    print("-" * 80)
    print()

    # Summary
    passed_count = sum(1 for _, passed, _ in results if passed)
    total_count = len(results)

    print(f"Total: {total_count} files")
    print(f"Passed: {passed_count}/{total_count}")
    print(f"Failed: {total_count - passed_count}")
    print()

    if all_passed:
        print("🎉 " * 20)
        print("ALL TEST FILES HAVE VALID SYNTAX!")
        print("🎉 " * 20)
        print()
        print("Next steps:")
        print("  1. Set up Fortran compiler (gfortran or ifort)")
        print("  2. Configure build: cp config/defaults/config.*.mk config/config.mk")
        print("  3. Compile ADflow: make clean && make")
        print("  4. Run tests: cd tests && python run_all_torus_tests.py --quick")
        return 0
    else:
        print("⚠ " * 20)
        print("SOME TEST FILES HAVE SYNTAX ERRORS!")
        print("⚠ " * 20)
        print()
        print("Please fix the errors above before proceeding.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
