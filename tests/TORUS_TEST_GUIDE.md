# Torus Time Spectral - Test Guide

**Complete guide for testing the torus time spectral implementation**

---

## Prerequisites

### 1. Compiler Setup

ADflow requires a Fortran compiler. Choose one of:

**Option A: GNU Fortran (gfortran)**
```bash
# macOS
brew install gcc

# Linux (Ubuntu/Debian)
sudo apt-get install gfortran

# Linux (RHEL/CentOS)
sudo yum install gcc-gfortran
```

**Option B: Intel Fortran (ifort/ifx)**
```bash
# Install Intel oneAPI HPC Toolkit
# https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html
```

### 2. Configure Build

Copy appropriate config file:
```bash
cd /Users/sichenghe/Desktop/adflow

# For gfortran
cp config/defaults/config.LINUX_GFORTRAN.mk config/config.mk

# For Intel Fortran
cp config/defaults/config.LINUX_INTEL.mk config/config.mk
```

Edit `config/config.mk` to set:
- CGNS library path
- Compiler flags
- MPI settings

### 3. Compile ADflow

```bash
cd /Users/sichenghe/Desktop/adflow
make clean
make
```

This will compile all Fortran sources including the new torus time spectral code:
- `src/modules/inputParam.F90` (torus parameters)
- `src/inputParam/inputParamRoutines.F90` (initialization)
- `src/initFlow/initializeFlow.F90` (torusTimeSpectralMatrices)
- `src/preprocessing/preprocessingAPI.F90` (calling the torus setup)

---

## Test Suite Organization

```
tests/
├── run_all_torus_tests.py          # Master test runner
├── unit_tests/                      # Fast, focused unit tests
│   ├── test_torus_differentiation_matrices.py
│   ├── test_torus_jacobian_pattern.py
│   ├── test_torus_grid_velocity.py
│   └── test_torus_mesh_deformation.py
├── verification/                    # Verification against theory/unsteady
│   └── test_torus_unsteady_interpolation.py
└── reg_tests/                       # Full solver regression tests
    └── test_torus_time_spectral_naca64A010.py
```

---

## Running Tests

### Quick Test (Unit Tests Only)

Run only the fast unit tests (~1-2 minutes):

```bash
cd /Users/sichenghe/Desktop/adflow/tests
python run_all_torus_tests.py --quick
```

**What it tests:**
1. ✓ Differentiation matrices (`dscalar`, `dscalar_torus`, `dscalar2_torus`)
2. ✓ Momentum differentiation (`dvector` block-diagonal structure)
3. ✓ Jacobian sparsity pattern (cross coupling, not dense)
4. ✓ Grid velocity computation
5. ✓ Mesh deformation on torus grid

### Full Test Suite

Run all tests including expensive verification (~30-60 minutes):

```bash
cd /Users/sichenghe/Desktop/adflow/tests
python run_all_torus_tests.py
```

**Additional tests:**
6. ✓ Unsteady interpolation verification (gold standard)
7. ✓ Full regression test (NACA 64A010, 2-frequency pitching)

### Verbose Output

For detailed diagnostic output:

```bash
python run_all_torus_tests.py --verbose
```

---

## Individual Test Details

### Test 1: Differentiation Matrices

**File:** `unit_tests/test_torus_differentiation_matrices.py`

**What it verifies:**
- Matrix dimensions: `dscalar(nSections, n1*n2, n1*n2)`
- Constant derivative is zero: `d/dθ(1) = 0`
- Sine derivative accuracy: `d/dθ(sin(θ)) ≈ cos(θ)` within 1e-10
- Cosine derivative accuracy: `d/dθ(cos(θ)) ≈ -sin(θ)` within 1e-10
- Sparsity pattern: Each row has n₁+n₂-1 non-zeros (cross pattern)
- dvector construction: Block-diagonal structure from dscalar

**Run individually:**
```bash
cd tests
python -m unittest unit_tests.test_torus_differentiation_matrices
```

**Key assertion:**
```python
# Each instance couples to n1 + n2 - 1 others (including itself)
expected_nnz = n1 + n2 - 1  # For 3×3: 3+3-1 = 5
actual_nnz = np.count_nonzero(np.abs(dscalar[0, i, :]) > 1e-14)
self.assertEqual(actual_nnz, expected_nnz)
```

### Test 2: Jacobian Sparsity Pattern

**File:** `unit_tests/test_torus_jacobian_pattern.py`

**What it verifies:**
- Jacobian has correct block structure
- Each spectral instance couples to n₁+n₂-1 neighbors
- Coupling only when `j1==i1` (same θ₁) OR `j2==i2` (same θ₂)
- Block size matches number of DOFs per instance

**Run individually:**
```bash
python -m unittest unit_tests.test_torus_jacobian_pattern
```

**Expected pattern (3×3 grid):**
```
Instance (i1,i2) couples to:
(1,1) → (1,1), (2,1), (3,1), (1,2), (1,3)  [5 neighbors]
(2,2) → (2,1), (2,2), (2,3), (1,2), (3,2)  [5 neighbors]
etc.
```

### Test 3: Grid Velocity Computation

**File:** `unit_tests/test_torus_grid_velocity.py`

**What it verifies:**
- Grid velocity formula: `v_grid = Σⱼ dscalar(i,j) · x(j)`
- Velocity magnitude is reasonable
- Geometric Conservation Law (GCL) satisfied
- Torus grid deforms correctly with prescribed motion

**Run individually:**
```bash
python -m unittest unit_tests.test_torus_grid_velocity
```

**Prescribed motion example:**
```python
# 2-frequency pitch: α(θ₁,θ₂) = α₀ + A₁·sin(θ₁) + A₂·sin(θ₂)
alpha_0 = 0.0  # deg
A1 = 1.0       # deg
A2 = 0.5       # deg
```

### Test 4: Mesh Deformation

**File:** `unit_tests/test_torus_mesh_deformation.py`

**What it verifies:**
- 2D torus indexing: `(i1, i2) → sps = (i2-1)*n1 + i1`
- Reverse mapping: `sps → (i1, i2)` is correct
- Prescribed motion applied to all instances
- Deformation pattern follows 2-frequency forcing
- Mesh quality maintained (no negative volumes)

**Run individually:**
```bash
python -m unittest unit_tests.test_torus_mesh_deformation
```

**Indexing verification:**
```python
# Forward: (i1, i2) → sps
sps = (i2 - 1) * n1 + i1

# Reverse: sps → (i1, i2)
i1_recovered = (sps - 1) % n1 + 1
i2_recovered = (sps - 1) // n1 + 1
```

### Test 5: Unsteady Interpolation (GOLD STANDARD)

**File:** `verification/test_torus_unsteady_interpolation.py`

**What it verifies:**
- Converged unsteady solution can be reconstructed on torus grid
- Residual after interpolation is small: `||R|| < 1e-6`
- Torus spectral method correctly represents quasi-periodic flow

**Procedure:**
1. Run unsteady simulation to convergence
2. Extract solution at time instances corresponding to torus grid points
3. Perform 2D FFT to get Fourier coefficients
4. Reconstruct solution on torus grid using inverse FFT
5. Load reconstructed solution into ADflow with torus setup
6. Compute residual - should be ≈ machine precision

**Run individually:**
```bash
python -m unittest verification.test_torus_unsteady_interpolation
```

**Mathematical basis:**
```python
# Fourier decomposition
u_hat[k1, k2] = (1/(n1*n2)) * Σᵢ₁Σᵢ₂ u[i1, i2] * exp(-2πi(k1*i1/n1 + k2*i2/n2))

# Reconstruction (should satisfy governing equations)
u[i1, i2] = Σₖ₁Σₖ₂ u_hat[k1, k2] * exp(2πi(k1*i1/n1 + k2*i2/n2))
```

### Test 6: Full Regression Test

**File:** `reg_tests/test_torus_time_spectral_naca64A010.py`

**What it verifies:**
- Complete torus time spectral solve converges
- Solution quality matches expected aerodynamic behavior
- Force/moment coefficients are physical
- Multi-frequency response is captured

**Test case:**
- Airfoil: NACA 64A010
- Grid: 3×3 torus (9 spectral instances)
- Frequency 1: ω₁ = 100 rad/s
- Frequency 2: ω₂ = 100√2 rad/s (incommensurate)
- Motion: α(θ₁,θ₂) = 0° + 1°·sin(θ₁) + 0.5°·sin(θ₂)
- Flow: M∞ = 0.5, Re = 1e6

**Run individually:**
```bash
python -m unittest reg_tests.test_torus_time_spectral_naca64A010
```

---

## Expected Results

### Unit Tests (Quick Mode)

```
================================================================================
Running: Unit Test 1: Differentiation Matrices
File: unit_tests/test_torus_differentiation_matrices.py
================================================================================
.....
----------------------------------------------------------------------
Ran 5 tests in 0.123s

OK

✓ PASS - Unit Test 1: Differentiation Matrices
       Tests: 5, Failures: 0, Errors: 0

[Similar for tests 2-4...]

================================================================================
TEST SUITE SUMMARY
================================================================================
Total: 4 test suites, 18 tests
Passed: 4/4 suites
Failed: 0 tests
Errors: 0 tests
--------------------------------------------------------------------------------

🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉
ALL TESTS PASSED! Torus time spectral implementation verified!
🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉 🎉
```

### Full Test Suite

Includes above plus:

```
================================================================================
Running: Verification Test: Unsteady Interpolation
File: verification/test_torus_unsteady_interpolation.py
================================================================================
NOTE: This test requires a converged unsteady reference solution.
      Skipping for now - requires manual setup.
...
----------------------------------------------------------------------
Ran 3 tests in 45.678s

OK

================================================================================
Running: Regression Test: Full Solver
File: reg_tests/test_torus_time_spectral_naca64A010.py
================================================================================
...
----------------------------------------------------------------------
Ran 3 tests in 156.234s

OK
```

---

## Troubleshooting

### Issue: ModuleNotFoundError: No module named 'adflow'

**Cause:** ADflow Python module not installed

**Fix:**
```bash
cd /Users/sichenghe/Desktop/adflow
pip install -e .
```

### Issue: ImportError: cannot import name 'adflow'

**Cause:** Fortran library not compiled or not in path

**Fix:**
```bash
make clean
make
# Check that libadflow.so exists in src/build/
```

### Issue: Test fails with "dscalar_torus not found"

**Cause:** ADflow compiled without torus modifications

**Fix:**
Ensure you compiled AFTER adding torus code to:
- `src/modules/inputParam.F90`
- `src/initFlow/initializeFlow.F90`

### Issue: Differentiation matrix row sums ≠ 0

**Cause:** Bug in torusTimeSpectralMatrices implementation

**Debug:**
```fortran
! In torusTimeSpectralMatrices, add:
do i = 1, ntot
    row_sum = sum(dscalar(1, i, :))
    if (abs(row_sum) > 1.0e-12_realType) then
        print *, 'WARNING: Row', i, 'sum =', row_sum
    end if
end do
```

### Issue: Sparsity pattern is dense instead of cross

**Cause:** Coupling logic error in matrix construction

**Check:**
```fortran
! Should only add non-zero when:
if (j1 == i1 .and. j2 /= i2) then
    ! Derivative w.r.t. theta2
else if (j2 == i2 .and. j1 /= i1) then
    ! Derivative w.r.t. theta1
end if
```

---

## Performance Benchmarks

### Expected Test Runtimes

| Test | Grid Size | Time | Description |
|------|-----------|------|-------------|
| Differentiation matrices | 3×3 | <1s | Matrix construction only |
| Jacobian pattern | 3×3 | <1s | Pattern analysis |
| Grid velocity | 3×3 | ~5s | Requires mesh |
| Mesh deformation | 3×3 | ~5s | Prescribed motion |
| **Quick suite total** | 3×3 | **~20s** | All unit tests |
| Unsteady interpolation | 3×3 | ~30-60min | Requires unsteady reference |
| Full regression | 3×3 | ~2-5min | Full ANK solve |
| **Full suite total** | 3×3 | **~1hr** | All tests |

### Memory Requirements

| Configuration | Memory per Process | Total (4 procs) |
|---------------|-------------------|-----------------|
| 3×3 torus, coarse mesh | ~500 MB | ~2 GB |
| 3×3 torus, medium mesh | ~2 GB | ~8 GB |
| 5×5 torus, coarse mesh | ~1.4 GB | ~5.6 GB |

---

## Continuous Integration

### Recommended CI Workflow

```yaml
# .github/workflows/torus_tests.yml
name: Torus Time Spectral Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2

      - name: Install dependencies
        run: |
          sudo apt-get install gfortran libopenmpi-dev
          pip install numpy mpi4py

      - name: Configure build
        run: cp config/defaults/config.LINUX_GFORTRAN.mk config/config.mk

      - name: Compile ADflow
        run: make

      - name: Run quick tests
        run: cd tests && python run_all_torus_tests.py --quick

      - name: Run full tests (nightly only)
        if: github.event_name == 'schedule'
        run: cd tests && python run_all_torus_tests.py
```

---

## Next Steps After Tests Pass

1. **Optimize for larger grids** (5×5, 7×7)
2. **Implement adjoint derivatives** for optimization
3. **Add phase optimization** for trajectory alignment
4. **Multi-frequency pitching validation** against experiments
5. **Rotor blade case** (3D, moving mesh)

---

## Related Documentation

- **[TORUS_IMPLEMENTATION_COMPLETE.md](../TORUS_IMPLEMENTATION_COMPLETE.md)** - Implementation overview
- **[TORUS_DATA_STORAGE.md](../TORUS_DATA_STORAGE.md)** - Data layout and indexing
- **[DVECTOR_THEORY.md](../DVECTOR_THEORY.md)** - Mathematical theory
- **[TORUS_SOLVER_ANALYSIS.md](../TORUS_SOLVER_ANALYSIS.md)** - ANK/NK solver details

---

**Status:** ⏳ Awaiting Fortran compiler setup and compilation

**Last Updated:** 2025-11-12

**Author:** Torus Time Spectral Implementation Team
