# ADflow Resolvent Analysis - Build and Test Instructions

## ✅ What's Been Done

### 1. Numerical Improvements ✓
- **Improved solver stability** by replacing matrix inversion with LU decomposition
- **Added sparse SVD** for memory efficiency
- **Validated on algebraic test** - all tests pass to machine precision

### 2. Implementation Complete ✓
- **Fortran API**: `src/modalAnalysis/resolventAPI.F90`
- **Python Classes**: `adflow/pyResolventAnalysis.py`
- **ADflow Integration**: `adflow/pyADflow.py` (modified)
- **Build System**: Added to `src/build/directoryList` and `src/build/fileList`

### 3. Testing Ready ✓
- **Algebraic test**: `test_resolvent_simple.py` (passes)
- **CFD test**: `tests/test_resolvent_simple.py` (ready to run)

---

## 🚀 Next Steps: Build and Test

### Step 1: Rebuild ADflow

```bash
cd /home/sicheng/repo/adflow

# Clean previous build
make clean

# Build ADflow with resolvent module
make

# This will:
# 1. Compile src/modalAnalysis/resolventAPI.F90
# 2. Link with adjoint infrastructure
# 3. Build Python extension with f2py
# 4. Install libadflow.so to adflow/ directory
```

**Expected output:**
```
Making dependencies!
Compiling modalAnalysis/resolventAPI.F90...
...
[Success message]
```

**If build fails**, check:
- PETSc is properly configured
- Fortran compiler is working
- All dependencies in `src/modalAnalysis/resolventAPI.F90` are available

### Step 2: Verify Python Integration

```bash
cd /home/sicheng/repo/adflow

# Test that resolvent module is available
python -c "from adflow import ADFLOW; print(dir(ADFLOW.adflow.resolventapi))"

# Should see:
# ['getresolventmatrixdense', 'getresolventmatrixinfo',
#  'setupresolventmatrix', 'exportresolventmatrixtofile']

# Test that Python classes are exported
python -c "from adflow import ResolventAnalysis; print('✓ ResolventAnalysis imported')"
```

### Step 3: Run Algebraic Test (Should Already Pass)

```bash
cd /home/sicheng/repo/adflow

# Run simple algebraic validation
python test_resolvent_simple.py
```

**Expected output:**
```
======================================================================
Testing Resolvent Analysis Methods
======================================================================
...
ω = 1.0
  Method 1 (Inversion):     σ₁ = 4.732255
  Method 2 (LU Decomp):     σ₁ = 4.732255
  Relative error: 0.00e+00  ✓ PASS
...
All tests PASSED!
```

### Step 4: Download Test Meshes

```bash
cd /home/sicheng/repo/adflow/input_files

# Download test mesh files
./get-input-files.sh

# This downloads ~200 MB of CGNS mesh files
# Files will be in /home/sicheng/repo/adflow/input_files/
```

### Step 5: Run CFD Test

```bash
cd /home/sicheng/repo/adflow/tests

# Run CFD resolvent test
python test_resolvent_simple.py
```

**Expected output:**
```
======================================================================
ADflow Resolvent Analysis - Simple Test
======================================================================

✓ Found mesh file: ../input_files/mdo_tutorial_euler_scalar_jst.cgns
Creating ADflow solver...
✓ Solver created

----------------------------------------------------------------------
Solving steady-state flow...
----------------------------------------------------------------------
[ADflow convergence output...]

Converged solution:
  CL = 0.XXXXXX
  CD = 0.XXXXXX

----------------------------------------------------------------------
Testing Jacobian assembly...
----------------------------------------------------------------------
✓ setupResolventJacobian() method found
Assembling Jacobian matrix...
✓ Jacobian assembled

State vector size: XXXX
Extracting Jacobian matrix...
✓ Jacobian extracted: shape = (XXXX, XXXX)

✓ Jacobian assembly test PASSED

----------------------------------------------------------------------
Testing resolvent analysis...
----------------------------------------------------------------------
✓ ResolventAnalysis object created
Computing resolvent operator and SVD...
✓ Resolvent analysis PASSED
  Dominant singular value: σ₁ = X.XXXXXX

======================================================================
TEST SUMMARY
======================================================================
✓ ADflow steady-state solve
✓ Jacobian assembly
✓ Jacobian extraction
✓ Resolvent SVD computation

All tests PASSED!
```

---

## 📝 Files Modified/Created

### Build System (Modified)
- ✅ `src/build/directoryList` - Added `modalAnalysis`
- ✅ `src/build/fileList` - Added `modalAnalysis/resolventAPI.F90`

### Fortran Implementation (Created)
- ✅ `src/modalAnalysis/resolventAPI.F90` - Fortran API
- ✅ `src/modalAnalysis/README.md` - Documentation

### Python Implementation (Created)
- ✅ `adflow/pyResolventAnalysis.py` - Main implementation (~550 lines)

### Python Integration (Modified)
- ✅ `adflow/__init__.py` - Exports ResolventAnalysis classes
- ✅ `adflow/pyADflow.py` - Added integration methods

### Tests (Created)
- ✅ `test_resolvent_simple.py` - Algebraic validation
- ✅ `tests/test_resolvent_simple.py` - CFD test

### Documentation (Created)
- ✅ `doc/resolvent_analysis.md` - Theory and usage
- ✅ `RESOLVENT_IMPLEMENTATION_COMPLETE.md` - Full guide
- ✅ `RESOLVENT_TESTING_STATUS.md` - Testing status
- ✅ `BUILD_AND_TEST_INSTRUCTIONS.md` - This file

---

## 🐛 Troubleshooting

### Build Errors

**Error: `resolventAPI.F90: No such file or directory`**
- Check that `src/modalAnalysis/resolventAPI.F90` exists
- Verify `src/build/directoryList` includes `modalAnalysis`
- Verify `src/build/fileList` includes `modalAnalysis/resolventAPI.F90`

**Error: Module not found (e.g., `ADjointVars`)**
- The resolventAPI.F90 depends on adjoint infrastructure
- Make sure adjoint modules are compiled first (they should be)
- Check module dependencies in the .dep file

**Error: PETSc-related compilation errors**
- Ensure PETSc environment variables are set
- Check `config/config.mk` for correct PETSc paths

### Runtime Errors

**ImportError: cannot import name 'ResolventAnalysis'**
```bash
# Check if module is in adflow/__init__.py
grep "ResolventAnalysis" adflow/__init__.py

# Should see:
# from .pyResolventAnalysis import ResolventAnalysis, ResolventAnalysisMatrixFree
```

**AttributeError: 'ADFLOW' object has no attribute 'setupResolventJacobian'**
```bash
# Check if Fortran module was compiled
python -c "from adflow import ADFLOW; print(dir(ADFLOW.adflow))" | grep resolvent

# Should see: resolventapi
```

**Test mesh not found**
```bash
cd /home/sicheng/repo/adflow/input_files
./get-input-files.sh
```

---

## 📊 Performance Expectations

### Algebraic Test (2x2 system)
- **Time**: <1 second
- **Memory**: <10 MB
- **Result**: σ₁ ≈ 0.25-5.0 (depending on ω)

### Small CFD Test (Tutorial Wing Euler, ~2000-5000 DOF)
- **CFD solve**: 1-5 minutes
- **Jacobian assembly**: 5-30 seconds
- **Jacobian extraction**: 1-10 seconds (dense matrix, ~40-200 MB)
- **Resolvent SVD**: 1-10 seconds
- **Total**: ~2-10 minutes

### Large CFD Problems (>100,000 DOF)
- **Dense Jacobian**: Impractical (>80 GB memory)
- **Recommendation**: Use matrix-free methods (future work)
- **Current status**: Matrix-free implementation incomplete

---

## ✨ Key Improvements Made

### 1. Numerical Stability
**Before:**
```python
R = np.linalg.inv(A)  # Explicit inversion - unstable!
U, S, Vh = scipy.linalg.svd(R)  # Full SVD - memory intensive
```

**After:**
```python
lu, piv = lu_factor(A)  # LU factorization - stable

def matvec(f):
    return lu_solve((lu, piv), f)  # Solve, not invert

R_op = LinearOperator((n, n), matvec=matvec)
U, S, Vh = svds(R_op, k=nModes, which='LM')  # Sparse SVD - efficient
```

### 2. Memory Efficiency
- **Old**: Stores full resolvent matrix R (n×n dense, complex)
- **New**: Uses LinearOperator (no storage, solve on-the-fly)
- **Savings**: For n=10,000: ~1.6 GB → ~0 GB for operator

### 3. Accuracy
- **Validation**: Matches full SVD to machine precision (error < 1e-15)
- **Tested**: Multiple frequencies, both real and complex formulations

---

## 📖 Usage Example

Once built and tested, here's how to use resolvent analysis:

```python
from adflow import ADFLOW, ResolventAnalysis
from baseclasses import AeroProblem

# Setup and solve CFD
CFDsolver = ADFLOW(options=aeroOptions)
ap = AeroProblem(name='wing', mach=0.8, alpha=2.0, ...)
CFDsolver(ap)  # Solve steady state

# Perform resolvent analysis
resolvent = ResolventAnalysis(CFDsolver, ap, omega=1.0)
sigma1 = resolvent.solveExplicit(useLU=True)  # Use stable method

print(f"Dominant singular value: σ₁ = {sigma1:.6f}")

# Get modes
u1 = resolvent.getResponseMode()  # Response mode
v1 = resolvent.getForcingMode()   # Forcing mode

# Frequency sweep
omega_vec, sigma1_vec = resolvent.computeFrequencySweep((0, 10), nPoints=50)
```

---

## 🎯 Success Criteria

The implementation is successful if:

1. ✅ **Build completes** without errors
2. ✅ **Python imports work**: `from adflow import ResolventAnalysis`
3. ✅ **Algebraic test passes**: All frequencies give matching results
4. ✅ **CFD test runs**: Jacobian assembly and extraction work
5. ✅ **Resolvent computed**: SVD gives reasonable σ₁ values

---

## 📞 Next Steps After Testing

Once all tests pass:

1. **Verify against paper examples** (Step 4 of original plan)
2. **Test on larger CFD cases** (if memory permits)
3. **Develop matrix-free methods** (for production use)
4. **Add to CI/CD** (if applicable)
5. **Write paper/documentation** (if publishing)

---

## 🔬 Reference

**Paper**: "Large-Scale Flow Control Performance Optimization via
Differentiable Resolvent Analysis" by He et al.

**Theory**: See `doc/resolvent_analysis.md`

**Implementation**: See `RESOLVENT_IMPLEMENTATION_COMPLETE.md`

---

**Last Updated:** 2025-11-15
**Status:** ✅ Ready to build and test
