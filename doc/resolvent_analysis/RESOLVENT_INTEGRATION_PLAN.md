# Resolvent Analysis: What's Needed for Full Integration

## Your Question: "Do we need to touch the source code? Need to add SLEPc?"

### Short Answer

**For basic functionality (small-medium problems):**
- ✅ **No SLEPc needed** - can use NumPy/SciPy SVD
- ⚠️ **Yes, minimal Fortran modification** - need to expose the dRdW matrix
- ~50-100 lines of Fortran code

**For large-scale CFD applications:**
- ✅ **Yes, SLEPc recommended** - for iterative SVD
- ⚠️ **More Fortran integration** - matrix-free operators
- ~200-500 lines total

## What We Have Now (Pure Python Framework)

```
adflow/
├── pyResolventAnalysis.py       ← New Python class (✓ Complete)
├── __init__.py                   ← Modified to export classes
└── ...

examples/
└── resolvent_analysis_example.py ← Usage example (✓ Complete)

doc/
└── resolvent_analysis.md         ← Documentation (✓ Complete)
```

**Status:** Framework is complete, but cannot run on actual ADflow solutions yet.

## What's Missing: Three Levels of Integration

### Level 1: Minimum Viable (Small Problems)
**Goal:** Make it work for small test cases

**What's needed:**
1. Add one Fortran subroutine to export dRdW matrix:

```fortran
! In src/adjoint/adjointAPI.F90
subroutine getdRdWMatrix(dRdW_ptr) bind(C, name='getdrdwmatrix')
    ! Return PETSc Mat object for assembled dRdW
    use ADjointPETSc, only: dRdW
    type(c_ptr), intent(out) :: dRdW_ptr

    ! Assemble the matrix if not already done
    if (.not. dRdWAllocated) call setupdRdW()

    ! Return pointer to PETSc matrix
    dRdW_ptr = c_loc(dRdW)
end subroutine getdRdWMatrix
```

2. Add Python wrapper in `pyADflow.py`:

```python
def getJacobianMatrix(self):
    """Get the assembled dRdW (Jacobian) matrix"""
    # Call Fortran routine to get PETSc Mat
    mat_ptr = self.adflow.adjointapi.getdrdwmatrix()

    # Convert PETSc Mat to NumPy (for small problems)
    # or keep as PETSc Mat for large problems
    return mat_ptr
```

**Then resolvent analysis works:**
```python
# In pyResolventAnalysis.py
def solveExplicit(self):
    J = self.CFDsolver.getJacobianMatrix()  # ← This now works!
    A = 1j*self.omega*np.eye(n) - J
    R = np.linalg.inv(A)
    U, S, Vh = np.linalg.svd(R)
    self.sigma1 = S[0]
    return self.sigma1
```

**Pros:**
- Minimal code changes (~50 lines Fortran, ~20 lines Python)
- No new dependencies
- Works for problems up to ~10,000 DOF

**Cons:**
- Memory intensive (full matrix)
- Not scalable to large CFD meshes
- O(N³) complexity for SVD

### Level 2: Matrix-Free (Medium-Large Problems)
**Goal:** Use existing Jacobian-vector products

**What's needed:**
ADflow already has these! Just need to wrap them properly:

```python
# Already exists in ADflow!
def _jacobianVectorProduct(self, v):
    """Compute J·v using existing ADflow routine"""
    Jv = self.CFDsolver.computeJacobianVectorProductFwd(
        wDot=v,
        residualDeriv=True
    )
    return Jv
```

**What's missing:**
- Linear solver for (jω·I - J)·u = f using J·v products
- This can be done entirely in Python using scipy.sparse.linalg.gmres!

```python
from scipy.sparse.linalg import gmres, LinearOperator

def _solveLinearSystem(self, rhs):
    """Solve (jω·I - J)·u = rhs using GMRES"""
    n = self.stateSize
    omega = self.omega

    def matvec(v):
        Jv = self._jacobianVectorProduct(v)
        return 1j*omega*v - Jv

    A = LinearOperator((n,n), matvec=matvec, dtype=complex)
    u, info = gmres(A, rhs, tol=1e-8)
    return u
```

**Pros:**
- No full matrix needed
- Works with existing ADflow infrastructure
- Can handle larger problems (~100,000 DOF)

**Cons:**
- Still need full SVD for singular values
- SVD of implicit operator is tricky

### Level 3: Full Large-Scale (SLEPc Integration)
**Goal:** Iterative SVD for massive problems

**What's needed:**
1. SLEPc interface for SVD
2. Adjoint operator for two-sided Arnoldi

```python
# Requires SLEPc
from slepc4py import SLEPc

def solve_slepc(self):
    """Compute dominant singular values using SLEPc"""
    # Create SLEPc SVD object
    svd = SLEPc.SVD()
    svd.create()

    # Define operators
    def mult(mat, x, y):
        # y = R·x = (jω·I - J)^{-1}·x
        x_array = x.getArray()
        u = self._solveLinearSystem(x_array)
        y.setArray(u)

    # Set up and solve
    # ... (SLEPc boilerplate)
```

**Pros:**
- Handles millions of DOF
- Memory efficient
- Parallel scalability

**Cons:**
- Requires SLEPc installation
- More complex implementation
- ~500 lines of code

## Recommended Path Forward

### Phase 1: Get Something Working (1-2 days)
Focus on **Level 1** - just expose the Jacobian matrix:

1. Add `getdRdWMatrix()` to Fortran
2. Add `getJacobianMatrix()` to Python
3. Test on small 2D cases

**Deliverable:** Working resolvent analysis on small problems

### Phase 2: Make it Practical (1 week)
Implement **Level 2** - matrix-free operators:

1. Wrap existing Jacobian-vector products
2. Implement GMRES-based linear solver
3. Add randomized SVD for top modes

**Deliverable:** Can analyze realistic 2D problems

### Phase 3: Production Ready (2-4 weeks)
Add **Level 3** - SLEPc integration:

1. Install/configure SLEPc
2. Implement iterative SVD
3. Add parallel communication
4. Optimize performance

**Deliverable:** Production-ready for 3D CFD

## What About the Existing Code in `article_resolvent_opt/`?

The code in `article_resolvent_opt/code/` is **fully working** but for algebraic systems, not CFD:

```
article_resolvent_opt/code/
├── resolvent.py          ← Works! (for algebraic systems)
├── nonlinear_eqn.py      ← Works!
├── example.py            ← Works! Run it now
└── reactor_*.py          ← Works! Brusselator example
```

**Key difference:**
- Algebraic systems: User provides explicit Jacobian function
- ADflow: Need to extract Jacobian from Fortran solver

**The algorithms are identical** - we just need to bridge to ADflow's Jacobian!

## Minimal Working Example

Here's what we need to add to make a minimal working version:

### 1. Fortran (src/adjoint/adjointAPI.F90)

```fortran
subroutine getDenseJacobian(J_array, n) bind(C)
    ! Export dRdW as dense array (small problems only!)
    use ADjointPETSc, only: dRdW
    integer(c_int), intent(in) :: n
    real(c_double), dimension(n,n), intent(out) :: J_array

    ! Convert PETSc Mat to dense array
    call MatDenseGetArrayF90(dRdW, J_array, ierr)
end subroutine getDenseJacobian
```

### 2. Python (adflow/pyADflow.py)

```python
def getJacobianDense(self):
    """Get Jacobian as dense NumPy array (small problems only!)"""
    n = self.getStateSize()
    J = numpy.zeros((n, n))
    self.adflow.adjointapi.getdensejacobian(J, n)
    return J
```

### 3. Update resolvent class (adflow/pyResolventAnalysis.py)

```python
def solveExplicit(self):
    """Solve using dense Jacobian"""
    J = self.CFDsolver.getJacobianDense()  # ← Now works!

    omega = self.omega
    n = J.shape[0]

    # Form resolvent: R = (jω·I - J)^{-1}
    A = 1j * omega * np.eye(n) - J
    R = np.linalg.inv(A)

    # SVD of resolvent
    U, S, Vh = np.linalg.svd(R)

    self.sigma1 = S[0]
    self.u1 = U[:, 0]
    self.v1 = Vh[0, :].conj()

    return self.sigma1
```

**Total code to add:** ~30 lines Fortran, ~10 lines Python

## Summary

| Level | Fortran Changes | Python Changes | Dependencies | Problem Size |
|-------|----------------|----------------|--------------|--------------|
| **1. Dense** | ~30 lines | ~50 lines | NumPy/SciPy | < 10k DOF |
| **2. Matrix-free** | None! | ~200 lines | NumPy/SciPy | < 100k DOF |
| **3. SLEPc** | ~50 lines | ~500 lines | SLEPc | Unlimited |

## Answer to Your Question

> "Do we need to touch the source code? Need to add SLEPc?"

**For basic functionality:**
- Yes, minimal Fortran changes (~30 lines)
- No SLEPc needed (use NumPy SVD)
- Can test immediately on small problems

**For production use:**
- More Fortran integration (~50-100 lines)
- SLEPc highly recommended
- Required for large 3D CFD cases

**But:** The Python framework we created is complete and ready to use once the Jacobian interface is added!
