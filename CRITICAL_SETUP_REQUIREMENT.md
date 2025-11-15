# CRITICAL SETUP REQUIREMENT for Resolvent Analysis

## Problem

When running resolvent analysis with explicit Jacobian extraction, you may encounter this error:

```
[0]PETSC ERROR: No support for this operation for this object type
[0]PETSC ERROR: No method zeroentries for Mat of type shell
```

## Root Cause

ADflow by default uses **matrix-free (shell)** Jacobian matrices for adjoint computations to save memory. These shell matrices only support matrix-vector products and **cannot** be explicitly extracted or zeroed.

The option controlling this behavior is:
- **`useMatrixFreedrdw = True`** (default): Uses shell matrix → Cannot extract
- **`useMatrixFreedrdw = False`** (required): Uses assembled AIJ matrix → Can extract

## Solution

**For explicit resolvent analysis (`ResolventAnalysis`), you MUST set:**

```python
aeroOptions = {
    # ... other options ...

    # CRITICAL: Force assembled Jacobian matrix (not matrix-free shell)
    'useMatrixFreedrdw': False,
}
```

## When is this required?

### ✓ Required (must set `useMatrixFreedrdw = False`)

1. **Explicit Resolvent Analysis** (`ResolventAnalysis` class)
   - Needs explicit Jacobian matrix extraction
   - Calls `getJacobianMatrix()` to extract J as dense/sparse array

2. **ILU Preconditioning for Matrix-Free** (`ResolventAnalysisMatrixFree` with ILU)
   - Even though matrix-free, ILU needs explicit Jacobian
   - Calls `getJacobianMatrix()` to build ILU factorization

### ✗ Not Required (can use default `useMatrixFreedrdw = True`)

1. **Pure Matrix-Free** (`ResolventAnalysisMatrixFree` without preconditioning)
   - Only uses Jacobian-vector products
   - Never calls `getJacobianMatrix()`
   - Can use shell matrix for memory efficiency

## Implementation Details

### What happens with `useMatrixFreedrdw = False`?

ADflow will:
1. Allocate sparse AIJ matrix for `dRdWT` (transpose of Jacobian)
2. Assemble matrix entries via `setupStateResidualMatrix()`
3. Support `MatZeroEntries()`, `MatConvert()`, extraction, etc.

**Memory cost**: O(nnz) where nnz ≈ 170 × n for typical 3D viscous cases
- n = 7,680: ~450 MB
- n = 100,000: ~6 GB
- n = 1,000,000: ~60 GB

### What happens with `useMatrixFreedrdw = True` (default)?

ADflow will:
1. Create shell matrix for `dRdWT`
2. Only define matrix-vector product operations
3. **Cannot** call `MatZeroEntries()`, `MatConvert()`, or extract values

**Memory cost**: O(1) - just function pointers
- Any size n: ~10 KB

## Complete Example

```python
from adflow import ADFLOW, ResolventAnalysis
from baseclasses import AeroProblem

# =========================================================================
# CRITICAL: Set useMatrixFreedrdw = False for resolvent analysis
# =========================================================================
aeroOptions = {
    'gridFile': 'naca64A010.cgns',
    'equationType': 'Euler',

    # ... standard flow solver options ...

    # REQUIRED for explicit Jacobian extraction
    'useMatrixFreedrdw': False,  # <-- CRITICAL!
}

# Create solver
CFDsolver = ADFLOW(options=aeroOptions)

# Setup and solve
ap = AeroProblem(
    name='naca64a010',
    mach=0.5,
    altitude=10000,
    alpha=0.0,
    areaRef=1.0,
    chordRef=1.0,
)
CFDsolver(ap)

# =========================================================================
# Now resolvent analysis will work
# =========================================================================
resolvent = ResolventAnalysis(CFDsolver, ap, omega=1.0)
sigma = resolvent.solve(nModes=3, method='scipy', form='complex')
```

## Error Messages

### If you forget to set `useMatrixFreedrdw = False`:

```
[0]PETSC ERROR: No support for this operation for this object type
[0]PETSC ERROR: No method zeroentries for Mat of type shell
[0]PETSC ERROR: #1 MatZeroEntries() at .../src/mat/interface/matrix.c:6059
```

**Fix**: Add `'useMatrixFreedrdw': False` to `aeroOptions`

### If you have memory issues with large problems:

```
ERROR: Cannot allocate 56 GB for Jacobian matrix
```

**Fix**: Use `ResolventAnalysisMatrixFree` instead of `ResolventAnalysis`

## Summary Table

| Method | Class | Requires `useMatrixFreedrdw=False`? | Memory Scaling |
|--------|-------|-------------------------------------|----------------|
| Explicit | `ResolventAnalysis` | **YES** | O(n²) dense or O(nnz) sparse |
| Matrix-Free + ILU | `ResolventAnalysisMatrixFree` with ILU | **YES** | O(nnz) for ILU |
| Pure Matrix-Free | `ResolventAnalysisMatrixFree` no precond | No | O(n) |

## Recommendation

**Default setup for most users:**

```python
aeroOptions = {
    # ... other options ...
    'useMatrixFreedrdw': False,  # Always set this for resolvent analysis
}
```

Only use `useMatrixFreedrdw = True` (default) if:
- You have very large problems (n > 1M)
- You're using pure matrix-free (no ILU)
- You're not extracting the Jacobian

---

**Last Updated**: 2025-11-15
**Status**: Critical requirement for all explicit Jacobian extraction
