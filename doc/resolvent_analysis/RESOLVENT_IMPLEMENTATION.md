# Resolvent Analysis Implementation in ADflow

## Summary

This implementation adds **resolvent analysis** capabilities to ADflow for frequency-domain flow stability analysis. The implementation is based on the paper "Large-Scale Flow Control Performance Optimization via Differentiable Resolvent Analysis" by He et al.

## What is Resolvent Analysis?

Resolvent analysis is a frequency-domain technique that extends linear stability analysis. It computes the frequency response of a linearized system:

```
R(ω) = (jω·I - J)^{-1}
```

The maximum amplification factor `σ₁` (dominant singular value of R) indicates how much the flow amplifies external forcing at frequency `ω`.

## Files Added

### Core Implementation

1. **`adflow/pyResolventAnalysis.py`**
   - Main resolvent analysis classes:
     - `ResolventAnalysis`: Base class using explicit Jacobian
     - `ResolventAnalysisMatrixFree`: Matrix-free implementation for large problems
   - Handles both complex and real doubled formulations
   - ~450 lines of Python code

2. **`adflow/__init__.py`** (modified)
   - Exports resolvent analysis classes
   - Integration with ADflow package structure

### Documentation

3. **`doc/resolvent_analysis.md`**
   - Comprehensive documentation
   - Theory background
   - Usage examples
   - Implementation details
   - References

4. **`examples/resolvent_analysis_example.py`**
   - Complete usage example
   - Shows integration with ADflow workflow
   - Demonstrates frequency sweep
   - Both explicit and matrix-free approaches

## Key Features

### ✓ Implemented

- **Class structure** following ADflow conventions
- **Complex/Real formulations** for compatibility
  - Complex arithmetic (default)
  - Real doubled system (for PETSc integration)
- **API design** for:
  - Single frequency analysis
  - Frequency sweep
  - Mode extraction (response and forcing modes)
- **Documentation** with theory and usage
- **Integration hooks** with ADflow infrastructure

### 🚧 Requires Fortran Integration

The following require deep integration with ADflow's Fortran core:

1. **Jacobian Access**
   ```fortran
   ! Need access to dRdW matrix assembly
   ! Currently handled in ADflow's Fortran routines
   ```

2. **Linear Solvers**
   - Integration with PETSc KSP solvers
   - Matrix-free Jacobian-vector products
   - Adjoint Jacobian-vector products

3. **Iterative SVD**
   - SLEPc integration for large-scale SVD
   - Arnoldi/Lanczos methods

## Complex vs. Real Formulation

### Question: Should we cast to real/imaginary parts?

**Answer:** Yes, and the implementation supports both!

The resolvent operator `(jω·I - J)` is inherently complex. We provide two formulations:

### 1. Complex Arithmetic (Straightforward)
```python
# Direct complex arithmetic
A = 1j*omega*np.eye(n) - J
R = np.linalg.inv(A)
U, S, Vh = np.linalg.svd(R)
```

**Pros:**
- Mathematically natural
- Smaller system size (N × N)
- Cleaner code

**Cons:**
- Requires complex linear solver
- Not all PETSc configurations support complex

### 2. Real Doubled Form (For PETSc)
```python
# Equivalent real system (2N × 2N)
# [ -J   -ωI ] [ u_r ]   [ f_r ]
# [  ωI  -J  ] [ u_i ] = [ f_i ]
```

**Pros:**
- Works with real-valued PETSc
- Compatible with existing ADflow infrastructure
- Standard in CFD applications

**Cons:**
- Doubles system size
- More complex bookkeeping
- Higher memory usage

**The implementation includes helper methods:**
```python
def _complexToRealForm(self, A_complex):
    """Convert complex matrix to real doubled form"""

def _realToComplexForm(self, v_real):
    """Convert real doubled vector to complex"""
```

## Usage Example

```python
from adflow import ADFLOW, ResolventAnalysis
from baseclasses import AeroProblem

# 1. Solve steady-state CFD
CFDsolver = ADFLOW(options=aeroOptions)
ap = AeroProblem(name='wing', mach=0.8, alpha=2.0, ...)
CFDsolver(ap)

# 2. Create resolvent analysis object
omega = 10.0  # Excitation frequency
resolvent = ResolventAnalysis(CFDsolver, ap, omega=omega)

# 3. Perform analysis (when Jacobian interface is complete)
sigma1 = resolvent.solveExplicit(useRealForm=True)
print(f"Maximum amplification: σ₁ = {sigma1}")

# 4. Get modes
u1 = resolvent.getResponseMode()  # Where flow responds
v1 = resolvent.getForcingMode()   # Where to apply forcing

# 5. Frequency sweep
omega_vec, sigma1_vec = resolvent.computeFrequencySweep((0, 100), nPoints=50)
```

## Working Examples

### Algebraic System (Complete Implementation)

The `article_resolvent_opt/` directory contains fully working examples:

```
article_resolvent_opt/code/
├── resolvent.py              # Resolvent class (complete)
├── nonlinear_eqn.py          # Base nonlinear equation class
├── example.py                # 2D dynamical system
├── reactor_optimization_2D.py # Brusselator example
└── README.md                 # Usage instructions
```

**These examples are fully functional** and demonstrate:
- Resolvent analysis on algebraic systems
- Adjoint-based sensitivity analysis
- Gradient-based optimization
- Frequency sweeps

Run them to verify the theory:
```bash
cd article_resolvent_opt/code
python example.py
```

## Integration Roadmap

To complete the ADflow integration:

### Phase 1: Basic Jacobian Access (Minimal)
```python
# Add method to ADflow class
def getJacobianMatrix(self):
    """Return assembled dRdW matrix"""
    # Interface to Fortran: call assemble_drdw()
    return J_petsc
```

This alone would enable `ResolventAnalysis.solveExplicit()`.

### Phase 2: Matrix-Free Operations (Recommended)
```python
# Already exists in ADflow!
computeJacobianVectorProductFwd(wDot, residualDeriv=True)  # J·v
computeJacobianVectorProductBwd(resBar, wDeriv=True)       # J^T·v
```

These can be used for matrix-free resolvent operator application.

### Phase 3: Iterative SVD (For Large Problems)
- Interface with SLEPc
- Arnoldi/Lanczos iterations
- Requires both forward and adjoint Jacobian-vector products

### Phase 4: Adjoint Derivatives (For Optimization)
- Derivatives of σ₁ with respect to design variables
- Integration with ADflow's adjoint infrastructure
- Enables resolvent-based shape optimization

## Testing

### Unit Tests (To Be Added)
```python
# tests/test_resolvent.py
def test_complex_to_real_conversion():
    """Test that complex and real forms are equivalent"""

def test_small_system():
    """Test resolvent on a small manufactured problem"""

def test_jacobian_products():
    """Verify J·v and J^T·v using finite differences"""
```

### Integration Tests
Compare with:
- Eigenvalue analysis (ω = 0)
- Known analytical solutions
- The algebraic examples in `article_resolvent_opt/`

## References

See `doc/resolvent_analysis.md` for complete references and theory.

Key paper:
- He, S., et al., "Large-Scale Flow Control Performance Optimization via
  Differentiable Resolvent Analysis," (submitted to AIAA Journal)

## Next Steps

1. **Implement Jacobian interface** (`getJacobianMatrix()`)
2. **Test on small problems** (2D cylinder, etc.)
3. **Add matrix-free linear solver** integration
4. **Implement iterative SVD** for large problems
5. **Add adjoint derivatives** for optimization
6. **Validate against** known results and eigenvalue analysis

## Notes

- The algebraic examples in `article_resolvent_opt/` are **fully working** and can be used as reference
- The ADflow infrastructure **already has** Jacobian-vector products - we just need to interface with them
- The real/complex formulation choice is **abstracted** in the API
- Documentation is **comprehensive** and based on the paper

---

**Implementation Status:** Framework complete, awaiting Fortran integration for Jacobian access.

**Author:** Based on resolvent analysis paper by He et al.
**Date:** 2025-11-15
