# Modal Analysis Module

This directory contains modal analysis tools for studying flow dynamics and stability.

## Overview

Modal analysis techniques decompose the flow into modes or structures to understand:
- Flow stability and instability mechanisms
- Frequency response characteristics
- Amplification mechanisms
- Optimal forcing and response patterns

## Current Implementation

### Resolvent Analysis (`resolventAPI.F90`)

Resolvent analysis is a frequency-domain stability analysis technique that computes the frequency response of a linearized system.

**Theory:**
The resolvent operator maps external forcing to flow response:
```
R(ω) = (jω·I - J)^{-1}
```
where:
- `ω` is the excitation frequency
- `J = ∂R/∂w` is the Jacobian at steady state
- `R(ω)` is the resolvent operator

**Key Features:**
- Assembles the Jacobian matrix from ADflow's adjoint infrastructure
- Exports matrix in dense format or to PETSc binary files
- Supports both frozen and coupled turbulence

**Usage:**
```fortran
! Setup the Jacobian for resolvent analysis
call setupResolventMatrix(frozenTurb=.true.)

! Get matrix information
call getResolventMatrixInfo(nRows, nCols, nnz)

! Export as dense (small problems only!)
call getResolventMatrixDense(J_array, n)

! Export to file for large problems
call exportResolventMatrixToFile('jacobian.bin')
```

**Python Interface:**
See `adflow.ResolventAnalysis` class in `/adflow/pyResolventAnalysis.py`

## Future Modal Analysis Tools

This directory is designed to accommodate additional modal analysis methods:

### Eigenvalue Analysis
- Global linear stability analysis
- Computing eigenvalues and eigenmodes of J
- Identifies unstable modes and growth rates

### Dynamic Mode Decomposition (DMD)
- Data-driven modal decomposition
- Extracts coherent structures from snapshots
- Applicable to unsteady flows

### Proper Orthogonal Decomposition (POD)
- Energy-based modal decomposition
- Identifies dominant flow structures
- Useful for reduced-order modeling

### Koopman Analysis
- Nonlinear extension of eigenvalue analysis
- Data-driven approach
- Useful for complex flows

## Organization

```
modalAnalysis/
├── README.md              # This file
├── resolventAPI.F90       # Resolvent analysis routines
├── eigenvalueAPI.F90      # (Future) Eigenvalue analysis
├── dmdAPI.F90             # (Future) DMD analysis
└── ...
```

## Integration with ADflow

Modal analysis tools leverage ADflow's existing infrastructure:

1. **Jacobian Assembly**: Uses `adjointUtils::setupStateResidualMatrix`
2. **PETSc Matrices**: Accesses `dRdWT` from `ADjointPETSc` module
3. **State Variables**: Works with converged flow solutions

## References

### Resolvent Analysis
1. He, S., et al., "Large-Scale Flow Control Performance Optimization via Differentiable Resolvent Analysis"
2. McKeon, B. J., and Sharma, A. S., "A critical-layer framework for turbulent pipe flow," JFM, 2010
3. Taira, K., et al., "Modal analysis of fluid flows: Applications and outlook," AIAA Journal, 2020

### General Modal Analysis
4. Schmid, P. J., "Dynamic mode decomposition of numerical and experimental data," JFM, 2010
5. Rowley, C. W., et al., "Spectral analysis of nonlinear flows," JFM, 2009

## Notes

- Modal analysis is **distinct from adjoint methods** (though it uses the same Jacobian)
- Adjoint methods: Compute gradients for optimization (`∂f/∂x`)
- Modal analysis: Study dynamics and stability (`eigenvalues, singular values`)
- Both use `J = ∂R/∂w`, but for different purposes!

## Contact

For questions about modal analysis implementation in ADflow, see:
- `/doc/resolvent_analysis.md` for resolvent-specific documentation
- `/examples/resolvent_analysis_example.py` for usage examples
