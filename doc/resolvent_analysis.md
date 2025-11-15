# Resolvent Analysis in ADflow

## Overview

Resolvent analysis is a frequency-domain stability analysis technique that extends classical linear stability analysis. It is used to study the frequency response of a linearized system in the vicinity of a steady-state or time-averaged solution.

This implementation is based on the paper:
> **"Large-Scale Flow Control Performance Optimization via Differentiable Resolvent Analysis"**
> by Sicheng He, Shugo Kaneko, Daning Huang, Chi-An Yeh, and Joaquim R. R. A. Martins

## Theory

### Resolvent Operator

The resolvent operator maps external forcing to flow response at a given frequency ω:

```
R(ω) = (jω·I - J)^(-1)
```

where:
- `ω` is the excitation frequency
- `J = ∂r/∂w` is the Jacobian matrix at the steady-state solution
- `j` is the imaginary unit
- `I` is the identity matrix

### Maximum Amplification

The maximum amplification factor is given by the dominant singular value of R:

```
σ₁ = max_δf (||δu||₂ / ||δf||₂)
```

where:
- `δf` is the external forcing
- `δu` is the flow response
- `σ₁` is obtained via Singular Value Decomposition (SVD) of R

### SVD of Resolvent Operator

```
R = U·Σ·V*
```

where:
- `U` contains left singular vectors (response modes)
- `V` contains right singular vectors (forcing modes)
- `Σ` contains singular values in decreasing order
- `V*` denotes complex conjugate transpose

## Real vs. Complex Formulation

The resolvent operator is inherently complex-valued. There are two approaches:

### 1. Complex Arithmetic (Default)

Work directly with complex matrices and vectors:
```
(jω·I - J)·u = f
```

### 2. Real Doubled Form

Convert to an equivalent real system of doubled size:

```
[ -J   -ωI ] [ u_r ]   [ f_r ]
[  ωI  -J  ] [ u_i ] = [ f_i ]
```

where `u = u_r + j·u_i` and `f = f_r + j·f_i`.

**Advantages of real form:**
- Works with existing real-valued PETSc solvers
- No complex arithmetic overhead
- Easier integration with ADflow's infrastructure

**Disadvantages:**
- Doubles the system size (2N instead of N)
- More memory usage

## Implementation in ADflow

### Classes

#### `ResolventAnalysis`
Basic resolvent analysis class using explicit Jacobian formation.

```python
from adflow import ADFLOW, ResolventAnalysis

# Create and solve CFD problem
CFDsolver = ADFLOW(options=aeroOptions)
CFDsolver(aeroProblem)

# Create resolvent analysis object
resolvent = ResolventAnalysis(CFDsolver, aeroProblem, omega=10.0)

# Solve for dominant singular value
sigma1 = resolvent.solveExplicit()

# Get modes
u1 = resolvent.getResponseMode()  # Response mode
v1 = resolvent.getForcingMode()   # Forcing mode
```

#### `ResolventAnalysisMatrixFree`
Matrix-free implementation using iterative methods (recommended for large problems).

```python
from adflow import ResolventAnalysisMatrixFree

# Create matrix-free resolvent analysis
resolvent = ResolventAnalysisMatrixFree(CFDsolver, aeroProblem, omega=10.0)
resolvent.setNModes(3)  # Compute top 3 modes

# Solve using iterative methods
sigma1 = resolvent.solve()
```

### Frequency Sweep

Analyze system response over a range of frequencies:

```python
omega_range = (0.0, 100.0)
nPoints = 50

omega_vec, sigma1_vec = resolvent.computeFrequencySweep(omega_range, nPoints)

# Find peak response
idx_max = np.argmax(sigma1_vec)
omega_peak = omega_vec[idx_max]
sigma1_peak = sigma1_vec[idx_max]
```

## Implementation Status

### ✓ Completed

- [x] Class structure and API design
- [x] Documentation based on resolvent theory paper
- [x] Integration framework with ADflow
- [x] Complex/real formulation handling
- [x] Frequency sweep interface

### 🚧 In Progress / Required

The following components require Fortran-level integration with ADflow:

1. **Jacobian Matrix Assembly**
   - Interface to ADflow's `dRdW` (∂R/∂w) matrix assembly
   - Access to PETSc Mat object for the Jacobian

2. **Matrix-Free Linear Solvers**
   - Integration with ADflow's PETSc-based linear solver infrastructure
   - Jacobian-vector products: `J·v`
   - Adjoint Jacobian-vector products: `J^T·v`

3. **Iterative SVD**
   - Arnoldi/Lanczos iteration for large-scale problems
   - Interface with SLEPc for eigenvalue/SVD computations

4. **Adjoint Derivatives** (for optimization)
   - Derivatives of singular values with respect to design variables
   - Integration with ADflow's adjoint infrastructure

## Example: Algebraic Model

A complete working example is available in `article_resolvent_opt/code/`:

```python
from nonlinear_eqn import nonlinear_eqn
from resolvent import resolvent

# Define the dynamical system residuals and derivatives
nonlinear_eqn_obj = nonlinear_eqn(
    ndof, x0,
    f_res, f_pres_pw, f_pres_px,
    f_p2res_pw2=f_p2res_pw2,
    f_p2res_pwpx=f_p2res_pwpx
)

# Create resolvent object
resolvent_obj = resolvent(nonlinear_eqn_obj, omega=1.0)

# Solve
s1 = resolvent_obj.solve()
print(f"Dominant singular value: {s1}")

# Compute derivatives (for optimization)
ds1_dx = resolvent_obj.compute_total_der()
```

See `article_resolvent_opt/code/example.py` for the complete 2D dynamical system example.

## Applications

Resolvent analysis can be used for:

1. **Flow Stability Analysis**
   - Identify amplification mechanisms
   - Locate critical frequencies
   - Study transition and turbulence

2. **Flow Control Design**
   - Determine optimal forcing locations and frequencies
   - Design passive flow control devices
   - Optimize active control systems

3. **Sensitivity Analysis**
   - Study how design changes affect stability
   - Gradient-based optimization with stability constraints

4. **Hopf Bifurcation Analysis**
   - Near a Hopf bifurcation, the resolvent becomes singular at the bifurcation frequency
   - Can be used to track bifurcation points

## Mathematical Details

### Relationship to Eigenvalue Analysis

Linear stability analysis solves:
```
J·q = λ·q
```

Resolvent analysis generalizes this by allowing non-zero forcing:
```
(jω·I - J)·u = f
```

When `ω = Im(λ)` for an eigenvalue λ, the resolvent becomes singular.

### Adjoint Sensitivity Analysis

For optimization with respect to design variables `x`, we need:
```
dσ₁/dx = ∂σ₁/∂x - ψ^T · ∂r/∂x
```

where `ψ` is the adjoint variable satisfying:
```
(∂r/∂w)^T · ψ = ∂σ₁/∂w
```

This allows computing derivatives with cost independent of the number of design variables.

## References

1. He, S., Kaneko, S., Huang, D., Yeh, C-A., and Martins, J. R. R. A.,
   "Large-Scale Flow Control Performance Optimization via Differentiable Resolvent Analysis,"
   AIAA Journal (submitted).

2. McKeon, B. J., and Sharma, A. S.,
   "A critical-layer framework for turbulent pipe flow,"
   Journal of Fluid Mechanics, Vol. 658, 2010, pp. 336-382.

3. Taira, K., et al.,
   "Modal analysis of fluid flows: Applications and outlook,"
   AIAA Journal, Vol. 58, No. 3, 2020, pp. 998-1022.

## Contact

For questions about the resolvent analysis implementation:
- See `article_resolvent_opt/code/` for working algebraic examples
- Refer to the paper for theoretical background
- Check `examples/resolvent_analysis_example.py` for ADflow usage
