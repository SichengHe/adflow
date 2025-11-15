# Resolvent Mode Extraction Guide

Quick reference for extracting and using resolvent mode shapes.

## What Are Resolvent Modes?

Resolvent analysis provides the SVD of the input-output operator:

```
R(ω) = (iω·I - J)^{-1} = U · Σ · V^H
```

Where:
- **V columns** = **Forcing modes** (optimal input shapes)
- **U columns** = **Response modes** (output responses)
- **Σ diagonal** = **Singular values** (amplification factors)

## Quick Start

### 1. Basic Mode Extraction

```python
from adflow import ADFLOW, ResolventAnalysisMatrixFree

# Setup and solve CFD (with useMatrixFreedrdw=False)
CFDsolver = ADFLOW(options=aeroOptions)
CFDsolver(ap)

# Compute resolvent modes
resolvent = ResolventAnalysisMatrixFree(CFDsolver, ap, omega=1.0)
resolvent.enablePreconditioner(precond_type='ilu')
resolvent.solve(nModes=5)

# Extract modes
v1 = resolvent.getForcingMode(0)    # Dominant forcing mode
u1 = resolvent.getResponseMode(0)   # Dominant response mode
sigma1 = resolvent.getSingularValue(0)  # Amplification factor

print(f"Amplification: {sigma1:.2f}x")
print(f"Forcing shape: {v1.shape}")
print(f"Response shape: {u1.shape}")
```

### 2. Extract Multiple Modes

```python
nModes = 5
for i in range(nModes):
    v_i = resolvent.getForcingMode(i)
    u_i = resolvent.getResponseMode(i)
    sigma_i = resolvent.getSingularValue(i)

    print(f"Mode {i+1}: σ = {sigma_i:.6f}, energy = {sigma_i**2:.2f}")
```

### 3. Save Modes for Later Use

```python
import numpy as np

# Save all modes
np.savez(
    'resolvent_modes.npz',
    forcing_modes=[resolvent.getForcingMode(i) for i in range(nModes)],
    response_modes=[resolvent.getResponseMode(i) for i in range(nModes)],
    singular_values=[resolvent.getSingularValue(i) for i in range(nModes)],
    omega=resolvent.omega
)

# Load later
data = np.load('resolvent_modes.npz')
v1 = data['forcing_modes'][0]
```

### 4. Save as CGNS for Visualization

```python
# Method 1: Use provided script
# See examples/save_resolvent_modes_cgns.py

# Method 2: Manual save
for i in range(nModes):
    v_i = resolvent.getForcingMode(i)
    u_i = resolvent.getResponseMode(i)

    # Set state and write
    CFDsolver.setStates(np.real(v_i))
    CFDsolver.writeSolution(f'mode{i+1}_forcing_real.cgns')

    CFDsolver.setStates(np.imag(v_i))
    CFDsolver.writeSolution(f'mode{i+1}_forcing_imag.cgns')
```

## Mode Interpretation

### Singular Values (Amplification)

```python
sigma1 = resolvent.getSingularValue(0)
```

- **σ > 1**: Input is amplified (unstable direction)
- **σ < 1**: Input is damped
- **Larger σ**: More effective for control

**Example**: σ₁ = 10 means forcing with shape v₁ produces response 10× larger

### Forcing Modes (Where to Apply Control)

```python
v1 = resolvent.getForcingMode(0)  # Optimal forcing shape
```

- **Physical meaning**: Spatial distribution of optimal input/forcing
- **For control**: Where to place actuators
- **High |v₁|**: Critical regions for forcing
- **Complex values**: Include phase information

### Response Modes (Expected Flow Response)

```python
u1 = resolvent.getResponseMode(0)  # Expected response shape
```

- **Physical meaning**: Flow response to forcing v₁
- **For control**: Where to place sensors
- **High |u₁|**: Regions with largest response
- **Relationship**: u₁ = σ₁ · R(ω) · v₁

## Practical Applications

### 1. Flow Control Design

```python
# Find optimal actuator placement
v1 = resolvent.getForcingMode(0)
forcing_magnitude = np.abs(v1)

# Regions where forcing_magnitude is high → place actuators
actuator_locations = np.where(forcing_magnitude > 0.8 * forcing_magnitude.max())[0]
```

### 2. Sensor Placement

```python
# Find optimal sensor locations
u1 = resolvent.getResponseMode(0)
response_magnitude = np.abs(u1)

# Regions where response_magnitude is high → place sensors
sensor_locations = np.where(response_magnitude > 0.8 * response_magnitude.max())[0]
```

### 3. Energy Analysis

```python
# Energy distribution across modes
nModes = 5
total_energy = sum(resolvent.getSingularValue(i)**2 for i in range(nModes))

for i in range(nModes):
    sigma_i = resolvent.getSingularValue(i)
    energy_fraction = (sigma_i**2 / total_energy) * 100
    print(f"Mode {i+1}: {energy_fraction:.1f}% of energy")
```

### 4. Modal Control Law

```python
# Design controller using dominant mode
v1 = resolvent.getForcingMode(0)
u1 = resolvent.getResponseMode(0)
sigma1 = resolvent.getSingularValue(0)

# Control input: f(t) = f_magnitude * Re(v1 * e^(iωt))
# Expected output: u(t) ≈ sigma1 * f_magnitude * Re(u1 * e^(iωt))

print(f"Expected amplification: {sigma1:.2f}x")
print(f"Energy amplification: {sigma1**2:.2f}x")
```

## API Reference

### Both Classes (ResolventAnalysis & ResolventAnalysisMatrixFree)

```python
# Get singular value
sigma = resolvent.getSingularValue(idx=0)

# Get forcing mode (optimal input)
v = resolvent.getForcingMode(idx=0)

# Get response mode (output)
u = resolvent.getResponseMode(idx=0)
```

**Parameters**:
- `idx` : int, Index of mode (0 = dominant/largest)

**Returns**:
- `sigma` : float, Singular value (amplification factor)
- `v` : numpy array (complex), Forcing mode shape
- `u` : numpy array (complex), Response mode shape

## Mode Properties

### Complex Values

Modes are complex-valued to capture phase information:

```python
v1 = resolvent.getForcingMode(0)

# Separate real and imaginary parts
v1_real = np.real(v1)
v1_imag = np.imag(v1)

# Magnitude and phase
v1_magnitude = np.abs(v1)
v1_phase = np.angle(v1)

# Time-varying signal: v(t) = Re(v1 * e^(iωt))
import matplotlib.pyplot as plt
t = np.linspace(0, 10, 100)
v_t = np.real(v1[0] * np.exp(1j * resolvent.omega * t))
plt.plot(t, v_t)
```

### Normalization

Modes are orthonormal:

```python
# Response modes are orthonormal: U^H @ U = I
# Forcing modes are orthonormal: V^H @ V = I

v1 = resolvent.getForcingMode(0)
v2 = resolvent.getForcingMode(1)

# Should be 1 for same mode, 0 for different modes
inner_product = np.vdot(v1, v2)
print(f"<v1, v2> = {inner_product:.6f}")
```

## Examples

See complete examples in:
- `examples/resolvent_mode_extraction.py` - Full tutorial
- `examples/save_resolvent_modes_cgns.py` - CGNS export for visualization

## Visualization Workflow

### 1. Extract Modes (Python)
```python
resolvent.solve(nModes=5)
for i in range(5):
    v_i = resolvent.getForcingMode(i)
    # Save as CGNS
```

### 2. Visualize (Tecplot/ParaView)
- Load CGNS files
- Plot forcing magnitude contours
- Identify actuator placement regions
- Plot response magnitude for sensor placement

### 3. Implement Control
- Design actuators at high-forcing regions
- Design sensors at high-response regions
- Use σ values to predict control effectiveness

## Troubleshooting

### "No resolvent analysis performed yet"
```python
# Must call solve() first
resolvent.solve(nModes=5)
# Then extract modes
v1 = resolvent.getForcingMode(0)
```

### Index out of bounds
```python
# Don't exceed number of computed modes
resolvent.solve(nModes=3)  # Only computed 3 modes
v4 = resolvent.getForcingMode(3)  # ERROR: only have indices 0,1,2
```

### Modes have wrong shape
```python
# Mode shape = (stateSize,)
# For Euler: stateSize = nCells * 5  (ρ, ρu, ρv, ρw, ρE)
# For RANS: stateSize = nCells * 6  (+ turbulence)

v1 = resolvent.getForcingMode(0)
print(f"State size: {v1.shape}")  # e.g., (7680,) for 1536 cells × 5 vars
```

## References

- Implementation: [IMPLEMENTATION_COMPLETE.md](IMPLEMENTATION_COMPLETE.md)
- Setup: [CRITICAL_SETUP_REQUIREMENT.md](CRITICAL_SETUP_REQUIREMENT.md)
- Theory: "Differentiable Resolvent Analysis" by He et al.

---

**Last Updated**: 2025-11-15
**Status**: Production Ready ✓
