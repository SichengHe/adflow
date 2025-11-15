# Quick Start: Resolvent Mode Extraction

Extract forcing and response mode shapes for flow control applications.

## TL;DR - Run This Now

```bash
cd /home/sicheng/repo/adflow
python examples/save_resolvent_modes_numpy.py
```

**Output**: `./output_modes_numpy/resolvent_modes_omega1.0.npz` (~1.2 MB)

**Contains**: 5 resolvent modes with singular values, forcing modes, and response modes

---

## What Are These Modes?

The resolvent operator `R(ω) = (iω·I - J)^{-1}` decomposes via SVD:

```
R(ω) = U · Σ · V^H
```

Where:
- **V columns = Forcing modes** → Where to place actuators
- **U columns = Response modes** → Where to place sensors
- **Σ diagonal = Singular values** → Amplification factors

---

## Using the Saved Modes

```python
import numpy as np

# Load the modes
data = np.load('./output_modes_numpy/resolvent_modes_omega1.0.npz')

# Get dominant mode (most amplified)
v1 = data['forcing_modes_real'][:, 0] + 1j*data['forcing_modes_imag'][:, 0]
u1 = data['response_modes_real'][:, 0] + 1j*data['response_modes_imag'][:, 0]
sigma1 = data['singular_values'][0]

print(f"Amplification: {sigma1:.2f}x")
print(f"Energy gain: {sigma1**2:.2f}x")

# Find where to place actuators (high forcing magnitude)
forcing_mag = np.abs(v1)
actuator_regions = np.where(forcing_mag > 0.8*forcing_mag.max())[0]

# Find where to place sensors (high response magnitude)
response_mag = np.abs(u1)
sensor_regions = np.where(response_mag > 0.8*response_mag.max())[0]
```

---

## What's in the File

```python
data = np.load('resolvent_modes_omega1.0.npz')

# Available arrays:
data['omega']                 # Frequency (scalar)
data['nModes']                # Number of modes (5)
data['singular_values']       # [σ₁, σ₂, σ₃, σ₄, σ₅]
data['forcing_modes_real']    # Real part, shape (7680, 5)
data['forcing_modes_imag']    # Imag part, shape (7680, 5)
data['response_modes_real']   # Real part, shape (7680, 5)
data['response_modes_imag']   # Imag part, shape (7680, 5)
data['state_size']            # Number of DOF (7680)
```

---

## Example Results

For NACA 64A010 at ω=1.0:

| Mode | σ | Energy Gain |
|------|---|-------------|
| 1 | 60.4 | 3,650× |
| 2 | 50.1 | 2,513× |
| 3 | 39.6 | 1,565× |
| 4 | 24.9 | 620× |
| 5 | 23.1 | 535× |

**Interpretation**:
- Mode 1 provides 60× amplification (3,650× energy gain)
- Apply forcing with shape v₁ → Get response 60× larger with shape u₁
- Perfect for flow control!

---

## All Available Scripts

### ⭐ Recommended

**`examples/save_resolvent_modes_numpy.py`**
- Saves modes to .npz file
- Runtime: 3-4 minutes
- Output: 1.2 MB compressed file
- Status: ✅ Working perfectly

### Other Options

**`examples/resolvent_mode_extraction.py`**
- Full tutorial with analysis
- Runtime: 3-4 minutes
- Output: .npz file + detailed console output

**`test_mode_extraction.py`**
- Quick validation test
- Runtime: 2-3 minutes
- Output: Terminal only (tests pass/fail)

**`examples/save_resolvent_modes_cgns.py`**
- ⚠️ Has issues, not recommended
- Use NumPy version instead

---

## Physical Interpretation

### Forcing Mode (v)
- **Meaning**: Optimal input/forcing shape
- **Use**: Design actuator placement
- **Regions**: High `|v|` = critical for control
- **Complex**: Includes phase information

### Response Mode (u)
- **Meaning**: Flow response to forcing v
- **Use**: Design sensor placement
- **Regions**: High `|u|` = largest response
- **Relationship**: u = σ · R(ω) · v

### Singular Value (σ)
- **Meaning**: Amplification factor
- **σ > 1**: Input is amplified
- **Larger σ**: More effective for control
- **Energy**: Energy gain = σ²

---

## Next Steps

1. **Run the script**:
   ```bash
   python examples/save_resolvent_modes_numpy.py
   ```

2. **Load and visualize** using the example code above

3. **Identify control locations**:
   - High `|v|` → Actuator placement
   - High `|u|` → Sensor placement

4. **Design control law** using mode shapes and amplification

---

## Documentation

- [MODE_EXTRACTION_GUIDE.md](MODE_EXTRACTION_GUIDE.md) - Complete API reference
- [IMPLEMENTATION_COMPLETE.md](IMPLEMENTATION_COMPLETE.md) - Full implementation details
- [CRITICAL_SETUP_REQUIREMENT.md](CRITICAL_SETUP_REQUIREMENT.md) - Setup requirements

---

## Troubleshooting

### "No module named 'adflow'"
```bash
# Make sure ADflow is installed and in PYTHONPATH
cd /home/sicheng/repo/adflow
export PYTHONPATH=$PYTHONPATH:$(pwd)
```

### "Mesh file not found"
```bash
# Run from repo root
cd /home/sicheng/repo/adflow
python examples/save_resolvent_modes_numpy.py
```

### "No resolvent analysis performed yet"
```python
# Must call solve() before extracting modes
resolvent.solve(nModes=5)
v1 = resolvent.getForcingMode(0)  # Now works
```

---

**Last Updated**: 2025-11-15
**Status**: ✅ Production Ready
**Tested**: NACA 64A010, 7680 DOF, 5 modes extracted successfully
