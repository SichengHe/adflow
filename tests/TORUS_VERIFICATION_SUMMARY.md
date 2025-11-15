# Torus Time Spectral Verification Summary

## Overview
This document summarizes the verification tests for the Torus Time Spectral method in ADflow, specifically focusing on the degenerate case where ω₁ = ω₂.

## Test Files

### 1. Unit Test: Spectral Operator Verification
**File**: `unit_tests/test_torus_degenerate_matches_ts.py`

**Purpose**: Verify that the torus spectral differentiation operator with n₂=1 matches the classical 1D time spectral operator.

**Setup**:
- Torus: n₁=4, n₂=1, ω₁=ω₂=95 rad/s
- Classic 1D TS: n=4, ω=95 rad/s

**Result**: ✅ **PASSED**
- Operators match with sign flip (convention difference)
- Maximum difference: < 1e-12

**Key Finding**: The torus with n₂=1 is mathematically equivalent to classical 1D time spectral.

---

### 2. Operator Comparison: Matrix Analysis
**File**: `compare_torus_n2_1_vs_1d.py`

**Purpose**: Compare the spectral derivative matrices (dscalar) between torus n₂=1 and 1D TS.

**Setup**:
- Torus: n₁=3, n₂=1
- 1D TS: n=3
- Both with ω=100 rad/s

**Result**: ✅ **Matrices match perfectly** (with sign convention difference)
- Max normalized difference: < 1e-16

---

### 3. CFD Comparison: Diagonal Instances vs 1D TS
**File**: `compare_torus_diagonal_vs_1d_ts.py`

**Purpose**: Verify that the diagonal instances of a 3×3 degenerate torus match 3-point 1D TS.

**Setup**:
- Torus: 3×3 grid, ω₁=ω₂=100 rad/s, A₁=A₂=1.0°
- Extract diagonal instances: (0,0), (1,1), (2,2)
- 1D TS: 3 instances, ω=100 rad/s, A=2.0°

**Phase Mapping for Diagonal**:
- Instance (0,0): θ = 0°, α = 0°
- Instance (1,1): θ = 120°, α = 1.732°
- Instance (2,2): θ = 240°, α = -1.732°

**Result**: ✅ **PERFECT MATCH**
- Instance (0,0): ΔCL = 8.013e-11
- Instance (1,1): ΔCL = 2.939e-10
- Instance (2,2): ΔCL = 3.787e-11
- **Maximum CL difference: 2.939e-10** (machine precision)

**Key Finding**: For degenerate torus (ω₁=ω₂, A₁=A₂):
- **Diagonal instances** (i₁=i₂) correspond to physical time instants
- **Phase mapping**: θ = θ₁ = θ₂ (NOT (θ₁+θ₂)/2)
- **Off-diagonal instances** do NOT correspond to physical time - they are coupling artifacts

---

### 4. Mesh and Grid Velocity Verification
**File**: `verify_mesh_gridvel_match.py`

**Purpose**: Verify that mesh positions match at corresponding phases.

**Result**: ✅ **Confirmed**
- At θ=0°: Mesh positions match exactly (diff = 0.000e+00)
- At other matching phases: Mesh positions differ because α values differ
- This confirms that only diagonal instances represent the same physical state

---

## Key Conclusions

### 1. Mathematical Equivalence
- **Torus with n₂=1 ≡ Classical 1D Time Spectral**
- The spectral derivative operators are identical (up to sign convention)

### 2. Degenerate Torus Structure
For a degenerate torus with ω₁=ω₂ and A₁=A₂:
- **Diagonal instances (i₁=i₂)**: True physical time instances
  - Phase: θ = θ₁ = θ₂
  - Motion: α(θ) = 2A·sin(θ)
  - Match 1D TS with amplitude 2A

- **Off-diagonal instances (i₁≠i₂)**: Spectral coupling artifacts
  - Do NOT correspond to physical time
  - Different α at same projected phase
  - Required for spectral accuracy but not physical states

### 3. Phase Mapping
**CORRECT mapping for diagonal**:
```
θ_physical = θ₁ = θ₂
```

**INCORRECT mapping** (previously used):
```
θ_proj = (θ₁ + θ₂) / 2  ❌
```

### 4. Sign Convention
- Torus and classical 1D TS use **opposite sign conventions** in dscalar
- This is a convention difference, not an error
- Both give identical results

---

## Test Coverage

| Test | Aspect | Status |
|------|--------|--------|
| Unit test | Operator equivalence (n₂=1) | ✅ PASSED |
| Matrix comparison | Operator matching | ✅ PASSED |
| CFD diagonal comparison | Physical solution matching | ✅ PASSED |
| Mesh verification | Geometric consistency | ✅ PASSED |

---

## Running the Tests

### Unit Test (fast - no CFD solve)
```bash
cd /home/sicheng/repo/adflow/adflow/tests
testflo unit_tests/test_torus_degenerate_matches_ts.py -v
```

### Operator Comparison (fast - no CFD solve)
```bash
python compare_torus_n2_1_vs_1d.py
```

### Diagonal CFD Comparison (slow - full CFD solve)
```bash
mpirun -np 4 python compare_torus_diagonal_vs_1d_ts.py
```

### Mesh Verification (medium - setup only)
```bash
mpirun -np 4 python verify_mesh_gridvel_match.py
```

---

## References

- **Torus Time Spectral Method**: Treats time as a 2D torus parameterized by (θ₁, θ₂)
- **Classical Time Spectral**: Treats time as 1D circle parameterized by θ
- **Degenerate Case**: When ω₁=ω₂, the 2D problem collapses to 1D on the diagonal

---

**Date**: 2025-11-14
**Verified by**: Claude Code + User Testing
**Status**: All tests passing ✅
