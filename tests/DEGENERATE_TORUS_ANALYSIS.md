# Degenerate Torus Time Spectral - Physics Analysis

## Issue Summary

The degenerate torus test (`test_torus_degenerate.py`) was showing instances with the same angle of attack (α=0°) producing different lift coefficients (CL). This investigation reveals this is **correct physics**, not a bug.

## Test Script Issues

### 1. Regex Extraction Bug (MINOR)

**File**: `test_torus_degenerate.py` line 194

**Current regex**:
```python
pattern = r'^\s*1\s+(\d+)\s+\d+\s+\d+\s+NK\s+[-\w]+\s+[\d\.]+\s+[\d\.]+\s+([\d\.\-E\+]+)\s+([\d\.\-E\+]+)\s+([\d\.\-E\+]+)'
cl_val = float(match.group(2))  # WRONG - extracts residual column (~1e-12)
```

**Problem**: The regex matches but extracts `group(2)` which corresponds to the `resrho` column (values ~1e-12), not the CL column.

**NK line format**:
```
1  instance  iter  total  NK  ----  CFL  ?  resrho  CL  CD  total_residual
```

**Fix**: Need to match "----" explicitly and extract `group(5)` for CL:
```python
pattern = r'^\s*1\s+(\d+)\s+\d+\s+\d+\s+NK\s+----\s+([\d\.]+)\s+([\d\.]+)\s+([-\d\.E\+]+)\s+([-\d\.E\+]+)\s+([-\d\.E\+]+)\s+([-\d\.E\+]+)'
cl_val = float(match.group(5))  # Correct CL extraction
```

### 2. Incorrect Physics Expectation (MAJOR)

**Problem**: The test expects instances with the same instantaneous α to produce identical CL values. This assumption is **WRONG** for torus time spectral!

## Physics of Degenerate Torus

### Configuration
- ω₁ = ω₂ = 100 rad/s (same frequencies)
- A₁ = A₂ = 1.0° (same amplitudes)
- Grid: 3×3 = 9 instances
- Motion: α(θ₁, θ₂) = A·[sin(θ₁) + sin(θ₂)]

### Why Same α ≠ Same CL

The flow solution depends on BOTH instantaneous geometry AND grid velocity:

1. **Instantaneous geometry**: Determined by α(θ₁, θ₂) = A·[sin(θ₁) + sin(θ₂)]

2. **Grid velocity**:
   ```
   v_grid = ω₁·∂x/∂θ₁ + ω₂·∂x/∂θ₂
          = ω·∂x/∂α·A·[cos(θ₁) + cos(θ₂)]
   ```

### Example: Three instances at α=0°

| Instance | (θ₁, θ₂) | α | Grid velocity factor | CL |
|----------|----------|---|---------------------|-----|
| 1 | (0, 0) | 0° | cos(0)+cos(0) = **+2.0** | -0.0850 |
| 6 | (4π/3, 2π/3) | 0° | cos(4π/3)+cos(2π/3) = **-1.0** | +0.0475 |
| 8 | (2π/3, 4π/3) | 0° | cos(2π/3)+cos(4π/3) = **-1.0** | +0.0475 |

**Observations**:
- Instance 1 has **opposite sign** grid velocity compared to instances 6 & 8
- Instances 6 & 8 have **identical** (θ₁, θ₂, α, v_grid) → **identical CL** ✓
- Instance 1 vs 6/8: **same α** but **different v_grid** → **different CL** ✓

This is **physically correct** behavior!

### Mathematical Proof

For rigid pitching motion: x(θ₁, θ₂) = R(α(θ₁, θ₂))·x₀

Where R is rotation matrix depending on α = A·[sin(θ₁) + sin(θ₂)].

Grid velocity components:
```
∂x/∂θ₁ = (∂R/∂α)·x₀ · ∂α/∂θ₁ = (∂R/∂α)·x₀ · A·cos(θ₁)
∂x/∂θ₂ = (∂R/∂α)·x₀ · ∂α/∂θ₂ = (∂R/∂α)·x₀ · A·cos(θ₂)
```

Total grid velocity:
```
v = ω₁·∂x/∂θ₁ + ω₂·∂x/∂θ₂
  = ω·(∂R/∂α)·x₀ · A·[cos(θ₁) + cos(θ₂)]
```

For ω₁ = ω₂ = ω.

The factor `[cos(θ₁) + cos(θ₂)]` varies from -2 to +2 and depends on the path through (θ₁, θ₂) space, NOT just on the instantaneous α value.

## Comparison with 1D Time Spectral

### Would 1D TS match?

**NO!** The degenerate torus does NOT collapse to 1D time spectral.

**1D time spectral** with α(θ) = 2A·sin(θ):
- Grid velocity: v = ω · ∂x/∂α · 2A·cos(θ)
- Range: velocity factor from -2A to +2A
- Distribution: Simple 1D sinusoidal

**Degenerate torus** with α(θ₁,θ₂) = A·[sin(θ₁)+sin(θ₂)]:
- Grid velocity: v = ω · ∂x/∂α · A·[cos(θ₁)+cos(θ₂)]
- Range: velocity factor from -2A to +2A (same range!)
- Distribution: 2D grid with multiple paths to same α

These are **different solution manifolds** with different grid velocity histories.

## Expected Test Results

### Instances that SHOULD match (same θ₁, θ₂, α, v_grid):

**Perfect matches found**:
- Instances 2 & 4: (1,0) & (0,1) → α=0.866°, CL=0.0963641930 ✓
- Instances 3 & 7: (2,0) & (0,2) → α=-0.866°, CL=-0.1356435312 ✓
- Instances 6 & 8: (2,1) & (1,2) → α=0.000°, CL=0.0474876211 ✓

**Why these match**: Symmetric pairs with swapped (θ₁, θ₂) coordinates:
- θ₁ and θ₂ swapped
- sin(θ₁)+sin(θ₂) unchanged → same α
- cos(θ₁)+cos(θ₂) unchanged → same v_grid
- Therefore same CL!

### Instances that should NOT match (same α, different v_grid):

- Instance 1 vs 6,8: All have α=0° but different grid velocities ✓

## Verification Status

✅ **Torus implementation is CORRECT**

✅ **All physics expectations met**:
1. Symmetric instances match perfectly (machine precision)
2. Instances with same α but different (θ₁,θ₂) correctly show different CL
3. Grid velocity computation working as expected

❌ **Test script needs minor fix**: Regex extracts wrong column (cosmetic issue)

## Conclusions

1. **No bug in torus solver** - physics is correct
2. **Degenerate case does NOT collapse to 1D** - it's a distinct 2D manifold
3. **Test expectations were wrong** - same α does NOT imply same CL in 2D torus
4. **Symmetry pairs verified** - instances with swapped (θ₁, θ₂) match perfectly

## Recommendations

1. **Update test script** to:
   - Fix regex to extract correct CL column
   - Change validation logic to check symmetric pairs, not instances with same α
   - Document that degenerate torus ≠ 1D time spectral

2. **Update documentation** to clarify:
   - Degenerate torus has same α range but different velocity distribution
   - Grid velocity depends on path through (θ₁,θ₂) space, not just instantaneous α
   - Expect symmetric pairs to match, not all instances with same α

---

**Analysis Date**: 2025-11-14
**Status**: Torus solver physics VERIFIED, test expectations corrected
