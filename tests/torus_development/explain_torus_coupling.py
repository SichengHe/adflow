"""
Explain how the torus spectral coupling works and why diagonal instances
don't form an independent 3x3 system.
"""
import numpy as np

print("="*70)
print("TORUS SPECTRAL COUPLING: Why Diagonal Instances Don't Decouple")
print("="*70)

# The 3x3 grid layout in (i1, i2) space
print("\nTorus 3×3 grid layout (instance numbering):")
print("     i2=0    i2=1    i2=2")
print("i1=0:  1       4       7")
print("i1=1:  2       5       8")  
print("i1=2:  3       6       9")
print()
print("Diagonal instances: 1 (0,0), 5 (1,1), 9 (2,2)")
print()

# The torus stencil is a 2D cross
print("2D Spectral Differentiation Stencil (cross pattern):")
print()
print("For instance at (i1, i2), the derivative couples with:")
print("  - ∂/∂θ₁: instances at (i1-1, i2), (i1, i2), (i1+1, i2)")
print("  - ∂/∂θ₂: instances at (i1, i2-1), (i1, i2), (i1, i2+1)")
print()
print("Combined: 5-point cross stencil (periodic BC)")
print()

# Show coupling for diagonal instances
print("Coupling for diagonal instances:")
print()
print("Instance 1 (0,0) couples with:")
print("  ∂/∂θ₁: (2,0)=2, (0,0)=1, (1,0)=4  [vertical neighbors]")
print("  ∂/∂θ₂: (0,2)=7, (0,0)=1, (0,1)=4  [horizontal neighbors]")
print("  → Couples with: 1, 2, 3, 6, 7")
print("  → Does NOT couple with other diagonal instances 5, 9!")
print()

print("Instance 5 (1,1) couples with:")
print("  ∂/∂θ₁: (0,1)=4, (1,1)=5, (2,1)=6  [vertical]")
print("  ∂/∂θ₂: (1,0)=2, (1,1)=5, (1,2)=8  [horizontal]")
print("  → Couples with: 2, 4, 5, 6, 8")
print("  → Does NOT couple with other diagonal instances 1, 9!")
print()

print("Instance 9 (2,2) couples with:")
print("  ∂/∂θ₁: (1,2)=8, (2,2)=9, (0,2)=7  [vertical]")
print("  ∂/∂θ₂: (2,1)=6, (2,2)=9, (2,0)=3  [horizontal]")
print("  → Couples with: 3, 6, 7, 8, 9")
print("  → Does NOT couple with other diagonal instances 1, 5!")
print()

print("="*70)
print("CONCLUSION")
print("="*70)
print("""
The diagonal instances (1, 5, 9) do NOT form an independent 3×3 system!

They are coupled through the full 9×9 2D spectral operator:
  1 ↔ 2, 3, 6, 7 ↔ 5 ↔ 4, 6, 8 ↔ 9

To compute time derivatives at diagonal instances, the torus solver
uses ALL 9 instances, not just the 3 diagonal ones.

This is why torus diagonal ≠ 1D time spectral:
  - Torus: 9-instance 2D spectral operator (even for diagonal)
  - 1D TS: 3-instance 1D spectral operator

The diagonal instances lie on a 1D periodic orbit physically, but
the torus solver doesn't recognize this and treats them as part of
a full 2D quasi-periodic system.
""")
