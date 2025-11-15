"""
Analyze what periodic solution the rational torus (omega1=omega2) represents
and how the 3x3 torus grid maps to time samples.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

print("="*70)
print("RATIONAL TORUS PERIODICITY ANALYSIS")
print("="*70)

# Torus parameters
n1 = 3
n2 = 3
omega1 = 100.0
omega2 = 100.0
A1 = 1.0
A2 = 1.0

# Period of the combined motion
period = 2.0 * np.pi / omega1  # Since omega1 = omega2

print(f"\nTorus grid: {n1} × {n2}")
print(f"ω₁ = {omega1} rad/s")
print(f"ω₂ = {omega2} rad/s")
print(f"Ratio ω₁/ω₂ = {omega1/omega2}")
print(f"\nMotion period T = 2π/ω = {period:.6f} s")

# Generate theta grid
theta1 = np.linspace(0.0, 2.0 * np.pi, n1, endpoint=False)
theta2 = np.linspace(0.0, 2.0 * np.pi, n2, endpoint=False)

print("\n" + "="*70)
print("TORUS GRID → TIME MAPPING")
print("="*70)

print("\nFor rational ratio ω₁/ω₂ = 1, each (θ₁, θ₂) corresponds to time:")
print("  θ₁ = ω₁·t  →  t = θ₁/ω₁")
print("  θ₂ = ω₂·t  →  t = θ₂/ω₂")
print("\nBut these must be CONSISTENT for a periodic solution!")
print("Inconsistency means the grid point is NOT on the periodic trajectory.\n")

# For each torus grid point, check if it's on the periodic orbit
print("Instance  (i1,i2)    θ₁       θ₂        t₁(s)     t₂(s)     Match?   α(deg)")
print("-" * 90)

consistent_points = []
inconsistent_points = []

for i2 in range(n2):
    for i1 in range(n1):
        idx = i2 * n1 + i1

        t1 = theta1[i1] / omega1
        t2 = theta2[i2] / omega2

        # Check if times are consistent (modulo period)
        t_diff = abs((t1 - t2) % period)
        is_consistent = (t_diff < 1e-10) or (abs(t_diff - period) < 1e-10)

        alpha_deg = A1 * np.sin(theta1[i1]) + A2 * np.sin(theta2[i2])

        match_str = "✅ YES" if is_consistent else "❌ NO "

        print(f"   {idx+1:2d}     ({i1},{i2})   {theta1[i1]:6.3f}  {theta2[i2]:6.3f}   "
              f"{t1:8.5f}  {t2:8.5f}   {match_str}   {alpha_deg:7.3f}°")

        if is_consistent:
            consistent_points.append((idx+1, i1, i2, t1, alpha_deg))
        else:
            inconsistent_points.append((idx+1, i1, i2, t1, t2, alpha_deg))

print("\n" + "="*70)
print("CONSISTENT POINTS (On Periodic Trajectory)")
print("="*70)

if consistent_points:
    print(f"\n{len(consistent_points)} points lie on the 1D periodic orbit:")
    print("\nInstance  (i1,i2)    Time(s)     α(deg)      Phase(rad)")
    print("-" * 60)
    for inst, i1, i2, t, alpha in consistent_points:
        phase = omega1 * t
        print(f"   {inst:2d}     ({i1},{i2})   {t:8.5f}    {alpha:7.3f}°    {phase:6.3f}")

    # Compare to standard 1D time spectral sampling
    print("\n" + "="*70)
    print("COMPARISON TO 1D TIME SPECTRAL")
    print("="*70)

    n_1d = len(consistent_points)
    theta_1d = np.linspace(0, 2*np.pi, n_1d, endpoint=False)

    print(f"\nStandard 1D TS with {n_1d} instances would sample at:")
    print("Instance    θ        t(s)       α(deg)")
    print("-" * 50)
    for i, theta in enumerate(theta_1d):
        t = theta / omega1
        alpha = 2 * A1 * np.sin(theta)  # Combined amplitude
        print(f"   {i+1:2d}     {theta:6.3f}   {t:8.5f}    {alpha:7.3f}°")
else:
    print("\n⚠️  NO points lie on the periodic orbit!")
    print("   The 3×3 torus grid doesn't sample the diagonal θ₁=θ₂")

print("\n" + "="*70)
print("INCONSISTENT POINTS (Off Periodic Trajectory)")
print("="*70)

if inconsistent_points:
    print(f"\n{len(inconsistent_points)} points are OFF the 1D periodic orbit:")
    print("\nThese points represent quasi-periodic motion, not simple periodic!")
    print("\nInstance  (i1,i2)    t₁(s)     t₂(s)     Δt(s)     α(deg)")
    print("-" * 70)
    for inst, i1, i2, t1, t2, alpha in inconsistent_points:
        dt = abs(t1 - t2)
        print(f"   {inst:2d}     ({i1},{i2})   {t1:8.5f}  {t2:8.5f}  {dt:8.5f}   {alpha:7.3f}°")

# Visualize the trajectory
print("\n" + "="*70)
print("GENERATING VISUALIZATION")
print("="*70)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Plot 1: (θ₁, θ₂) space with torus grid
theta1_fine = np.linspace(0, 2*np.pi, 100)
theta2_fine = np.linspace(0, 2*np.pi, 100)

# Periodic diagonal: θ₁ = θ₂
ax1.plot(theta1_fine, theta1_fine, 'r-', linewidth=2, label='Periodic orbit (θ₁=θ₂)')
ax1.plot(theta1_fine, theta1_fine + 2*np.pi, 'r--', linewidth=1, alpha=0.5)
ax1.plot(theta1_fine, theta1_fine - 2*np.pi, 'r--', linewidth=1, alpha=0.5)

# Plot torus grid points
for i2 in range(n2):
    for i1 in range(n1):
        idx = i2 * n1 + i1 + 1
        is_on_orbit = any(inst == idx for inst, _, _, _, _ in consistent_points)

        color = 'green' if is_on_orbit else 'blue'
        marker = 'o' if is_on_orbit else 'x'
        size = 100 if is_on_orbit else 50

        ax1.scatter(theta1[i1], theta2[i2], c=color, marker=marker, s=size,
                   edgecolors='black', linewidth=1.5, zorder=10)
        ax1.text(theta1[i1]+0.1, theta2[i2]+0.1, str(idx), fontsize=8)

ax1.set_xlabel('θ₁ (rad)', fontsize=12)
ax1.set_ylabel('θ₂ (rad)', fontsize=12)
ax1.set_title('Torus Grid in (θ₁, θ₂) Space\n(Green = on periodic orbit, Blue = off orbit)', fontsize=12)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(-0.5, 2*np.pi + 0.5)
ax1.set_ylim(-0.5, 2*np.pi + 0.5)
ax1.set_aspect('equal')
ax1.legend()

# Plot 2: α(t) trajectory
t_fine = np.linspace(0, period, 1000)
alpha_periodic = 2 * A1 * np.sin(omega1 * t_fine)

ax2.plot(t_fine, alpha_periodic, 'r-', linewidth=2, label='True periodic: α=2A·sin(ωt)')

# Plot torus grid points
for inst, i1, i2, t, alpha in consistent_points:
    ax2.scatter(t, alpha, c='green', marker='o', s=100, edgecolors='black',
               linewidth=1.5, zorder=10, label='On orbit' if inst == consistent_points[0][0] else '')

for inst, i1, i2, t1, t2, alpha in inconsistent_points:
    # Plot at both t1 and t2 to show inconsistency
    alpha_at_t1 = 2 * A1 * np.sin(omega1 * t1)
    ax2.scatter(t1, alpha, c='blue', marker='x', s=100, linewidth=2,
               zorder=10, label='Off orbit' if inst == inconsistent_points[0][0] else '')
    ax2.plot([t1, t1], [alpha, alpha_at_t1], 'b--', alpha=0.5, linewidth=1)

ax2.set_xlabel('Time (s)', fontsize=12)
ax2.set_ylabel('α (degrees)', fontsize=12)
ax2.set_title('Angle of Attack vs Time\n(Torus samples vs True periodic)', fontsize=12)
ax2.grid(True, alpha=0.3)
ax2.legend()

plt.tight_layout()
plt.savefig('torus_periodicity_analysis.png', dpi=150, bbox_inches='tight')
print("\n✅ Saved visualization to: torus_periodicity_analysis.png")

# Summary
print("\n" + "="*70)
print("CONCLUSION")
print("="*70)

if len(consistent_points) == 0:
    print("""
⚠️  The 3×3 torus grid does NOT sample the periodic orbit!

For ω₁ = ω₂, the periodic trajectory is θ₁ = θ₂ (the diagonal).
The 3×3 uniform grid only hits this diagonal at:
  - (0, 0)
  - (2π/3, 2π/3)  ← NOT in a 3-point grid with endpoint=False!
  - (4π/3, 4π/3)  ← NOT in a 3-point grid with endpoint=False!

So NONE of the 9 grid points (except possibly (0,0)) lie on the
periodic orbit!

This is why the torus solution doesn't match 1D time spectral -
it's sampling OFF the periodic trajectory into the quasi-periodic
region of the 2D torus.
""")
else:
    print(f"""
✅ {len(consistent_points)} out of 9 torus grid points lie on the periodic orbit.

These points should match a {len(consistent_points)}-instance 1D time spectral
solution with amplitude 2A = {2*A1}° at the same time samples.

The remaining {len(inconsistent_points)} points sample the quasi-periodic
region and will have different aerodynamics.
""")

print("="*70)
