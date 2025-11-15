"""
Visualize torus time spectral vs time accurate comparison results

This script demonstrates spectral interpolation of torus TS solution
into the time domain for direct comparison with time-accurate results.
"""
import numpy as np
import matplotlib.pyplot as plt
import sys

# Check if ta_results.npz exists (from actual run, not manual)
try:
    ta_data = np.load('ta_results.npz')
    print("Loaded TA results from ta_results.npz")
except FileNotFoundError:
    print("Error: ta_results.npz not found. Run time-accurate solver first:")
    print("  mpirun -np 16 python compare_torus_ts_vs_ta.py ta")
    sys.exit(1)

time_ta = ta_data['time_ta']
cl_ta = ta_data['cl_ta']
omega1 = ta_data['omega1']
omega2 = ta_data['omega2']

print("="*70)
print("VISUALIZATION: Torus TS vs Time Accurate")
print("="*70)

print(f"\nTime-Accurate Results:")
print(f"  Frequencies: ω₁={omega1:.4f}, ω₂={omega2:.4f}")
print(f"  Ratio: ω₁/ω₂ = {omega1/omega2:.6f} = 1/√2")
print(f"  Time points: {len(time_ta)}")
print(f"  Time range: [{time_ta[0]:.4f}, {time_ta[-1]:.4f}]")
print(f"  CL range: [{np.min(cl_ta):.6f}, {np.max(cl_ta):.6f}]")

# Skip initial transient (first 20% of time)
transient_idx = int(0.2 * len(time_ta))
time_ta_steady = time_ta[transient_idx:]
cl_ta_steady = cl_ta[transient_idx:]
print(f"\nAfter removing transient (first 20%):")
print(f"  Steady-state points: {len(time_ta_steady)}")
print(f"  Time range: [{time_ta_steady[0]:.4f}, {time_ta_steady[-1]:.4f}]")

# Load torus results
try:
    torus_data = np.load('torus_results.npz')
except FileNotFoundError:
    print("\nError: torus_results.npz not found. Run torus solver first:")
    print("  mpirun -np 4 python compare_torus_ts_vs_ta.py torus")
    sys.exit(1)

print("\nTorus results loaded! Creating comparison plot...")

cl_grid = torus_data['cl_grid']
theta1_torus = torus_data['theta1']
theta2_torus = torus_data['theta2']
n1, n2 = cl_grid.shape

print(f"  Torus grid: {n1}×{n2} = {n1*n2} instances")
print(f"  CL range: [{np.min(cl_grid):.6f}, {np.max(cl_grid):.6f}]")

# Import spectral interpolation and phase alignment
from compare_torus_ts_vs_ta import spectral_interp, phase_alignment

# Convert TA time points to phase points (use steady state only)
# theta1 = (omega1 * t) mod 2π, theta2 = (omega2 * t) mod 2π
# This maps time domain to torus phase space
print("\n" + "="*70)
print("SPECTRAL INTERPOLATION & PHASE ALIGNMENT")
print("="*70)
print(f"\nMapping time to torus phase space:")
print(f"  θ₁(t) = (ω₁ · t) mod 2π,  where ω₁ = {omega1:.4f} rad/s")
print(f"  θ₂(t) = (ω₂ · t) mod 2π,  where ω₂ = {omega2:.4f} rad/s")

# Compute phase values for steady-state TA solution
theta1_ta_steady = (omega1 * time_ta_steady) % (2 * np.pi)
theta2_ta_steady = (omega2 * time_ta_steady) % (2 * np.pi)

# Find optimal phase shifts to align TA with torus solution
print("\nComputing optimal phase shifts...")
phi1, phi2 = phase_alignment(cl_grid, theta1_ta_steady, theta2_ta_steady, cl_ta_steady)
print(f"  Optimal phase shifts: φ₁ = {phi1:.6f} rad ({phi1*180/np.pi:.2f}°)")
print(f"                        φ₂ = {phi2:.6f} rad ({phi2*180/np.pi:.2f}°)")

# Apply phase shifts
theta1_aligned = (theta1_ta_steady + phi1) % (2 * np.pi)
theta2_aligned = (theta2_ta_steady + phi2) % (2 * np.pi)

# Generate dense time grid for smooth torus curve
n_interp = 1000  # Many points for smooth visualization
t_min, t_max = time_ta_steady[0], time_ta_steady[-1]
time_smooth = np.linspace(t_min, t_max, n_interp)
theta1_smooth = (omega1 * time_smooth + phi1) % (2 * np.pi)
theta2_smooth = (omega2 * time_smooth + phi2) % (2 * np.pi)

# Interpolate torus solution on dense grid for smooth curve
cl_torus_smooth_raw = spectral_interp(cl_grid, theta1_smooth, theta2_smooth)

# Also interpolate at TA time points for direct comparison
cl_torus_at_ta_raw = spectral_interp(cl_grid, theta1_aligned, theta2_aligned)

# Add constant offset to align mean values (like paper example)
# This accounts for any constant shift between solutions
C_star = np.mean(cl_ta_steady - cl_torus_at_ta_raw)
cl_torus_smooth = cl_torus_smooth_raw + C_star
cl_torus_at_ta = cl_torus_at_ta_raw + C_star

print(f"\nSpectral interpolation completed:")
print(f"  Smooth curve: {n_interp} points for visualization")
print(f"  TA points: {len(time_ta_steady)} points for comparison (all plotted)")
print(f"  Mean offset C* = {C_star:.6e} (applied to align means)")
print(f"  CL_torus range: [{np.min(cl_torus_smooth):.6f}, {np.max(cl_torus_smooth):.6f}]")

# Create comparison figure
fig, axes = plt.subplots(2, 1, figsize=(14, 10))

# Plot 1: CL comparison with ALL TA points and smooth torus curve
ax = axes[0]
# Plot ALL TA points (no decimation)
ax.plot(time_ta_steady, cl_ta_steady, '.', markersize=4, label=f'Time Accurate ({len(time_ta_steady)} points)', alpha=0.7, zorder=3)
ax.plot(time_smooth, cl_torus_smooth, '-', linewidth=2, label='Torus TS (spectral interp, phase-aligned)', alpha=0.9, zorder=2)
ax.set_xlabel('Time (s)', fontsize=12)
ax.set_ylabel('CL', fontsize=12)
ax.set_title('Comparison: Torus Time Spectral vs Time Accurate (Phase-Aligned)', fontsize=14)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=11)

# Plot 2: Error for all points
ax = axes[1]
error = np.abs(cl_ta_steady - cl_torus_at_ta)
ax.semilogy(time_ta_steady, error, '.', markersize=4, color='red', alpha=0.7)
ax.set_xlabel('Time (s)', fontsize=12)
ax.set_ylabel('|CL_TA - CL_Torus|', fontsize=12)
ax.set_title('Absolute Error (After Phase Alignment)', fontsize=14)
ax.grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.savefig('torus_vs_ta_comparison.png', dpi=300, bbox_inches='tight')
print("\nSaved comparison to torus_vs_ta_comparison.png")

print(f"\nError Statistics (Phase-Aligned):")
print(f"  Max error: {np.max(error):.6e}")
print(f"  Mean error: {np.mean(error):.6e}")
print(f"  RMS error: {np.sqrt(np.mean(error**2)):.6e}")
print(f"  Relative L2 error: {np.linalg.norm(error)/np.linalg.norm(cl_ta_steady):.6e}")

# ========== NEW FIGURE: Longer time period with PSD ==========
print("\n" + "="*70)
print("QUASI-PERIODICITY DEMONSTRATION")
print("="*70)

# Generate longer time series from torus solution
n_long = 5000
t_long_duration = 100.0  # Much longer time period
time_long = np.linspace(0, t_long_duration, n_long)
theta1_long = (omega1 * time_long) % (2 * np.pi)
theta2_long = (omega2 * time_long) % (2 * np.pi)

# Interpolate torus solution for long time series
cl_torus_long = spectral_interp(cl_grid, theta1_long, theta2_long)

print(f"\nLong-duration TS signal:")
print(f"  Duration: {t_long_duration:.1f} s")
print(f"  Time points: {n_long}")
print(f"  Sampling rate: {n_long/t_long_duration:.2f} Hz")
print(f"  CL range: [{np.min(cl_torus_long):.6f}, {np.max(cl_torus_long):.6f}]")

# Compute Power Spectral Density
print("\nComputing Power Spectral Density...")
from scipy import signal

# Detrend and window the signal
cl_detrended = signal.detrend(cl_torus_long)
window = signal.windows.hann(n_long)
cl_windowed = cl_detrended * window

# Compute PSD using Welch's method for smoother spectrum
dt = time_long[1] - time_long[0]
freqs, psd = signal.welch(cl_torus_long, fs=1.0/dt, nperseg=min(1024, n_long//4),
                          scaling='density', detrend='linear')

print(f"  Frequency resolution: {freqs[1] - freqs[0]:.6f} Hz")
print(f"  Max frequency: {freqs[-1]:.6f} Hz")

# Fundamental frequencies in Hz
f1 = omega1 / (2 * np.pi)
f2 = omega2 / (2 * np.pi)

print(f"\nFundamental frequencies:")
print(f"  ω₁ = {omega1:.6f} rad/s → f₁ = ω₁/(2π) = {f1:.6f} Hz")
print(f"  ω₂ = {omega2:.6f} rad/s → f₂ = ω₂/(2π) = {f2:.6f} Hz")

# Find peaks in PSD
peak_indices, properties = signal.find_peaks(psd, height=np.max(psd)*0.01, distance=5)
peak_freqs = freqs[peak_indices]
peak_heights = psd[peak_indices]

# Sort peaks by height
sorted_idx = np.argsort(peak_heights)[::-1]
top_peaks = sorted_idx[:min(15, len(sorted_idx))]

print(f"\nTop frequency peaks detected in PSD:")

# Generate expected combination frequencies (linear combinations of f1 and f2)
tolerance = 0.005  # Hz
expected_combos = []
for m in range(-3, 4):
    for n in range(-3, 4):
        if m == 0 and n == 0:
            continue
        f_combo = abs(m * f1 + n * f2)
        if 0 < f_combo < max(3*f1, 3*f2):  # Only frequencies in our plot range
            expected_combos.append((m, n, f_combo))

for i, idx in enumerate(top_peaks[:10]):
    freq_val = peak_freqs[idx]
    height_val = peak_heights[idx]
    print(f"  Peak {i+1}: f = {freq_val:.6f} Hz, Power = {height_val:.2e}")

    # Try to match with expected combinations
    matched = False
    for m, n, f_combo in expected_combos:
        if abs(freq_val - f_combo) < tolerance:
            if m == 1 and n == 0:
                print(f"          → f₁ = {f1:.6f} Hz")
            elif m == 0 and n == 1:
                print(f"          → f₂ = {f2:.6f} Hz")
            elif m == -1 and n == 0:
                print(f"          → -f₁ = {f1:.6f} Hz")
            elif m == 0 and n == -1:
                print(f"          → -f₂ = {f2:.6f} Hz")
            else:
                sign_m = "+" if m >= 0 else ""
                sign_n = "+" if n >= 0 else ""
                print(f"          → {sign_m}{m}·f₁ {sign_n}{n}·f₂ = {f_combo:.6f} Hz")
            matched = True
            break
    if not matched:
        print(f"          → (unidentified combination)")

# Create two-panel figure: Time series (left) and PSD (right)
fig2, axes2 = plt.subplots(1, 2, figsize=(16, 6))

# Left panel: Long time series
ax_time = axes2[0]
ax_time.plot(time_long, cl_torus_long, '-', linewidth=0.8, color='#1f77b4', alpha=0.8)
ax_time.set_xlabel('Time (s)', fontsize=12)
ax_time.set_ylabel('CL', fontsize=12)
ax_time.set_title(f'Torus TS Solution (Extended Duration: {t_long_duration:.0f}s)', fontsize=13)
ax_time.grid(True, alpha=0.3)
ax_time.text(0.02, 0.98, f'ω₁ = {omega1:.4f} rad/s\nω₂ = {omega2:.4f} rad/s\nω₁/ω₂ = 1/√2 (irrational)\n→ Quasi-periodic',
             transform=ax_time.transAxes, fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

# Right panel: Power Spectral Density
ax_psd = axes2[1]
ax_psd.semilogy(freqs, psd, '-', linewidth=1.5, color='#ff7f0e', alpha=0.8, label='PSD')

# Mark the fundamental frequencies
ax_psd.axvline(f1, color='red', linestyle='--', linewidth=2.0, alpha=0.8,
               label=f'f₁ (ω₁=1.0 rad/s)')
ax_psd.axvline(f2, color='green', linestyle='--', linewidth=2.0, alpha=0.8,
               label=f'f₂ (ω₂=√2 rad/s)')

# Add text annotations near the peaks with ω values
y_text = ax_psd.get_ylim()[1] * 0.3
ax_psd.text(f1, y_text, f'ω₁ = 1.0\nf₁ = {f1:.4f} Hz',
            fontsize=10, ha='center', va='bottom', color='red',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='red', alpha=0.8))
ax_psd.text(f2, y_text, f'ω₂ = √2\nf₂ = {f2:.4f} Hz',
            fontsize=10, ha='center', va='bottom', color='green',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='green', alpha=0.8))

# Mark some combination frequencies
combo_freqs = [
    (2*f1, '2f₁', 'blue'),
    (2*f2, '2f₂', 'cyan'),
    (f1+f2, 'f₁+f₂', 'magenta'),
    (abs(f2-f1), '|f₂-f₁|', 'orange'),
]
for freq, label, color in combo_freqs:
    if freq < max(3*f1, 3*f2):
        ax_psd.axvline(freq, color=color, linestyle=':', linewidth=1.2, alpha=0.5)
        # Add text annotation at the top
        y_pos = ax_psd.get_ylim()[1] * 0.5
        ax_psd.text(freq, y_pos, label, rotation=90, fontsize=8,
                   verticalalignment='bottom', color=color, alpha=0.7)

# Mark detected peaks
ax_psd.plot(peak_freqs, peak_heights, 'x', color='purple', markersize=8, markeredgewidth=2,
            label=f'{len(peak_freqs)} peaks', alpha=0.8)

ax_psd.set_xlabel('Frequency (Hz)', fontsize=12)
ax_psd.set_ylabel('Power Spectral Density', fontsize=12)
ax_psd.set_title('Power Spectral Density (Quasi-Periodic Signal)', fontsize=13)
ax_psd.grid(True, alpha=0.3, which='both')
ax_psd.legend(fontsize=9, loc='upper right')
ax_psd.set_xlim([0, max(3*f1, 3*f2)])  # Focus on relevant frequency range

# Add text box explaining quasi-periodicity
textstr = f'Quasi-periodic signal:\nf₁/f₂ = ω₁/ω₂ = 1/√2 (irrational)\n→ Dense spectrum with peaks\n   at m·f₁ + n·f₂'
ax_psd.text(0.98, 0.05, textstr, transform=ax_psd.transAxes, fontsize=9,
            verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

plt.tight_layout()
plt.savefig('torus_quasiperiodic_psd.png', dpi=300, bbox_inches='tight')
print("\nSaved quasi-periodicity demonstration to torus_quasiperiodic_psd.png")

print("\n" + "="*70)
print("QUASI-PERIODICITY EXPLANATION")
print("="*70)
print(f"The frequency ratio ω₁/ω₂ = 1/√2 ≈ 0.707107 is IRRATIONAL.")
print(f"This means the signal NEVER exactly repeats (no fundamental period).")
print(f"The PSD shows multiple peaks at frequencies f₁, f₂, and their combinations:")
print(f"  - Linear combinations: m·f₁ + n·f₂ (m, n = integers)")
print(f"  - These peaks create a dense spectrum characteristic of quasi-periodicity")
print(f"  - Unlike a periodic signal (which has discrete harmonics of ONE frequency)")
print("="*70)

plt.show()
