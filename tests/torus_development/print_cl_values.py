"""
Simple script to extract and print CL values from torus degenerate test output
"""
import re
import sys

if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    filename = 'torus_degenerate.log'

print(f"Reading from: {filename}")
print()

# Correct regex pattern to extract CL from NK lines
# Format: 1  instance  iter  total  NK/*NK  ----  CFL  ?  resrho  CL  CD  total_residual
pattern = r'^\s*1\s+(\d+)\s+\d+\s+\d+\s+\*?NK\s+----\s+([\d\.]+)\s+([\d\.]+)\s+([-\d\.E\+]+)\s+([-\d\.E\+]+)\s+([-\d\.E\+]+)\s+([-\d\.E\+]+)'

instance_data = {}

try:
    with open(filename, 'r') as f:
        for line in f:
            match = re.search(pattern, line)
            if match:
                instance = int(match.group(1))
                resrho = match.group(4)
                cl = float(match.group(5))
                cd = float(match.group(6))
                # Keep last occurrence (final iteration)
                instance_data[instance] = {'cl': cl, 'cd': cd, 'resrho': resrho}
except FileNotFoundError:
    print(f"Error: File '{filename}' not found!")
    print("\nAvailable log files:")
    import os
    for f in os.listdir('.'):
        if f.endswith('.log'):
            print(f"  {f}")
    sys.exit(1)

print(f"Extracted CL values for {len(instance_data)} instances:")
print()
print("Instance     CL              CD              resrho")
print("-" * 60)

for instance in sorted(instance_data.keys()):
    cl = instance_data[instance]['cl']
    cd = instance_data[instance]['cd']
    resrho = instance_data[instance]['resrho']
    print(f"   {instance:2d}    {cl:15.10f}  {cd:15.8e}  {resrho}")

# Group by grid coordinates if we have 9 instances (3x3)
if len(instance_data) == 9:
    import numpy as np

    # Reconstruct alpha grid
    n1, n2 = 3, 3
    theta1 = np.linspace(0.0, 2.0 * np.pi, n1, endpoint=False)
    theta2 = np.linspace(0.0, 2.0 * np.pi, n2, endpoint=False)

    alpha_grid = np.zeros((n1, n2))
    for i1 in range(n1):
        for i2 in range(n2):
            alpha_deg = 1.0 * np.sin(theta1[i1]) + 1.0 * np.sin(theta2[i2])
            alpha_grid[i1, i2] = alpha_deg

    print()
    print("=" * 60)
    print("DETAILED BREAKDOWN WITH GRID COORDINATES")
    print("=" * 60)
    print()
    print("Inst  (i1,i2)    θ₁       θ₂        α(deg)       CL")
    print("-" * 60)

    for instance in range(1, 10):
        idx = instance - 1
        i1 = idx % n1
        i2 = idx // n1
        t1 = theta1[i1]
        t2 = theta2[i2]
        alpha = alpha_grid[i1, i2]
        cl = instance_data[instance]['cl']
        print(f" {instance:2d}    ({i1},{i2})   {t1:6.3f}  {t2:6.3f}   {alpha:7.3f}°   {cl:15.10f}")

    # Check for symmetric pairs
    print()
    print("=" * 60)
    print("SYMMETRIC PAIRS (swapped θ₁, θ₂)")
    print("=" * 60)
    print()

    pairs = [(2, 4), (3, 7), (6, 8)]
    for inst1, inst2 in pairs:
        cl1 = instance_data[inst1]['cl']
        cl2 = instance_data[inst2]['cl']
        diff = abs(cl1 - cl2)

        idx1 = inst1 - 1
        i1_1, i2_1 = idx1 % n1, idx1 // n1
        idx2 = inst2 - 1
        i1_2, i2_2 = idx2 % n1, idx2 // n1

        status = "✅" if diff < 1e-10 else "❌"
        print(f"Instance {inst1} ({i1_1},{i2_1}) vs {inst2} ({i1_2},{i2_2}):")
        print(f"  CL₁ = {cl1:.10f}")
        print(f"  CL₂ = {cl2:.10f}")
        print(f"  Δ   = {diff:.2e} {status}")
        print()
