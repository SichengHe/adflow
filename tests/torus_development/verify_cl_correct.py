import re

# Read the log file
with open('torus_verify.log', 'r') as f:
    lines = f.readlines()

# Extract final NK iteration for each instance
final_values = {}
pattern = r'^\s*1\s+(\d+)\s+\d+\s+\d+\s+NK\s+----\s+([\d\.]+)\s+([\d\.]+)\s+([-\d\.E\+]+)\s+([-\d\.E\+]+)\s+([-\d\.E\+]+)\s+([-\d\.E\+]+)'

for line in lines:
    match = re.search(pattern, line)
    if match:
        instance = int(match.group(1))
        resrho = float(match.group(4))
        cl = float(match.group(5))
        cd = float(match.group(6))
        residual = float(match.group(7))
        final_values[instance] = {'cl': cl, 'cd': cd, 'resrho': resrho, 'residual': residual}

# Alpha values for each instance (from the test)
import numpy as np
n1, n2 = 3, 3
A1, A2 = 1.0, 1.0
theta1 = np.linspace(0.0, 2.0 * np.pi, n1, endpoint=False)
theta2 = np.linspace(0.0, 2.0 * np.pi, n2, endpoint=False)

alpha_grid = {}
for i2 in range(n2):
    for i1 in range(n1):
        instance = i2 * n1 + i1 + 1
        alpha_deg = A1 * np.sin(theta1[i1]) + A2 * np.sin(theta2[i2])
        alpha_grid[instance] = (i1, i2, alpha_deg)

# Print results
print("="*90)
print("TORUS TIME SPECTRAL - CONVERGED CL VALUES")
print("="*90)
print(f"\nFinal converged values ({len(final_values)} instances):\n")
print(f"{'Inst':<6} {'(i1,i2)':<8} {'α (deg)':<10} {'CL':<18} {'CD':<15} {'Residual':<12}")
print("-"*90)

for inst in sorted(final_values.keys()):
    i1, i2, alpha = alpha_grid[inst]
    cl = final_values[inst]['cl']
    cd = final_values[inst]['cd']
    res = final_values[inst]['residual']
    print(f"{inst:<6} ({i1},{i2}){'':<5} {alpha:>8.3f}   {cl:>16.10f}  {cd:>13.6e}  {res:>10.2e}")

# Group by alpha and check matching
print("\n" + "="*90)
print("DEGENERATE CASE VERIFICATION: Instances with Same α Should Have Identical CL")
print("="*90)

alpha_groups = {}
for inst, (i1, i2, alpha) in alpha_grid.items():
    alpha_key = round(alpha, 4)
    if alpha_key not in alpha_groups:
        alpha_groups[alpha_key] = []
    alpha_groups[alpha_key].append(inst)

print()
all_match = True
for alpha_key in sorted(alpha_groups.keys()):
    instances = alpha_groups[alpha_key]
    print(f"α = {alpha_key:>7.3f}°  →  Instances: {instances}")
    
    if len(instances) > 1:
        cl_values = [final_values[inst]['cl'] for inst in instances]
        cl_diff = max(cl_values) - min(cl_values)
        print(f"{'':20}  CL values: {[f'{cl:.10f}' for cl in cl_values]}")
        print(f"{'':20}  CL range:  {cl_diff:.2e}", end="")
        
        if cl_diff < 1e-10:  # Machine precision
            print("  ✅ MATCH (machine precision)")
        elif cl_diff < 1e-6:
            print("  ✅ MATCH")
        else:
            print("  ❌ MISMATCH")
            all_match = False
    print()

print("="*90)
if all_match:
    print("✅ SUCCESS: All instances with same α have identical CL!")
    print("   Torus correctly collapses to 1D behavior when A1=A2 and omega1=omega2")
else:
    print("❌ FAILURE: Some instances have mismatched CL")
print("="*90)

# Show mean CL
mean_cl = sum(final_values[i]['cl'] for i in final_values) / len(final_values)
print(f"\nMean CL across all instances: {mean_cl:.8f}")
print("="*90)
