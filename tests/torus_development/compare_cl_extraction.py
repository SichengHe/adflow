import re

# Read the log file
with open('torus_full_output.log', 'r') as f:
    content = f.read()

print("="*90)
print("CL EXTRACTION COMPARISON")
print("="*90)

# Method 1: Test script's regex (what the test uses)
print("\n1. TEST SCRIPT REGEX (from test_torus_degenerate.py line 191):")
print("   Pattern: r'^\\s*1\\s+(\\d+)\\s+\\d+\\s+\\d+\\s+NK\\s+[-\\w]+\\s+[\\d\\.]+\\s+[\\d\\.]+\\s+([\\d\\.\\-E\\+]+)\\s+([\\d\\.\\-E\\+]+)\\s+([\\d\\.\\-E\\+]+)'")
print()

pattern1 = r'^\s*1\s+(\d+)\s+\d+\s+\d+\s+NK\s+[-\w]+\s+[\d\.]+\s+[\d\.]+\s+([\d\.\-E\+]+)\s+([\d\.\-E\+]+)\s+([\d\.\-E\+]+)'
test_cls = {}
for line in content.split('\n'):
    match = re.search(pattern1, line)
    if match:
        instance = int(match.group(1))
        cl_val = float(match.group(2))
        cd_val = float(match.group(3))
        test_cls[instance] = cl_val

print(f"   Matched {len(test_cls)} instances")
if test_cls:
    for inst in sorted(test_cls.keys())[:5]:
        print(f"   Instance {inst}: CL = {test_cls[inst]}")
else:
    print("   NO MATCHES FOUND!")

# Method 2: Correct regex (actual column positions)
print("\n2. CORRECT REGEX (based on actual column mapping):")
print("   Pattern: Accounts for NK line format with proper column positions")
print()

pattern2 = r'^\s*1\s+(\d+)\s+\d+\s+\d+\s+NK\s+----\s+([\d\.]+)\s+([\d\.]+)\s+([-\d\.E\+]+)\s+([-\d\.E\+]+)\s+([-\d\.E\+]+)\s+([-\d\.E\+]+)'
correct_cls = {}
for line in content.split('\n'):
    match = re.search(pattern2, line)
    if match:
        instance = int(match.group(1))
        # Column order: instance, CFL, ?, resrho, CL, CD, total_residual
        cl_val = float(match.group(5))
        cd_val = float(match.group(6))
        correct_cls[instance] = cl_val

print(f"   Matched {len(correct_cls)} instances")
if correct_cls:
    for inst in sorted(correct_cls.keys()):
        print(f"   Instance {inst}: CL = {correct_cls[inst]:.10f}")

# Show sample lines for debugging
print("\n" + "="*90)
print("SAMPLE NK LINES FROM OUTPUT:")
print("="*90)
nk_lines = [line for line in content.split('\n') if re.search(r'^\s*1\s+\d+\s+\d+\s+\d+\s+NK', line)]
for i, line in enumerate(nk_lines[-9:]):  # Show last 9 (one per instance)
    print(f"{i+1}. {line}")

# Compare what the test script sees vs reality
print("\n" + "="*90)
print("COMPARISON:")
print("="*90)
print()
if not test_cls:
    print("⚠️  TEST SCRIPT REGEX DOESN'T MATCH NK LINES!")
    print("    This explains why the test shows [NOT FOUND] for all instances")
    print()

if correct_cls:
    print("✓  Correct regex successfully extracts CL values")
    print(f"   Found CL values for all {len(correct_cls)} instances")
    
    # Group by alpha and check matching
    import numpy as np
    n1, n2 = 3, 3
    theta1 = np.linspace(0, 2*np.pi, n1, endpoint=False)
    theta2 = np.linspace(0, 2*np.pi, n2, endpoint=False)
    
    print("\n" + "="*90)
    print("GROUPING BY ALPHA:")
    print("="*90)
    
    alpha_groups = {}
    for i2 in range(n2):
        for i1 in range(n1):
            inst = i2 * n1 + i1 + 1
            alpha_deg = np.sin(theta1[i1]) + np.sin(theta2[i2])
            alpha_key = round(alpha_deg, 4)
            if alpha_key not in alpha_groups:
                alpha_groups[alpha_key] = []
            alpha_groups[alpha_key].append(inst)
    
    for alpha_key in sorted(alpha_groups.keys()):
        instances = alpha_groups[alpha_key]
        print(f"\nα = {alpha_key:7.4f}°  →  Instances: {instances}")
        if len(instances) > 1:
            cls = [correct_cls[inst] for inst in instances]
            print(f"   CL values: {[f'{cl:.10f}' for cl in cls]}")
            cl_diff = max(cls) - min(cls)
            print(f"   CL range:  {cl_diff:.2e}", end="")
            if cl_diff < 1e-10:
                print("  ✅ MATCH (machine precision)")
            elif cl_diff < 1e-6:
                print("  ✅ MATCH")
            else:
                print("  ❌ MISMATCH")

