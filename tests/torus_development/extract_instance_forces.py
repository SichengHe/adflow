"""
Utility to extract instance-specific force coefficients from ADflow time spectral output
"""
import re
import numpy as np


def extract_instance_forces(output_text, n_instances):
    """
    Extract force coefficients for each time spectral instance from solver output.

    Parameters
    ----------
    output_text : str
        Captured stdout from CFDSolver call
    n_instances : int
        Number of time spectral instances

    Returns
    -------
    forces : dict
        Dictionary with instance number as key, containing:
        - 'cl': Lift coefficient
        - 'cd': Drag coefficient
        - 'cmx': Moment coefficient about x-axis
        - 'cmy': Moment coefficient about y-axis
        - 'cmz': Moment coefficient about z-axis
        - 'fx': Force coefficient in x
        - 'fy': Force coefficient in y
        - 'fz': Force coefficient in z
        - 'cfx': Friction force coefficient in x
        - 'cfy': Friction force coefficient in y
        - 'cfz': Friction force coefficient in z
        - 'resrho': Density residual

    cl_array : ndarray
        Array of CL values indexed by instance (1-based indexing matches Fortran)
    cd_array : ndarray
        Array of CD values
    cm_array : ndarray
        Array of CMz values (pitching moment)
    """

    # Pattern to match NK iteration lines with force coefficients
    # Format: 1  instance  iter  total  NK/*NK  ----  CFL  ?  resrho  CL  CD  residual
    pattern = r'^\s*1\s+(\d+)\s+\d+\s+\d+\s+\*?NK\s+----\s+([\d\.]+)\s+([\d\.]+)\s+([-\d\.E\+]+)\s+([-\d\.E\+]+)\s+([-\d\.E\+]+)\s+([-\d\.E\+]+)'

    instance_forces = {}

    for line in output_text.split('\n'):
        match = re.search(pattern, line)
        if match:
            instance = int(match.group(1))
            resrho = float(match.group(4))
            cl = float(match.group(5))
            cd = float(match.group(6))

            # Store last occurrence (final iteration for this instance)
            instance_forces[instance] = {
                'cl': cl,
                'cd': cd,
                'resrho': resrho
            }

    # Also try to extract from force summary lines if available
    # These contain more complete force information
    force_summary_pattern = r'Inst\s+(\d+):\s+CL\s*=\s*([-\d\.E\+]+)\s+CD\s*=\s*([-\d\.E\+]+)\s+CMz\s*=\s*([-\d\.E\+]+)'

    for line in output_text.split('\n'):
        match = re.search(force_summary_pattern, line)
        if match:
            instance = int(match.group(1))
            cl = float(match.group(2))
            cd = float(match.group(3))
            cmz = float(match.group(4))

            if instance in instance_forces:
                instance_forces[instance]['cmz'] = cmz
            else:
                instance_forces[instance] = {'cl': cl, 'cd': cd, 'cmz': cmz}

    # Convert to arrays for easier use
    cl_array = np.full(n_instances + 1, np.nan)  # 1-indexed to match Fortran
    cd_array = np.full(n_instances + 1, np.nan)
    cm_array = np.full(n_instances + 1, np.nan)

    for inst, forces in instance_forces.items():
        if 1 <= inst <= n_instances:
            cl_array[inst] = forces['cl']
            cd_array[inst] = forces['cd']
            cm_array[inst] = forces.get('cmz', np.nan)

    return instance_forces, cl_array, cd_array, cm_array


def print_instance_forces(instance_forces, n_instances):
    """
    Pretty print force coefficients for each instance.

    Parameters
    ----------
    instance_forces : dict
        Dictionary returned by extract_instance_forces
    n_instances : int
        Number of time spectral instances
    """

    print("\n" + "="*80)
    print("INSTANCE-SPECIFIC FORCE COEFFICIENTS")
    print("="*80)
    print()
    print("Instance      CL              CD              CMz             resrho")
    print("-"*80)

    for inst in range(1, n_instances + 1):
        if inst in instance_forces:
            forces = instance_forces[inst]
            cl = forces['cl']
            cd = forces['cd']
            cmz = forces.get('cmz', np.nan)
            resrho = forces.get('resrho', np.nan)

            if np.isnan(cmz):
                print(f"  {inst:3d}     {cl:15.10f}  {cd:15.10e}      N/A         {resrho:.6e}")
            else:
                print(f"  {inst:3d}     {cl:15.10f}  {cd:15.10e}  {cmz:15.10e}  {resrho:.6e}")
        else:
            print(f"  {inst:3d}     [NOT FOUND]")

    print()


def save_instance_forces(filename, instance_forces, n_instances, alpha_values=None):
    """
    Save instance force coefficients to a text file.

    Parameters
    ----------
    filename : str
        Output filename (e.g., 'instance_forces.dat')
    instance_forces : dict
        Dictionary returned by extract_instance_forces
    n_instances : int
        Number of time spectral instances
    alpha_values : array-like, optional
        Angle of attack for each instance (degrees)
    """

    with open(filename, 'w') as f:
        f.write("# Instance-specific force coefficients from ADflow time spectral\n")
        if alpha_values is not None:
            f.write("# Instance  Alpha(deg)     CL              CD              CMz\n")
        else:
            f.write("# Instance     CL              CD              CMz\n")

        for inst in range(1, n_instances + 1):
            if inst in instance_forces:
                forces = instance_forces[inst]
                cl = forces['cl']
                cd = forces['cd']
                cmz = forces.get('cmz', 0.0)

                if alpha_values is not None:
                    alpha = alpha_values[inst-1] if inst-1 < len(alpha_values) else 0.0
                    f.write(f"{inst:4d}  {alpha:10.5f}  {cl:15.10f}  {cd:15.10e}  {cmz:15.10e}\n")
                else:
                    f.write(f"{inst:4d}  {cl:15.10f}  {cd:15.10e}  {cmz:15.10e}\n")


if __name__ == "__main__":
    """
    Example usage: python extract_instance_forces.py <logfile>
    """
    import sys

    if len(sys.argv) < 2:
        print("Usage: python extract_instance_forces.py <logfile>")
        print("\nExample: python extract_instance_forces.py torus_degenerate.log")
        sys.exit(1)

    filename = sys.argv[1]

    try:
        with open(filename, 'r') as f:
            output_text = f.read()
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found!")
        sys.exit(1)

    # Try to detect number of instances from output
    pattern = r'nTimeIntervalsSpectral.*=\s*(\d+)'
    match = re.search(pattern, output_text)
    if match:
        n_instances = int(match.group(1))
    else:
        # Default to 9 for torus 3x3
        n_instances = 9
        print(f"Could not detect number of instances, assuming {n_instances}")

    instance_forces, cl_array, cd_array, cm_array = extract_instance_forces(output_text, n_instances)

    print_instance_forces(instance_forces, n_instances)

    print("\nArrays (1-indexed, element 0 unused):")
    print(f"CL: {cl_array[1:]}")
    print(f"CD: {cd_array[1:]}")
    print(f"CM: {cm_array[1:]}")

    # Save to file
    output_file = filename.replace('.log', '_forces.dat')
    save_instance_forces(output_file, instance_forces, n_instances)
    print(f"\nForce coefficients saved to: {output_file}")
