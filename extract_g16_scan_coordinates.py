import re

def extract_geometries_from_log(log_content):
    """Extracts all geometries from Gaussian log content."""
    geometries = []
    geometry = None
    capture = False
    skip_lines = 4

    for line in log_content:
        if "Standard orientation:" in line:
            # Start capturing after 'Standard orientation:'
            capture = True
            geometry = []
            skip_lines = 4  # Reset skip counter for each new geometry
            continue

        if capture:
            if skip_lines > 0:
                skip_lines -= 1  # Skip the next 5 lines
                continue
            if "-----" in line:
                # Stop capturing when encountering dashed lines again
                geometries.append(geometry)
                capture = False
            else:
                parts = line.split()
                atom_number = int(parts[1])
                x, y, z = map(float, parts[3:6])
                geometry.append((atom_number, x, y, z))

    return geometries

def convert_to_xyz(geometry):
    """Converts a single geometry to XYZ format string."""
    atom_symbols = {1: "H", 6: "C", 7: "N", 8: "O"}
    xyz_content = f"%s\n%s\n" % (len(geometry), "test")
    for atom_number, x, y, z in geometry:
        symbol = atom_symbols.get(atom_number, 'X')  # Default to 'X' if atom number is not in the dictionary
        xyz_content += f"{symbol} {x:.6f} {y:.6f} {z:.6f}\n"
    return xyz_content

def save_geometries_to_files(geometries, base_filename):
    """Saves each geometry to a separate XYZ file."""
    for i, geometry in enumerate(geometries, start=1):
        xyz_content = convert_to_xyz(geometry)
        filename = f"{base_filename}{i}.xyz"
        with open(filename, 'w') as file:
            file.write(xyz_content)

# Example usage:
log_file_path = 'H2O_angle_scan.log'  # Replace with your log file path
base_filename = "H2O_angle_scan"
with open(log_file_path, 'r') as log_file:
    log_content = log_file.readlines()

geometries = extract_geometries_from_log(log_content)
save_geometries_to_files(geometries, base_filename)

