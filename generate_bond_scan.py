# Open the original file and read its content
with open('HF_bond_scan.xyz', 'r') as file:
    lines = file.readlines()

# Extract the first two lines (which are constant for all files)
header = lines[:2]
coordinates = lines[2:]

# Initialize starting z-coordinate for H atom
start_z = 0.500
end_z = 2.000
increment = 0.0001

# Generate the new coordinates and write them to new files
file_count = 1

# Loop to create modified files
z_coord = start_z
while z_coord <= end_z:
    # Update the z-coordinate of H atom in the coordinates list
    new_coordinates = coordinates[:]
    new_coordinates[1] = f"H       0.000   0.000   {z_coord:.3f}\n"

    # Create a new filename for each increment
    new_filename = f'HF_bond_scan{file_count}.xyz'

    # Write the header and updated coordinates to the new file
    with open(new_filename, 'w') as new_file:
        new_file.writelines(header)
        new_file.writelines(new_coordinates)

    # Increment z-coordinate and file count
    z_coord += increment
    file_count += 1

print("Bond scan files generated successfully.")

