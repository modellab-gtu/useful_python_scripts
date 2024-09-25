import os
import pandas as pd
import argparse
import subprocess

def obrms_rmsd_calc_2(conf, ref):
    try:
        result = subprocess.run(
            ['obrms', ref, conf, '-m'],  # Open Babel command to calculate RMSD
            capture_output=True,
            text=True,
            check=True
        )
        # Check if there's any output
        output_lines = result.stdout.strip().splitlines()
        
        if len(output_lines) == 0:
            print(f"Error: No output from obrms for {ref} and {conf}")
            return None
        
        # Parse the last line for the RMSD value
        rmsd_line = output_lines[-1]
        rmsd_value = rmsd_line.split()[-1]

        try:
            rmsd = float(rmsd_value)  # Extract and convert the RMSD value
        except ValueError:
            print(f"Error: Unable to parse RMSD value from {rmsd_value}")
            return None

        return rmsd

    except subprocess.CalledProcessError as e:
        print(f"Error running obrms for {ref} and {conf}: {e}")
        return None

def rmsd_matrix_ob(struct_dir):
    all_files = [f for f in os.listdir(struct_dir) if f.endswith(".sdf")]
    num_files = len(all_files)

    # Create an empty DataFrame for the RMSD values
    d = pd.DataFrame(0, index=all_files, columns=all_files).astype(float)

    # Fill the DataFrame with RMSD values using obrms
    for i in range(num_files):
        for j in range(i, num_files):
            mol1 = os.path.join(struct_dir, all_files[i])
            mol2 = os.path.join(struct_dir, all_files[j])
            rmsd = obrms_rmsd_calc_2(mol1, mol2)
            if rmsd is not None:
                d.iloc[i, j] = rmsd
                d.iloc[j, i] = rmsd  # Ensure the matrix is symmetric
            print(f"RMSD between {all_files[i]} and {all_files[j]}: {rmsd if rmsd is not None else 'N/A'}")

    # Save the results to a CSV file
    d.to_csv(f"Results_{os.path.basename(struct_dir)}_obrms.csv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute RMSD matrix using Open Babel and cluster molecules.')
    parser.add_argument('struct_dir', type=str, help='Directory containing structure files')

    args = parser.parse_args()

    rmsd_matrix_ob(args.struct_dir)

