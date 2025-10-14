#!/usr/bin/env python3
import pymol2
import os
import subprocess

# ----------------- User Settings -----------------
PROTEIN_FOLDER = "6CDL"  # folder containing index.ndx and pdb subfolder
PDB_SUBFOLDER = "pdb"              # subfolder containing PDB frames
BASENAME = "trjmol"                 # PDB basename
GROUP_NUM = 13                      # 0-based ligand group number
START = 0
END = 1000
# -------------------------------------------------

INDEX_FILE = os.path.join(PROTEIN_FOLDER, "index.ndx")
PDB_FOLDER = os.path.join(PROTEIN_FOLDER, PDB_SUBFOLDER)
PROTEIN_NAME = os.path.basename(PROTEIN_FOLDER.rstrip("/"))
OUTPUT_FILE = f"residues_{PROTEIN_NAME}.txt"

def extract_nearby_residues_allframes(start=START, end=END):
    """
    Extract combined unique residues near ligand from all frames.
    Saves combined list to OUTPUT_FILE and returns it as a sorted list.
    """
    with open(INDEX_FILE) as f:
        lines = [line.strip() for line in f if line.strip().startswith('[')]
    if GROUP_NUM >= len(lines):
        raise ValueError(f"Group number {GROUP_NUM} exceeds number of groups in index file")
    ligand_group = lines[GROUP_NUM].strip('[] ').strip()
    print(f"Ligand group name: {ligand_group}")

    all_residues = set()
    
    with pymol2.PyMOL() as pymol:
        cmd = pymol.cmd

        for i in range(start, end+1):
            pdb_file = os.path.join(PDB_FOLDER, f"{BASENAME}{i}.pdb")
            if not os.path.isfile(pdb_file):
                print(f"Skipping missing file: {pdb_file}")
                continue

            cmd.reinitialize()
            cmd.load(pdb_file, "complex")
            cmd.select("ligand", f"resn {ligand_group}")
            cmd.select("nearby", "byres (ligand around 4 and polymer)")

            # Collect unique residue indices (CA only)
            res_indices = set()
            cmd.iterate("nearby and name CA", "res_indices.add(int(resi))", space={'res_indices': res_indices})

            sorted_res = sorted(res_indices)
            print(f"{i}: {sorted_res}")
            all_residues.update(res_indices)

    # Save all unique residues
    with open(OUTPUT_FILE, "w") as f:
        for res in sorted(all_residues):
            f.write(f"{res}\n")

    print(f"All unique residues saved in {OUTPUT_FILE}")
    return sorted(all_residues)


def make_ndx_with_residues(residues, frame_pdb=None):
    """
    Calls GROMACS gmx make_ndx to create a new index group including the residues.
    residues: list of integers
    frame_pdb: path to PDB file to use (frame1 recommended)
    """
    if frame_pdb is None:
        frame_pdb = os.path.join(PDB_FOLDER, f"{BASENAME}1.pdb")
    if not os.path.isfile(frame_pdb):
        raise FileNotFoundError(f"PDB file not found: {frame_pdb}")

    # Each residue prefixed with 'r' and joined by '|' for GROMACS
    residue_str = "|".join(f"r{r}" for r in residues)
    ndx_input = f"{residue_str}\nq\n"

    cmd_args = ["gmx", "make_ndx", "-f", frame_pdb, "-n", INDEX_FILE, "-o", INDEX_FILE]
    print(f"Running: {' '.join(cmd_args)}")
    
    process = subprocess.run(cmd_args, input=ndx_input.encode(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if process.returncode != 0:
        print("Error running gmx make_ndx:")
        print(process.stderr.decode())
    else:
        print("gmx make_ndx completed successfully.")
        print(process.stdout.decode())


def rename_ndx_group_by_residues(index_file, residues, new_name="nearbylig"):
    """
    Rename the GROMACS index group created from a specific residue list.
    residues: list of integers used to create the group
    """
    # Auto-generated name from gmx make_ndx uses '_' instead of '|'
    auto_name = "_".join(f"r_{r}" for r in residues)

    with open(index_file, "r") as f:
        lines = f.readlines()

    # Find the line starting with '[' and containing the auto-generated name
    for i, line in enumerate(lines):
        if line.strip().startswith('[') and auto_name in line:
            lines[i] = f"[ {new_name} ]\n"
            print(f"Renamed group '{auto_name}' to '{new_name}' in index file.")
            break
    else:
        print(f"Warning: Could not find group with residues {auto_name} in {index_file}")

    # Write back
    with open(index_file, "w") as f:
        f.writelines(lines)


if __name__ == "__main__":
    # Step 1: Extract combined unique residues from all frames
    combined_residues = extract_nearby_residues_allframes(start=START, end=END)
    
    # Step 2: Generate new index group in GROMACS (using frame 1 PDB)
    make_ndx_with_residues(combined_residues)
    
    # Step 3: Rename the auto-generated group to 'nearbylig'
    rename_ndx_group_by_residues(INDEX_FILE, combined_residues, "nearbylig")

