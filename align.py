import os
import argparse
from rdkit.Chem import SDMolSupplier, SDWriter, RemoveHs, AddHs
from rdkit.Chem.rdMolAlign import GetBestRMS, AlignMol
import rdkit
import rdkit.Chem
import rdkit.Chem.rdMolAlign

def alignTwoMol(mol1, mol2):
    """Align mol1 to mol2 and return the aligned molecule."""
    if mol1 is None or mol2 is None:
        raise ValueError("One or both of the molecules are not loaded correctly.")
    # Align mol1 to mol2
    success = AlignMol(mol1, mol2)
    if not success:
        raise RuntimeError("Alignment failed. Ensure that the molecules have compatible structures.")
    return mol1

def calc_BestRMS(mol1, mol2):
    """Calculate and return the best RMSD between mol1 and mol2."""
    if mol1 is None or mol2 is None:
        raise ValueError("One or both of the molecules are not loaded correctly.")
    return GetBestRMS(mol1, mol2)

def mainAlingTwoMol():
    parser = argparse.ArgumentParser(description="Align molecules and save as SDF.")
    parser.add_argument("-i", "--mol_path", type=str, required=True, help="Path to the input SDF file.")
    parser.add_argument("-r", "--mol_ref_path", type=str, required=True, help="Path to the reference SDF file.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Path to the output SDF file.")

    args = parser.parse_args()
    mol_path = args.mol_path
    mol_ref_path = args.mol_ref_path
    output = args.output

    mol_supplier = SDMolSupplier(mol_path, removeHs=True)
    mol_ref_supplier = SDMolSupplier(mol_ref_path, removeHs=True)

    # Assuming we are aligning the first molecule in the files
    mol = mol_supplier[0]
    mol_ref = mol_ref_supplier[0]

    if mol is None or mol_ref is None:
        raise ValueError("Failed to load one or both molecules.")

    # Align the molecule (without hydrogens)
    mol_aligned = alignTwoMol(mol, mol_ref)

    # Calculate and print RMSD without hydrogens
    rmsd = calc_BestRMS(mol_aligned, mol_ref)
    print("RMSD-->", rmsd)

    # Add hydrogens after alignment
    mol_aligned_with_h = AddHs(mol_aligned, addCoords=True)

    # Save the aligned molecule with added hydrogens to the output SDF file
    writer = SDWriter(output)
    writer.write(mol_aligned_with_h)
    writer.close()

mainAlingTwoMol()

