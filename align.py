import os
import pandas as pd
import argparse
from scipy.cluster.hierarchy import ward, fcluster
from scipy.spatial.distance import pdist
from rdkit.Chem import SDMolSupplier, SDWriter
from rdkit.Chem.rdMolAlign import GetBestRMS, AlignMol
import rdkit
import rdkit.Chem
import rdkit.Chem.rdMolAlign

def alignTwoMol(mol1, mol2):
    return AlignMol(mol1, mol2)  # Align mol1 to mol2

def calc_BestRMS(mol1,mol2):
    rmsd = GetBestRMS(mol1, mol2)
    return rmsd

def mainAlingTwoMol():
    parser = argparse.ArgumentParser(description="Give something ...")
    parser.add_argument("-i", "--mol_path", type=str, required=True, help="")
    parser.add_argument("-r", "--mol_ref_path", type=str, required=True, help="")
    parser.add_argument("-o", "--output", type=str, required=True, help="")

    args = parser.parse_args()
    mol_path = args.mol_path
    mol_ref_path = args.mol_ref_path
    output = args.output

    mol_supplier = SDMolSupplier(mol_path, removeHs=True)
    mol_ref_supplier = SDMolSupplier(mol_ref_path, removeHs=True)

    # Assuming we are aligning the first molecule in the files
    mol = mol_supplier[0]
    mol_ref = mol_ref_supplier[0]

    aligned = alignTwoMol(mol, mol_ref)
    writer = SDWriter(output)
    writer.write(mol)
    writer.close()

    rmsd = GetBestRMS(mol, mol_ref)
    print("RMSD-->", rmsd)

mainAlingTwoMol()
