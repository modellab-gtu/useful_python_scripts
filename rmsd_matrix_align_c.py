import os
import pandas as pd
import argparse
from scipy.cluster.hierarchy import ward, fcluster
from scipy.spatial.distance import pdist
from rdkit.Chem import SDMolSupplier
from rdkit.Chem.rdMolAlign import GetBestRMS
from rdkit.Chem import AssignAtomChiralTagsFromStructure, AssignStereochemistry

def alignTwoMol_2(mol1, mol2):

    try:
        rmsd = GetBestRMS(mol1, mol2)
    except RuntimeError:
        print(f"Warning: No substructure match found between molecules {mol1.GetProp('_Name')} and {mol2.GetProp('_Name')}")
        rmsd = None
    return rmsd

def rmsd_calc_2(conf, ref, remove_hs=False):
    suppl = SDMolSupplier(conf, removeHs=remove_hs)
    suppl_ref = SDMolSupplier(ref, removeHs=remove_hs)

    rmsd = alignTwoMol_2(suppl_ref[0], suppl[0])
    return rmsd

def rmsd_matrix(struct_dir, remove_hs=False):
    all_files = [f for f in os.listdir(struct_dir) if f.endswith(".sdf")]
    num_files = len(all_files)

    # Create an empty DataFrame for the RMSD values
    d = pd.DataFrame(0, index=all_files, columns=all_files).astype(float)

    # Fill the DataFrame with RMSD values
    for i in range(num_files):
        for j in range(i, num_files):
            mol1 = os.path.join(struct_dir, all_files[i])
            mol2 = os.path.join(struct_dir, all_files[j])
            rmsd = rmsd_calc_2(mol1, mol2, remove_hs)
            if rmsd is not None:
                d.iloc[i, j] = rmsd
                d.iloc[j, i] = rmsd  # Ensure the matrix is symmetric
            print(f"RMSD between {all_files[i]} and {all_files[j]}: {rmsd if rmsd is not None else 'N/A'}")

    # Continue with clustering as before
    condensed_dist_matrix = pdist(d.values)
    Z = ward(condensed_dist_matrix)
    clustering = fcluster(Z, t=1.0, criterion='distance')

    d['cluster'] = clustering
    d.to_csv(f"Results_{os.path.basename(struct_dir)}.csv")

#def rmsd_matrix(struct_dir, remove_hs=False):
#    #all_files = os.listdir(struct_dir)
#    all_files = [f for f in os.listdir(struct_dir) if f.endswith(".sdf")]
#    num_files = len(all_files)
#    
#    # Create an empty DataFrame for the RMSD values
#    d = pd.DataFrame(0, index=all_files, columns=all_files).astype(float)
#    
#    # Fill the DataFrame with RMSD values
#    for i in range(num_files):
#        for j in range(i, num_files):
#            mol1 = os.path.join(struct_dir, all_files[i])
#            mol2 = os.path.join(struct_dir, all_files[j])
#            rmsd = rmsd_calc_2(mol1, mol2, remove_hs)
#            d.iloc[i, j] = rmsd
#            d.iloc[j, i] = rmsd  # Ensure the matrix is symmetric
#            print(f"RMSD between {all_files[i]} and {all_files[j]}: {rmsd}")
#    
#    # Convert the DataFrame to a condensed distance matrix for clustering
#    condensed_dist_matrix = pdist(d.values)
#    
#    # Perform hierarchical clustering
#    Z = ward(condensed_dist_matrix)
#    clustering = fcluster(Z, t=0.5, criterion='distance')
#    
#    # Add cluster labels to the DataFrame
#    d['cluster'] = clustering
#    
#    # Save the results to a CSV file
#    d.to_csv(f"Results_{os.path.basename(struct_dir)}.csv")
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute RMSD matrix and cluster molecules.')
    parser.add_argument('struct_dir', type=str, help='Directory containing structure files')
    parser.add_argument('--remove_hs', action='store_true', help='Remove hydrogens for RMSD calculation')
    
    args = parser.parse_args()
    
    rmsd_matrix(args.struct_dir, args.remove_hs)

