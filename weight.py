from Bio.PDB import PDBParser, PDBIO
import os
import sys
import pandas as pd 
# 1 Dalton is approximately equal to the mass of one hydrogen atom, which is about 1.00784 atomic mass units (u).
# Script is close by not exact...
# Protein Weight (Da) = Protein Weight (kDa) * 1000

atomic_weights = {
    'H': 1.008,
    'He': 4.003,
    'Li': 6.941,
    'Be': 9.012,
    'B': 10.81,
    'C': 12.01,
    'N': 14.01,
    'O': 16.00,
    'F': 19.00,
    'Ne': 20.18,
    'Na': 22.99,
    'Mg': 24.31,
    'Al': 26.982,
    'Si': 28.085,
    'P': 30.974,
    'S': 32.06,
    'Cl': 35.45,
    'K': 39.098,
    'Ar': 39.95,
    'Ca': 40.08,
    'Sc': 44.956,
    'Ti': 47.87,
    'V': 50.942,
    'Cr': 51.996,
    'Mn': 54.938,
    'Fe': 55.845,
    'Ni': 58.693,
    'Co': 58.933,
    'Cu': 63.546,
    'Zn': 65.38,
    'Ga': 69.723,
    'Ge': 72.630,
    'As': 74.922,
    'Se': 78.971,
    'Br': 79.904,
    'Kr': 83.798,
    'Rb': 85.468,
    'Sr': 87.62,
    'Y': 88.906,
    'Zr': 91.224,
    'Nb': 92.906,
    'Mo': 95.95,
    'Tc': 98,
    'Ru': 101.1,
    'Rh': 102.9,
    'Pd': 106.4,
    'Ag': 107.9,
    'Cd': 112.4,
    'In': 114.8,
    'Sn': 118.7,
    'Sb': 121.8,
    'I': 126.9,
    'Te': 127.6,
    'Xe': 131.3,
    'Cs': 132.9,
    'Ba': 137.3,
    'La': 138.9,
    'Ce': 140.1,
    'Pr': 140.9,
    'Nd': 144.2,
    'Pm': 145,
    'Sm': 150.4,
    'Eu': 152.0,
    'Gd': 157.3,
    'Tb': 158.9,
    'Dy': 162.5,
    'Ho': 164.9,
    'Er': 167.3,
    'Tm': 168.9,
    'Yb': 173.0,
    'Lu': 175.0,
    'Hf': 178.5,
    'Ta': 180.9,
    'W': 183.8,
    'Re': 186.2,
    'Os': 190.2,
    'Ir': 192.2,
    'Pt': 195.1,
    'Au': 197.0,
    'Hg': 200.6,
    'Tl': 204.4,
    'Pb': 207.2,
    'Bi': 208.9,
    'Th': 232.0,
    'Pa': 231.0,
    'U': 238.0,
    'Np': 237.0,
    'Pu': 244.0,
    'Am': 243.0,
    'Cm': 247.0,
    'Bk': 247.0,
    'Cf': 251.0,
    'Es': 252.0,
    'Fm': 257.0,
    'Md': 258.0,
    'No': 259.0,
    'Lr': 262.0,
}

def calculate_structure_weight(structure):
    weight = 0.0
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    element = atom.element
                    if element in atomic_weights:
                        weight += atomic_weights[element]
    return weight/1000.0

def main(spacegroup):
    # Get the current working directory
    current_directory = os.getcwd()

    # Construct the path to the target directory
    target_directory = os.path.join(current_directory, f"run_sfall/pdb/pdb_{spacegroup}")

    # Check if the target directory exists
    if not os.path.exists(target_directory):
        print(f"Directory not found: {target_directory}")
        return
    # Create empty lists to store results
    pdb_ids = []
    structure_weights = []

    # Iterate over files in the target directory
    for filename in os.listdir(target_directory):
        # Check if the file has a .pdb extension
        if filename.endswith(".pdb"):
            pdb_file_path = os.path.join(target_directory, filename)

            # Create a PDB parser
            parser = PDBParser(QUIET=True)

            try:
                # Load the PDB structure from the file
                structure = parser.get_structure("protein", pdb_file_path)
    
                #calculate the weight of the structure
                weight = calculate_structure_weight(structure)
                
                pdb_ids.append(filename)
                structure_weights.append(weight)

                # Print the results every 10th iteration
                if len(pdb_ids) % 50 == 0:
                    print('... calculating structure weights ...', '\n')

            except Exception as e:
                print(f"Error processing PDB file {filename}: {e}")
                continue
    
    # Create a DataFrame from the lists
    structure_weights_df = pd.DataFrame({
        'PDB_ID': pdb_ids,
        'Calculated Structure Weight (kDa)': structure_weights
    })
    # Print the final dataframe
    # print(structure_weights_df)
    return structure_weights_df

# if __name__ == "__main__":
#     if len(sys.argv) < 2:
#         print("Please provide a spacegroup as an argument.")
#     else:
#         main(spacegroup=sys.argv[1])