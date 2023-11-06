import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels as sm
import os, sys
from Bio.PDB import PDBParser

from plot import plot as p
from plot import dist as d
from format_scripts import write as w
from scrape_scripts import download_files as dl
from scrape_scripts import scrape as s
import weight as wt
# Research Question: whether there is a linear dependence of X-ray diffraction intensities on the number of unit cells exposed to the X-rays.

# input_data
# h, k, l = orientation of the crystal lattice planes in three-dimensional reciprocal space
# FREE = which reflections were used for the refinement
# FP = the observed structure factor amplitudes, (observed/measured amplitudes of the structure factors)
# SIGFP = the estimated standard deviations of the observed structure factor amplitudes
# FC = model-based structure factor amplitudes
# PHIC = calculated phase angles associated with the FC column
# FC_ALL, PHIC_ALL = additional structure factor amplitudes and phase angles (useful for refinement/model building)
# DELFWT, PHDELWT = difference in Fourier coefficients (also used for model building/refinement)
# FOM = figure of merit (FOM), quantifies the reliability of point estimates
# FC_ALL_LS, PHIC_ALL_LS = represents the least square Fourier coefficients (map calculations)

# output_data
# h, k, l = orientation of the crystal lattice planes in three-dimensional reciprocal space
# I (observed) = observed structure factor amplitudes (experimentally measured values of diffraction patterns)
        # provide information about the strength of the reflections in the diffraction pattern
# SIGFP (Sigma-Fobserved) = estimated standard deviations of the observed amplitudes (FP), indicate uncertainty in the observed amplitudes  
# FC (Fcalc) = calculated structure factor amplitudes based on the model (predicted by a structure model compared to observed amplitudes (FP) during refinement process)
# PHIC (PHIcalc) = calculated phase angles associated with the FC column (used in calculating electron density maps)
 
# FP and FC contain information about the observed and calculated structure factor amplitudes.
def analyze_crystals(space_group):
    results_list = []
    intensity_dict = {}
    phase_dict = {}
    data_dir = f"/Users/adamkurth/Documents/vscode/CXFEL_Image_Analysis/CXFEL/unitcell_project/run_sfall/data/data_{space_group}"
    
    # Iterate over .hkl files in the directory
    for file_name in os.listdir(data_dir):
        if file_name.endswith(".hkl"):
            pdb_id = file_name[:4]  # Extract the PDB ID from the first 4 characters of the file name
            file_path = os.path.join(data_dir, file_name)
            try:
                # Check if the file is empty
                if os.path.getsize(file_path) == 0:
                    print(f"Warning: {file_path} is empty.")
                    continue
                df_in = pd.read_table(file_path, delim_whitespace=True, header=None)
                
            except FileNotFoundError:
                print(f"Error: {file_path} not found.")
                continue
            except pd.errors.EmptyDataError:
                print(f"Error: {file_path} is empty or contains only whitespace.")
                continue
            except pd.errors.ParserError:
                print(f"Error parsing {file_path}. Skipping.")
                continue
        
            intensity_dict[pdb_id] = df_in.iloc[:, 3]
            phase_dict[pdb_id] = df_in.iloc[:, 4]
            
            if df_in.shape[1] >= 3:
                intensity_mean = df_in.iloc[:, 3].mean()
                intensity_max = df_in.iloc[:, 3].max()
                intensity_min = df_in.iloc[:, 3].min()
                max_min_diff_intensity = intensity_max - intensity_min
                phase_mean = df_in.iloc[:, 4].mean()
                phase_max = df_in.iloc[:, 4].max()
                phase_min = df_in.iloc[:, 4].min()
                max_min_diff_phase = phase_max - phase_min
            else:
                intensity_mean = None
                intensity_max = None
                intensity_min = None
                max_min_diff_intensity = None
                phase_mean = None
                phase_max = None
                phase_min = None
                max_min_diff_phase = None
                
            results_list.append({'PDB_ID': pdb_id, 
                                 'Spacegroup': space_group, 
                                 'Mean Intensity': intensity_mean,
                                 'Max Intensity': intensity_max, 
                                 'Min Intensity': intensity_min,
                                 'Max-Min Intensity Difference': max_min_diff_intensity, 
                                 'Mean Phase': phase_mean,
                                 'Max Phase': phase_max, 
                                 'Min Phase': phase_min,
                                 'Max-Min Phase Difference': max_min_diff_phase
                                 })
    
    intensity_df = pd.DataFrame(intensity_dict)
    phase_df = pd.DataFrame(phase_dict)
    for i in (intensity_df, phase_df):
        intensity = intensity_df.dropna() #desired output
        phase = phase_df.dropna() #desired output

    unit_cell_df = extract_unit_cell_attributes(space_group) 
    intensity_df = pd.DataFrame(results_list)
    intensity_df = pd.concat([unit_cell_df, intensity_df], axis=1, join="outer")
    # Drop the duplicate 'PDB_ID' column
    intensity_df = intensity_df.loc[:, ~intensity_df.columns.duplicated()]
    # gives unitcell df with all unit cell attributes and intensity df with all intensity attributes
    final_df = calculateVolume(unit_cell_df['alpha'], unit_cell_df['beta'], unit_cell_df['gamma'],unit_cell_df['a'], unit_cell_df['b'], unit_cell_df['c'], intensity_df)
    # print("Before merge: \n", final_df)
    #retrieves weights from weight.py
    weights = wt.main(space_group)
    weights['PDB_ID'] = weights['PDB_ID'].str.replace('.pdb', '', regex=False)
    #combine
    final_df = pd.merge(final_df, weights, on='PDB_ID')

    # Finally reorder columns
    final_df = final_df[['PDB_ID', 'Spacegroup', 'Calculated Structure Weight (kDa)', 'a', 'b', 'c', 'alpha', 'beta', 'gamma', 'UnitCellVolume', 'CrystalVolume', 'VolumeToUnitCellVolRatio', 'Mean Intensity', 'Max Intensity', 'Min Intensity', 'Max-Min Intensity Difference', 'Mean Phase', 'Max Phase', 'Min Phase', 'Max-Min Phase Difference']]
    final_df, intensity, phase = [df.dropna() for df in [final_df, intensity, phase]]
    return final_df, intensity, phase

def extract_unit_cell_attributes(space_group):
    unit_cell_attributes_list = []
    pdb_dir = f"/Users/adamkurth/Documents/vscode/CXFEL_Image_Analysis/CXFEL/unitcell_project/run_sfall/pdb/pdb_{space_group}"
    
    # structure_weights = wt.main(space_group)
    # structure_weight_df = pd.DataFrame(structure_weights)
    
    for pdb_file_name in os.listdir(pdb_dir):
        if pdb_file_name.endswith(".pdb"):
            pdb_name = os.path.splitext(pdb_file_name)[0]
            pdb_path = os.path.join(pdb_dir, pdb_file_name)

            try:
                with open(pdb_path, "r", encoding="utf-8") as pdb_file:
                    pdb_lines = pdb_file.readlines()

                cryst_line = next((line for line in pdb_lines if line.startswith("CRYST1")), None)
                if cryst_line:
                    cryst_components = cryst_line.split()
                    unit_cell_dict = {
                        'PDB_ID': pdb_name,
                        'a': float(cryst_components[1]),
                        'b': float(cryst_components[2]),
                        'c': float(cryst_components[3]),
                        'alpha': float(cryst_components[4]),
                        'beta': float(cryst_components[5]),
                        'gamma': float(cryst_components[6])
                    }
                    unit_cell_attributes_list.append(unit_cell_dict)
                else:
                    print(f"Warning: CRYST1 line not found in {pdb_path}. Skipping.")

            except FileNotFoundError:
                print(f"Error: {pdb_path} not found.")
                continue
    
    unit_cell_attributes_df = pd.DataFrame(unit_cell_attributes_list)
    return unit_cell_attributes_df



def calculateVolume(alpha, beta, gamma, a, b, c, df_in):
    # Calculate unit cell volume (monoclinic)
    unitcell_vol = a * b * c
    # Calculate the volume of the parallelepiped
    crystal_vol = a * b * c * np.sqrt(1 - np.cos(np.deg2rad(alpha))**2 - np.cos(np.deg2rad(beta))**2 - np.cos(np.deg2rad(gamma))**2 + 2 * np.cos(np.deg2rad(alpha)) * np.cos(np.deg2rad(beta)) * np.cos(np.deg2rad(gamma)))
    
    # Append the calculated values to the DataFrame
    df_in['UnitCellVolume'] = unitcell_vol
    df_in['CrystalVolume'] = crystal_vol
    df_in['VolumeToUnitCellVolRatio'] = crystal_vol / unitcell_vol
    return df_in


# def extract_unit_cell_attributes(space_group, pdb_files):
#     # Create an empty list to store dictionaries of unit cell attributes
#     unit_cell_attributes = []

#     # Loop through each pdb file and calculate the atomic weight
#     for pdb_file in pdb_files:
#         atomic_weight = calculate_atomic_weight(pdb_file)

#         # Append the calculated values to the list of unit cell attributes
#         unit_cell_attributes.append({
#             'SpaceGroup': space_group,
#             'AtomicWeight': atomic_weight,
#         })

#     return unit_cell_attributes




if __name__ == "__main__":

    """Get pdb ids for each spacegroup"""
    # from scrape_scripts dir import s
    # MONOCLINIC
    # P 1 2 1
    # P 1 21 1
    # C 1 2 1
    # P 1 21/c 1
    # C 1 2/c 1
    # # C 1 c 1
    
    ids_P121 = s.get_pdb_ids("P 1 2 1", "monoclinic", limit=100)
    ids_P1211 = s.get_pdb_ids("P 1 21 1", "monoclinic", limit=100)
    ids_C121 = s.get_pdb_ids("C 1 2 1", "monoclinic", limit=100)
    # smaller sample size
    # ids_P121 = s.get_pdb_ids("P 1 21/c 1", "monoclinic", limit=100)
    # s.get_pdb_ids("P 1 21/c 1", "monoclinic", limit=50)
    # s.get_pdb_ids("C 1 2/c 1", "monoclinic", limit=50)
    
    """Download files from each spacegroup"""

    # dl = dl.ProteinDownloader()
    # dl.download_files(ids_P121, "P121")
    # dl.download_files(ids_P1211, "P1211")
    # dl.download_files(ids_C121, "C121")

    """Analyze each spacegroup""" 
    
    
    # weights in kDa (kilo Daltons)
    
    wd = os.getcwd()
    P1211_df, P1211_intensities, P1211_phases = analyze_crystals("P1211")    
    print('P1211 Intensity DF: \n', P1211_intensities)
    print('P1211 Phase DF: \n', P1211_phases)
    print('P1211 Main DF: \n', P1211_df)  
    # # w.write(P1211_df, P1211_intensities, P1211_phases, 'P1211', wd)
    P121_df, P121_intensities, P121_phases = analyze_crystals("P121")
    print('P121 Intensity DF: \n', P121_intensities)
    print('P121 Phase DF: \n', P121_phases)
    print('P121 Main DF: \n', P121_df)  
    # # w.write(P121_df, P121_intensities, P121_phases, 'P121', wd)
    C121_df, C121_intensities, C121_phases  = analyze_crystals("C121")
    print('C121 Intensity DF: \n', C121_intensities)
    print('C121 Phase DF: \n', C121_phases)
    print('C121 Main DF: \n', C121_df)  
    # # # w.write(C121_df, C121_intensities, C121_phases, 'C121', wd)

    for df in [P1211_intensities, P121_intensities, C121_intensities]:
        print(df.describe())
        print("\n")

    """"""
    d.plot_hist(P1211_intensities)
    d.plot_hist_phases(P1211_phases)
    
    d.contour(P1211_intensities)
    d.contour(P1211_phases)
    
    d.plot_density(P1211_intensities)
    d.plot_density(P1211_phases)
    
    # d.plot_boxplot(P1211_intensities)   
    
    # lm = p.linear_model(['VolumeToUnitCellVolRatio'], 'Mean Intensity', P1211_df, P1211_df)
    # lm = p.linear_model(['VolumeToUnitCellVolRatio'], 'Mean Phase', P1211_df, P1211_df)
    # lm = p.linear_model(['VolumeToUnitCellVolRatio'], 'Calculated Structure Weight (kDa)', P1211_df, P1211_df)

    # lm = p.linear_model(['VolumeToUnitCellVolRatio'], 'Mean Intensity', P121_df, P121_df)
    # lm = p.linear_model(['VolumeToUnitCellVolRatio'], 'Mean Phase', P121_df, P121_df)
    # lm = p.linear_model(['VolumeToUnitCellVolRatio'], 'Calculated Structure Weight (kDa)', P121_df, P121_df)

    # lm = p.linear_model(['VolumeToUnitCellVolRatio'], 'Mean Intensity', C121_df, C121_df)
    # lm = p.linear_model(['VolumeToUnitCellVolRatio'], 'Mean Phase', C121_df, C121_df)
    # lm = p.linear_model(['VolumeToUnitCellVolRatio'], 'Calculated Structure Weight (kDa)', C121_df, C121_df)

    # residuals = p.plot_residuals(lm, 'VolumeToUnitCellVolRatio', 'Max-Min Difference', P1211_df, P1211_df)
    
    # linear_model(['VolumeToUnitCellVolRatio'], 'Mean FP', results_df, results_df)
    # conditition number is large, indicates stable model
    # lm = p.linear_model(['VolumeToUnitCellVolRatio'], 'Mean FP', results_df, results_df)
    # p.plot_residuals(lm, 'VolumeToUnitCellVolRatio', 'Mean FP', results_df, results_df)

    # for spacegroup P121 
    
    # def analyze_spacegroup(space_group):
    #     crystal_data = []
    #     if space_group == "P1211":
    #         crystal_data = [("104m", "P1211"), ("137l", "P1211"), ("169l", "P1211"), ("1a28", "P1211"), ("1a2a", "P1211"), ("105m", "P1211"), ("153l", "P1211"), ("154l", "P1211"), ("157d", "P1211"), ("176l", "P1211"), ("180l", "P1211"), ("1a3a", "P1211"), ("1a01", "P1211"), ("1a02", "P1211"), ("1a0o", "P1211"), ("1a2a", "P1211"), ("1a2j", "P1211"), ("135l", "P1211")]
    #     elif space_group == "P121":
    #         crystal_data = [("1gh4", "P1211"), ("1gxb", "P1211"), ("1hjc", "P1211"), ("1hkn", "P1211"), ("1ic1", "P1211"), ("1ijy", "P1211"), ("1jde", "P1211"), ("1kbl", "P1211"), ("1kc7", "P1211"), ("1kog", "P1211"), ("1kwy", "P1211"), ("1kyz", "P1211"), ("1lf3", "P1211"), ("1lf4", "P1211"), ("1vgl", "P1211"), ("1vjh", "P1211"), ("1vl9", "P1211"), ("1xpp", "P1211"), ("1xtf", "P1211"), ("1y5x", "P1211"), ("1y99", "P1211")]
    #     else:
    #         raise ValueError("Invalid space group name")
        
    #     dfs_in = []
    #     dfs_out = []
    #     unit_cells = []
        
    #     for name, space_group in crystal_data:
    #         try:
    #             df_in, df_out, unit_cell = analyze_crystal(name, space_group)
    #             dfs_in.append(df_in)
    #             dfs_out.append(df_out)
    #             unit_cells.append(unit_cell)
    #         except FileNotFoundError:
    #             print(f"Error: File not found for crystal {name}")
    #         except Exception as e:
    #             print(f"Error analyzing crystal {name}: {str(e)}")
        
    #     crystal_names = [name[0] for name in crystal_data]
    #     results_df = create_crystal_df(crystal_names, space_group, weights, dfs_in, unit_cells, dfs_out)
    #     intensity_df = create_intensity_df(crystal_names, dfs_out)
        
    #     return crystal_data, dfs_in, dfs_out, unit_cells, results_df, intensity_df
    

    
    

    
    
    
# sm.graphics.plot_fit(lm, 'VolumeToUnitCellVolRatio', vlines=True)

# Adding the total structure weights
