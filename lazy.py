import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels as sm
import os, sys

from plot import plot as p
from plot import dist as d
from format_scripts import write as w
from scrape_scripts import download_files as dl

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
# FP (observed) = observed structure factor amplitudes (experimentally measured values of diffraction patterns)
        # provide information about the strength of the reflections in the diffraction pattern
# SIGFP (Sigma-Fobserved) = estimated standard deviations of the observed amplitudes (FP), indicate uncertainty in the observed amplitudes  
# FC (Fcalc) = calculated structure factor amplitudes based on the model (predicted by a structure model compared to observed amplitudes (FP) during refinement process)
# PHIC (PHIcalc) = calculated phase angles associated with the FC column (used in calculating electron density maps)
 
# FP and FC contain information about the observed and calculated structure factor amplitudes.
def analyze_crystal(name, space_group, working_dir):
    os.chdir(working_dir)

    # Set the working directory
    pdb_name = f"{name}.pdb"
    input_name = f"{name}_input_formatted.txt"
    output_name = f"{name}_output_formatted.txt"
    # # Check if the name contains '_hkl_out_formatted'
    # if "_hkl_out_formatted" in name or "_output_formatted" in name:
    #     output_name = f"{name}.txt"
    # else:
    #     output_name = f"{name}_hkl_out_formatted.txt"
    # working_dir = "/Users/adamkurth/Documents/vscode/CXFEL_Image_Analysis/CXFEL/unitcell_project"

    # Retrieve the pdb data
    pdb_path = f"{working_dir}/pdb_files/{pdb_name}"
    try:
        with open(pdb_path, "r", encoding="utf-8") as pdb_file:
            pdb_lines = pdb_file.read().splitlines()
    except FileNotFoundError:
        print(f"Error: {pdb_path} not found.")
        return

    # Extract CRYST1 line
    cryst_line = next(line for line in pdb_lines if line.startswith("CRYST1"))
    cryst_components = cryst_line.split()

    # Dimension lengths
    a = float(cryst_components[1])
    b = float(cryst_components[2])
    c = float(cryst_components[3])

    # Angles (in degrees)
    alpha = float(cryst_components[4])
    beta = float(cryst_components[5])
    gamma = float(cryst_components[6])

    global unit_cell_attributes
    unit_cell_attributes = {'a': a, 'b': b, 'c': c, 'alpha': alpha, 'beta': beta, 'gamma': gamma}

    # Retrieve the input data
    input_path = f"{working_dir}/input_data/{input_name}"
    if not os.path.exists(input_path):
        print(f"Error: {input_path} not found.")
        return

    with open(input_path, "r", encoding="utf-8") as input_file:
        input_data = input_file.read().splitlines()

    # Skip any blank lines
    input_data = [line for line in input_data if line.strip() != ""]
    # Split the data and extract columns
    data_rows = [line.split() for line in input_data[1:]]

    h = [int(line[0]) for line in data_rows]
    k = [int(line[1]) for line in data_rows]
    l = [int(line[2]) for line in data_rows]
    free = [int(line[3]) for line in data_rows]
    fp = [float(line[4]) for line in data_rows]
    sigfp = [float(line[5]) for line in data_rows]
    fc = [float(line[6]) for line in data_rows]
    phic = [float(line[7]) for line in data_rows]
    fc_all = [float(line[8]) for line in data_rows]
    phic_all = [float(line[9]) for line in data_rows]
    fwt = [float(line[10]) for line in data_rows]
    phwt = [float(line[11]) for line in data_rows]
    delfwt = [float(line[12]) for line in data_rows]
    phdelfwt = [float(line[13]) for line in data_rows]
    fom = [float(line[14]) for line in data_rows]
    fc_all_ls = [float(line[15]) for line in data_rows]
    phic_all_ls = [float(line[16]) for line in data_rows]
    
    column_names = ['h', 'k', 'l', 'FREE', 'FP', 'SIGFP', 'FC', 'PHIC', 'FC_ALL', 'PHIC_ALL', 'FWT', 'PHWT', 'DELFWT', 'PHDELWT', 'FOM', 'FC_ALL_LS', 'PHIC_ALL_LS']
    data_dict = {'h': h, 'k': k, 'l': l, 'FREE': free, 'FP': fp, 'SIGFP': sigfp, 'FC': fc, 'PHIC': phic, 'FC_ALL': fc_all, 'PHIC_ALL': phic_all, 'FWT': fwt, 'PHWT': phwt, 'DELFWT': delfwt, 'PHDELWT': phdelfwt, 'FOM': fom, 'FC_ALL_LS': fc_all_ls, 'PHIC_ALL_LS': phic_all_ls}
    df_in = pd.DataFrame(data_dict)
    
    calculateVolume(alpha, beta, gamma, a, b, c, df_in)
    
    # Retrieve the output data
    output_path = f"{working_dir}/output_data/{output_name}"
    with open(output_path, "r", encoding="utf-8") as output_file:
        output_data = output_file.read().splitlines()

    # Skip any blank lines
    output_data = [line for line in output_data if line.strip() != ""]

    # Split the data and extract columns
    data_rows = [line.split() for line in output_data[1:]]
    h = [int(line[0]) for line in data_rows]
    k = [int(line[1]) for line in data_rows]
    l = [int(line[2]) for line in data_rows]
    free = [int(line[3]) for line in data_rows]
    fp = [float(line[4]) for line in data_rows]
    sigfp = [float(line[5]) for line in data_rows]
    fcalc = [float(line[6]) for line in data_rows]
    phicalc = [float(line[7]) for line in data_rows]
    
    data_dict = {'h': h, 'k': k, 'l': l, 'FREE': free, 'FP': fp, 'SIGFP': sigfp, 'FCalc': fcalc, 'PHICalc': phicalc}
    df_out = pd.DataFrame(data_dict)
    
    return df_in, df_out, unit_cell_attributes

def calculateVolume(alpha, beta, gamma, a, b, c, df_in):
    # Calculate unit cell volume (monoclinic)
    unitcell_vol = a * b * c
    # Calculate the volume of the parallelepiped
    crystal_vol = a * b * c * np.sqrt(1 - np.cos(np.deg2rad(alpha))**2 - np.cos(np.deg2rad(beta))**2 - np.cos(np.deg2rad(gamma))**2 + 2 * np.cos(np.deg2rad(alpha)) * np.cos(np.deg2rad(beta)) * np.cos(np.deg2rad(gamma)))
    df_in['UnitCellVolume'] = unitcell_vol
    df_in['CrystalVolume'] = crystal_vol
    df_in['VolumeToUnitCellVolRatio'] = df_in['CrystalVolume'] / df_in['UnitCellVolume']
    return df_in

def create_crystal_df(crystal_names, spacegroup, weights, dfs_in, unit_cells, dfs_out):

    results_df = pd.DataFrame(columns=['Spacegroup', 'Crystal', 'a', 'b', 'c', 'alpha', 'beta', 'gamma', 'Unit Cell Volume', 'Crystal Volume', 'VolumeToUnitCellVolRatio', 'Mean FP', 'Structure Weight'])

    # Loop over the crystal names and dataframes
    for i, name in enumerate(crystal_names):
        # Calculate the unit cell volume and crystal volume
        a, b, c, alpha, beta, gamma = unit_cells[i].values()
        unit_cell_volume = dfs_in[i]['UnitCellVolume'].iloc[0]
        crystal_volume = dfs_in[i]['CrystalVolume'].iloc[0]
        ratio = dfs_in[i]['VolumeToUnitCellVolRatio'].iloc[0]
        # Calculate the mean FP value
        mean_fp = dfs_out[i]['FP'].mean()
        weight = weights[name]
        # Create a new row for the crystal
        row = [spacegroup, name, a, b, c, alpha, beta, gamma, unit_cell_volume, crystal_volume, ratio, mean_fp, weight]
        # Append the row to the results DataFrame
        results_df.loc[len(results_df)] = row
    results_df = results_df.round({'Unit Cell Volume': 2, 'Crystal Volume': 2, 'VolumeToUnitCellVolRatio': 3, 'Mean FP': 3})
    return results_df


def create_intensity_df(crystal_names, dfs_out):
    intensity_dfs = [df_out[["FP"]].rename(columns={"FP": name}) for name, df_out in zip(crystal_names, dfs_out)]
    # Concatenate the list of dataframes into a single dataframe
    intensity_df = pd.concat(intensity_dfs, axis=1)
    # Remove any NaN values
    # intensity_df = intensity_df.dropna()
    return intensity_df

def analyze_spacegroup(space_group, working_dir):
    os.chdir(working_dir)  # change to the specified working directory
    print(os.getcwd(), "\n")
    crystal_data = []
    if space_group == "P1211":
        crystal_data = [("104m", "P1211"), ("137l", "P1211"), ("169l", "P1211"), ("1a28", "P1211"), ("1a2a", "P1211"), ("105m", "P1211"), ("153l", "P1211"), ("154l", "P1211"), ("157d", "P1211"), ("176l", "P1211"), ("180l", "P1211"), ("1a3a", "P1211"), ("1a01", "P1211"), ("1a02", "P1211"), ("1a0o", "P1211"), ("1a2a", "P1211"), ("1a2j", "P1211"), ("135l", "P1211")]
    elif space_group == "P121":
        crystal_data = [("1gh4", "P121"), ("1gxb", "P121"), ("1hkn", "P121"), ("1kwy", "P121"), ("1hjc", "P121"), ("1ic1", "P121"), ("1ijy", "P121"), ("1jde", "P121"), ("1kbl", "P121"), ("1kc7", "P121"), ("1kog", "P121"), ("1kyz", "P121"), ("1lf3", "P121"), ("1lf4", "P121"), ("1vgl", "P121"), ("1vjh", "P121"), ("1vl9", "P121"), ("1xpp", "P121"), ("1xtf", "P121"), ("1y5x", "P121"), ("1y99", "P121")]
    else:
        raise ValueError("Invalid space group name")
    
    dfs_in = []
    dfs_out = []
    unit_cells = []
    
    for name, space_group in crystal_data:
        try:
            df_in, df_out, unit_cell = analyze_crystal(name, space_group, working_dir)
            dfs_in.append(df_in)
            dfs_out.append(df_out)
            unit_cells.append(unit_cell)
        except FileNotFoundError:
            print(f"Error: File not found for crystal {name}")
        except Exception as e:
            print(f"Error analyzing crystal {name}: {str(e)}")
    
    crystal_names = [name[0] for name in crystal_data]
    results_df = create_crystal_df(crystal_names, space_group, weights, dfs_in, unit_cells, dfs_out)
    intensity_df = create_intensity_df(crystal_names, dfs_out)
    print(results_df)
    print(intensity_df)
    return crystal_data, dfs_in, dfs_out, unit_cells, results_df, intensity_df, working_dir
    
if __name__ == "__main__":
    
    """Download files from each spacegroup"""
    # for spacegroup P1211
    
    # dl_P1211 = dl.ProteinDownloader()
    # p1211_spacegroup = ["104m", "137l", "169l", "1a28", "1a2a", "105m", "153l", "154l", "157d", "176l", "180l", "1a3a", "1a01", "1a02", "1a0o", "1a2a", "1a2j", "135l"]
    # dl_P1211.download_files(p1211_spacegroup, "P1211")
    
    # dl_p121 = dl.ProteinDownloader()
    # p121_spacegroup = ["1gh4", "1gxb", "1hjc", "1hkn", "1ic1", "1ijy", "1jde", "1kbl", "1kc7", "1kog", "1kwy", "1kyz", "1lf3", "1lf4", "1vgl", "1vjh", "1vl9", "1xpp", "1xtf", "1y5x", "1y99"]
    # dl_p121.download_files(p121_spacegroup, "P121")  

    #  FOR NEXT SPACEGROUP
    # dl_p121 = dl.ProteinDownloader()
    # p121_spacegroup = [""]
    # dl_p121.download_files(p121_spacegroup, "P121")  
    
    """Analyze each spacegroup""" 
    
    
    # weights in kDa (kilo Daltons)
    weights = {
        '104m': 18.03,
        '137l': 37.38, 
        '169l': 90.9, 
        '1a28': 59.74,
        '1a2a': 111.11,
        '19hc': 73.97,
        '105m': 18.03, 
        '153l': 20.41,
        '154l': 21.03,
        '157d' : 7.7,
        '176l' : 37.24,
        '180l' : 37.31,
        '1a01' : 64.38,
        '1a02' : 59.85, 
        '1a0o' : 114.17,
        '1a2j' : 21.16, 
        '1a3a' : 65.39,
        '1a3n' : 64.55,
        '12e8' :  94.84,
        '135l' : 14.23, 
        
        '1gh4' : 14.01,
        '1gxb' : 151.25, 
        '1hcj' : 107.55,
        '1hkn' : 94.71,
        '1ic1' : 44.22,
        '1ijy' : 29.89,
        '1jde' : 96.77,
        '1kbl' : 97.14,
        '1kc7' : 97.22, 
        '1kog' : 472.58,
        '1kyw' : 120.95, 
        '1kyz' : 121.71, 
        '1lf3' : 37.93, 
        '1lf4' : 37.12,
        '1o17' : 150.48,
        '1vfg' : 140.76, 
        '1vgl' : 48.57, 
        '1vjh' : 27.75,   
        '1vl9' : 14.64, 
        '1xpp' : 54.26,
        '1xtf' : 98.08,  
        '1y5x' : 86.39,
        '1y99' : 14.98
        
        
        }
    
    working_dir_p1211 = "/Users/adamkurth/Documents/vscode/CXFEL_Image_Analysis/CXFEL/unitcell_project/spacegroup/P1211"
    working_dir_p121 = "/Users/adamkurth/Documents/vscode/CXFEL_Image_Analysis/CXFEL/unitcell_project/spacegroup/P121"
    analyze_spacegroup("P1211", working_dir_p1211)
    analyze_spacegroup("P121", working_dir_p121)
    
    # top scoring crystals selected from the 20 crystals in the P1211 spacegroup
    # crystal_data = [("104m", "P1211"), ("137l", "P1211"), ("169l", "P1211"), ("1a28", "P1211"), ("1a2a", "P1211"), ("105m", "P1211"), ("153l", "P1211"), ("154l", "P1211"), ("157d", "P1211"), ("176l", "P1211"), ("180l", "P1211"), ("1a3a", "P1211"), ("1a01", "P1211"), ("1a02", "P1211"), ("1a0o", "P1211"), ("1a2a", "P1211"), ("1a2j", "P1211"), ("135l", "P1211")]
    # dfs_in = []
    # dfs_out = []
    # unit_cells = []
    
    # for name, space_group in crystal_data:
    #     try:
    #         df_in, df_out, unit_cell = analyze_crystal(name, space_group)
    #         dfs_in.append(df_in)
    #         dfs_out.append(df_out)
    #         unit_cells.append(unit_cell)
    #     except FileNotFoundError:
    #         print(f"Error: File not found for crystal {name}")
    #     except Exception as e:
    #         print(f"Error analyzing crystal {name}: {str(e)}")
    
    # crystal_names = [name[0] for name in crystal_data]
    # results_df = create_crystal_df(crystal_names, 'P1211', weights, dfs_in, unit_cells, dfs_out)
    # print(results_df)
    # intensity_df = create_intensity_df(crystal_names, dfs_out)
    # print(intensity_df)

    # w.write(results_df, intensity_df, 'P1211')
    
    # PLOTS
    # plot the distribution of the FP variable
    # d.plot_hist(intensity_df)
    # d.contour(intensity_df)
    # d.plot_density(intensity_df)
    # d.plot_boxplot(intensity_df)   
    
    # lm = p.linear_model(['VolumeToUnitCellVolRatio'], 'Mean FP', results_df, results_df)
    # residuals = p.plot_residuals(lm, 'VolumeToUnitCellVolRatio', 'Mean FP', results_df, results_df)

    # shows slight linear relationship between FP and Volume to Unit Cell Volume Ratio
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
