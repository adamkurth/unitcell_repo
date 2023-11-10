rm(list=ls())
library(dplyr)
library(readr)

setwd("/Users/adamkurth/Documents/vscode/CXFEL_Image_Analysis/CXFEL/unitcell_repo/r_scripts")
getwd()

### CRYSTAL DF
# -------------------------------------------------------------------------
createdf = function(spacegroup_ids) {
  # Initialize an empty list to store results
  results = list()
  
  for(spacegroup_id in spacegroup_ids) {
    # Read in the crystal data
    df_crystal = read_table2(paste0("../", spacegroup_id, "_output/", spacegroup_id, "_crystal_data.txt"), 
                             col_types = cols())
    # Read in the intensity data
    df_intensity = read_table2(paste0("../", spacegroup_id, "_output/", spacegroup_id, "_intensity_data.txt"), 
                               col_types = cols())
    # Read in the phases data
    df_phases = read_table2(paste0("../", spacegroup_id, "_output/", spacegroup_id, "_phases_data.txt"), 
                            col_types = cols())
    
    # Remove empty NA columns from df_crystal
    df_crystal = df_crystal[, sapply(df_crystal, function(col) !all(is.na(col)))]
    
    # Remove an extra column if present
    # Assuming you want to remove the first column (which is usually the row names)
    df_crystal = df_crystal[, -1]
    
    # Rename columns appropriately
    colnames(df_crystal) = c("PDB.ID", "Spacegroup", "Calc.Structure.weight", "a", "b", "c", "Alpha", "Beta", 
                             "Gamma", "Unit.Vol", "Crystal.Vol", "Vol.Unit.Ratio", "Mean.Intensity", 
                             "Max.Intensity", "Min.Intensity", "MaxMin.Intensity.Dif", "Mean.Phase", 
                             "Max.Phase", "Min.Phase", "MaxMin.Phase.Difference")
    colnames(df_intensity) = sub("^X", "", colnames(df_intensity))
    colnames(df_phases) = sub("^X", "", colnames(df_phases))
    
    # Store the formatted data frames in a list
    formatted_dfs = list(df_crystal, df_intensity, df_phases)
    
    # Initialize data frames with the same structure but filled with NA
    initialized_dfs = lapply(formatted_dfs, function(df) {
      data.frame(matrix(NA, nrow = nrow(df), ncol = ncol(df), 
                        dimnames = list(NULL, colnames(df))))
    })
    
    # Combine the data frames with data and the initialized data frames
    results[[spacegroup_id]] = list('formatted' = formatted_dfs, 'initialized' = initialized_dfs)
  }
  
  return(results)
}

# Example usage:
spacegroup_ids <- c("P1211", "P121", "C121")
dfs = createdf(spacegroup_ids)

dfs_initialized <- lapply(dfs, function(spacegroup_list) {
  lapply(spacegroup_list, function(df_list) {
    lapply(df_list, function(df) {
      # Assuming the structure (number of rows and columns) is to be preserved
      array(NA, dim = dim(df))
    })
  })
})

dfs = dfs_initialized









c1.sg1 = crystaldf("../P1211_output/P1211_crystal_data.txt")


# Initialize an empty data frame to store the combined data
combined_df <- data.frame()

# Spacegroup identifiers
spacegroup_ids <- c("P1211", "P121", "C121")




c1.sg1 = read.table("../P1211_output/P1211_crystal_data.txt", header = TRUE, fill = TRUE)
i1.sg1 = read.table("../P1211_output/P1211_intensity_data.txt", header=TRUE, fill = TRUE)
p1.sg1 = read.table("../P1211_output/P1211_phase_data.txt", header=TRUE, fill = TRUE)


# Loop through each spacegroup ID to read and combine data
for (spacegroup_id in spacegroup_ids) {
  # Read spacegroup, intensity, and phase data
  spacegroup_data <- read_and_clean(paste0("../", spacegroup_id, "_output/", spacegroup_id, "_crystal_data.txt"))
  intensity_data <- read_and_clean(paste0("../", spacegroup_id, "_output/", spacegroup_id, "_intensity_data.txt"))
  phases_data <- read_and_clean(paste0("../", spacegroup_id, "_output/", spacegroup_id, "_phases_data.txt"))
  
  # Combine the data for the current spacegroup
  spacegroup_df <- cbind(spacegroup_data, intensity_data, phases_data)
  spacegroup_df
  # Add a column for the spacegroup ID
  # P1211 = 0, P121 = 1, C121 = 2
  spacegroup_df$SpaceGroupID <- match(spacegroup_id, spacegroup_ids) - 1  
  
  # Combine with the main dataframe
  combined_df <- rbind(combined_df, spacegroup_df)
}

combined_df

# Now combined_df contains all the data with a SpaceGroupID to distinguish them

