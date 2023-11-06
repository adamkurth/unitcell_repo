rm(list=ls())
# libraries
library(ggplot2)
setwd("/Users/adamkurth/Documents/vscode/CXFEL_Image_Analysis/CXFEL/unitcell_project/")
getwd()

spacegroup_P1211 <- read.table("P1211_output/P1211_crystal_data.txt", header = TRUE, fill = TRUE)
# For some reason there's empty NA columns when reading.
spacegroup_P1211 <- spacegroup_P1211[, sapply(spacegroup_P1211, function(col) !all(is.na(col)))]
intensity_P1211 <- read.table("P1211_output/P1211_intensity_data.txt", header = TRUE)
phases_P1211  <- read.table("P1211_output/P1211_phases_data.txt", header = TRUE)  
colnames(spacegroup_P1211) = c("PDB_ID", "Spacegroup", "a", "b", "c", "alpha", "beta", "gamma", "UnitCellVolume", "CrystalVolume", "VolumeToUnitCellVolRatio", "Mean_Intensity", "Max_Intensity", "Min_Intensity", "Max-Min_Difference")
colnames(intensity_P1211) = sub("^X", "", colnames(intensity_P1211))
colnames(phases_P1211) = sub("^X", "", colnames(phases_P1211))

head(spacegroup_P1211)
head(intensity_P1211)
head(phases_P1211)

# intensity_stnd = scale(intensity_P1211)
spacegroup_P1211$Ratio_Stnd = scale(spacegroup_P1211$VolumeToUnitCellVolRatio)
spacegroup_P1211$Mean_Int_Stnd = scale(spacegroup_P1211$Mean_Intensity)
# ploynomial looking??
plot(Ratio_Stnd ~ Mean_Int_Stnd , data = spacegroup_P1211)



