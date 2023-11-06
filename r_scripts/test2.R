rm(list=ls())
setwd("/Users/adamkurth/Documents/vscode/CXFEL Image Analysis/CXFEL/unitcell_project/")

source("analyzeCrystal.R")

crystal1 = analyzeCrystal("104m.pdb","104m_input_formatted.txt", "104m_hkl_out_formatted.txt")
crystal2 = analyzeCrystal("137l.pdb","137l_input_formatted.txt", "137l_hkl_out_formatted.txt")
crystal3 = analyzeCrystal("1a2a.pdb","1a2a_input_formatted.txt", "1a2a_hkl_out_formatted.txt")
crystal4 = analyzeCrystal("1a28.pdb","1a28_input_formatted.txt", "1a28_hkl_out_formatted.txt")
crystal5 = analyzeCrystal("169l.pdb","169l_input_formatted.txt", "169l_hkl_out_formatted.txt")

df_104m_in = crystal1$df_in
df_104m_out = crystal1$df_out

df_137l_in = crystal2$df_in
df_137l_out = crystal2$df_out

df_1a2a_in = crystal3$df_in
df_1a2a_out = crystal3$df_out

df_1a28_in = crystal4$df_in
df_1a28_out = crystal4$df_out

df_169l_in = crystal5$df_in
df_169l_out = crystal5$df_out

# lm_104m = lm(df_104m_in$FC ~ df_104m_in$VolumeToUnitCellVolRatio)
mean(df_104m_in$FC)
mean(df_104m_in$VolumeToUnitCellVolRatio)
# create a dataframe of the mean FC and volume ratio
plot(df_104m_in$FP, residuals(lm_104m)) 
# Observed intensities vs residuals of model.

# CrystalVolume_104m = df_104m_in$CrystalVolume[1]
# UnitCellVolume_104m = df_104m_in$UnitCellVolume[1]
# Ratio_104m = CrystalVolume_104m / UnitCellVolume_104m



