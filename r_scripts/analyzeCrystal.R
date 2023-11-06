analyzeCrystal = function(pdb_name, input_name, output_name){
  setwd( "/Users/adamkurth/Documents/vscode/CXFEL Image Analysis/CXFEL/unitcell_project")
  # retrieve the pdb data
  pdb_path = file.path(getwd(), "formatted_data/pdb_files", pdb_name)
  pdb_lines = readLines(pdb_path, encoding = "UTF-8")  
  cryst_line = grep("CRYST1", pdb_lines, value = TRUE)
  cryst_components = strsplit(cryst_line, "\\s+")
  #dimension lengths
  a = as.numeric(cryst_components[[1]][2])
  b = as.numeric(cryst_components[[1]][3])
  c = as.numeric(cryst_components[[1]][4])
  #angles (in degrees)
  alpha = as.numeric(cryst_components[[1]][5])
  beta = as.numeric(cryst_components[[1]][6])
  gamma = as.numeric(cryst_components[[1]][7])
  
  #retrieve the input data
  input_path = file.path(getwd(), "formatted_data/input_data/", input_name)
  input_data = read.table(input_path)
  input_data = na.omit(input_data)
  
  # check if all columns are not empty
  if (all(sapply(input_data, function(x) !is.null(x) && !all(is.na(x))))){
    h = input_data$h  
    k = input_data$k 
    l = input_data$l
    FREE = input_data$FREE
    FP = input_data$FP
    SIGFP = input_data$SIGFP
    FC = input_data$FC
    PHIC = input_data$PHIC
    FC_ALL = input_data$FC_ALL
    PHIC_ALL = input_data$PHIC_ALL
    FWT = input_data$FWT
    PHWT = input_data$PHWT
    DELFWT = input_data$DELFWT
    PHDELWT = input_data$PHDELWT
    FOM = input_data$FOM
    FC_ALL_LS = input_data$FC_ALL_LS
    PHIC_ALL_LS = input_data$PHIC_ALL_LS  
    nrow(PHIC_ALL_LS)
    df_in = data.frame(h,k,l,FREE, FP, SIGFP, FC, PHIC, FC_ALL, PHIC_ALL, FWT, PHWT, DELFWT, PHDELWT , FOM, FC_ALL_LS, PHIC_ALL_LS)
    
    # calculate unit cell volume (monoclinic)
    unitcell_vol = a*b*c
    
    # volume of parallelepiped based on unit cell parameters (in degrees) 
    crystal_vol = a * b * c * sqrt(1 - cos(alpha)^2 - cos(beta)^2 - cos(gamma)^2 + 2 * cos(alpha) * cos(beta) * cos(gamma))
    
    df_in$UnitCellVolume = unitcell_vol
    df_in$CrystalVolume = crystal_vol
    df_in$VolumeToUnitCellVolRatio = df_in$CrystalVolume / df_in$UnitCellVolume
  } else {
    df_in = NULL
    message("Input data contains empty columns. Skipping creation of df_in.")
  }
  
  #retrieve the output data
  output_path = file.path(getwd(), "formatted_data/output_data/", output_name)
  output_data = read.table(output_path)
  output_data = na.omit(output_data)
  
  # check if all columns are not empty
  if (all(sapply(output_data, function(x) !is.null(x) && !all(is.na(x))))){
    h = output_data$h 
    k = output_data$k 
    l = output_data$l
    FREE = output_data$FREE
    FP = output_data$FP 
    SIGFP = output_data$SIGFP
    FCalc = output_data$FCalc 
    PHICalc = output_data$PHICalc
    
    df_out = data.frame(h,k,l, FREE, FP, SIGFP, FCalc, PHICalc)
  } else {
    df_out = NULL
    message("Output data contains empty columns. Skipping creation of df_out.")
  }
  
  return(list(df_in = df_in, df_out = df_out))
}
