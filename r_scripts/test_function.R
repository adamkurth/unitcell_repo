rm(list=ls())
# libraries
library(ggplot2)
library(dplyr)
library(ggplot2) 
library(GGally) 
library(car)
library(tidyr)
library(reshape2)

######## Custom Functions
# -------------------------------------------------------------------------
custom.residual.plot <- function(model) {
  # Check if arg is a linear model
  if (class(model) != "lm") {
    stop("The input must be a linear model object of class 'lm'.")
  }
  fitted_values <- fitted(model)
  residuals <- residuals(model)
  
  # plotting the residuals
  plot(fitted_values, residuals,
       xlab = "Fitted Values",
       ylab = "Residuals",
       main = "Residual Plot",
       pch = 20)
  
  # indicate the baseline
  abline(h = 0, col = "red", lty = 2)
  invisible(list(fitted_values = fitted_values, residuals = residuals))
}

assumption.1.check.linearity <- function(m) {
  cat("Assumption 1: Checking for linearity...\n")
  custom.residual.plot(m)
}

# Assumption 2: Check for excessive outiers.
# Cooks distance: aggregates measure of influence of the ith case on all fitted values.
# does not require repeatedly fitting the regression model omitting the ith case.
# "product" of leverage and outliers.
assumption.2.check.outliers <- function(m) {
  cat("Assumption 2: Checking for no excessive outliers...\n")
  
  cooks_dist <- cooks.distance(m)
  plot(cooks_dist, type = "h",  
       xlab = "Observation Index", 
       ylab = "Cook's Distance",
       main = "Cook's Distance Plot for Outlier Detection")
  
  cutoff <- 4 / length(cooks_dist)
  significant_outliers <- which(cooks_dist > cutoff)
  if(length(significant_outliers) > 0) {
    cat("Observations with a Cook's distance greater than", cutoff, "may be considered influential:\n")
    print(significant_outliers)
  } else {
    cat("No individual observations are significantly influential.\n")
  }
  return(cooks_dist)
}

# Assumption 3: Constant Variance (homoskedasticity)
assumption.3.check.homoskedasticity <- function(m) {
  cat("Assumption 3: Checking for constant variance (homoskedasticity)...\n")
  
  plot(m$fitted.values, sqrt(abs(m$residuals)),
       xlab = "Fitted Values",
       ylab = "Sqrt(|Residuals|)",
       main = "Scale-Location Plot")
  
  ncv <- car::ncvTest(m)
  cat(sprintf("Non-Constant Variance Score Test Statistic: %f\n", ncv$statistic))
  cat(sprintf("p-value: %f\n", ncv$p.value))
}


# Assumption 4: Normally Distributed Errors
assumption.4.check.normality <- function(m) {
  cat("Assumption 4: Checking for normally distributed errors...\n")
  qqnorm(m$residuals)
  qqline(m$residuals)
}

# Assumption 5: Independence of Errors
assumption.5.check.independence <- function(m) {
  cat("Assumption 5: Checking for independence of errors...\n")
  
  dw_test <- car::durbinWatsonTest(m)
  cat(sprintf("Durbin-Watson Test Statistic: %f\n", dw_test$statistic))
  cat(sprintf("p-value: %f\n", dw_test$p.value))
  return(dw_test)
}

# Assumption 6: Multicollinearity
assumption.6.check.multicollinearity <- function(m) {
  cat("Assumption 6: Checking for multicollinearity...\n")
  
  vif_values <- car::vif(m)
  print(vif_values)
}
# -------------------------------------------------------------------------

check.assumptions <- function(m) {
  if (class(m) != "lm") {
    stop("Input must be a linear model with class 'lm'.")
  }
  
  assumption.1.check.linearity(m)
  readline(prompt="Press [Enter] to continue...\n")
  
  assumption.2.check.outliers(m)
  readline(prompt="Press [Enter] to continue...\n")
  
  assumption.3.check.homoskedasticity(m)
  readline(prompt="Press [Enter] to continue...\n")
  
  assumption.4.check.normality(m)
  readline(prompt="Press [Enter] to continue...\n")
  
  assumption.5.check.independence(m)
  readline(prompt="Press [Enter] to continue...\n")
  
  assumption.6.check.multicollinearity(m)
}

fix_header_df1 <- function(df) {
  # Reorder the columns with "Calculated Structure Weight" as the 3rd column
  colnames(df) <- c("Spacegroup.ID", "PDB.ID", "Calculated.Structure.Weight", "a", "b", "c", "Alpha", "Beta", "Gamma", "Unit.Vol", "Crystal.Vol", "Vol.Unit.Ratio", "Mean.Intensity", "Max.Intensity", "Min.Intensity", "MaxMin.Intensity.Dif", "Mean.Phase", "Max.Phase", "Min.Phase", "MaxMin.Phase.Difference")
  return(df)
}
####### SETUP
# -------------------------------------------------------------------------
setwd("/Users/adamkurth/Documents/vscode/CXFEL_Image_Analysis/CXFEL/unitcell_repo/r_scripts")
# getwd()
# -------------------------------------------------------------------------
df1_P1211 = read.table("../P1211_output/P1211_crystal_data.txt", header = TRUE, fill = TRUE)
# Remove empty NA columns when reading 
df1_P1211 = df1_P1211[, sapply(df1_P1211, function(col) !all(is.na(col)))]
df2_P1211 = read.table("../P1211_output/P1211_intensity_data.txt", header = TRUE)
df3_P1211 = read.table("../P1211_output/P1211_phases_data.txt", header = TRUE)  
# -------------------------------------------------------------------------
df1_P121 = df1_P1211[, sapply(df1_P1211, function(col) !all(is.na(col)))]
df2_P121 = read.table("../P121_output/P121_intensity_data.txt", header = TRUE)
df3_P121 = read.table("../P121_output/P121_phases_data.txt", header = TRUE)
# -------------------------------------------------------------------------
df1_C121 = df1_P1211[, sapply(df1_P1211, function(col) !all(is.na(col)))]
df2_C121 = read.table("../C121_output/C121_intensity_data.txt", header = TRUE)
df3_C121 = read.table("../C121_output/C121_phases_data.txt", header = TRUE) 

df1_P1211 = fix_header_df1(df1_P1211)

# extra column in df1
df1_P1211 = df1_P1211[, -1]
df1_P121 = df1_P121[,-1]
df1_C121 = df1_C121[,-1]
# -------------------------------------------------------------------------
########## FORMATTING 
# remove the improper colnames
colnames(df1_P1211) = c("PDB.ID", "Spacegroup", "Calc.Structure.weight", "a", "b", "c", "Alpha", "Beta", "Gamma", "Unit.Vol", "Crystal.Vol", "Vol.Unit.Ratio", "Mean.Intensity", "Max.Intensity", "Min.Intensity", "MaxMin.Intensity.Dif", "Mean.Phase", "Max.Phase", "Min.Phase", "MaxMin.Phase.Difference")
colnames(df1_P1211) = sub("^X", "", colnames(df1_P1211))
colnames(df2_P1211) = sub("^X", "", colnames(df2_P1211))
colnames(df2_P1211) <- sub("^X", "", colnames(df2_P1211))
colnames(df2_P121) <- sub("^X", "", colnames(df2_P121))
colnames(df2_C121) <- sub("^X", "", colnames(df2_C121))
colnames(df3_P1211) = sub("^X", "", colnames(df3_P1211))
colnames(df3_P121) = sub("^X", "", colnames(df3_P121))
colnames(df3_C121) = sub("^X", "", colnames(df3_C121))

df1_P1211 = as.data.frame(df1_P1211)
df2_P1211 = as.data.frame(df2_P1211)
df3_P1211 = as.data.frame(df3_P1211)

df1_P121 = as.data.frame(df1_P121)
df2_P121 = as.data.frame(df2_P121)
df3_P121 = as.data.frame(df3_P121)

df1_C121 = as.data.frame(df1_C121)
df2_C121 = as.data.frame(df2_C121)
df3_C121 = as.data.frame(df3_C121)

new_colnames = c("PDB.ID", "Spacegroup.ID", "Calc.Structure.weight", "a", "b", "c", "Alpha",
                 "Beta", "Gamma", "Unit.Vol", "Crystal.Vol", "Vol.Unit.Ratio",
                 "Mean.Intensity", "Max.Intensity", "Min.Intensity", "MaxMin.Intensity.Dif",
                 "Mean.Phase", "Max.Phase", "Min.Phase", "MaxMin.Phase.Difference")
colnames(df1_P1211) = new_colnames
colnames(df1_P121) = new_colnames
colnames(df1_C121) = new_colnames
# -------------------------------------------------------------------------
#### FOR CONVERTING DF2, DF3

transform_df <- function(df) {
  # Add a row_id to track the original row order
  df$row_id <- seq_len(nrow(df))
  
  # Melt the dataframe by 'row_id' to keep track of the original rows
  long_df <- melt(df, id.vars = 'row_id')
  
  # Cast the melted dataframe into a wide format with each PDB ID as a row
  # The values will be spread across columns named with the original row numbers
  new_df <- dcast(long_df, variable ~ row_id, value.var = "value")
  
  # Return the transformed dataframe
  return(new_df)
}

df2_P1211 = transform_df(df2_P1211)
df2_P121 = transform_df(df2_P121)
df2_C121 = transform_df(df2_C121)

df3_P1211 = transform_df(df3_P1211)
df3_P121 = transform_df(df3_P121)
df3_C121 = transform_df(df3_C121)

df1_P1211$Spacegroup.ID = 0
df2_P1211$Spacegroup.ID = 0
df3_P1211$Spacegroup.ID = 0
# Initialize all values in P121 to 1
df1_P121$Spacegroup.ID = 1
df2_P121$Spacegroup.ID = 1
df3_P121$Spacegroup.ID = 1
# Initialize all values in C121 to 2
df1_C121$Spacegroup.ID = 2
df2_C121$Spacegroup.ID = 2
df3_C121$Spacegroup.ID = 2
# -------------------------------------------------------------------------
# Reorder df1, df2, df3 and append spacegroup.id to 1st column
df1_P1211 = data.frame(Spacegroup.ID = 0, df1_P1211)
df1_P1211 = df1_P1211 %>% select(-Spacegroup.ID.1)
df1_P121 = data.frame(Spacegroup.ID = 1, df1_P121)
df1_P121 = df1_P121 %>% select(-Spacegroup.ID.1)
df1_C121 = data.frame(Spacegroup.ID = 2, df1_C121)
df1_C121 = df1_C121 %>% select(-Spacegroup.ID.1)

df2_P1211 = data.frame(Spacegroup.ID = 0, df2_P1211)
df2_P121 = data.frame(Spacegroup.ID = 1, df2_P121)
df2_C121 = data.frame(Spacegroup.ID = 2, df2_C121)

df3_P1211 = data.frame(Spacegroup.ID = 0, df3_P1211)
df3_P121 = data.frame(Spacegroup.ID = 1, df2_P121)
df3_C121 = data.frame(Spacegroup.ID = 2, df3_C121)
# -------------------------------------------------------------------------
df1_list <- list(df1_P1211, df1_P121, df1_C121)
df2_list <- list(df2_P1211, df2_P121, df2_C121)
df3_list <- list(df3_P1211, df3_P121, df3_C121)
all_df1 <- data.frame(); all_df2 <- data.frame(); all_df3 <- data.frame()
all_df1 <- bind_rows(df1_P1211, df1_P121, df1_C121)
all_df2 <- bind_rows(df2_P1211, df2_P121, df2_C121)
all_df3 <- bind_rows(df3_P1211, df3_P121, df3_C121)
head(all_df1)
head(all_df2)
head(all_df3)
# -------------------------------------------------------------------------