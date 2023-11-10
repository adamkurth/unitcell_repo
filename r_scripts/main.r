rm(list=ls())
# libraries
library(ggplot2)
library(dplyr)
library(ggplot2) 
library(GGally) 
library(car)


####### SETUP SPACEGROUP P1211
# -------------------------------------------------------------------------

setwd("/Users/adamkurth/Documents/vscode/CXFEL_Image_Analysis/CXFEL/unitcell_repo/r_scripts")
getwd()

spacegroup_P1211 = read.table("../P1211_output/P1211_crystal_data.txt", header = TRUE, fill = TRUE)
# For some reason there's empty NA columns when reading.
spacegroup_P1211 = spacegroup_P1211[, sapply(spacegroup_P1211, function(col) !all(is.na(col)))]
intensity_P1211 = read.table("../P1211_output/P1211_intensity_data.txt", header = TRUE)
phases_P1211 = read.table("../P1211_output/P1211_phases_data.txt", header = TRUE)  

# extra column in spacegroup_P1211
spacegroup_P1211 = spacegroup_P1211[, -1]

attach(spacegroup_P1211)
spacegroup_P1211
########## FORMATTING P1211

# remove the improper colnames
colnames(spacegroup_P1211) = c("PDB.ID", "Spacegroup", "Calc.Structure.weight", "a", "b", "c", "Alpha", "Beta", "Gamma", "Unit.Vol", "Crystal.Vol", "Vol.Unit.Ratio", "Mean.Intensity", "Max.Intensity", "Min.Intensity", "MaxMin.Intensity.Dif", "Mean.Phase", "Max.Phase", "Min.Phase", "MaxMin.Phase.Difference")
colnames(intensity_P1211) = sub("^X", "", colnames(intensity_P1211))
colnames(phases_P1211) = sub("^X", "", colnames(phases_P1211))

head(spacegroup_P1211)
head(intensity_P1211)
head(phases_P1211)

spacegroup_P1211_df = as.data.frame(spacegroup_P1211)

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

########## MODELS P1211
# -------------------------------------------------------------------------

# model P1211
m_P1211 = lm(Mean.Intensity ~ Calc.Structure.weight + Vol.Unit.Ratio, spacegroup_P1211)
# m_P1211 = lm(scale(Mean.Intensity) ~ scale(Calc.Structure.weight) + scale(Vol.Unit.Ratio), spacegroup_P1211)
plot(m_P1211)

# m_P1211_R = lm(Mean.Intensity ~  Vol.Unit.Ratio, spacegroup_P1211)
# plot(m_P1211_R)
# 
# m_P1211_R = lm(Mean.Intensity ~  Vol.Unit.Ratio, spacegroup_P1211)
# plot(m_P1211_R)

################ Diagnostics
# -------------------------------------------------------------------------

residual.plot = custom.residual.plot(m_P1211)


# check assumptions
# check.assumptions(m_P1211) # all in one function.
# individually
assumption.1.check.linearity(m_P1211)
assumption.2.check.outliers(m_P1211) # most definitely influence of outliers.
assumption.3.check.homoskedasticity(m_P1211)
assumption.4.check.normality(m_P1211)
assumption.5.check.independence(m_P1211)
assumption.6.check.multicollinearity(m_P1211)
# check for multicolinearity again.
vif(m_P1211)

anova(m_P1211)

# Visualize outliers, leverage points, and influential points
car::influencePlot(m_P1211,id=list(labels=row.names(spacegroup_P1211)))

# Outliers detected rows: 4, 53, 59, 5

cooks = cooks.distance(m_P1211)
p = 2
n = nrow(spacegroup_P1211)
percentile = 100*pf(q=cooks, df1=p, df2=n-p)
cooks_df = data.frame(spacegroup_P1211, cooks=round(cooks,4),percentile=round(percentile,1))
plot(cooks)

# retrieve rows that are above threshold
which(cooks_df$cooks>0.04)

spacegroup_P1211_numeric <- spacegroup_P1211_df %>% select(
    -a, -b, -c, -Alpha, -Beta, -Gamma, 
    -Max.Intensity, -Min.Intensity, -MaxMin.Intensity.Dif, 
    -Max.Phase, -Min.Phase, -MaxMin.Phase.Difference
  )
ggpairs(spacegroup_P1211_numeric) # check for linear relationships

# check 1
m1 = lm(Vol.Unit.Ratio ~ Mean.Intensity, data=spacegroup_P1211_numeric)
plot(m1)
vif(m1)
car::influencePlot(m1,id=list(labels=row.names(spacegroup_P1211_numeric)))



m2 = lm(Vol.Unit.Ratio ~ Mean.Intensity + Mean.Phase + Mean.Intensity:Mean.Phase, data=spacegroup_P1211_numeric)
plot(m2)
vif(m2) # higher VIF

spacegroup_P1211_numeric$Vol.Unit.Ratio.c = spacegroup_P1211_numeric$Vol.Unit.Ratio - mean(spacegroup_P1211_numeric$Vol.Unit.Ratio)
spacegroup_P1211_numeric$Mean.Intensity.c = spacegroup_P1211_numeric$Mean.Intensity - mean(spacegroup_P1211_numeric$Mean.Intensity) 
spacegroup_P1211_numeric$Mean.Phase.c = spacegroup_P1211_numeric$Mean.Phase - mean(spacegroup_P1211_numeric$Mean.Phase)

# centered 
m3 = lm(Vol.Unit.Ratio.c ~ Mean.Intensity.c + Mean.Phase.c + Mean.Intensity.c:Mean.Phase.c, data = spacegroup_P1211_numeric)
plot(m3)
vif(m2) # lower VIF

anova(m3,m2)
cooks = cooks.distance(m3)
spacegroup_P1211_numeric$cooks = cooks
which(spacegroup_P1211_numeric$cooks>0.04) # why this threshold?






#   Mean.Phase.c 
# + Mean.Intensity.c:Mean.Phase.c



plot(m_P1211$fitted.values, sqrt(abs(m_P1211$residuals)),
     xlab="Fitted Values",
     ylab= "sqrt(|residuals|)")




#residual plot
plot(m_P1211$residuals)
abline(0,0)

intensity_P1211 = as.matrix(intensity_P1211)
phases_P1211 = as.matrix(phases_P1211)
plot(intensity_P1211, phases_P1211)




detach(spacegroup_P1211)




plot(intensity_P1211 ~ phases_P1211)
 
plot(Mean.Intensity ~ Calc.Structure.weight + Mean.Intensity, data = spacegroup_P1211)
plot(scale(Mean.Intensity) ~ scale(Calc.Structure.weight) + scale(Mean.Intensity), data = spacegroup_P1211)

intensity_stnd = scale(intensity_P1211)
phases_stnd = scale(phases_P1211)

# Calculate mean standardized mean for each pdb column
intensity_stnd_mean_P1211 <- apply(intensity_stnd, 2, function(x) mean(x))
phase_stnd_mean_P1211 <- apply(phases_stnd, 2, function(x) mean(x))

# Append mean standardized mean to spacegroup_P1211
spacegroup_P1211$intensity_stnd_mean_P1211 <- intensity_stnd_mean_P1211
spacegroup_P1211$phase_stnd_mean_P1211 <- phase_stnd_mean_P1211

# ploynomial looking??
plot(Ratio_Stnd ~ Mean_Int_Stnd , data = spacegroup_P1211)




# pairs()
# cor()


##########
# CHANGE TO ADD SPACEGROUP AS 0, 1, 2, ...

