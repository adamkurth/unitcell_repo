
# Unitcell Repository

## Overview
This document explains the motivations and details of a project centered around the field of crystallography, specifically focusing on protein crystallography. The project aims to understand the three-dimensional structures and interactions of proteins using X-ray diffraction, a technique with significant implications in drug and vaccine development.

## Key Components

### Potential Predictors
- **Mean.Intensity:** Average intensity of diffracted X-rays, providing insights into the crystal structure's overall electron density.
- **Mean.Phase:** Average measure of phase angles in diffracted waves, essential for determining electron density maps.
- **Unit Cell Parameters (a, b, c, alpha, beta, gamma):** Fundamental for determining the system and size of the crystal lattice.
- **Unitcell.Volume:** Indicates the crystal structure's density.
- **Packing.Density:** Ratio of space within the crystal occupied by atoms.
- **Space Group:** Indicates the crystal's group symmetry, key to understanding its structure.
- **Calc.Structure.Weight:** Total mass of atoms in a molecule, particularly important in protein crystals.

### Structure Factor Calculation
Utilizing CCP4, a crystallography program, for calculating the structure factor, which is critical for interpreting protein structure characteristics and the diffraction pattern.

### Research Question
Investigating the linear dependence of X-ray diffraction intensities on the number of unit cells exposed to X-rays and the effect of packing density on diffracted intensity.

### Model
A linear model considering factors like calculated structure weight and packing density to predict mean intensity.

### Statistical Analysis
- **Cook’s Distance:** Used for outlier detection in multiple regression, ensuring independence of residuals.
- **Durbin-Watson (D-W) Test Statistic:** Used to check for autocorrelation among residuals, ensuring no significant correlation is present.
- **Influence Plots:** Visualized influence points in data using `car::influencePlots()`.
- **Residual Plots:** Visualized the residual fit of model using `plot()` as well as `ggplot()`.
- **Leverage Plots:** Used leverage plots to see outlier points in the data with creating a model, then using `plot()` function.

### Confidence/Prediction Interval Plots
Analyzing the relationship between calculated structure weight, packing density, and mean intensity. The analysis indicates a wide variance in the linear relationship, with the possibility of a zero slope.

## Motivation
The primary motivation is to advance the understanding of protein structures through X-ray crystallography, contributing to the field of drug discovery and vaccine development. The choice of predictors and the modeling approach are designed to provide deeper insights into the intricate details of crystal structures.

## Conclusion
This README serves as a guide to understanding the objectives, methods, and findings of the project, highlighting its significance in the realm of protein crystallography and its potential impact on medical research and development.