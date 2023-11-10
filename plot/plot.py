import numpy as np
import matplotlib.pyplot as plt
import statsmodels as sm
import statsmodels.api as sm
import scipy.stats as stats
import sys
import os
import math


def plot(x_col, y_col, df, x_label=None, y_label=None, title=None):
    plt.scatter(df[x_col], df[y_col])
    if x_label:
        plt.xlabel(x_label)
    else:
        plt.xlabel(x_col)
    if y_label:
        plt.ylabel(y_label)
    else:
        plt.ylabel(y_col)
    if title:
        plt.title(title)
    else:
        plt.title(f"{y_col} vs {x_col}")
    plt.show()

def linear_model(x_cols, y_col, df_in, df_out):
    # Create a linear regression model
    X = df_in[x_cols]
    y = df_out[y_col]
    X = sm.add_constant(X)
    model = sm.OLS(y, X).fit()

    # Print the model summary
    print(model.summary())

    # Creating a 3x2 grid of subplots
    fig, ax = plt.subplots(3, 2, figsize=(15, 15))  # Adjust the figsize accordingly
    plt.subplots_adjust(hspace=0.4, wspace=0.4)

    for i, col in enumerate(x_cols):
        if i < 2:  # Since we only have two spots in the first row
            ax[0, i].scatter(df_in[col], y)
            ax[0, i].set_xlabel(col)
            ax[0, i].set_ylabel(y_col)
            ax[0, i].set_title(f'{col} vs {y_col}')
        
    # Q-Q plot of the residuals
    sm.graphics.qqplot(model.resid, line='45', fit=True, ax=ax[1, 0])
    ax[1, 0].set_title('Q-Q Plot of the Residuals')

    # Residuals vs. Fitted plot
    ax[1, 1].scatter(model.fittedvalues, model.resid)
    ax[1, 1].axhline(y=0, color='red', linestyle='dashed')
    ax[1, 1].set_xlabel('Fitted values')
    ax[1, 1].set_ylabel('Residuals')
    ax[1, 1].set_title('Residuals vs. Fitted')

    # Leverage vs. Residuals squared plot
    sm.graphics.plot_leverage_resid2(model, ax=ax[2, 0])
    ax[2, 0].set_title('Leverage vs. Residuals squared')

    # Influence plot
    sm.graphics.influence_plot(model, ax=ax[2, 1], criterion="cooks")
    ax[2, 1].set_title('Influence Plot')
    plt.show()

    return model

def plot_residuals(model, x_cols, y_col, df_in, df_out):
    """
    Plots the residuals of a linear regression model.

    Parameters:
    model (statsmodels.regression.linear_model.RegressionResultsWrapper): The linear regression model.
    x_col (str): The name of the independent variable column.
    y_col (str): The name of the dependent variable column.
    df (pandas.DataFrame): The DataFrame containing the data.

    Returns:
    None
    """
    fig, ax = plt.subplots()
    ax.scatter(df_in[x_cols], model.resid)
    ax.axhline(y=0, color='r', linestyle='-')
    ax.set_xlabel(x_cols)
    ax.set_ylabel("Residuals")
    ax.set_title(f"Residuals vs {x_cols}")
    err = np.std(model.resid)
    ax.errorbar(df_in[x_cols], model.resid, yerr=err, fmt='o', color='blue', ecolor='green', capsize=5, capthick=2)
    for i, val in enumerate(model.resid):
        ax.annotate(round(val, 2), (df_in[x_cols][i], val), textcoords="offset points", xytext=(0,10), ha='center')
    plt.show()
    

# shows slight linear relationship between FP and Volume to Unit Cell Volume Ratio
# linear_model(['VolumeToUnitCellVolRatio'], 'Mean FP', results_df, results_df)

# conditition number is large, indicates stable model
# lm = linear_model(['VolumeToUnitCellVolRatio'], 'Mean FP', results_df, results_df)

# plot_residuals(lm, 'VolumeToUnitCellVolRatio', 'Mean FP', results_df, results_df)