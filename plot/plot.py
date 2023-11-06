import numpy as np
import matplotlib.pyplot as plt
import statsmodels as sm
import statsmodels.api as sm
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

    # Plot the data and the regression line
    fig, ax = plt.subplots()
    ax.scatter(df_in[x_cols[0]], df_out[y_col])
    ax.plot(df_in[x_cols[0]], model.predict(X), color='red')
    ax.set_xlabel(x_cols[0])
    ax.set_ylabel(y_col)
    for i in range(1, len(x_cols)):
        ax.scatter(df_in[x_cols[i]], df_out[y_col])
        ax.plot(df_in[x_cols[i]], model.predict(X), color='red')
        ax.set_xlabel(x_cols[i])
    ax.set_title(f"{y_col} vs {', '.join(x_cols)}")

    # Plot the error bars
    y_pred = model.predict(X)
    y_err = y - y_pred
    # ax.errorbar(df_in[x_cols[0]], y_pred, yerr=abs(y_err), fmt='o', color='red', ecolor='red', capsize=5, capthick=2, linestyle='dotted')
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