import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from scipy.stats import gaussian_kde
from sklearn.decomposition import PCA
import numpy as np
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D

def plot_hist(values_df):
    # Scale the values to be between 0 and 1
    scaler = StandardScaler()
    
    values_df_stnd = pd.DataFrame(scaler.fit_transform(values_df), columns=values_df.columns)

    # Get the column names
    column_names = values_df_stnd.columns[:-1]

    # Create a figure with a larger size
    fig, ax = plt.subplots(figsize=(12, 10))

    # Plot all columns on one plot
    for column_name in column_names:
        ax.hist(values_df_stnd[column_name], bins=200, alpha=0.5, label=column_name)

        # Add a vertical line at the mean value of the column
        mean_value = values_df_stnd[column_name].mean()
        ax.axvline(mean_value, color='red')

    ax.legend(loc='upper right')
    
    # Set x and y limits
    ax.set_xlim(-6, 6)
    # ax.set_xlim(-2, 5)
    ax.set_ylim(0, 300)
    plt.show()
    return None

def plot_hist_phases(values_df):
    # Scale the values to be between 0 and 1
    scaler = StandardScaler()
    
    values_df_stnd = pd.DataFrame(scaler.fit_transform(values_df), columns=values_df.columns)

    # Get the column names
    column_names = values_df_stnd.columns[:-1]

    # Create a figure with a larger size
    fig, ax = plt.subplots(figsize=(12, 10))

    # Plot all columns on one plot
    for column_name in column_names:
        ax.hist(values_df_stnd[column_name], bins=200, alpha=0.5, label=column_name)

        # Add a vertical line at the mean value of the column
        mean_value = values_df_stnd[column_name].mean()
        ax.axvline(mean_value, color='red')

    ax.legend(loc='upper right')
    
    # Set x and y limits
    ax.set_xlim(-4, 4)
    # ax.set_xlim(-2, 5)
    ax.set_ylim(0, 300)
    plt.show()
    return None

    
def contour(values_df):
    # Scale the values to be between 0 and 1
    scaler = MinMaxScaler()
    values_df_stnd = pd.DataFrame(scaler.fit_transform(values_df), columns=values_df.columns)

    # Perform PCA on the data
    pca = PCA(n_components=2)
    values_df_pca = pd.DataFrame(pca.fit_transform(values_df_stnd), columns=['PC1', 'PC2'])

    # Create a figure with a larger size
    fig, ax = plt.subplots(figsize=(12, 10))

    # Define a color map for the points
    cmap = plt.get_cmap('tab10')

    # Plot the points
    ax.scatter(values_df_pca['PC1'], values_df_pca['PC2'], c=cmap(1))

    # Compute the kernel density estimate
    kde = gaussian_kde(values_df_pca.T)

    # Define a grid of points to evaluate the density at
    x, y = np.mgrid[values_df_pca['PC1'].min():values_df_pca['PC1'].max():200j,
                    values_df_pca['PC2'].min():values_df_pca['PC2'].max():200j]
    positions = np.vstack([x.ravel(), y.ravel()])

    # Evaluate the density at the grid points
    z = np.reshape(kde(positions).T, x.shape)

    # Plot the density
    ax.contourf(x, y, z, cmap='tab10')
    # Set x and y limits
    ax.set_xlim(values_df_pca['PC1'].min(), values_df_pca['PC1'].max())
    ax.set_ylim(values_df_pca['PC2'].min(), values_df_pca['PC2'].max())
    plt.show()
    return None

def plot_density(values_df):
    # Scale the values to be between 0 and 1
    scaler = StandardScaler()
    values_df_stnd = pd.DataFrame(scaler.fit_transform(values_df), columns=values_df.columns)

    # Get the column names
    column_names = values_df_stnd.columns[:-1]

    # Create a figure with a larger size
    _, ax = plt.subplots(figsize=(12, 10))

    # Create density plots for each column
    for column_name in column_names:
        sns.kdeplot(values_df_stnd[column_name], label=column_name)

    ax.legend(loc='upper right')
    
    # Set x and y limits
    ax.set_xlim(-6, 6)
    
    # Set x and y labels
    ax.set_xlabel('Values')
    ax.set_ylabel('Estimated Density')

    plt.show()
    return None

def plot_boxplot(values_df):
    # Create a figure with a larger size
    _, ax = plt.subplots(figsize=(12, 10))

    # Create a boxplot for each column
    ax.boxplot(values_df)

    # Set x and y labels
    ax.set_xlabel('Column')
    ax.set_ylabel('Values')

    plt.show()
    return None
