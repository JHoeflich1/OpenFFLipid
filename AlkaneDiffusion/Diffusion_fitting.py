#This file takes the Ds and Sterr multiindex dataframe and fits it to the function:
#   D(1/L) = A*(1/L) + D_inf 
# and returns D_inf. 
# the function is fit using WLS 

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Load your Ds_df from pickle file
Ds_df = pd.read_pickle("DS_final.pkl")

molecules = ['TIP3P', 'hexane', 'heptane', 'octane', 'decane', 'pentadecane']
sizes = {0:512, 0:1024, 0:2048, 0:4096}  # change to a dictionary for box sizes 

def model(x, A, D_inf):
    return A * x + D_inf

fit_results = {}

for molecule in molecules:
    Size_data = []
    DS_data = []
    DS_data_err = []
    
    # get data for the current molecule and sizes
    for box_size, mol_size in sizes.items():
        matrix = Ds_df.loc[(molecule, mol_size)]
        Size_data.append(box_size)
        DS_data.append(matrix['DS'])
        DS_data_err.append(matrix['Stderr'])
        print(f"Molecule: {molecule}, Size: {box_size}, DS: {matrix['DS']}, Stderr: {matrix['Stderr']}")

    
    # Perform curve fitting for the current molecule
    x_data = np.array(Size_data)
    y_data = np.array(DS_data)
    y_err = np.array(DS_data_err)
    
    initial_guess = [1.0, 0.0]  # Initial guess for parameters A and D_inf
    
    params, covariance = curve_fit(model, 1 / x_data, y_data, sigma=y_err, p0=initial_guess, absolute_sigma=True)
    perr = np.sqrt(np.diag(covariance))
    print(f"Standard errors in parameters for {molecule}: {perr}")
    
    fit_results[molecule] = {'params': params,'covariance': covariance}
    
    # plot
    plt.errorbar(x_data, y_data, yerr=y_err, label=f'{molecule} Data')
    x_fit = np.linspace(min(x_data), max(x_data), 100)
    y_fit = model(1.0 / x_fit, *params)
    plt.plot(x_fit, y_fit, label=f'{molecule} Fit')
    plt.xlabel('Size (nm)')
    plt.ylabel('DS')
    plt.legend()
    plt.savefig('')

# Print fitted parameters and covariance for each molecule
for molecule, result in fit_results.items():
    print(f"Molecule: {molecule}, Parameters (A, D_inf): {result['params']}, Covariance: {result['covariance']}")
