#This file takes the Ds and Sterr multiindex dataframe and fits it to the function:
#   D(1/L) = A*(1/L) + D_inf 
# and returns D_inf. 
# the function is fit using WLS 

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

Ds_df = pd.read_pickle("DS_final.pkl")

molecules = ['TIP3P','hexane', 'heptane', 'octane', 'decane', 'pentadecane']
sizes = [512, 1024, 2048, 4096]  #These need to be changed to size of the boxes is nm

# for molecule in molecules:
#     for size in sizes:
#         matrix = msd_df.loc[(molecule, size)]
#         DS = matrix[0]
#         Sterr = matrix[1]
#         print(f"Molecule: {molecule}, Size: {size}, DS: {DS}, Stderr: {Stderr}")

def model(params, x): # where x is 1/L
    A, D_inf = params
    return A * x + D_inf

def weighted_residuals(params, x, y, weights):
    return (y - model(params, x)) / weights

fig, axs = plt.subplots(len(molecules), 1, figsize=(8, 6), sharex=True)

for i, molecule in enumerate(molecules):
    ax = axs[i]
    for size in sizes:
        matrix = Ds_df.loc[(molecule, size)]
        DS = matrix['Ds']
        Stderr = matrix['Stderr']
        L_inverse = 1 / size  # Assuming size corresponds to 1/L

        # plot DS with error bars
        ax.errorbar(L_inverse, DS, yerr=Stderr, fmt='o', label=f'Size {size}')

        # Perform weighted least squares fitting
        weights = 1 / Stderr**2
        initial_params = np.array([1.0, 0.0])  # Initial guess for A and D_inf
        result = minimize(fun=lambda params: np.sum(weighted_residuals(params, L_inverse, DS, weights)**2),
                          x0=initial_params)

        fitted_params = result.x
        A_fit, D_inf_fit = fitted_params

        # Plot the fitted line
        x_fit = np.linspace(min(L_inverse), max(L_inverse), 100)
        y_fit = model(fitted_params, x_fit)
        ax.plot(x_fit, y_fit, label=f'Fit: D_inf = {D_inf_fit:.2f}', linestyle='--')

    ax.set_title(f'Molecule: {molecule}')
    ax.set_xlabel('1/L')
    ax.set_ylabel('DS')
    ax.legend()

plt.tight_layout()
plt.show()

