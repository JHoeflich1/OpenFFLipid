#This file takes the Ds and Sterr multiindex dataframe and fits it to the function:
#   D(1/L) = A*(1/L) + D_inf 
# and returns D_inf. 
# the function is fit using WLS 

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import subprocess
from scipy import stats

# load data for DS and Box size
Ds_df = pd.read_pickle("diffusion_coefficients.pkl")
box = pd.read_csv('box_sizes.csv')


diffusion_df = pd.DataFrame(columns = ['molecule', 'average diffusion curvefit', 'stdev diffusion curvefit', 'average diffusion linear fit'])
# fit a model to the function D(1/L) = A*(1/L) + D_inf
def model(x, A, D_inf): #where x is 1/L (independent variable_), A and D-inf (parameters to fit)
    return A*x + D_inf

fit_results = {}

for molecule in box['molecule'].unique():
    # Filter for current molecule
    molecule_data = box[box['molecule'] == molecule]
    print(molecule_data, "molecule data")
    
    Size_data = []
    DS_data = []
    DS_data_err = []
    
    for _, row in molecule_data.iterrows(): # iteracte over dataframe rows, returns a series for each row
        box_size = row['box_length_avg']
        print(box_size, 'box size')
        mol_size = row['size']
        
        # Map the (molecule, mol_size) combination to Ds_df
        matrix = Ds_df.loc[(molecule, mol_size)]
        Size_data.append(box_size)
        DS_data.append(matrix['Ds'])  # save Ds
        DS_data_err.append(matrix['Stderr'])  # save Stderr column
        # print(f"Molecule: {molecule}, Size: {box_size}, DS: {matrix['Ds']}, Stderr: {matrix['Stderr']}")
    
    # Convert to numpy arrays
    x_data = np.array(Size_data[:3])
    print(x_data, "size x data")
    y_data = np.array(DS_data[:3])
    print(y_data, "Ds y data")
    y_err = np.array(DS_data_err[:3])
    
    # Perform curve fitting
    guess = [1.0, 0.0] #initial guess for model
    params, covariance = curve_fit(model, 1 / x_data, y_data, sigma=y_err, p0=guess, absolute_sigma=False) #weightes determined by sigma
    # params = popt, minimized values
    # covariance = pcov
    perr = np.sqrt(np.diag(covariance))

    fit_results[molecule] = {'params': params, 'covariance': covariance, 'errors': perr}
    
    # Plot
    plt.errorbar(x_data, y_data, yerr=y_err, label=f'{molecule} Data')
    x_fit = np.linspace(min(x_data), max(x_data), 100)
    y_fit = model(1.0 / x_fit, *params)
    plt.plot(x_fit, y_fit, label=f'{molecule} Fit')
    plt.xlabel('Size (nm)')
    plt.ylabel('DS')
    plt.title(f'Diffusion Fitting for {molecule}')
    plt.legend()
    plt.savefig(f'{molecule}_Diffusion_Fitting.png')
    plt.clf()  


    params_cm2_s = params / 10000 #convert to m^2/s to cm^2/s
    perr_cm2_s = perr / 10000

    new_row = pd.DataFrame([{
        'molecule': molecule,
        'average diffusion curvefit': params_cm2_s[1],
        'stdev diffusion curvefit': perr_cm2_s[1],
        'average diffusion linear fit': None
    }])

    diffusion_df = pd.concat([diffusion_df, new_row], ignore_index=True)

### another way to estimate diffusion using a liner fit
molecules = ['pentane', 'hexane', 'heptane', 'octane', 'decane', 'pentadecane']
sizes = [512, 1024, 2048]

for molecule in molecules:
    DS_data_line = []
    DS_data_err_line = []
    for size in sizes:
        matrix = Ds_df.loc[(molecule, size)]
        DS_data_line.append(matrix['Ds'])  # save Ds
        DS_data_err_line.append(matrix['Stderr'])  # save Stderr column

    x_data = np.array(sizes)
    x_data_scaled = 1/x_data**(1/3)
    y_data = np.array(DS_data_line)
    y_err = np.array(DS_data_err_line)

    res = stats.linregress(x_data_scaled, y_data)

    plt.plot(x_data_scaled, y_data, 'o', label='original data')
    plt.plot(x_data_scaled, res.intercept + res.slope*x_data_scaled, 'r', label=f'fitted line\n int = {res.intercept}')
    plt.xlabel('1 / Size^(1/3)')
    plt.ylabel('DS')
    plt.title(f'Linear Diffusion Fitting for {molecule}')
    plt.legend()
    plt.savefig(f'{molecule}_Linear_Fitting.png')
    plt.clf()
    
    diffusion_df.loc[diffusion_df['molecule'] == molecule, 'average diffusion linear fit'] = res.intercept


diffusion_df.to_csv('diffusion_results.csv', index=False)

# Print fitted parameters and covariance for each molecule
for molecule, result in fit_results.items():
    params_cm2_s = result['params'] / 10000  # convert to cm^2/s
    perr_cm2_s = result['errors'] / 10000
    covariance_cm2_s = result['covariance'] / 10000**2  # covariance must be divided by 10000 squared
    
    print(f"Molecule: {molecule}, Parameters (A, D_inf): {params_cm2_s}, Standard Errors: {perr_cm2_s}, Covariance: {covariance_cm2_s}")
