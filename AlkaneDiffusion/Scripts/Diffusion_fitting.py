import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats

# Load data for DS and box size
Ds_df = pd.read_pickle("diffusion_final.pkl")
box = pd.read_csv('box_sizes.csv')

diffusion_df = pd.DataFrame(columns=['molecule', 'average diffusion curvefit', 'stdev diffusion curvefit', 'average diffusion linear fit'])

# Define the model function for curve fitting
def model(x, A, D_inf):  # x is 1/L, A and D_inf are parameters to fit
    return A * x + D_inf

fit_results = {}

# Fit the model for each molecule
for molecule in box['molecule'].unique():
    # Filter for the current molecule
    molecule_data = box[box['molecule'] == molecule]
    
    Size_data = []
    DS_data = []
    DS_data_err = []
    
    for _, row in molecule_data.iterrows():
        box_size = row['box_length_avg']
        mol_size = row['size']
        
        # Map the (molecule, mol_size) combination to Ds_df
        matrix = Ds_df.loc[(molecule, mol_size)]
        Size_data.append(box_size)
        DS_data.append(matrix['Ds'])  # save Ds
        DS_data_err.append(matrix['Stderr'])  # save Stderr column
    
    # Convert to numpy arrays
    x_data = np.array(Size_data)
    y_data = np.array(DS_data)
    y_err = np.array(DS_data_err)
    
    # Perform curve fitting
    guess = [1.0, 0.0]  # initial guess for model
    params, covariance = curve_fit(
        model, 1 / x_data, y_data, sigma=y_err, p0=guess, absolute_sigma=True
    )
    perr = np.sqrt(np.diag(covariance))

    fit_results[molecule] = {'params': params, 'covariance': covariance, 'errors': perr}
    
    # Plot
    plt.errorbar(1 / x_data, y_data, yerr=y_err, fmt='o', label=f'{molecule} Data')
    x_fit = np.linspace(min(1 / x_data), max(1 / x_data), 100)
    y_fit = model(x_fit, *params)
    plt.plot(x_fit, y_fit, label=f'{molecule} Fit')
    plt.xlabel('1 / Size (nm)')
    plt.ylabel('DS')
    plt.title(f'Diffusion Fitting for {molecule}')
    plt.legend()
    plt.savefig(f'{molecule}_Diffusion_Fitting.png')
    plt.clf()
    
    params_cm2_s = params / 10000  # convert to cm^2/s
    perr_cm2_s = perr / 10000

    new_row_d = pd.DataFrame({
        'molecule': [molecule],
        'average diffusion curvefit': [params_cm2_s[1]],
        'stdev diffusion curvefit': [perr_cm2_s[1]],
        'average diffusion linear fit': [None]
    })
    diffusion_df = pd.concat([diffusion_df, new_row_d], ignore_index=True)

# Perform linear fit as an alternative method
molecules = ['water', 'pentane', 'hexane', 'heptane', 'octane', 'decane', 'pentadecane']
sizes = [512, 1024, 2048]

for molecule in molecules:
    DS_data_line = []
    DS_data_err_line = []
    for size in sizes:
        matrix = Ds_df.loc[(molecule, size)]
        DS_data_line.append(matrix['Ds'])  # save Ds
        DS_data_err_line.append(matrix['Stderr'])  # save Stderr column

    x_data = np.array(sizes)
    x_data_scaled = 1 / x_data**(1 / 3)
    y_data = np.array(DS_data_line)
    y_err = np.array(DS_data_err_line)

    res = stats.linregress(x_data_scaled, y_data)

    plt.plot(x_data_scaled, y_data, 'o', label='Original Data')
    plt.plot(x_data_scaled, res.intercept + res.slope * x_data_scaled, 'r', label=f'Fitted Line\nint = {res.intercept:.4f}')
    plt.xlabel('1 / Size^(1/3)')
    plt.ylabel('DS')
    plt.title(f'Linear Diffusion Fitting for {molecule}')
    plt.legend()
    plt.savefig(f'{molecule}_Linear_Fitting.png')
    plt.clf()
    
    diffusion_df.loc[diffusion_df['molecule'] == molecule, 'average diffusion linear fit'] = res.intercept * 1e-4  # convert to m²/s from cm²/s

diffusion_df.to_csv('diffusion_fitted.csv', index=False)

# Print fitted parameters and covariance for each molecule
for molecule, result in fit_results.items():
    params_m2_s = result['params'] * 1e-4  # convert to m²/s
    perr_m2_s = result['errors'] * 1e-4
    covariance_m2_s = result['covariance'] * 1e-8  # covariance must be divided by 10000 squared
    
    print(f"Molecule: {molecule}")
    print(f"Parameters (A, D_inf) in m²/s: {params_m2_s}")
    print(f"Standard Errors: {perr_m2_s}")
    print(f"Covariance Matrix: {covariance_m2_s}\n")
