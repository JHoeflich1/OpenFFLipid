#This script is used to determine average density of production run
#
#

import os
import numpy as np
import pdb
import subprocess
import matplotlib.pyplot as plt
import pandas as pd

sizes = [512, 1024, 2048]
alkaneMWs = {
    'pentane': 72.15, # (g/mol) 
    'hexane': 86.17, 
    'heptane':  100.21, 
    'octane':  114.23, 
    'decane': 142.29, 
    'pentadecane':  212.42, 
    'water':  18.01
}
densities =[]


df_densities = pd.DataFrame(columns=['molecule', 'size', 'density_avg', 'density_std'])

for molecule, MW in alkaneMWs.items():
    for size in sizes:
        command_V = f"echo Volume| gmx energy -f npt2_{molecule}_{size}.edr -o {molecule}_{size}_V2.xvg"
        subprocess.run(command_V, shell=True, check=True)

        #From the volume, calculate the length of the box 
        V_len = np.loadtxt(f'{molecule}_{size}_V2.xvg', comments=['#','@'])
        V_ps = V_len[:,0]
        V_nm3 = V_len[:,1]
        Vavg= np.mean(V_nm3)
        Vstd = np.std(V_nm3)
        x = Vavg**(1/3)
        x_std = Vstd * (1/3)*Vavg**(-2/3) #Note that you need error prop for nonlinear scaling nm^3 to nm
        x_2 = round(x,2)
        x_std_2 = round(x_std,3)



        #calcualte the total density of the simulation as well
        density = size*MW/6.02e23 /Vavg *1e21 #molec * moles/molec * g/mol /nm^3 * nm3/ml
        density_std = size * MW / 6.02e23 * (1/1e-21) * Vstd / Vavg**2
        density_2 = round(density, 3)
        density_std_2 = round(density_std, 3)
        densities.append(f'the final average density for {molecule},{size} is {density_2} +/- {density_std_2} g/ml')

        new_row_d = pd.DataFrame({
            'molecule': [molecule],
            'size': [size],
            'density_avg': [density],
            'density_std': [density_std]
        })
        df_densities = pd.concat([df_densities, new_row_d], ignore_index=True)


for density in densities:
    print(density)

df_densities.to_csv('densities_npt2.csv')
