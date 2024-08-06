#This script is used after min, nvt and npt equilibration scripts are run on gromacs, prior to running nvt production.
#
#Input is the nvt1_{}_{}.xtc, .tpr file used to calculate the average volume of the box
#Then the npt.gro file is modified to account for the new calcualted volume of the box/box dimensions 
#

import os
import numpy as np
import pdb
import subprocess
import matplotlib.pyplot as plt
import pandas as pd

sizes = [512,1024,2048,4096]
molecules = ['water']
# # all simulations were run for 2 ns of NPT. You only want to average box size after it has relaxed. For each simulation, plot the total volume over the full trajectory 
# for molecule in molecules:
#     for size in sizes:
#         command_V = f"echo Volume | gmx energy -f npt1_{molecule}_{size}.edr -o {molecule}_{size}_V.xvg"
#         subprocess.run(command_V, shell=True, check=True)

#         #plot volume 
#         V_len = np.loadtxt(f'{molecule}_{size}_V.xvg', comments=['#','@'])
#         V_ps = V_len[:,0]
#         V_nm3 = V_len[:,1]

#         plt.figure()
#         plt.plot(V_ps, V_nm3, label=f"{molecule} {size}")
#         plt.xlabel('Time (ps)')
#         plt.ylabel('Volume (nm^3)')
#         plt.title(f'Volume vs Time for {molecule} {size}')
#         plt.legend()
#         plt.savefig(f"{molecule}_{size}.png")
#         plt.close()

#####
#run above code, determine the time after volume is fully relaxed/equilibrated, and run the rest of the code

length_b = 4000  #before running, check how long you run npt for and only average last part of simulation for the average. 
length_e = 7000  # the npt simulation was run for 10 ns, average last 5 nm, probably dont need to run that long next time


alkaneMWs = {
    'water':  18.01
}
densities =[]
avg_box = []

molecules = ['water']


df_boxsize = pd.DataFrame(columns=['molecule', 'size', 'box_length_avg', 'box_length_std'])
df_densities = pd.DataFrame(columns=['molecule', 'size', 'density_avg', 'density_std'])

for molecule, MW in alkaneMWs.items():
    for size in sizes:
        #what gmx command gives volume?
        command_V = f"echo Volume| gmx energy -f npt1_{molecule}_{size}.edr -b {length_b} -e {length_e} -o {molecule}_{size}_V.xvg"
        subprocess.run(command_V, shell=True, check=True)

        #From the volume, calculate the length of the box 
        V_len = np.loadtxt(f'{molecule}_{size}_V.xvg', comments=['#','@'])
        V_ps = V_len[:,0]
        V_nm3 = V_len[:,1]
        Vavg= np.mean(V_nm3)
        Vstd = np.std(V_nm3)
        x = Vavg**(1/3)
        x_std = Vstd * (1/3)*Vavg**(-2/3) #Note that you need error prop for nonlinear scaling nm^3 to nm
        x_2 = round(x,2)
        x_std_2 = round(x_std,3)
        avg_box.append(f'the final average box size for {molecule},{size} is {x_2} +/- {x_std_2} nm')

        new_row_b = pd.DataFrame({
            'molecule': [molecule],
            'size': [size],
            'box_length_avg': [x],
            'box_length_std': [x_std]
        })
        df_boxsize = pd.concat([df_boxsize, new_row_b], ignore_index=True)

        #Use editconf to change the box dimensions of npt.gro file to the average 
        command = f"gmx editconf -f npt1_{molecule}_{size}.gro -box {x} {x} {x} -o npt1_box_{molecule}_{size}.gro" 
        subprocess.run(command, shell=True, check=True)

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

for size in avg_box:
    print(size)

df_boxsize.to_csv('box_sizes.csv')
df_densities.to_csv('densities_npt1.csv')
