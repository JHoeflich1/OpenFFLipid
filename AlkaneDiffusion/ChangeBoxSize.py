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

sizes = [512, 1024, 2048]
molecules = ['pentane','hexane','heptane','octane','decane','pentadecane', 'water']
# all simulations were run for 2 ns of NPT. You only want to average box size after it has relaxed. For each simulation, plot the total volume over the full trajectory 
for molecule in molecules:
    for size in sizes:
        command_V = f"echo Volume | gmx energy -f npt1_{molecule}_{size}.edr -o {molecule}_{size}_V.xvg"
        subprocess.run(command_V, shell=True, check=True)

        #plot volume 
        V_len = np.loadtxt(f'{molecule}_{size}_V.xvg', comments=['#','@'])
        V_ps = V_len[:,0]
        V_nm3 = V_len[:,1]

        plt.figure()
        plt.plot(V_ps, V_nm3, label=f"{molecule} {size}")
        plt.xlabel('Time (ps)')
        plt.ylabel('Volume (nm^3)')
        plt.title(f'Volume vs Time for {molecule} {size}')
        plt.legend()
        plt.savefig(f"{molecule}_{size}.png")
        plt.close()

#####
#run above code, determine the time after volume is fully relaxed/equilibrated, and run the rest of the code

# length_b = 5000  #before running, check how long you run npt for and only average last part of simulation for the average. 
# length_e = 10000  # the npt simulation was run for 10 ns, average last 5 nm, probably dont need to run that long next time

# for molecule in molecules:
#     for size in sizes:
#         #what gmx command gives volume?
#         command_V = f"echo 18 | gmx energy -f npt1_{molecule}_{size}.edr -b {length_b} -e {length_e} -o {molecule}_{size}_V.xvg"
#         subprocess.run(command_V, shell=True, check=True)

#         #From the volume, calculate the length of the box 
#         V_len = np.loadtxt(f'{molecule}_{size}_V.xvg', comments=['#','@'])
#         V_ps = V_len[:,0]
#         V_nm3 = V_len[:,1]
#         Vavg= np.mean(V_nm3)
#         x = Vavg**(1/3)
#         print(x)


#         #Use editconf to change the box dimensions of npt.gro file to the average 
#         command = f"gmx editconf -f npt1_{molecule}_{size}.gro -box {x} {x} {x} -o npt1_box_{molecule}_{size}.gro" 
#         subprocess.run(command, shell=True, check=True)