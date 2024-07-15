#This script is used after min and npt scripts are run on gromacs, prior to running nvt.
#
#Input is the nvt.xtc, .tpr file used to calculate the average volume of the box
#Then the npt.gro file is modified to account for the new calcualted volume of the box/box dimensions 
#

import os
import numpy as np
import pdb
import subprocess

sizes = [512, 1024, 2048, 4096]
molecules = ['TIP3P','hexane','heptane','octane','decane','pentadecane']
length_b = 5000  #before running, check how long you run npt for and only average last part of simulation for the average. 
length_e = 10000  # the npt simulation was run for 10 ns, average last 5 nm, probably dont need to run that long next time

for molecule in molecules:
    for size in sizes:

        #Note that 
        #Note that npt was run semiisotropic (if you repeat simulation change npt pressure coupling to isotropic), so 
        # you need dimensions in x and z. If isotropic you can switch to total volume fluctuations and extract box size there
        command_x = f"echo 18 | gmx energy -f npt_{molecule}_{size}.edr -b {length_b} -e {length_e} -o {molecule}_{size}_x.xvg"
        subprocess.run(command_x, shell=True, check=True)

        command_z = f"echo 20 | gmx energy -f npt_{molecule}_{size}.edr -b {length_b} -e {length_e} -o {molecule}_{size}_z.xvg"
        subprocess.run(command_z, shell=True, check=True)

        #Calculate the average x and z box lengths over the trajectory
        x_len = np.loadtxt(f'{molecule}_{size}_x.xvg', comments=['#','@'])
        x_ps = x_len[:,0]
        x_nm = x_len[:,1]
        x= np.mean(x_nm)
        # print(x)

        z_len = np.loadtxt(f'{molecule}_{size}_z.xvg', comments=['#','@'])
        z_ps = z_len[:,0]
        z_nm = z_len[:,1]
        z= np.mean(z_nm)
        #  print(z)

        #Use editconf to change the box dimensions of npt.gro file to the average 
        command = f"gmx editconf -f npt_{molecule}_{size}.gro -box {x} {x} {z} -o npt_box_{molecule}_{size}.gro" 
        subprocess.run(command, shell=True, check=True)