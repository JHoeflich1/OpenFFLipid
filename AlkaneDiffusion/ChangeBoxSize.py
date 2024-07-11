#This script is used after min and npt scripts are run on gromacs, prior to running nvt.
#
#Input is the nvt.xtc, .tpr file used to calculate the average volume of the box
#Then the npt.gro file is modified to account for the new calcualted volume of the box/box dimensions 
#

import os
import numpy as np
import pdb
import subprocess

molecules = ['hexane','heptane','octane','decane','pentadecane']
sizes = [512,1024,2048,4096]
length = 1  #before running, check how long you run npt for and only average last part of simulation for the average. 

for molecule in molecules:
    for size in sizes:

        #Note that the x and y dimension should be the same 
        command_x = f"echo 18 | gmx energy -f {molecule}_{size}.edr -o {molecule}_{size}_x.xvg"
        subprocess.run(command_x, shell=True, check=True)

        command_z = f"echo 20 | gmx energy -f {molecule}_{size}.edr -o {molecule}_{size}_z.xvg"
        subprocess.run(command_z, shell=True, check=True)

        #Calculate the average x and z box lengths over the entire trajcctory 
        x_len = np.loadtxt(f'{molecule}_{size}_x.xvg', comments=['#','@'])
        x_ps = x_len[:,0]
        x_nm = x_len[:,1]
        x= np.mean(x_nm)

        z_len = np.loadtxt(f'{molecule}_{size}_z.xvg', comments=['#','@'])
        z_ps = z_len[:,0]
        z_nm = z_len[:,1]
        z= np.mean(z_nm)

        #Use editconf to change the box dimensions of npt.gro file to the average 
        command = f"gmx editconf -f npt_{molecule}_{size}.gro -box {x} {x} {z}"
        subprocess.run(command, shell=True, check=True)