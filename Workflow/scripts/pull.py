import numpy as np
import pandas as pd
import os
import multiprocessing
import glob

def pullLipid(Lipid_name,lipid_gro, lipid_top,pull_atom):
    '''
    Runs the gmx pulling command to straighten out the lipid and cleans up generated files.
    
    Parameters:
    lipid_gro : str - The GRO file for the lipid.
    lipid_top : str - The TOP file for the lipid.
    pull_atom : int/str - The atom index or identifier from which the lipid will be pulled.
    '''
    command_ndx = f"gmx make_ndx -f {lipid_gro} -o {pull_atom}.ndx"
    os.system(f"echo 'a {pull_atom}\nname 3 pull_atom\nq\n' | {command_ndx}")

    command_grompp = f"gmx grompp -f scripts/pull.mdp -c {lipid_gro} -p {lipid_top} -n {pull_atom}.ndx -o pull.tpr -maxwarn 1" 
    os.system(command_grompp)

    command_pull = "gmx mdrun -deffnm pull" 
    os.system(command_pull)

    command_pdb = f"gmx editconf -f pull.gro -o {Lipid_name}.pdb"
    os.system(command_pdb)

    # Clean up by deleting all files starting with "pull."
    files_to_delete = glob.glob('pull*')
    for file in files_to_delete:
        os.remove(file)



