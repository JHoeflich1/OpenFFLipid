import numpy as np
import pandas as pd
import os

import multiprocessing
from multiprocessing.pool import Pool
from functools import partial
import time

def calculate_msd(alkane,alkane_atoms, alkane_size):
    '''
    Calculates the mean square displacement for a moleucle
    inputs:
        alkane : your molecule
        alkane_atoms : the total number of atoms in your moleucle. example pentane has 17 atoms/molecule
        alkane_size : how many total molecules you have in your simulation box
    '''
    # first shorten .xtc trajectory from 7001 frames to 701 with 10 ps timestep 
    command = f"gmx trjconv -f nvt2_{alkane}_{alkane_size}.xtc -skip 10 -o nvt2_short_{alkane}_{alkane_size}.xtc"
    os.system(command)

    for i in range(alkane_size):
        # Skip calculations that have already been done
        msd_file = f"msds/msd_{alkane}_{alkane_size}_{i}.xvg"
        if os.path.exists(msd_file):
            print(f"{msd_file} already exists, skipping this calculation.")
            continue
        with open(f"ndxs/ndxs_{alkane}_{alkane_size}_{i}.ndx", "w") as f:
            f.write(f"[ p{i} ]\n")
            n0 = alkane_atoms * i + 1
            f.write(f"{n0}")
            #print(f'ndxs_test/ndxs_{alkane}_{alkane_size}_{i}.ndx is written')

        # Run gmx msd with each index file
        #print('moving on to gmx')
        command = f"echo 0 | gmx msd -f nvt2_short_{alkane}_{alkane_size}.xtc -s nvt2_{alkane}_{alkane_size}.tpr -o msds/msd_{alkane}_{alkane_size}_{i}.xvg -n ndxs/ndxs_{alkane}_{alkane_size}_{i}.ndx -rmpbc -pbc"
        os.system(command)
        # -rmpbc means that molecules are made whole for each frame
        # -pbc means to se periodic boundary conditions for distance calculation
        # Read the MSD values from the output file
        #with open(f"msds/msd_{alkane}_{alkane_size}_{i}.xvg") as f:
        #    lines = f.readlines()
        
        #itv = 0
        #for l in lines:
        #    if l[0] != '#' and l[0] != '@':
        #        vals = l.split()
        #        msd_df.loc[(alkane, alkane_size, i), itv] = float(vals[1])
        #        itv += 1


if __name__ == '__main__':

    molecules = {26:'octane',32:'decane',47:'pentadecane'}#{3: 'water', 17:'pentane',20:'hexane',23:'heptane',26 :'octane',32:'decane',47:'pentadecane'}
    sizes = [512,1024, 2048]
    tlen = 701  # trajectory length

    if not os.path.exists('msds'):
        os.makedirs('msds')

    if not os.path.exists('ndxs'):
        os.makedirs('ndxs')

    # store dataframes for each moelecule/ size combination
    dfs = []

    for mol_key, mol_name in molecules.items():
        for size in sizes:
            index = pd.MultiIndex.from_product(
                [[mol_name], [size], range(size)],
                names=['molecule', 'sizes', 'particle']
            )
            columns = range(tlen)
            df = pd.DataFrame(np.zeros((len(index), tlen)), index=index, columns=columns)
            dfs.append(df)

    # print(dfs)
    # concat to a single dataframe
    msd_df = pd.concat(dfs)
    
    # nprocs = multiprocessing.cpu_count()
    # print(f'Running on {nprocs} CPUs:')

    for mol_key, mol_name in molecules.items():
        for size in sizes:
            calculate_msd(mol_name,mol_key,size)
            
    msd_df.to_pickle("msd_data.pkl")




# Wrote this script for parallel processing, Worked but the msd data was all wrong..
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import itertools
# import multiprocessing 
# import os

# def calculate_msd(args):
#     '''
#     Calculates the mean square displacement for a moleucle
#     inputs: args, a single touple containing three things
#         alkane : your molecule (string)
#         alkane_atoms : the total number of atoms in your moleucle. example pentane has 17 atoms/molecule (int)
#         alkane_size : how many total molecules you have in your simulation box (int)
#     '''
#     alkane,alkane_atoms, alkane_size = args
#     # first shorten .xtc trajectory from 6001 frames to 601
#     # command = f"gmx trjconv -f nvt2_{alkane}_{alkane_size}.xtc -skip 10 -o nvt2_short_{alkane}_{alkane_size}.xtc"
#     # os.system(command)

#     for i in range(alkane_size):
#         # Skip calculations that have already been done
#         msd_file = f"msds_test/msd_{alkane}_{alkane_size}_{i}.xvg"
#         if os.path.exists(msd_file):
#             print(f"{msd_file} already exists, skipping this calculation.")
#             continue
#         with open(f"ndxs_test/ndxs_{alkane}_{alkane_size}_{i}.ndx", "w") as f:
#             f.write(f"[ p{i} ]\n")
#             n0 = alkane_atoms * i + 1
#             f.write(f"{n0}")
#             #print(f'ndxs_test/ndxs_{alkane}_{alkane_size}_{i}.ndx is written')

#         # Run gmx msd with each index file
#         #print('moving on to gmx')
#         command = f"echo 0 | gmx msd -f nvt2_short_{alkane}_{alkane_size}.xtc -s nvt2_{alkane}_{alkane_size}.tpr -o msds_test/msd_{alkane}_{alkane_size}_{i}.xvg -n ndxs_test/ndxs_{alkane}_{alkane_size}_{i}.ndx -rmpbc -pbc"
#         os.system(command)
#         # -rmpbc means that molecules are made whole for each frame
#         # -pbc means to se periodic boundary conditions for distance calculation
#         # Read the MSD values from the output file
#         #with open(f"msds/msd_{alkane}_{alkane_size}_{i}.xvg") as f:
#         #    lines = f.readlines()
        
#         #itv = 0
#         #for l in lines:
#         #    if l[0] != '#' and l[0] != '@':
#         #        vals = l.split()
#         #        msd_df.loc[(alkane, alkane_size, i), itv] = float(vals[1])
#         #        itv += 1

# if __name__ == '__main__':
#     molecules = ['water/3','pentane/17','hexane/20','heptane/23','octane/26','decane/32','pentadecane/47']
#     sizes = [512, 1024, 2048]

#     # Create the list of tuples
#     combinations = list(itertools.product(molecules, sizes))
#     new_list = []

#     for entry in combinations:
#         molecule, number = entry[0].split('/')
#         number = int(number)  # Convert number to integer value
#         new_list.append([molecule, number, entry[1]])
#     # print(new_list)  #this creates a touple for 

#     if not os.path.exists('msds_test'):
#         os.makedirs('msds_test')

#     if not os.path.exists('ndxs_test'):
#         os.makedirs('ndxs_test')

#     pool = multiprocessing.Pool()
#     pool.map(calculate_msd,new_list)
#     pool.close()
#     pool.join()
