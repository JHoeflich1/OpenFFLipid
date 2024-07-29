import numpy as np
import pandas as pd
import os

import multiprocessing

def calculate_msd(alkane,alkane_atoms, alkane_size):
    '''
    Calculates the mean square displacement for a moleucle
    inputs:
        alkane : your molecule
        alkane_atoms : the total number of atoms in your moleucle. example pentane has 17 atoms/molecule
        alkane_size : how many total molecules you have in your simulation box
    '''
    # first shorten .xtc trajectory from 6001 frames to 601
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
    
    nprocs = multiprocessing.cpu_count()
    print(f'Running on {nprocs} CPUs:')

    molecules_process = ['water','pentane','hexane','heptane','octane','decane','pentadecane']
    processes = len(molecules_process) # each molecule will be one process
    pool = multiprocessing.Pool(process=process)

 
    results = pool.map(calculate_msd, args=(molecules_process,mol_key,size))
    
    pool.close()
    #msd_df.to_pickle("msd_data.pkl")


###################################
#before making into a function for paralell processing

# # populate df with msd 
# for key, mol in molecules.items():
#     for num in sizes:
#         for i in range(num):
#             # Create index file per particle
#             with open(f"ndxs_{mol}_{num}_{i}.ndx", "w") as f:
#                 f.write(f"[ p{i} ]\n")
#                 n0 = key * i + 1
#                 f.write(f"{n0}")

#             # Run gmx msd with each index file
#             command = f"echo 0 | gmx msd -f nvt2_{mol}_{num}.xtc -s nvt2_{mol}_{num}.tpr -o msds/msd_{mol}_{num}_{i}.xvg -n ndxs_{mol}_{num}_{i}.ndx -rmpbc -pbc"
#             os.system(command)
#             # -rmpbc means that molecules are made whole for each frame
#             # -pbc means to se periodic boundary conditions for distance calculation
#             # Read the MSD values from the output file
#             with open(f"msds/msd_{mol}_{num}_{i}.xvg") as f:
#                 lines = f.readlines()
            
#             itv = 0
#             for l in lines:
#                 if l[0] != '#' and l[0] != '@':
#                     vals = l.split()
#                     msd_df.loc[(mol, num, i), itv] = float(vals[1])
#                     itv += 1

# # Save the DataFrame to a file
# msd_df.to_pickle("msd_data.pkl")




################ TRYING AGAIN
# import numpy as np
# import pandas as pd
# import os
# import multiprocessing
# from multiprocessing.pool import Pool
# from functools import partial

# def prep_msd(alkane, alkane_atoms, alkane_size):
#     '''
#     Shortens full trajectory so 6001 frames decreases to 601, writes all index files 
#     inputs:
#         alkane : your molecule
#         alkane_atoms : the total number of atoms in your molecule. example pentane has 17 atoms/molecule
#         alkane_size : how many total molecules you have in your simulation box
#         i : index of the molecule in the simulation box
#         tlen : trajectory length
#     '''
#     command = f"gmx trjconv -f nvt2_{alkane}_{alkane_size}.xtc -skip 10 -o nvt2_short_{alkane}_{alkane_size}.xtc"
#     os.system(command)
#     for i in range(alkane_size):
#     with open(f"ndxs_test/ndxs_{alkane}_{alkane_size}_{i}.ndx", "w") as f:
#         f.write(f"[ p{i} ]\n")
#         n0 = alkane_atoms * i + 1
#         f.write(f"{n0}")

# def calculate_msd(alkane, alkane_atoms, alkane_size, i, tlen):
#     '''
#     Calculates the mean square displacement for a molecule
#     inputs:
#         alkane : your molecule
#         alkane_atoms : the total number of atoms in your molecule. example pentane has 17 atoms/molecule
#         alkane_size : how many total molecules you have in your simulation box
#         i : index of the molecule in the simulation box
#         tlen : trajectory length
#     '''

#     # Run gmx msd with each index file
#     command = f"echo 0 | gmx msd -f nvt2_short_{alkane}_{alkane_size}.xtc -s nvt2_{alkane}_{alkane_size}.tpr -o msds_test/msd_{alkane}_{alkane_size}_{i}.xvg -n ndxs_test/ndxs_{alkane}_{alkane_size}_{i}.ndx -rmpbc -pbc"
#     os.system(command)

#     # Read the MSD values from the output file
#     with open(f"msds_test/msd_{alkane}_{alkane_size}_{i}.xvg") as f:
#         lines = f.readlines()

#     itv = 0
#     for l in lines:
#         if l[0] != '#' and l[0] != '@':
#             vals = l.split()
#             msd_values[itv] = float(vals[1])
#             itv += 1
#     print(f'msd values complete for {alkane} {alkane_size}')
#     return (alkane, alkane_size, i, msd_values)


# if __name__ == '__main__':

#     molecules = {3: 'water', 17: 'pentane', 20: 'hexane', 23: 'heptane', 26: 'octane', 32: 'decane', 47: 'pentadecane'}
#     sizes = [512, 1024, 2048]
#     tlen = 6001  # frame length

#     if not os.path.exists('msds_test'):
#         os.makedirs('msds_test')

#     if not os.path.exists('ndxs_test'):
#         os.makedirs('ndxs_test')

#     # store dataframes for each molecule/ size combination
#     dfs = []

#     for mol_key, mol_name in molecules.items():
#         for size in sizes:
#             index = pd.MultiIndex.from_product(
#                 [[mol_name], [size], range(size)],
#                 names=['molecule', 'sizes', 'particle']
#             )
#             columns = range(tlen)
#             df = pd.DataFrame(np.zeros((len(index), tlen)), index=index, columns=columns)
#             dfs.append(df)

#     # concat to a single dataframe
#     msd_df = pd.concat(dfs)

#     # Create a partial function to pass the common arguments
#     calculate_msd_partial = partial(calculate_msd, tlen=tlen)

#     # Prepare the arguments for the multiprocessing pool
#     tasks = []
#     for mol_key, mol_name in molecules.items():
#         for size in sizes:
#             for i in range(size):
#                 tasks.append((mol_name, mol_key, size, i))

#     # Run the tasks in parallel
#     with Pool(multiprocessing.cpu_count()) as pool:
#         results = pool.starmap(calculate_msd_partial, tasks)

#     # Populate the msd_df with the results
#     for alkane, alkane_size, i, msd_values in results:
#         msd_df.loc[(alkane, alkane_size, i)] = msd_values

#     msd_df.to_pickle("msd_data.pkl")
