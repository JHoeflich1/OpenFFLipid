import numpy as np
import pandas as pd
import os


def populateDataframe(molecules, sizes, dataframe):
    '''
    Poopulates MSD dataframe with data from MSD
    inputs:
        alkane : list of molecules
        alkane_size : how many total molecules you have in your simulation box
    '''

    # Read the MSD values from the output file
    for mol_name in molecules:
        for size in sizes:
            for i in range(size):
                with open(f"msds/msd_{mol_name}_{size}_{i}.xvg") as f:
                    lines = f.readlines()
                #read lines and save msds to dataframe    
                itv = 0
                for l in lines:
                    if l[0] != '#' and l[0] != '@':
                        vals = l.split()
                        dataframe.loc[(mol_name, size, i), itv] = float(vals[1])
                        itv += 1
            print(f'{mol_name}, {size} is finished populating df')



if __name__ == '__main__':
    
    molecules = ['water','pentane','hexane','heptane','octane','decane','pentadecane']
    sizes = [512,1024, 2048]
    tlen = 701     #in my msds I have 1001 frames at a 5 ps timestep

    # store dataframes for each moelecule/ size combination
    dfs = []

    for mol_name in molecules:
        for size in sizes:
            index = pd.MultiIndex.from_product(
                [[mol_name], [size], range(size)],
                names=['molecule', 'sizes', 'particle']
            )
            columns = range(tlen)
            df = pd.DataFrame(np.full((len(index), tlen), np.nan), index=index, columns=columns)
            dfs.append(df)


    msd_df = pd.concat(dfs)
 
    populateDataframe(molecules,sizes, msd_df)
    msd_df.to_pickle("msd_data.pkl")



#### tried parallel processing, worked byt the msd data was all out of order

# import numpy as np
# import pandas as pd
# import os
# import multiprocessing
# import itertools


# def populateDataframe(args):
#     molecule, size, df = args
#     '''
#     Populates MSD dataframe with data from MSD
#     inputs:
#         molecule : name of the molecule
#         size : how many total molecules you have in your simulation box
#         df : DataFrame to be populated
#     '''

#     for i in range(size):
#         with open(f"msds/msd_{mol_name}_{size}_{i}.xvg") as f:
#             lines = f.readlines()
#         #read lines and save msds to dataframe    
#         itv = 0
#         for l in lines:
#             if l[0] != '#' and l[0] != '@':
#                 vals = l.split()
#                 df.loc[(molecule, size, i), itv] = float(vals[1])
#                 itv += 1
#     print(f'{molecule}, {size} is finished populating df')
#     return df

# if __name__ == '__main__':
    
#     molecules = ['water','pentane','hexane','heptane','octane','decane','pentadecane']
#     sizes = [512,1024, 2048]
#     tlen = 701

#     # store dataframes for each moelecule/ size combination
#     dfs = []

#     for mol_name in molecules:
#         for size in sizes:
#             index = pd.MultiIndex.from_product(
#                 [[mol_name], [size], range(size)],
#                 names=['molecule', 'sizes', 'particle']
#             )
#             columns = range(tlen)
#             df = pd.DataFrame(np.full((len(index), tlen), np.nan), index=index, columns=columns)
#             dfs.append(df)

#     combinations = list(itertools.product(molecules, sizes))
#     new_list = [(mol_name, size, dfs[i]) for i, (mol_name, size) in enumerate(combinations)]
    
#     # Use multiprocessing to populate DataFrames
#     with multiprocessing.Pool() as pool:
#         dfs_filled = pool.map(populateDataframe, new_list)

#     if not os.path.exists('msd_data'):
#         os.makedirs('msd_data')

#     for (mol_name, size), df in zip(combinations, dfs_filled):
#         df.to_csv(f'msd_data/msd_data_{mol_name}_{size}.csv')
