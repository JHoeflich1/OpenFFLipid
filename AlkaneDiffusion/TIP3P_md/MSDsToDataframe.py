import numpy as np
import pandas as pd
import os
from multiprocessing import Pool


def read_and_populate(mol_name, size, tlen):
    '''
    Reads MSD values for a given molecule and size and populates a DataFrame.
    '''
    index = pd.MultiIndex.from_product(
        [[mol_name], [size], range(size)],
        names=['molecule', 'sizes', 'particle']
    )
    columns = range(tlen)
    df = pd.DataFrame(np.full((len(index), tlen), np.nan), index=index, columns=columns)
    
    for i in range(size):
        with open(f"msds/msd_{mol_name}_{size}_{i}.xvg") as f:
            lines = f.readlines()
        # Read lines and save MSDs to dataframe    
        itv = 0
        for l in lines:
            if l[0] != '#' and l[0] != '@':
                vals = l.split()
                df.loc[(mol_name, size, i), itv] = float(vals[1])
                itv += 1
    print(f'{mol_name}, {size} is finished populating df')
    return df


if __name__ == '__main__':
    molecules = ['water']#, 'pentane', 'hexane', 'heptane', 'octane', 'decane', 'pentadecane']
    sizes = [512, 1024, 2048]
    tlen = 5001 # In my MSDs I have 5001 frames at a 1 ps timestep

    pool = Pool()
    tasks = [(mol_name, size, tlen) for mol_name in molecules for size in sizes]
    print(tasks)
    results = pool.starmap(read_and_populate, tasks)
    pool.close()
    pool.join()

    msd_df = pd.concat(results)
    msd_df.to_pickle("msd_data_new.pkl")
