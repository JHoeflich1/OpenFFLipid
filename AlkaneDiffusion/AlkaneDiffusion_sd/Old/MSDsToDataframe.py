import os
import numpy as np
import pandas as pd

def populateDataframe(molecules, sizes, dataframe):
    '''
    Populates MSD dataframe with data from MSD files.
    inputs:
        molecules : list of molecule names
        sizes : list of molecule sizes (number of molecules in the simulation box)
        dataframe : DataFrame to be populated with MSD data
    '''
    for mol_name in molecules:
        for size in sizes:
            for i in range(size):
                file_path = f"msds_PPtest/msd_{mol_name}_{size}_{i}.xvg"
                if not os.path.isfile(file_path):
                    print(f"File not found: {file_path}")
                    continue

                with open(file_path) as f:
                    lines = f.readlines()

                itv = 0
                for l in lines:
                    if l[0] != '#' and l[0] != '@':
                        try:
                            vals = l.split()
                            dataframe.loc[(mol_name, size, i), itv] = float(vals[1])
                            itv += 1
                        except (IndexError, ValueError) as e:
                            print(f"Error parsing line in {file_path}: {l}")
                            print(f"Error: {e}")
            print(f'{mol_name}, {size} is finished populating df')
            print(dataframe.head())

if __name__ == '__main__':
    molecules = ['water', 'pentane', 'hexane', 'heptane', 'octane', 'decane', 'pentadecane']
    sizes = [512, 1024, 2048]
    tlen = 1001  

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
    populateDataframe(molecules, sizes, msd_df)
    msd_df.to_pickle("msd_data_new.pkl")

    # Check the last few entries in the DataFrame
    print(msd_df.tail())
