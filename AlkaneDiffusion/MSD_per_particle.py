import numpy as np
import pandas as pd
import os

molecules = {3: 'TIP3P', 6:'hexane',7:'heptane',8:'octane',10:'decane',15:'pentadecane'}
sizes = [512,1024, 2048, 4096]
tlen = 5  # trajectory length

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

# concat to a single dataframe
msd_df = pd.concat(dfs)

# print(msd_df)

# populate df with msd 
for key, mol in molecules.items():
    for num in sizes:
        for i in range(num):
            # Create index file per particle
            with open(f"ndxs_{mol}_{num}_{i}.ndx", "w") as f:
                f.write(f"[ p{i} ]\n")
                n0 = key * i
                f.write(f"{n0}")

            # Run gmx msd with each index file
            command = f"gmx msd -f nvt_{mol}_{num}.xtc -s nvt_{mol}_{num}.tpr -o msds/msd_{mol}_{num}_{i}.xvg -n ndxs_{mol}_{num}_{i}.ndx"
            os.system(command)

            # Read the MSD values from the output file
            with open(f"msds/msd_{mol}_{num}_{i}.xvg") as f:
                lines = f.readlines()
            
            itv = 0
            for l in lines:
                if l[0] != '#' and l[0] != '@':
                    vals = l.split()
                    msd_df.loc[(mol, num, i), itv] = float(vals[1])
                    itv += 1

# Save the DataFrame to a file
msd_df.to_pickle("msd_data.pkl")