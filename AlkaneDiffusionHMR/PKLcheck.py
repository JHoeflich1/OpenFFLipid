import pandas as pd
import numpy as np

file_path = "DS_final.pkl"
df = pd.read_pickle(file_path)
# Uncomment these lines to inspect the DataFrame if needed
print(df.head())
print(df.tail())
print(np.shape(df))

molecules = ['water', 'pentane', 'hexane', 'heptane', 'octane', 'decane', 'pentadecane']
sizes = [512, 1024, 2048]

for mol in molecules:
    for size in sizes:
        ds = df.loc[(mol, size), 'Ds']
        sd = df.loc[(mol, size), 'Stderr']
        print(f"Diffusion for {mol} {size} is {ds:.4g} +/- {sd:.4g} cm^2/s")