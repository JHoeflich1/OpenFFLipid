import pandas as pd

file_path = "msd_data.pkl"
msd_df = pd.read_pickle(file_path)
# Uncomment these lines to inspect the DataFrame if needed
print(msd_df)
msd_df.to_csv('msd_data.csv')

# molecules = ['water', 'pentane', 'hexane', 'heptane', 'octane', 'decane', 'pentadecane']
# sizes = [512, 1024, 2048]

# for mol in molecules:
#     for size in sizes:
#         ds = msd_df.loc[(mol, size), 'Ds']
#         sd = msd_df.loc[(mol, size), 'Stderr']
#         print(f"Diffusion for {mol} {size} is {ds:.4g} +/- {sd:.4g} cm^2/s")