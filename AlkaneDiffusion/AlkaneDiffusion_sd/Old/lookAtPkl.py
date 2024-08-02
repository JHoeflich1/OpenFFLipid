import pandas as pd

file_path = "msd_data.pkl"
msd_df = pd.read_pickle(file_path)
print(msd_df.head())

