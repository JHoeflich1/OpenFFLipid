import pandas as pd

file_path = "msd_data.pkl"
msd_df = pd.read_pickle(file_path)
print(msd_df.head())
print(tail -n 40 msd_df)
