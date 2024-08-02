import pandas as pd

file_path = "msd_data_new.pkl"
msd_df = pd.read_pickle(file_path)
# Uncomment these lines to inspect the DataFrame if needed
print(msd_df)
msd_df.to_csv('msd_data_new.csv')

molecules = ['water']
sizes = [512, 1024]#, 2048]

# for mol in molecules:
#     for size in sizes:
#         ds = msd_df.loc[(mol, size), 'Ds']
#         sd = msd_df.loc[(mol, size), 'Stderr']
#         print(f"Diffusion for {mol} {size} is {ds:.4g} +/- {sd:.4g} cm^2/s")
# Check if all rows have 1001 columns
import pandas as pd

# Load the DataFrame
msd_data = pd.read_pickle('msd_data_new.pkl')

# Check if all rows have 1001 columns and no NaN values
correct_length = msd_data.shape[1] == 1001

if correct_length:
    print("All rows have 1001 entries.")
else:
    print("Some rows do not have 1001 entries.")

# Additional check for each row
all_numeric_and_no_nan = True

for idx, row in msd_data.iterrows():
    if len(row) != 1001:
        print(f"Row {idx} does not have 1001 entries. It has {len(row)} entries.")
        all_numeric_and_no_nan = False
    elif row.isna().any():
        print(f"Row {idx} contains NaN values.")
        all_numeric_and_no_nan = False
    elif not pd.api.types.is_numeric_dtype(row):
        print(f"Row {idx} contains non-numeric values.")
        all_numeric_and_no_nan = False

if all_numeric_and_no_nan:
    print("All rows have 1001 numerical entries with no NaN values.")
else:
    print("Some rows do not meet the criteria.")
