import pandas as pd
import matplotlib.pyplot as plt
import os

# Define the folders containing the data
folders = {
    "AlkaneDiffusion_md": "Non-HMR",
    "AlkaneDiffusionHMR_md": "HMR"
}

# Initialize an empty DataFrame
combined_df = pd.DataFrame()

# Loop through each folder and read the CSV files
for folder, label in folders.items():
    file_path = os.path.join(folder, "diffusion_results.csv")
    if os.path.exists(file_path):
        df = pd.read_csv(file_path)
        df['label'] = label  # Add a column to indicate HMR or Non-HMR
        combined_df = pd.concat([combined_df, df], ignore_index=True)
    else:
        print(f"Warning: File {file_path} does not exist.")

# Add a column for the alkane chain length
chain_length_map = {
    "pentane": 5,
    "hexane": 6,
    "heptane": 7,
    "octane": 8,
    "decane": 10,
    "pentadecane": 15,
    "hexadecane": 16,  # Corrected to match experimental data
    "tetradecane": 14  # Added missing chain length
}

combined_df['chain_length'] = combined_df['molecule'].map(chain_length_map)

# Add experimental data and convert to match the scale
experimental_data = {
    "pentane": 5.45,
    "hexane": 4.18,
    "heptane": 3.12,
    "octane": 2.00,
    "decane": 1.31,
    "tetradecane": 0.56,  # Added experimental data for tetradecane
    "hexadecane": 0.39,
    "pentadecane": 0.56  # Added experimental data for pentadecane
}

# Create a DataFrame for experimental data
exp_df = pd.DataFrame({
    'chain_length': [chain_length_map[mol] for mol in experimental_data.keys()],
    'experimental': [val * 1e-9 for val in experimental_data.values()]  # Adjust if necessary
})

# Merge the simulation data with experimental data
combined_df = pd.merge(combined_df, exp_df, on='chain_length')

# Calculate the percentage difference (signed)
combined_df['percent_difference'] = ((combined_df['average diffusion curvefit'] - combined_df['experimental']) / combined_df['experimental']) * 100

# Filter out heptadecane and hexadecane
molecules_to_exclude = ["hexadecane", "heptadecane"]
filtered_df = combined_df[~combined_df['molecule'].isin(molecules_to_exclude)]

# Plotting the percentage difference
plt.figure(figsize=(10, 6))

for label in filtered_df['label'].unique():
    subset = filtered_df[filtered_df['label'] == label]
    plt.plot(
        subset['chain_length'],
        subset['percent_difference'],
        label=f'{label} Percent Difference',
        marker='o',
        linestyle='-'
    )

plt.xlabel('Alkane Chain Length')
plt.ylabel('Percent Difference ((Sim - Exp) / Exp) * 100')
plt.title('Chain-Length Dependent Deviations in Simulated vs. Experimental Diffusion')
plt.legend(title='Dataset')
plt.savefig('percent_difference_plot.png')
plt.show()
