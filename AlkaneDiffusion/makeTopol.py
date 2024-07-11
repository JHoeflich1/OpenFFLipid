import os

# Molecules and their respective sizes
molecules = ['hexane', 'heptane', 'octane', 'decane', 'pentadecane']
sizes = [512, 1024, 2048, 4096]

# Function to create new .top files
def create_new_top_files(molecule, original_top_file, sizes):
    # Read original .top file content
    with open(original_top_file, 'r') as f:
        original_content = f.readlines()

    # Create new .top files for each size
    for size in sizes:
        new_top_file = f"{molecule}_{size}.top"
        with open(new_top_file, 'w') as f:
            inside_molecules_section = False
            for line in original_content:
                if line.startswith('[ molecules ]'):
                    inside_molecules_section = True
                    f.write(line)
                    f.write(';name\tnumber\n')
                    f.write(f"ALK\t{size}\n")  # Add the modified line once
                elif inside_molecules_section and line.startswith('ALK'):
                    continue  # Skip the original ALK line
                else:
                    f.write(line)
                if inside_molecules_section and line.strip() == '':
                    inside_molecules_section = False  # End of [ molecules ] section

# Iterate over each molecule and process
for molecule in molecules:
    original_top_file = f"{molecule}.top"
    if os.path.exists(original_top_file):
        create_new_top_files(molecule, original_top_file, sizes)
        print(f"Processed {molecule}")
    else:
        print(f"Error: {original_top_file} does not exist.")


size_gmx = {13.0:512, 18.0: 1024,25.0: 2048, 35.0:4096} #note that packmol is in angstroms and gromacs in nm


# Convert pdb to gro for each molecule and size
for molecule in molecules:
    for box, num in size_gmx.items():
        pdb_file = f"{molecule}_{num}.pdb"
        gro_file = f"{molecule}_{num}.gro"
        if os.path.exists(pdb_file):
            os.system(f"gmx editconf -f {pdb_file} -o {gro_file} -bt cubic -box {box} {box} {box}")
            print(f"Converted {pdb_file} to {gro_file}")
        else:
            print(f"Error: {pdb_file} does not exist.")