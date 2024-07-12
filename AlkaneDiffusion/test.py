import numpy as np

molecules = {6:'hexane', 7:'heptane'}
num_molecules = [512, 1024]
tlen = 2001 #trajectory length
n = 3 #number of atoms in your molecule

msd_dict = {}

for key, value in molecules.items():
    msd_dict[value] = {}
    for num in num_molecules:
        msd_dict[value][num] = np.zeros([num, tlen])

for molecule, sizes in msd_dict.items():
    for size, msd_array in sizes.items():
        print(f"Molecule: {molecule}, Size: {size}, Shape: {msd_array.shape}")