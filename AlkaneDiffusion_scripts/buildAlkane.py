from openff.toolkit import Molecule, Topology, ForceField
from openff.interchange import Interchange
from openff.units import unit, Quantity 
import numpy as np

import openff.nagl
from openff.nagl import GNNModel 
from openff.nagl_models import list_available_nagl_models
from openff.interchange.components._packmol import pack_box


#make moleucles to run 
pentane = Molecule.from_smiles("C" * 5, allow_undefined_stereo=True)
hexane = Molecule.from_smiles("C" * 6, allow_undefined_stereo=True)
heptane = Molecule.from_smiles("C" * 7, allow_undefined_stereo=True)
octane = Molecule.from_smiles("C" * 8, allow_undefined_stereo=True)
decane = Molecule.from_smiles("C" * 10, allow_undefined_stereo=True)
pentadecane = Molecule.from_smiles("C" * 15, allow_undefined_stereo=True)
# water = Molecule.from_smiles("O", allow_undefined_stereo=True)

alkanes = {'pentane': pentane,'hexane':hexane, 'heptane':heptane, 'octane':octane,'decane':decane, 'pentadecane':pentadecane} #'water': water,


# Experimental densities in g/nm^3 grabbed from online. Multiply by molar mass to get mol/nm^3
densities = {
    'pentane': 6.26e-22 / 72.15, # g/nm^3 * 1/(g/mol) 
    'hexane': 6.59e-22 / 86.17, 
    'heptane': 6.8e-22 / 100.21, 
    'octane': 7.03e-22 / 114.23, 
    'decane': 7.3e-22 / 152.29, 
    'pentadecane': 7.69e-22 / 212.42, 
    # 'water': 1e-21 / 18.01
}

sizes = [512, 1024, 2048]  # number of molecules

# dictionay to store box sizes
box_sizes = {key: [] for key in densities.keys()}

# Calculate box sizes for each molecule and size
for molecule, density in densities.items():
    for size in sizes:
        box_size = ((size / (density* 6.02e23))*3)**(1/3) #multiply density by avogadros number to get molecules/nm^3, and multipy volume by 3 to give wiggle room with packing
        box_sizes[molecule].append(round(box_size,2))

for molecule, size_list in box_sizes.items():
    print(f"Box sizes for {molecule}: {size_list}")


for name, alkane in alkanes.items():

    alkane.generate_conformers(n_conformers=1)
    print(f"{name} conformers: {len(alkane.conformers)}")

    if not alkane.conformers:
        print(f"Error: No conformers generated for {name}")
        continue


    
    alkane.name = "ALK"
    for i, atom in enumerate(alkane.atoms, 3):
        atom.metadata["residue_name"] = 'ALK'
    alkane.generate_unique_atom_names()

    for size in sizes:
            cubic_box_size = box_sizes[name][sizes.index(size)]
            print(cubic_box_size,f'printing out box size for {name}, {size}')
            cubic_box = unit.Quantity(10 *cubic_box_size* np.eye(3), unit.angstrom) #multiply by 10 to convert nm to angstroms
            try:
                packed_topol = pack_box(molecules =[alkane], number_of_copies=[size], solute=None, tolerance= 2*unit.angstrom, box_vectors = cubic_box)
                packed_interchange = Interchange.from_smirnoff(
                     force_field = ForceField("openff-2.1.0.offxml"),
                     topology= packed_topol,
                     box = cubic_box)
                packed_interchange.to_gromacs(f"{name}_{size}")
            except Exception as e:
                print(f"Error packing box for {name} with size {size}: {e}")