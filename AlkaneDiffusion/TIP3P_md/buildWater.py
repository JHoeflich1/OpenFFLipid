#build TIP3P simulations 
from openff.toolkit import Molecule, Topology, ForceField
from openff.interchange import Interchange
from openff.units import unit, Quantity 
import numpy as np
from openff.interchange.components._packmol import pack_box



water = Molecule.from_smiles("O", allow_undefined_stereo=True)
molecules = {'water': water}

# Experimental densities in g/nm^3 grabbed from online. Multiply by molar mass to get mol/nm^3
densities = {
    'water': 1e-21 / 18.01
}

sizes = [512, 1024, 2048, 4096]  # number of molecules

# dictionay to store box sizes
box_sizes = {key: [] for key in densities.keys()}

# Calculate box sizes for each molecule and size
for molecule, density in densities.items():
    for size in sizes:
        box_size = ((size / (density* 6.02e23))*3.5)**(1/3) #multiply density by avogadros number to get molecules/nm^3, and multipy volume by 3 to give wiggle room with packing
        box_sizes[molecule].append(round(box_size,2))

for molecule, size_list in box_sizes.items():
    print(f"Box sizes for {molecule}: {size_list}")


for name, moleucle in molecules.items():

    moleucle.generate_conformers(n_conformers=1)
    print(f"{name} conformers: {len(moleucle.conformers)}")

    if not moleucle.conformers:
        print(f"Error: No conformers generated for {name}")
        continue
 
    
    moleucle.name = "TIP"
    for i, atom in enumerate(moleucle.atoms, 3):
        atom.metadata["residue_name"] = 'TIP'
    moleucle.generate_unique_atom_names()

    for size in sizes:
            cubic_box_size = box_sizes[name][sizes.index(size)]
            print(cubic_box_size,f'printing out box size for {name}, {size}')
            cubic_box = unit.Quantity(10 *cubic_box_size* np.eye(3), unit.angstrom) #multiply by 10 to convert nm to angstroms
            try:
                # Consider increasing the tolerance to 2 for waters because my sims are crashing
                packed_topol = pack_box(molecules =[moleucle], number_of_copies=[size], solute=None, tolerance= 2*unit.angstrom, box_vectors = cubic_box)
                packed_interchange = Interchange.from_smirnoff(
                     force_field = ForceField("openff-2.1.0.offxml"),
                     topology= packed_topol,
                     box = cubic_box)
                packed_interchange.to_gromacs(f"{name}_{size}")
            except Exception as e:
                print(f"Error packing box for {name} with size {size}: {e}")