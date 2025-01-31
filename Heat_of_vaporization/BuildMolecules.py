from openff.toolkit import Molecule, Topology, ForceField
from openff.interchange import Interchange
from openff.units import unit, Quantity 
import numpy as np

import openff.nagl
from openff.nagl import GNNModel 
from openff.nagl_models import list_available_nagl_models
from openff.interchange.components._packmol import pack_box


# Make molecules
pentane = Molecule.from_smiles("C" * 5, allow_undefined_stereo=True)
hexane = Molecule.from_smiles("C" * 6, allow_undefined_stereo=True)
heptane = Molecule.from_smiles("C" * 7, allow_undefined_stereo=True)
octane = Molecule.from_smiles("C" * 8, allow_undefined_stereo=True)
decane = Molecule.from_smiles("C" * 10, allow_undefined_stereo=True)
pentadecane = Molecule.from_smiles("C" * 15, allow_undefined_stereo=True)

alkanes = {
    'pentane': pentane, 'hexane': hexane, 'heptane': heptane, 
    'octane': octane, 'decane': decane, 'pentadecane': pentadecane
}

# Experimental densities in g/nm^3 for approximate packing
densities = {
    'pentane': 6.26e-22 / 72.15, 
    'hexane': 6.59e-22 / 86.17, 
    'heptane': 6.8e-22 / 100.21, 
    'octane': 7.03e-22 / 114.23, 
    'decane': 7.3e-22 / 152.29, 
    'pentadecane': 7.69e-22 / 212.42, 
}

size = 1000  # Number of molecules

# Dictionary to store single box size values for liquid simulations
box_sizes = {}

# Calculate box sizes for each molecule
for molecule, density in densities.items():
    box_size = ((size / (density * 6.02e23)) * 3) ** (1/3)  # Add packing room
    box_sizes[molecule] = round(box_size, 2)  # Store as single float

# for molecule, box_size in box_sizes.items():
#     print(f"Box size for {molecule}: {box_size} nm")

for name, alkane in alkanes.items():
    alkane.generate_conformers(n_conformers=1)
    # print(f"{name} conformers: {len(alkane.conformers)}")

    if not alkane.conformers:
        print(f"Error: No conformers generated for {name}")
        continue

    alkane.name = "ALK"
    for i, atom in enumerate(alkane.atoms, 3):
        atom.metadata["residue_name"] = 'ALK'
    alkane.generate_unique_atom_names()

    # Get the gas phase. Lets use a 10 nm box 
    gas_box = unit.Quantity(100 * np.eye(3), unit.angstrom)
    gas_top = Topology.from_molecules([alkane])
    ff = ForceField("openff-1.3.1.offxml")
    gas_interchange = Interchange.from_smirnoff(
        force_field=ff,
        topology=gas_top,
        box=gas_box,
    )
    gas_interchange.to_gromacs(f"{name}_gas")

    # Get the liquid gromacs inputs
    cubic_box_size = box_sizes[name]  # Directly access float
    
    cubic_box = unit.Quantity(10 * cubic_box_size * np.eye(3), unit.angstrom)  # Convert nm to Å


    try:
        packed_topol = pack_box(
            molecules=[alkane], number_of_copies=[size], solute=None,
            tolerance=2 * unit.angstrom, box_vectors=cubic_box
        )
        packed_interchange = Interchange.from_smirnoff(
            force_field=ff,
            topology=packed_topol,
            box=cubic_box
        )
        packed_interchange.to_gromacs(f"{name}_liquid")
    except Exception as e:
        print(f"Error packing liquid box for {name}: {e}")
