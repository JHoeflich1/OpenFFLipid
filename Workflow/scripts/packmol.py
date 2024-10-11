from openff.toolkit import Molecule, Quantity, RDKitToolkitWrapper, Topology, unit #type:ignore
from openff.units import Quantity, unit#type:ignore
import pandas as pd #type:ignore
import numpy as np
import mdtraj#type:ignore
import os
import shutil 
import subprocess
from typing import List
import json
from datetime import datetime
from openff.interchange import Interchange#type:ignore
from openff.toolkit import ForceField, Molecule, Topology#type:ignore

def init_json(file_name):
    ''' Initiates a json structure and writes it to a file. This will be carried throughout the 
    lipid build to document important inforaiton 
    
    Parameters:
        file_name: name of the json file
        '''
    
    current_date = datetime.now().isoformat()
    data = {
        "experiment": "Bilayer Build",
        "date": current_date,
        "parameters": {}
    }
    # Write the initial data to the JSON file
    with open(file_name, 'w') as json_file:
        json.dump(data, json_file, indent=4)
    
    print(f"JSON initialized and saved to {file_name}")
    return data

def append_and_save_json(file_name, data, new_params):
    """
    Appends new parameters to the JSON data structure and saves it to a file.
    
    Parameters:
    - file_name (str): The JSON file to save to.
    - data (dict): The current JSON data structure.
    - new_params (dict): New parameters to add to the 'parameters' section.
    """
    # Append new parameters
    data["parameters"].update(new_params)
    
    # Save the updated data back to the JSON file
    with open(file_name, 'w') as json_file:
        json.dump(data, json_file, indent=4)
    
    print(f"New parameters appended and JSON data successfully saved to {file_name}")


def load_top_pdb_files(molecule_names: List[str], source_dir: str, dest_dir: str) -> List[str]:
    """
    Copy the PDB and top file of a lipid from the source directory to the destination directory.

    Parameters
    ----------
    lipid_name : str
        The name of the lipid for which the PDB file needs to be copied (without extension).
    source_dir : str
        The path to the directory containing the lipid PDB files.
    dest_dir : str
        The path to the destination directory where the PDB file should be copied.

    Returns
    -------
    list of str
        A list of paths to the copied PDB  and top files in the destination directory.
    """
    copied_files = []

    for mol in molecule_names:
        # Define source and destination file paths for each lipid
        source_pdb_path = os.path.join(source_dir, f'{mol}.pdb')
        dest_pdb_path = os.path.join(dest_dir, f'{mol}.pdb')
        source_top_path = os.path.join(source_dir, f'{mol}.top')
        dest_top_path = os.path.join(dest_dir, f'{mol}.top')

        # Check if the source PDB file exists
        if not os.path.isfile(source_pdb_path):
            raise FileNotFoundError(f"The file '{source_pdb_path}' does not exist.")

        # Check if the source .top file exists
        if not os.path.isfile(source_top_path):
            raise FileNotFoundError(f"The file '{source_top_path}' does not exist.")

        # Copy the PDB file
        shutil.copy(source_pdb_path, dest_pdb_path)
        copied_files.append(dest_pdb_path)

        # Copy the top file
        shutil.copy(source_top_path, dest_top_path)
        copied_files.append(dest_top_path)

    return copied_files

def _build_input_file(
    lipid_file_names: List[str], 
    lipid_counts: List[int],
    solvent: str,
    solvent_count: int,
    tolerance: float,
    output: str,
    ) -> tuple[str, str, List[float]]:
    """
    Calculate box size.
    Construct the packmol input file.
    Initiate config.json file

    Parameters
    ----------
    lipid_file_names: list of str
        The paths to the molecule pdb files.
    molecule_counts: list of int
        The number of each lipid to add.
    solvent: str
        The path to the structure to solvate.
    solvent_count: int
        The number of solvent molecules.
    tolerance: float
        The packmol convergence tolerance.
    output: string
        The packmol input file name. 


    Returns
    -------
    str
        The path to the input file.
    str
        The path to the output file.
    list of float
        The calculated box dimensions [x, y, z].

    """
    if len(lipid_file_names) != len(lipid_counts):
        raise ValueError("The number of lipid names must match the number of lipid counts.")
    
    #save all data to a json file
    file_name = 'config.json'
    data = init_json(file_name)

    Library = 'PulledLipid.csv'
    lipid_library = pd.read_csv(f'Dictionary/{Library}') #load in the lipid library

    #move lipids pdb from dictionary to working dir
    cwd = os.getcwd()
    dest_dir = cwd

    # Ensure PDB files are copied for each lipid
    for lipid_name in lipid_file_names:
        source_dir = os.path.join(cwd,'Dictionary', 'lipids_parameterized', lipid_name)
        load_top_pdb_files([lipid_name], source_dir, dest_dir)  # Call function with list

    # Ensure PDB files are copied for solvent
    if solvent is not None:
        source_dir = os.path.join(cwd,'Dictionary', 'lipids_parameterized', solvent)
        load_top_pdb_files([solvent], source_dir, dest_dir) #solvent must be wrapped in  a list

    #initiate a variable that stores the largest hg_tg_distance to then pack the solvent 
    max_z_distance = 0

    # Add the global header options in packmol
    output_file = "packmol_output.pdb"
    input_lines = [
        f"tolerance {tolerance:f}",
        "filetype pdb",
        f"output {output_file}",
        "",]
    

  


    ############################
    # I am assuming that lipid area in x-y space is 5 angstoms by 5 agstoms to give
    # room with packing.
    n = 6

    for lipid_name, lipid_count in zip(lipid_file_names, lipid_counts):
        # Find the corresponding row in the lipid library
        lipid_info = lipid_library[lipid_library['Name'] == lipid_name]
        if lipid_info.empty:
            raise ValueError(f"Lipid '{lipid_name}' not found in the library at {Library}")

        # Load relevant lipid data
        hg_tg_distance = float(lipid_info['HG/TG distance'].values[0])
        tg_index = int(lipid_info['Tailgroup Atom Index'].values[0])
        hg_index = int(lipid_info['Headgroup Atom Index'].values[0])

        # Update max_z_distance
        if hg_tg_distance > max_z_distance:
            max_z_distance = hg_tg_distance

        # Top leaflet
        input_lines.extend([
            f"structure {lipid_name}.pdb",
            f"  number {int(lipid_count / 2)}",
            f"  inside box 0. 0. 0. {np.sqrt(sum(lipid_counts)) * n:.2f} {np.sqrt(sum(lipid_counts)) * n:.2f} {hg_tg_distance * 1.20:.2f}",
            f"  atoms {int(tg_index)}",
            f"    below plane 0. 0. 1. {hg_tg_distance * .20:.2f}",
            f"  end atoms",
            f"  atoms {int(hg_index)}",
            f"    over plane 0. 0. 1. {hg_tg_distance - hg_tg_distance * .20:.2f}",
            f"  end atoms",
            "end structure",
            "",
            # Bottom leaflet
            f"structure {lipid_name}.pdb",
            f"  number {int(lipid_count / 2)}",
            f"  inside box 0. 0. -{hg_tg_distance * 1.20:.2f} {np.sqrt(sum(lipid_counts)) * n:.2f} {np.sqrt(sum(lipid_counts)) * n:.2f} 0.",
            f"  atoms {int(tg_index)}",
            f"    over plane 0. 0. 1. -{hg_tg_distance * .20:.2f}",
            f"  end atoms",
            f"  atoms {int(hg_index)}",
            f"    below plane 0. 0. 1. -{hg_tg_distance - hg_tg_distance * .20:.2f}",
            f"  end atoms",
            "end structure",
            "",
        ])

    z_distance_water = 0

    if solvent is not None:
        # Add the section of the molecule to solvate if provided.
        density_water = 1e-21 / 18.01  # g/nm^3 / (g/mol) = mol/nm^3
        water_per_layer = solvent_count / 2
        water_volume_per_layer = water_per_layer / (density_water * 6.02e23) * 3 * 1000  # molecules/[(mol/nm^3)*(molecules/mol) *scale volume by 3 to allow for wiggle room * (angstrom^3/nm^3) = A^3
        z_distance_water = water_volume_per_layer / (np.sqrt(sum(lipid_counts)) * n) ** 2

        input_lines.extend([
            f"structure {solvent}.pdb",
            f"  number {int(solvent_count / 2)}",
            f"  inside box 0. 0. {max_z_distance * 1.20:.2f} {np.sqrt(sum(lipid_counts)) * n:.2f} {np.sqrt(sum(lipid_counts)) * n:.2f} {max_z_distance * 1.20 + 20:.2f}",
            "end structure",
            "",
            f"structure {solvent}.pdb",
            f"  number {int(solvent_count / 2)}",
            f"  inside box 0. 0. -{max_z_distance * 1.20 + 20:.2f} {np.sqrt(sum(lipid_counts)) * n:.2f} {np.sqrt(sum(lipid_counts)) * n:.2f} -{max_z_distance * 1.20:.2f}",
            "end structure",
            "",
        ])

    # Total box size
    x = np.sqrt(sum(lipid_counts)) * n
    y = np.sqrt(sum(lipid_counts)) * n
    z = 2 * max_z_distance * 1.1 + z_distance_water

    dims = [x, y, z]

    packmol_input = "\n".join(input_lines)

    # Write packmol input
    packmol_file_name = "packmol_input.inp"

    with open(packmol_file_name, "w") as file_handle:
        file_handle.write(packmol_input)

    #save all iumportant information to json file
    new_parameters_1 = {
        "lipid_file_names": lipid_file_names,
        "lipid counts": lipid_counts,
        "solvent": solvent,
        "solvent_count": solvent_count,
        "packmol input file": output,
        "box dimensions":dims,
        "packmol output file": output_file
    }

    append_and_save_json(file_name, data, new_parameters_1)


    return packmol_file_name, output_file, dims


def runPackmol(
    inp: str,
    parameterize : bool,
    force_field: str,
    hmr:bool): 
    """ Accepts a packmol input file, runs packmol, and parameterizes the output
    
    Parameters:
        inp = packmol input file
        force_field: OpenFF Forcefield to parameterize with"""

    command = f"packmol < {inp}"
    
    # Open the log file and redirect stdout and stderr
    log_file='packmol_output.log'
    with open(log_file, 'w') as f:
        subprocess.run(command, shell=True, check=True, stdout=f, stderr=subprocess.STDOUT)

    command_print = 'tail -n 28 packmol_output.log'
    subprocess.run(command_print, shell=True)

    #we can also sppecify if you want to parameterize your bilayer
    if parameterize == True:
        #first load in config data from experiment.json file
        with open('config.json', 'r') as config_file:
            config = json.load(config_file)

        dims = config['parameters']['box dimensions']
        lipids = config['parameters']['lipid_file_names']
        solvent = config['parameters']['solvent']
        number_of_lipids = config['parameters']['lipid counts']
        number_of_solvent = config['parameters']['solvent_count']
        output_packmol = config['parameters']['packmol output file']


        lipid_molecules = []
        for lipid in lipids:
            path = f'{lipid}.pdb'
            molecule = Molecule.from_file(path)  # Load the molecule from the PDB file
            lipid_molecules.append(molecule) # Load the molecule from the PDB file

        total_lipid_molecules = [lipid_molecules[i] for i, count in enumerate(number_of_lipids) for _ in range(count)]
        
        # Load solvent molecule
        solvent_path = f'{solvent}.pdb'
        solvent_molecule = Molecule.from_file(solvent_path)

        total_molecules = total_lipid_molecules + number_of_solvent * [solvent_molecule]

        topology = Topology.from_molecules(total_molecules)
        # Packmol bilayer to parametrize
        path = mdtraj.load(output_packmol)
        topology.set_positions(path.xyz[0] * unit.nanometer)
        topology.box_vectors =  np.array(dims)*0.1 * unit.nanometer #convert from angstoms (packmol default) to openff default nm

        forcefield = ForceField(force_field)

        interchange = Interchange.from_smirnoff(
            force_field=forcefield,
            topology=topology,
            charge_from_molecules = lipid_molecules
        )
        interchange  

        if hmr==True:
            interchange.to_gromacs(prefix = "bilayer", decimal = 3, hydrogen_mass = 3)
            print("Files parameterized and saved as bilayer.top and bilayer.gro. System configured for HMR")
        else:
            interchange.to_gromacs(prefix = "bilayer")
            print("Files parameterized and saved as bilayer.top and bilayer.gro. System not configured for HMR")
    else:
        print("System is not parameterized. Output= packmol.output.pdb")

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    #add both commands to
    parser.add_argument('command', choices=['_build_input_file', 'runPackmol'], help='Command to run')

    #add arguments for _build_input_file command
    parser.add_argument('-ln', '--lipidnames', nargs='+', type=str, help='Lipid names, e.g., POPC POPE')
    parser.add_argument('-lc', '--lipidcounts', nargs='+', type=int, help='Total lipid counts corresponding to lipid names, e.g., 64 64')
    parser.add_argument('-s', '--solvent', type=str, help='Structure to solvate (ex: TIP3P)')
    parser.add_argument('-sc', '--solventcount', type=int, help='Total number of solvent molecules (ex: 1000)')
    parser.add_argument('-tol', '--tolerance', type=float, default=2, help='Packmol convergence tolerance (default: 2)')
    parser.add_argument('-o', '--output', type=str, default='packmol_input.inp', help='Output file name (default: packmol_input.inp)')
    

    #add arguments for runPackmol command
    parser.add_argument('-i', '--fileinput', type=str, default = 'packmol_input.inp', help= 'Input file for Packmol (default: packmol_input.inp)')
    parser.add_argument('-p', '--parameterize', type=bool, default='True',help='Parameterize with OpenFF')
    parser.add_argument('-f', '--forcefield', type=str, default='openff-2.2.0.offxml',help='Parameterize with OpenFF force field, default openff-2.2.0.offxml')
    parser.add_argument('-hr', '--hmr', type=bool, default='True',help='Parameterize with hydrogen mass repartitioning. False=no HMR')

    args = parser.parse_args()

    if args.command == '_build_input_file':
        _build_input_file(lipid_file_names=args.lipidnames, 
            lipid_counts=args.lipidcounts, 
            solvent=args.solvent, 
            solvent_count=args.solventcount, 
            tolerance=args.tolerance, 
            output=args.output)
    elif args.command == 'runPackmol':
        runPackmol(inp=args.fileinput,
                   parameterize=args.parameterize,
                   force_field=args.forcefield,
                   hmr=args.hmr)
