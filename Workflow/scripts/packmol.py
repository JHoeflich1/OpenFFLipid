from openff.toolkit import Molecule, Quantity, RDKitToolkitWrapper, Topology, unit
import pandas as pd
import numpy as np
import os
import shutil 
import subprocess

def load_pdb_files(lipid_name: str, source_dir: str, dest_dir: str):
    """
    Copy the PDB file of a lipid from the source directory to the destination directory.

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
    str
        The path to the copied PDB file in the destination directory.
    """
    # Define source and destination file paths
    source_file_path = os.path.join(source_dir, f'{lipid_name}.pdb')
    dest_file_path = os.path.join(dest_dir, f'{lipid_name}.pdb')

    # Check if the source file exists
    if not os.path.isfile(source_file_path):
        raise FileNotFoundError(f"The file '{source_file_path}' does not exist.")

    # Copy the file
    shutil.copy(source_file_path, dest_file_path)

def _build_input_file(
    lipid_file_name: str,#list[str], at some point extend this to multiple lipids
    lipid_counts: int,# list[int],
    solvent: str,
    solvent_count: int,
    tolerance: int,
    Library: str, # this is the name of lipid library of pulled lipids csv
    ) -> tuple[str, str]:
    """
    Calculate box size.
    Construct the packmol input file.

    Parameters
    ----------
    lipid_file_names: list of str
        The paths to the molecule pdb files.
    molecule_counts: list of int
        The number of each molecule to add.
    structure_to_solvate: str, optional
        The path to the structure to solvate.
    box_size: openff.units.Quantity
        The lengths of each side of the box we want to fill. This is the box
        size of the rectangular brick representation of the simulation box; the
        packmol box will be shrunk by the tolerance.
    tolerance: openff.units.Quantity
        The packmol convergence tolerance.

    Returns
    -------
    str
        The path to the input file.
    str
        The path to the output file.

    """

    lipid_library = pd.read_csv(f'Dictionary/{Library}') #load in the lipid library

    #move lipids pdb from dictionary to working dir
    cwd = os.getcwd()
    source_dir = os.path.join(cwd, 'Dictionary', 'lipids_parameterized', lipid_file_name)
    dest_dir = cwd

    load_pdb_files(lipid_file_name, source_dir, dest_dir)

    # Add the global header options in packmol
    output_file_path = "packmol_output.pdb"
    input_lines = [
        f"tolerance {tolerance:f}",
        "filetype pdb",
        f"output {output_file_path}",
        "",]
    
    # # notice how this is only for symetrical bilayers
    # total_lipids_leaflet = (np.sum(lipid_counts))/2

    # for lipid_name, count in zip(lipid_file_names, lipid_counts):
        
    # Find the corresponding row in the lipid library
    lipid_info = lipid_library[lipid_library['Name'] == lipid_file_name]
    if lipid_info.empty:
        raise ValueError(f"Lipid '{lipid_file_name}' not found in the library at {Library}")

    #load in the atom indexes and distance
    hg_tg_distance = float(lipid_info['HG/TG distance'].values[0])
    tg_index = float(lipid_info['Tailgroup Atom Index'].values[0])
    hg_index = float(lipid_info['Headgroup Atom Index'].values[0])
        
    ###########
    # I am assuming that lipid area in x-y space is 5 angstoms by 5 agstoms to give
    # room with packing.
    n = 6

    input_lines.extend(
        [
            ### top leaflet
            f"structure {lipid_file_name}.pdb",
            f"  number {int(lipid_counts/2)}",
            f"  inside box 0. 0. 0. {np.sqrt(lipid_counts)*n:.2f} {np.sqrt(lipid_counts)*n:.2f} {hg_tg_distance*1.20:.2f}",
            f"  atoms {int(tg_index)}",
            f"    below plane 0. 0. 1. {hg_tg_distance*.20:.2f}",
            f"  end atoms",
            f"  atoms {int(hg_index)}",
            f"    over plane 0. 0. 1. {hg_tg_distance-hg_tg_distance*.20:.2f}",
            f"  end atoms",
            "end structure",
            "",
            ### bottom leaflet
            f"structure {lipid_file_name}.pdb",
            f"  number {int(lipid_counts/2)}",
            f"  inside box 0. 0. -{hg_tg_distance*1.20:.2f} {np.sqrt(lipid_counts)*n:.2f} {np.sqrt(lipid_counts)*n:.2f} 0.",
            f"  atoms {int(tg_index)}",
            f"    over plane 0. 0. 1. -{hg_tg_distance*.20:.2f}",
            f"  end atoms",
            f"  atoms {int(hg_index)}",
            f"    below plane 0. 0. 1. -{hg_tg_distance-hg_tg_distance*.20:.2f}",
            f"  end atoms",
            "end structure",
            "",
        ],
    )
    # Add the section of the molecule to solvate if provided.
    # Please note that I assume all tip3p will fit in z length 20 angstoms 
    if solvent is not None:
        input_lines.extend(
            [
                f"structure {solvent}.pdb",
                f"  number {int(solvent_count/2)}",
                f"  inside box 0. 0. {hg_tg_distance*1.20:.2f} {np.sqrt(lipid_counts)*n:.2f} {np.sqrt(lipid_counts)*n:.2f} {hg_tg_distance*1.20 + 20:.2f}",
                "end structure",
                "",
                f"structure {solvent}.pdb",
                f"  number {int(solvent_count/2)}",
                f"  inside box 0. 0. -{hg_tg_distance*1.20 + 20:.2f} {np.sqrt(lipid_counts)*n:.2f} {np.sqrt(lipid_counts)*n:.2f} -{hg_tg_distance*1.20:.2f}",
                "end structure",
                "",
            ],
        )
    
    #total box size

    x = np.sqrt(lipid_counts)*n
    y = np.sqrt(lipid_counts)*n
    z = 2 * hg_tg_distance*1.20 + 20

    dims = [x,y,z]

    packmol_input = "\n".join(input_lines)

    # Write packmol input
    packmol_file_name = "packmol_input.inp"

    with open(packmol_file_name, "w") as file_handle:
        file_handle.write(packmol_input)

    return packmol_file_name, output_file_path, dims


def runPackmol(inp):
    command= f"packmol < {inp}" 
    subprocess.run(command, shell=True, check=True)

