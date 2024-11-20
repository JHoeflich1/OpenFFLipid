import os
import pandas as pd
import sys
import numpy as np
import glob
import time
import shutil
import subprocess
from openff.toolkit import ForceField, Molecule, Topology
from openff.interchange import Interchange
from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
import MDAnalysis as mda
from MDAnalysis.analysis import distances


def print_pulled_lipids():
    '''Lists the available lipids in the PulledLipid.csv dataframe. These are pre-parameterized and 
    straightened out lipids that are optimized for packing in packmol'''
    cwd = os.getcwd()
    parent_directory = os.path.dirname(cwd)
    target_directory = os.path.join(parent_directory, 'Workflow', 'Dictionary', 'PulledLipid.csv')
    
    pulled_lib = pd.read_csv(target_directory)

    print("Available lipids: ")
    for name in pulled_lib['Name']:
        print(name)



def parameterize_new_lipid(lipid_name, lipid_smiles):
    ''' Parameterize a molecule object with openff Interchange
    Inputs:
    lipid: a smiles string 
    Lipid_name: a shorthand name for the lipid used for residue names ex.POPC
    Returns:
    TOP and GRO file
    '''
    lipid = Molecule.from_smiles(lipid_smiles,allow_undefined_stereo=True)
    lipid.generate_conformers()
    
    lipid.assign_partial_charges("openff-gnn-am1bcc-0.1.0-rc.3.pt", toolkit_registry=NAGLToolkitWrapper())
    lipid.partial_charges

    lipid.name = lipid_name
    for i, atom in enumerate(lipid.atoms, 0):
        atom.metadata["residue_name"] = lipid_name
    lipid.generate_unique_atom_names() 
    topology = Topology.from_molecules([lipid])

    # Specify forcefield
    forcefield = ForceField("openff-2.2.0.offxml")

    interchange = Interchange.from_smirnoff(
        force_field=forcefield,
        topology=topology,
        charge_from_molecules = [lipid]
    )
    interchange  

    interchange.to_top(f"{lipid_name}.top") #, decimal = 3, hydrogen_mass = 3) #for hygrogen mass repartitioning.. probably need another topology file filder for HMR
    # or is there from grompp option to change topology file that can be used downstream?
    interchange.to_gro(f"{lipid_name}.gro")


def pull_new_lipid(Lipid_name, pull_atom):
    '''Runs the gmx pulling command to straighten out the lipid and cleans up generated files. Print time taken to run
    Inputs:
    Lipid_name: a shorthand name for the lipid used for residue names ex.POPC
    pull_atom: the atom in the lipid headroup from which the lipid will be straightened out
    Returns:
    straightneed out lipid in pdb file format
    '''
    with open(os.devnull, 'w') as devnull:
        command_ndx = f"gmx make_ndx -f {Lipid_name}.gro -o {pull_atom}.ndx"
        os.system(f"echo 'a {pull_atom}\nname 3 pull_atom\nq\n' | {command_ndx}")

        command_grompp = f"gmx grompp -f scripts/runPull.mdp -c {Lipid_name}.gro -p {Lipid_name}.top -n {pull_atom}.ndx -o pull.tpr -maxwarn 2"
        subprocess.run(command_grompp, shell=True)

        command_pull = "gmx mdrun -deffnm pull"
        subprocess.run(command_pull, shell=True)

        command_pdb = f"gmx editconf -f pull.gro -o {Lipid_name}.pdb"
        subprocess.run(command_pdb, shell=True)
        

    # Clean up by deleting all files starting with "pull."
    files_to_delete = glob.glob('pull*')
    for file in files_to_delete:
        os.remove(file)
    # remove index file
    index_files = glob.glob('*.ndx')
    for file in index_files:
        os.remove(file)

def saveInterchange(lipid_name, file_paths, parent_directory):
    ''' Lipids that you want to save to the Lipids_parameterized dictionary after they are pulled
    inputs:
    lipid_name: the lipid name
    file_paths: the names of the PDB and TOP files that you want saved ex. ['POPC.pdb','POPC.top']
    parent_directory: ignore for now
    '''

    target_directory = os.path.join(parent_directory, 'Dictionary', 'lipids_parameterized', lipid_name)

    # Create the directory if it doesn't exist
    if not os.path.exists(target_directory):
        os.makedirs(target_directory)
    
    # copy each file to the target directory
    for file_path in file_paths:

        file_name = os.path.basename(file_path)
        destination_path = os.path.join(target_directory, file_name)
        shutil.copy(file_path, destination_path)
        print(f'{file_path} saved to {destination_path}')

class Lipid(object):
    """ Saving lipid parameters into an object"""

    def __init__(self, name, headgroup_atom, headgroup_atom_index, tailgroup_atom,tailgroup_atom_index, distance, experimental_density):

        self.name = name
        self.headgroup_atom = headgroup_atom
        self.headgroup_atom_index = headgroup_atom_index
        self.tailgroup_atom = tailgroup_atom
        self.tailgroup_atom_index = tailgroup_atom_index
        self.distance = distance
        self.experimental_density = experimental_density #this is in g/cm^3


def calcLipidLength(lipid, Lipid_name, parent_directory):
    '''Calculate the distance between headgroup and tailgroup atoms of a lipid,
       and update the lipid object with this information
       
       Inputs:
       lipid: An object representing the lipid with attributes headgroup_atom and tailgroup_atom
       Lipid_name: The name of the lipid used to locate the corresponding pdb
    '''
    pdb_path = os.path.join(parent_directory, 'Dictionary', 'lipids_parameterized', Lipid_name, f'{Lipid_name}.pdb')

    
    u_pdb = mda.Universe(pdb_path)
    hg_atom = lipid.headgroup_atom
    tg_atom = lipid.tailgroup_atom
    head_group = u_pdb.select_atoms(f'name {hg_atom}')
    tail_group = u_pdb.select_atoms(f'name {tg_atom}')
    
    # calculate the distance and save into object
    calc_distance = distances.distance_array(head_group.positions, tail_group.positions)
    lipid.distance = calc_distance[0][0] 
    
    # write in the atom index into object
    lipid.headgroup_atom_index = head_group[0].index 
    lipid.tailgroup_atom_index = tail_group[0].index 

def lipidToDict(lipid):
    return {
        'Name': lipid.name,
        'Headgroup Atom': lipid.headgroup_atom,
        'Headgroup Atom Index': lipid.headgroup_atom_index,
        'Tailgroup Atom': lipid.tailgroup_atom,
        'Tailgroup Atom Index': lipid.tailgroup_atom_index,
        'HG/TG distance': lipid.distance,
        'Experimental Density': lipid.experimental_density
    }

def loadExistingData(cwd):
    """Load existing data from the CSV file or create an empty DataFrame with specified columns."""
    csv_file_path = os.path.join(cwd, 'Dictionary', 'PulledLipid.csv')

    if os.path.exists(csv_file_path):
        df = pd.read_csv(csv_file_path)
    else:
        # Create an empty DataFrame with specified columns if the file does not exist
        df = pd.DataFrame(columns=['Name', 'Volume', 'Headgroup Atom', 'Tailgroup Atom', 'Experimental Density'])
    
    return df, csv_file_path  # Ensure both values are returned consistently


def saveLipidCsv(lipid, cwd):
    """Save a lipid to a CSV file, appending it if it doesn't already exist."""
    # Load existing data
    df, csv_file_path = loadExistingData(cwd)
    
    # Convert lipid to dictionary
    lipid_dict = lipidToDict(lipid)
    
    # Does lipid exist in dataframe
    if lipid_dict['Name'] not in df['Name'].values:
        new_df = pd.DataFrame([lipid_dict])
        
        df = pd.concat([df, new_df], ignore_index=True)
        
        # Save to csv if not already in there
        df.to_csv(csv_file_path, index=False)
        print(f"Lipid '{lipid.name}' saved to CSV location: {csv_file_path}.")
    else:
        print(f"Lipid '{lipid.name}' already exists in CSV location: {csv_file_path}")

def add_new_lipid(lipid_name, lipid_headgroup, lipid_tailgroup):
    """This function is meant to allow the user to add their own lipid smiles strings.
    It takes in a unique lipid, pulls it and saves it to the lipid dictionry"""
    start_time = time.time() 

    cwd = os.getcwd()
    # parent_directory = os.path.dirname(cwd)

    #create .top file and .gro file to be pulled
    # parameterize_new_lipid(lipid_smiles, lipid_name)
    pull_new_lipid(lipid_name, lipid_headgroup)
    
    # #save parameterized top file and pulled lipid coordinate file
    file_paths = [f'{lipid_name}.pdb', f'{lipid_name}.top'] 
    saveInterchange(lipid_name, file_paths, cwd)

    #Now I want to import a lipid object and calculate the distance and 
    lipid = Lipid(
        name=lipid_name,
        headgroup_atom=lipid_headgroup,  
        headgroup_atom_index=None, 
        tailgroup_atom=lipid_tailgroup,
        tailgroup_atom_index=None,
        distance=None,  
        experimental_density=None  # Example value in g/cm^3, not necessary but may be helpful in the future?
)

    # Calculate lipid length (distance between headgroup and tailgroup)
    calcLipidLength(lipid, lipid_name, cwd)

    # Save the lipid to the CSV file with all its corresponding information
    saveLipidCsv(lipid, cwd)

    #delete the files created in scripts folder
    # files_to_delete = glob.glob(f'{lipid_name}.*')
    # for file in files_to_delete:
    #     os.remove(file)


    end_time = time.time()  # End timing
    elapsed_time = end_time - start_time

    print(f"Total time taken to run add_new_lipid for {lipid_name}: {elapsed_time:.2f} seconds")


if __name__ == '__main__':
    #globals()[sys.argv[1]]()
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['print_pulled_lipids','parameterize_new_lipid', 'add_new_lipid'], help='Command to run')
    parser.add_argument('-ln', '--lipidname', type=str, help='Lipid name, e.g., POPC')
    parser.add_argument('-ls', '--lipidsmiles', type=str, help='Lipid isomeric smiles string')
    parser.add_argument('-hg', '--headgroup', type=str, help='Name of atom in head group')
    parser.add_argument('-tg', '--tailgroup', type=str, help='Name of terminal carbon atom in SN1 or SN2 tail')
    
    args = parser.parse_args()

    if args.command == 'print_pulled_lipids':
        print_pulled_lipids()
    elif args.command == 'parameterize_new_lipid':
        if args.lipidname and args.lipidsmiles:
            parameterize_new_lipid(args.lipidname, args.lipidsmiles)
    elif args.command == 'add_new_lipid':
        if args.lipidname and args.headgroup and args.tailgroup:
            add_new_lipid(args.lipidname, args.headgroup, args.tailgroup)
        else:
            print("Error: Missing required arguments for add_new_lipid.")