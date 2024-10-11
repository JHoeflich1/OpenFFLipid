import os
import shutil
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import distances


def saveInterchange(lipid_name, file_paths):
    '''This function take two inputs, the lipid name and the file paths that you want to save to the Lipids_parameterized dictionary 
    after pulling a lipid'''
    cwd = os.getcwd()
    target_directory = os.path.join(cwd, 'Dictionary', 'lipids_parameterized', lipid_name)

    # Create the directory if it doesn't exist
    if not os.path.exists(target_directory):
        os.makedirs(target_directory)
    
    # Copy each file to the target directory
    for file_path in file_paths:

        file_name = os.path.basename(file_path)
        destination_path = os.path.join(target_directory, file_name)
        shutil.copy(file_path, destination_path)
        print(f'{file_path} saved to {destination_path}')

class Lipid(object):
    """ Saving lipid parameters """

    def __init__(self, name, headgroup_atom, headgroup_atom_index, tailgroup_atom,tailgroup_atom_index, distance, experimental_density):

        self.name = name
        self.headgroup_atom = headgroup_atom
        self.headgroup_atom_index = headgroup_atom_index
        self.tailgroup_atom = tailgroup_atom
        self.tailgroup_atom_index = tailgroup_atom_index
        self.distance = distance
        self.experimental_density = experimental_density #this is in g/cm^3

def calcLipidLength(lipid, Lipid_name):
    '''Calculate the distance between headgroup and tailgroup atoms of a lipid,
       and update the lipid object with this information
       
       Parameters:
       lipid: An object representing the lipid with attributes headgroup_atom and tailgroup_atom
       Lipid_name: The name of the lipid used to locate the corresponding pdb
    '''
    cwd = os.getcwd()
    pdb_path = os.path.join(cwd, 'Dictionary', 'lipids_parameterized', Lipid_name, f'{Lipid_name}.pdb')
    print(pdb_path)

    ##getting an error with MDAnalysis 
    u_pdb = mda.Universe(pdb_path)
    # for atom in u_pdb.atoms:
    #     print(f"Atom: {atom.name}, Element: {atom.element}")

    hg_atom = lipid.headgroup_atom
    tg_atom = lipid.tailgroup_atom
    head_group = u_pdb.select_atoms(f'name {hg_atom}')
    tail_group = u_pdb.select_atoms(f'name {tg_atom}')
    
    # calculate the distance
    calc_distance = distances.distance_array(head_group.positions, tail_group.positions)
    lipid.distance = calc_distance[0][0] 
    
    # write in the atom index
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

folder_path = 'Dictionary'
csv_file_name = 'PulledLipid.csv'
csv_file_path = os.path.join(folder_path, csv_file_name)

def loadExistingData():
    """Load existing data from the CSV file or create an empty DataFrame with specified columns."""
    if os.path.exists(csv_file_path):
        return pd.read_csv(csv_file_path)
    else:
        return pd.DataFrame(columns=['Name', 'Volume', 'Headgroup Atom', 'Tailgroup Atom', 'Experimental Density'])


def saveLipidCsv(lipid):
    """Save a lipid to a CSV file, appending it if it doesn't already exist."""
    # Load existing data
    df = loadExistingData()
    
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

