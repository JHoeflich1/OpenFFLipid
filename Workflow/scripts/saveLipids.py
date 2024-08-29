import os
import shutil
import pandas as pd
import MDAnalysis as mda

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

class Lipid(object):
    """ Saving lipid parameters """

    def __init__(self, name, headgroup_atom, tailgroup_atom, distance, experimental_volume):

        self.name = name
        self.headgoup_atom = headgroup_atom
        self.tailgroup_atom = tailgroup_atom
        self.distance = distance
        self.experimental_volume = experimental_volume #this is in g/cm^3

def calcLipidLength(lipid, Lipid_name):
    '''This code accepts a lipid object as well as the pulled lipid pdb
    This code will use mdanalysis to calcualte the distance between the headgroup and tail group, and add them to the lipid object'''
    u = mda.Universe(f'/Dictionary/lipids_parameterized/{Lipid_name}/{Lipid_name}.pdb')
    hg_atom = lipid.headgoup_atom
    tg_atom = lipid.tailgroup_atom
    head_group = u.select_atoms(f'name {hg_atom}') #select atom in head group, in this case Nitrogen
    tail_group = u.select_atoms(f'name {tg_atom}') #select atom in tail group, in this case last carbon in unsaturated tail
    calc_distance = mda.distances.distance_array(head_group.positions,tail_group.positions)
    lipid.distance = calc_distance


def lipid_to_dict(lipid):
    return {
        'Name': lipid.name,
        'Volume': lipid.volume,
        'Headgroup Atom': lipid.headgroup_atom,
        'Tailgroup Atom': lipid.tailgroup_atom,
        'Experimental Density': lipid.experimental_density
    }

def load_existing_data(file_path):
    if os.path.exists(file_path):
        return pd.read_csv(file_path)
    else:
        return pd.DataFrame(columns=['Name', 'Volume', 'Headgroup Atom', 'Tailgroup Atom', 'Experimental Density'])

def save_lipid_to_csv(lipid, file_path):
    # Load existing data
    df = load_existing_data(file_path)
    
    # Convert lipid to dictionary
    lipid_dict = lipid_to_dict(lipid)
    
    # Check if the lipid already exists in the DataFrame
    if lipid_dict['Name'] not in df['Name'].values:
        # Append new lipid data
        df = df.append(lipid_dict, ignore_index=True)
        
        # Save updated DataFrame to CSV
        df.to_csv(file_path, index=False)
        print(f"Lipid '{lipid.name}' saved to CSV.")
    else:
        print(f"Lipid '{lipid.name}' already exists in CSV.")
