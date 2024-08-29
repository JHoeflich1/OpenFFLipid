import os
import shutil

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