# This code creates .inp files for molecules and sizes and then runs them in packmol and spits out the pdb file


import os

molecules = ['hexane','heptane','octane','decane','pentadecane']
sizes = {130:512, 180: 1024,250: 2048, 350:4096} #note that packmol is in angstroms 

for mol in molecules:
    #convert gro to pdb using gromacs
    gro = f"{mol}.gro"
    pdb = f"{mol}.pdb"
    os.system(f"gmx editconf -f {gro} -o {pdb}")

    for key, value in sizes.items():
        #create .inp file
        inp = f"{mol}_{value}.inp"
        with open(inp, 'w') as f:
            f.write(f"tolerance 2.0\n\n")
            f.write(f"filetype pdb\n\n")
            f.write(f"output {mol}_{value}.pdb\n\n")
            f.write(f"structure {mol}.pdb\n")
            f.write(f"  number {value}\n")
            f.write(f"  inside box 0. 0. 0. {key}. {key}. {key}.\n")
            f.write(f"end structure\n")

        os.system(f"packmol < {inp}")
