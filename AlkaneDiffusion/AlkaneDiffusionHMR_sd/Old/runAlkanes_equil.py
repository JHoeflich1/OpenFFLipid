import subprocess

# Define the molecule sizes
sizes = [512, 1024, 2048]
molecules = ['pentane','hexane','heptane','octane','decane','pentadecane','water']

# SLURM job template
job_template = """#!/bin/bash

#SBATCH --account=ucb500_asc1
#SBATCH --time=24:00:00
#SBATCH --partition=amilan
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --job-name={molecule}_{size}_equil1
#SBATCH --output={molecule}_{size}_equil1.%j.out


module purge
source /projects/nasc4134/pkgs/gromacs-2023.1/bin/GMXRC
module load gcc
module load openmpi/4.1.1
module load gromacs


# Energy minimization
gmx grompp -p {molecule}_{size}.top -f min.mdp -c {molecule}_{size}.gro -o min_{molecule}_{size}.tpr
gmx mdrun -deffnm min_{molecule}_{size}

# First NVT for equilibration
gmx grompp -p {molecule}_{size}.top -f nvt1.mdp -c min_{molecule}_{size}.gro -o nvt1_{molecule}_{size}.tpr
gmx mdrun -deffnm nvt1_{molecule}_{size}

# NPT
gmx grompp -p {molecule}_{size}.top -f npt1.mdp -c nvt1_{molecule}_{size}.gro -o npt1_{molecule}_{size}.tpr
gmx mdrun -deffnm npt1_{molecule}_{size}

"""

# Create and submit job scripts for each configuration
for molecule in molecules:
    for size in sizes:
        # Create the job script content
        job_script = job_template.format(molecule = molecule, size=size)

        # Write the job script to a file
        job_filename = f"{molecule}_{size}_equil1.sh"
        with open(job_filename, "w") as job_file:
            job_file.write(job_script)

        # Submit the job script using subprocess
        # subprocess.run(["sbatch", job_filename])  # for some reason alpine does not like this subprocess call and it wont let me load modules. hell = truw 
