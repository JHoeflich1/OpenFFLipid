import subprocess

# Define the molecule sizes
sizes = [512, 1024, 2048, 4096]
molecules = ['TIP3P','hexane','heptane','octane','decane','pentadecane']

# SLURM job template
job_template = """#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --account=ucb500_asc1
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=24:00:00
#SBATCH --job-name={molecule}_{size}_nvt
#SBATCH --output={molecule}_{size}_nvt.%j.out

module purge
source /projects/nasc4134/pkgs/gromacs-2023.1/bin/GMXRC
module load gcc
module load openmpi/4.1.1

# NVT
mpirun -np 1 gmx_mpi grompp -p {molecule}_{size}.top -f nvt.mdp -c npt_{molecule}_{size}.gro -o nvt_{molecule}_{size}.tpr
mpirun -np 64 gmx_mpi mdrun -deffnm npt_{molecule}_{size}
"""

# Create and submit job scripts for each configuration
for molecule in molecules:
    for size in sizes:
        # Create the job script content
        job_script = job_template.format(molecule = molecule, size=size)
        
        # Write the job script to a file
        job_filename = f"{molecule}_{size}_nvt.sh"
        with open(job_filename, "w") as job_file:
            job_file.write(job_script)
        
        # Submit the job script using subprocess
        subprocess.run(["sbatch", job_filename])