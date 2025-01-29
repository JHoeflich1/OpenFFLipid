import subprocess

molecules = ['pentane','hexane','heptane','octane','decane','pentadecane']

# SLURM job template for liquids
job_template_liquids = """#!/bin/bash

#SBATCH --account=ucb500_asc1
#SBATCH --partition=amilan
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=24:00:00
#SBATCH --job-name={molecule}_liq
#SBATCH --output={molecule}_liq.%j.out

module purge
source /projects/nasc4134/pkgs/gromacs-2023.1/bin/GMXRC
ml gcc/14.2.0
ml openmpi/5.0.6

# Energy minimization
mpirun -np 1 gmx_mpi grompp -p {molecule}_liquid.top -f min.mdp -c {molecule}_liquid.gro -o min_{molecule}_liquid.tpr
mpirun -np 64 gmx_mpi mdrun -deffnm min_{molecule}_liquid

# NVT 
mpirun -np 1 gmx_mpi grompp -p {molecule}_liquid.top -f nvt_liquid.mdp -c min_{molecule}_liquid.gro -o nvt_{molecule}_liquid.tpr
mpirun -np 64 gmx_mpi mdrun -deffnm nvt_{molecule}_liquid

# NPT 
mpirun -np 1 gmx_mpi grompp -p {molecule}_liquid.top -f npt_liquid.mdp -c nvt_{molecule}_liquid.gro -o npt_{molecule}_liquid.tpr
mpirun -np 64 gmx_mpi mdrun -deffnm npt_{molecule}_liquid

"""


job_template_gas = """#!/bin/bash

#SBATCH --account=ucb500_asc1
#SBATCH --partition=amilan
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=24:00:00
#SBATCH --job-name={molecule}_gas
#SBATCH --output={molecule}_gas.%j.out

module purge
source /projects/nasc4134/pkgs/gromacs-2023.1/bin/GMXRC
ml gcc/14.2.0
ml openmpi/5.0.6

# Energy minimization
mpirun -np 1 gmx_mpi grompp -p {molecule}_gas.top -f min.mdp -c {molecule}_gas.gro -o min_{molecule}_gas.tpr
mpirun -np 64 gmx_mpi mdrun -deffnm min_{molecule}_gas

# NVT 
mpirun -np 1 gmx_mpi grompp -p {molecule}_gas.top -f nvt_gas.mdp -c min_{molecule}_gas.gro -o nvt_{molecule}_gas.tpr
mpirun -np 64 gmx_mpi mdrun -deffnm nvt_{molecule}_gas

"""
# Create and submit job scripts for each configuration
for molecule in molecules:
    # Create the job script content
    job_script_liq = job_template_liquids.format(molecule = molecule)
    job_script_gas = job_template_gas.format(molecule = molecule)
    
    # Write the job script to a file
    job_filename_l = f"{molecule}_liq.sh"
    with open(job_filename_l, "w") as job_file:
        job_file.write(job_script_liq)
    
    job_filename_g = f"{molecule}_gas.sh"
    with open(job_filename_g, "w") as job_file:
        job_file.write(job_script_gas)
        
    # Submit the job script using subprocess
    #subprocess.run(["sbatch", job_filename])  instead submit in bash: "for script in *.sh; do sbatch "$script"; done"