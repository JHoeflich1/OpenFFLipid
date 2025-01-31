import subprocess

molecules = ['pentane','hexane','heptane','octane','decane','pentadecane']

# SLURM job template for liquids
job_template_liquids = """#!/bin/bash

#SBATCH --account=ucb500_asc1
#SBATCH --partition=amilan
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=24:00:00
#SBATCH --job-name=liq_{molecule}
#SBATCH --output=liq_{molecule}.%j.out

module purge
source /projects/nasc4134/pkgs/gromacs-2023.1/bin/GMXRC
ml gcc/14.2.0
ml openmpi/5.0.6

# Energy minimization
gmx grompp -p {molecule}_liquid.top -f min.mdp -c {molecule}_liquid.gro -o min_{molecule}_liquid.tpr
gmx mdrun -deffnm min_{molecule}_liquid

# NVT 
gmx grompp -p {molecule}_liquid.top -f nvt_liquid.mdp -c min_{molecule}_liquid.gro -o nvt_{molecule}_liquid.tpr
gmx mdrun -deffnm nvt_{molecule}_liquid

# NPT 
gmx grompp -p {molecule}_liquid.top -f npt2_liquid.mdp -c nvt_{molecule}_liquid.gro -o npt2_{molecule}_liquid.tpr
gmx mdrun -deffnm npt2_{molecule}_liquid

"""


job_template_gas = """#!/bin/bash

#SBATCH --account=ucb500_asc1
#SBATCH --partition=amilan
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=24:00:00
#SBATCH --job-name=gas_{molecule}
#SBATCH --output=gas_{molecule}.%j.out

module purge
source /projects/nasc4134/pkgs/gromacs-2023.1/bin/GMXRC
ml gcc/14.2.0
ml openmpi/5.0.6

# Energy minimization
gmx grompp -p {molecule}_gas.top -f min.mdp -c {molecule}_gas.gro -o min_{molecule}_gas.tpr
gmx mdrun -deffnm min_{molecule}_gas

# NVT 
gmx grompp -p {molecule}_gas.top -f nvt_gas.mdp -c min_{molecule}_gas.gro -o nvt_{molecule}_gas.tpr
gmx mdrun -deffnm nvt_{molecule}_gas

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