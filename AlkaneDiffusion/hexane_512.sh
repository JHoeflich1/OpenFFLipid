#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --account=ucb500_asc1
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=24:00:00
#SBATCH --job-name=hexane_512
#SBATCH --output=hexane_512.%j.out

module purge
source /projects/nasc4134/pkgs/gromacs-2023.1/bin/GMXRC
module load gcc
module load openmpi/4.1.1

# Energy minimization
mpirun -np 1 gmx_mpi grompp -p hexane_512.top -f min.mdp -c hexane_512.gro -o min_512.tpr
mpirun -np 64 gmx_mpi mdrun -deffnm min_512

# NPT
mpirun -np 1 gmx_mpi grompp -p hexane_512.top -f npt.mdp -c min_512.gro -o npt_512.tpr
mpirun -np 64 gmx_mpi mdrun -deffnm npt_512
