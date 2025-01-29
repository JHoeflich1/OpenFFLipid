#!/bin/bash

#SBATCH --account=ucb500_asc1
#SBATCH --partition=amilan
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=24:00:00
#SBATCH --job-name=pentane_gas
#SBATCH --output=pentane_gas.%j.out

module purge
source /projects/nasc4134/pkgs/gromacs-2023.1/bin/GMXRC
ml gcc/14.2.0
ml openmpi/5.0.6

# Energy minimization
mpirun -np 1 gmx_mpi grompp -p pentane_gas.top -f min.mdp -c pentane_gas.gro -o min_pentane_gas.tpr
mpirun -np 64 gmx_mpi mdrun -deffnm min_pentane_gas

# NVT 
mpirun -np 1 gmx_mpi grompp -p pentane_gas.top -f nvt_gas.mdp -c min_pentane_gas.gro -o nvt_pentane_gas.tpr
mpirun -np 64 gmx_mpi mdrun -deffnm nvt_pentane_gas

