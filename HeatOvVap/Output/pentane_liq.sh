#!/bin/bash

#SBATCH --account=ucb500_asc1
#SBATCH --partition=amilan
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=24:00:00
#SBATCH --job-name=pentane_liq
#SBATCH --output=pentane_liq.%j.out

module purge
source /projects/nasc4134/pkgs/gromacs-2023.1/bin/GMXRC
ml gcc/14.2.0
ml openmpi/5.0.6

# Energy minimization
mpirun -np 1 gmx_mpi grompp -p pentane_liquid.top -f min.mdp -c pentane_liquid.gro -o min_pentane_liquid.tpr
mpirun -np 64 gmx_mpi mdrun -deffnm min_pentane_liquid

# NVT 
mpirun -np 1 gmx_mpi grompp -p pentane_liquid.top -f nvt_liquid.mdp -c min_pentane_liquid.gro -o nvt_pentane_liquid.tpr
mpirun -np 64 gmx_mpi mdrun -deffnm nvt_pentane_liquid

# NPT 
mpirun -np 1 gmx_mpi grompp -p pentane_liquid.top -f npt_liquid.mdp -c nvt_pentane_liquid.gro -o npt_pentane_liquid.tpr
mpirun -np 64 gmx_mpi mdrun -deffnm npt_pentane_liquid

