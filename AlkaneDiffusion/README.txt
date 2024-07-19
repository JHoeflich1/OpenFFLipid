Validate diffusion parameters (coefficient, viscosity) for OpenFF Alkanes

Pentane,Hexane,Heptane,Octane,Decane,Pentadecane, water @ 300K, 1 bar


Pipeline:
1. Run buildAlkane.py, buildWater.py to create alkanes and water from smiles string, pack using openFF packmol toolkit, parameterize with interchange + use HMR
2. Run runAlkanes_equil.py to create .sh run scripts for each simulation in RC account 
    submit .sh with "for script in *equil1.sh; do sbatch "$script"; done " to run
3. Run ChangeBoxSize.py to modify box size before production run

######################################################################
6. Run runAlkances_nvt.py to create .sh run scripts
    submit .sh with " for script in *nvt_1.sh; do sbatch "$script"; done " to run NVT
7. Run MSD_per_particle.py for each molecule size and each moelcule. Will need to adjust for all molecules 
8. Run Bootstrapping_diffusion to get the average Ds and Sterr (computed from bootstrapping)
8. Run Diffusion_fitting.py to calcaulte and plot the diffusion coefficients for different box sizes
9. Create folders for each moelcule and separate files into their molecules 
