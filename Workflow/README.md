# Steps to build and parameterize a lipid bilayer
Workflow for building and parameterizing a lipid bilayer from scratch. 

## Python environment 
All steps in workflow were performed with the attached conda environment `environment.yml`.
This can be installed using `mamba env create -f test.yml`

To activate environemnt: `conda activate openff_workflow`

## Workflow 
Navigate to the Workflow directory. This contains folders scripts, dictionary, a Workflow jupyter notebook. This workflow is configured so that everything can be run on the command line. 

1. First step is to print the list of available lipids in the current library. These lipids are parameterized with openFF Sage 2.2.0 and NAGL model `openff-gnn-am1bcc-0.1.0-rc.3.pt`.

`python scripts/postPull.py print_pulled_lipids`


If the lipid you wish to simulate is not in the lipid library, look at section below `Adding Lipids`


2. After selecting your lipids and counts that you want to pack, generate the packmol input structure file as follows:


`python scripts/packmol.py _build_input_file -ln POPC POPE -lc 64 64 -s TIP3P -sc 1000 -t 2` 

use --help for detailed list on flags
This script returns a config.json file that stores relevant data, and a packmol input file (default packmol_input.inp) that contains the packing constraints for your system. 


3. Initiate packmol and parameterize your system:


`python scripts/packmol.py runPackmol -i packmol_input.inp`


Parameterize your system with HMR is default, if no HMR is desired, use flag `-h False`.  This produces files named `bilayer.gro` and `biayer.top`. It is recommended to visualize this structure in VMD or another molecular visualization software before running simulations.
For a 128 lipid system this shoudl take approx ~5 minutes


4. Runnuing molecular dynamics simulations:
It is not recommended to modify equilibration steps. If necessary, be aware that bilayer may phase separate or form pores. These steps are configured to run on my CU Boulder reserach computing account.


Energy minimization 
NVT 1fs
NPT 1fs
Repeat NVT and NPT steps 2x with larger timsteps 2 fs

## Adding Lipids

If the lipid you wish to simulate is not in the list of available lipids, you may add it by following the instructions below:


Run parameterize_new_lipid script, the two inputs are your lipid name `{Lipid}` and it`s corresponding smiles string


`python postPull.py parameterize_new_lipid --lipidname {Lipid} --lipidsmiles {smiles_string}`



ex. {Lipid} = POPC ; {smiles_string} = `O([P@@]([O-])(OCC[NH3+])(=O))C[C@H](COC(=O)CCCCCCCCCCCCCCC)(OC(=O)CCCCCCC\\C=C/CCCCCCCC)`



This outputs `{Lipid}.gro` and `{Lipid}.top`. Visualize the lipid in VMD or another preferred software program and select one atom in the lipid`s headgroup (preferably a terminal heavy atom) as well as one terminal heavy atom in the SN1 or SN2 tail. This will be used to straighten out the lipid and estimate the lipids total length.


Run add_new_lipid, with the {Lipid}, {smiles_string}, and the atom name of the head group atom {hg_name}, and tail atom {tg_name}


`python postPull.py add_new_lipid -ln {Lipid} -ls {smiles_string} -hg {hg_name} -tg {tg_name}`


ex. {hg_name} = N1x ; {tg_name} = C25x

This saves your lipid topology, and pulled coordinate file in `/Dictionary/lipids_parameterized/{Lipid}`, as well as adds an entry in the `/Dictionary/PulledLipid.csv` table. Please note that the default openFF force field version and NAGL version used for parameterization is: `openff-gnn-am1bcc-0.1.0-rc.3.pt`, `openff-2.2.0.offxml`









