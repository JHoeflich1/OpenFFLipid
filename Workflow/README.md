# Steps for how to build and parameterize a lipid bilayer

#for this I have been using my openff_clone environemnt. Will need to confugre a proper env 

Navigate to the Wofkflow directory. This contains folders scripts, dictionary, a Workflow jupyter notebook.

First step is to print the list of available lipids in the current library. These lipids are parameterized with openFF Sage 2.2.0 and NAGL.
> python scripts/postPull.py print_pulled_lipids

_______________________________________________________
If the lipid you wish to simulate is not in the list of available lipids, you may add it by following the instructions below:
Run parameterize_new_lipid script, the two inputs are your lipid name {Lipid} and it's corresponding smiles string
> python postPull.py parameterize_new_lipid --lipidname {Lipid} --lipidsmiles {smiles_string}
ex. {Lipid} = POPC ; {smiles_string} = 'O([P@@]([O-])(OCC[NH3+])(=O))C[C@H](COC(=O)CCCCCCCCCCCCCCC)(OC(=O)CCCCCCC\\C=C/CCCCCCCC)'
This outputs {Lipid}.gro and {Lipid}.top. Visualize the lipid in VMD or another preferred software program and select one atom in the lipid's headgroup (preferably a terminal heavy atom)
as well as one terminal heavy atom in the SN1 or SN2 tail. This will be used to straighten out the lipid and estimate the lipids total length
Run add_new_lipid, with the {Lipid}, {smiles_string}, and the atom name of the head group atom {hg_name}, and tail atom {tg_name}
python postPull.py add_new_lipid -ln {Lipid} -ls {smiles_string} -hg {hg_name} -tg {tg_name}
ex. {hg_name} = N1x ; {tg_name} = C25x
This saves your lipid topology, and pulled coordinate file in /Dictionary/lipids_parameterized/{Lipid}, as well as adds an entry in the /Dictionary/PulledLipid.csv table

Please note that the default openFF force field version and NAGL version used for parameterization is:
"openff-gnn-am1bcc-0.1.0-rc.3.pt"
"openff-2.2.0.offxml"
_______________________________________________________


After selecting your lipids and counts that you want to pack, generate the packmol input structure file as follows:
> python scripts/packmol.py _build_input_file -ln POPC POPE -lc 64 64 -s TIP3P -sc 1000 -t 2
This script also stores a config.json file that stores relevant data. This creates a file (default packmol_input.inp) that contains the packing constraints for your system. 

Initiate packmol and barameterize your system by calling 
> python scripts/packmol.py runPackmol -i packmol_input.inp
The default forcefield is -f openff-2.2.0.offxml, but can be changed
Parameterize your system with HMR is default, if no HMR is desired, use flag -h False 
This produces a files named bilayer.gro and biayer.top. It is recommended to visualize this structure in VMD or another molecular visualization software before running simulations

Runnuing molecular dynamics simulations:
It is not recommended to modify equilibration steps. If necessary, proceed with caution. 
These steps are configured to run on my CU Boulder reserach computing account.

Energy minimization 
NVT 1fs
NPT 1fs
Repeat NVT and NPT steps 2x with larger timsteps 2 fs







