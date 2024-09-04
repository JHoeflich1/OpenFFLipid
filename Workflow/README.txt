Steps for how to build and parameterize a lipid bilayer

First step is to print the list of available lipids. These lipids are parameterized with openFF Sage 2.2.0 and NAGL.
> python postPull.py print_pulled_lipids

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


Select a lipid (or a list of lipids) and load them into your directory using 
> load_




