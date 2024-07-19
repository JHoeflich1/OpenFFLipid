import pandas as pd
from rdkit.Chem.rdmolfiles import MolFromSmiles
import os
import sys

# read lipid dictionary
df= pd.read_csv("/Dictionary/LipidDictionary.csv", header=[0, 1] )
print(df)


## create smiles string
def makeLipidSmiles(Lipid):
    n = len(Lipid)  # length of the name
    Lipid_sn12 = Lipid[:2]
    Lipid_hg = Lipid[2:n]



    HG_smiles = df.loc[('HG', Lipid_hg, 'HG'), 'Structure']
    SN1_smiles = df.loc[('TG', Lipid_sn12, 'SN1'), 'Structure']
    SN2_smiles = df.loc[('TG', Lipid_sn12, 'SN2'), 'Structure']

    Lipid_smiles = f"{SN1_smiles}C(=O)OC[C@H](OC(=O){SN2_smiles})(C{HG_smiles})"
    print(f'Output smiles for {Lipid} is {"{SN1_smiles}C(=O)OC[C@H](OC(=O){SN2_smiles})(C{HG_smiles})"}')

    replacements = {
        '{SN1}': SN1_smiles,
        '{SN2}': SN2_smiles,
        '{HG}': HG_smiles
        }

    molecule_rdkit = MolFromSmiles("{SN1}C(=O)OC[C@H](OC(=O){SN2})(C{HG})",True, replacements = replacements)
    return molecule_rdkit

        # mol = Molecule.from_rdkit(molecule)
        # mol.generate_conformers()
        # mol.visualize()
if __name__=='__main__':
    if len(sys.argv) != 2:
        print("Usage: python3 makeDict.py <Lipid>")
        sys.exit(1)
    Lipid = sys.argv[1]
    makeLipidSmiles(Lipid)