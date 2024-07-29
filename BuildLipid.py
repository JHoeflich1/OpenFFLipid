import pandas as pd
import os

def load_dataframe():
    cwd = os.getcwd()
    dictionary_path = os.path.join(cwd, 'Dictionary', 'LipidLibrary.csv')
    df = pd.read_csv(dictionary_path, header=0)
    df.set_index(['layer', 'acronym', 'id'], inplace=True)
    return df

def makeLipidSmiles(Lipid, df):
    n = len(Lipid)  # Slipid Lipid string (first 2 letters is the Tail, last 2-3 is the HG)
    Lipid_sn12 = Lipid[:2]
    Lipid_hg = Lipid[2:n]

    HG_smiles = df.loc[('HG', Lipid_hg, 'HG'), 'Structure']
    SN1_smiles = df.loc[('TG', Lipid_sn12, 'SN1'), 'Structure']
    SN2_smiles = df.loc[('TG', Lipid_sn12, 'SN2'), 'Structure']

    Lipid_smiles = f"{SN1_smiles}C(=O)OC[C@H](OC(=O){SN2_smiles})(C{HG_smiles})"
    return Lipid_smiles

if __name__ == '__main__':
    df = load_dataframe()

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--lipid', type=str, default='POPC', help='Lipid to build, default POPC')
    args = parser.parse_args()
    lipid_smiles = makeLipidSmiles(args.lipid, df)
    print(lipid_smiles)



    # replacements = {
    #     '{SN1}': SN1_smiles,
    #     '{SN2}': SN2_smiles,
    #     '{HG}': HG_smiles
    #     }

    # molecule_rdkit = MolFromSmiles("{SN1}C(=O)OC[C@H](OC(=O){SN2})(C{HG})",True, replacements = replacements)
    # return molecule_rdkit

    #     # mol = Molecule.from_rdkit(molecule)
    #     # mol.generate_conformers()
    #     # mol.visualize()