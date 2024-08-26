import pandas as pd
import os

def load_dataframe():
    cwd = os.getcwd()
    dictionary_path = os.path.join(cwd, 'Dictionary', 'LipidLibrary.csv')
    df = pd.read_csv(dictionary_path, header=0)
    df.set_index(['layer', 'acronym', 'id'], inplace=True)
    return df

def list_available_lipids(df):
    # Filter the DataFrame based on the 'layer' column values
    hg_acronyms = df[df.index.get_level_values('layer') == 'HG'].index.get_level_values('acronym').unique().tolist()
    tg_acronyms = df[df.index.get_level_values('layer') == 'TG'].index.get_level_values('acronym').unique().tolist()
    sterol_acronyms = df[df.index.get_level_values('layer') == 'Sterol'].index.get_level_values('acronym').unique().tolist()
    
    print("Available headgroups: ", hg_acronyms)
    print("Available tailgroups: ", tg_acronyms)
    print("Available sterols: ", sterol_acronyms)
    
    # Create all possible lipid combinations by pairing each TG with each HG
    # available_lipids = [tg + hg for tg in tg_acronyms for hg in hg_acronyms]
    
    # return available_lipids

def makeLipidSmiles(Lipid, df):
    if len(Lipid) > 4:
        # Handle sterols
        if Lipid not in df[df.index.get_level_values('layer') == 'Sterol'].index.get_level_values('acronym'):
            raise ValueError(f"Sterol {Lipid} not found in the dictionary.")
        Sterol_smiles = df.loc[('Sterol', Lipid, 'st'), 'Structure']
        return Sterol_smiles
    
    Lipid_sn12 = Lipid[:2]
    Lipid_hg = Lipid[2:]
    
    # Check if the TG and HG exist in the DataFrame
    if not Lipid_sn12 in df[df.index.get_level_values('layer') == 'TG'].index.get_level_values('acronym'):
        raise ValueError(f"Tail group {Lipid_sn12} not found in the dictionary.")
    if not Lipid_hg in df[df.index.get_level_values('layer') == 'HG'].index.get_level_values('acronym'):
        raise ValueError(f"Head group {Lipid_hg} not found in the dictionary.")
    
    # Retrieve the SMILES strings
    HG_smiles = df.loc[('HG', Lipid_hg, 'hg'), 'Structure']
    HG_pull_atom = df.loc[('HG', Lipid_hg, 'hg'), 'pull']
    SN1_smiles = df.loc[('TG', Lipid_sn12, 'sn1'), 'Structure']
    SN2_smiles = df.loc[('TG', Lipid_sn12, 'sn2'), 'Structure']
    
    # Construct the final SMILES string
    Lipid_smiles = f"{HG_smiles}C[C@H](COC(=O){SN1_smiles})(OC(=O){SN2_smiles})"
    
    return Lipid_smiles, HG_pull_atom

if __name__ == '__main__':
    df = load_dataframe()
    list_available_lipids(df)  # Optionally list available lipids
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--lipid', type=str, default='POPC', help='Specify lipid to build, default POPC')
    args = parser.parse_args()
    lipid_smiles = makeLipidSmiles(args.lipid, df)
