import pandas as pd
import os

# # Sterols Add at some point
# #           'Cholesterol' : 'C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C',
# #           'Ergosterol' : 'C[C@H](/C=C/[C@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3C2=CC=C4[C@@]3(CC[C@@H](C4)O)C)C',
# #           'beta-Sitosterol' : 'CC[C@H](CC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C)C(C)C',
# #           'Stigmasterol' : 'CC[C@H](/C=C/[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C)C(C)C'}


#############################################################
# Build Dictrionary
Lipid_split = ['HG', 'TG']
head_acronyms = ['PA', 'PC', 'PE', 'PG', 'PS', 'PP']
tail_acronyms = ['DD', 'DC', 'DL', 'DM', 'DP', 'DS', 'PS', 'PO', 'PL', 'PE', 'SO', 'SL', 'DR', 'DO']
head_ids = ['HG']
tail_ids = ['SN1', 'SN2']


multi_index_hg = pd.MultiIndex.from_product([head_acronyms, head_ids], names=["acronym", "head_id"])
multi_index_tg = pd.MultiIndex.from_product([tail_acronyms, tail_ids], names=["acronym", "tail_id"])
index_hg = pd.MultiIndex.from_product([['HG'], head_acronyms, head_ids], names=["layer", "acronym", "id"])
index_tg = pd.MultiIndex.from_product([['TG'], tail_acronyms, tail_ids], names=["layer", "acronym", "id"])
full_index = index_hg.append(index_tg)
df = pd.DataFrame(index=full_index)#, columns=['Smiles'])

# add smiles strings
df.loc[('HG', 'PA', 'HG'), 'Structure'] = 'O[P@@]([O-])(O)(=O)'       # PA
df.loc[('HG', 'PC', 'HG'), 'Structure'] = 'O[P@@]([O-])(OCC[N+](C)(C)C)(=O)' # PC
df.loc[('HG', 'PE', 'HG'), 'Structure'] = 'O[P@@]([O-])(OCC[NH3+])(=O)'  # PE
df.loc[('HG', 'PG', 'HG'), 'Structure'] = 'O[P@@]([O-])(OC[C@@H](O)CO)(=O)' # PG
df.loc[('HG', 'PS', 'HG'), 'Structure'] = 'O[P@@]([O-])(OC[C@@H]([NH3+])(C(=O)([O-])))(=O)' # PS
df.loc[('HG', 'PP', 'HG'), 'Structure'] = 'O[P@@]([O-])(O[P@@]([O-])(O)(=O))(=O)' # PP
# Figure out smiles for PI 
df.loc[('TG', 'DD', 'SN1'), 'Structure'] = 'C'*9
df.loc[('TG', 'DD', 'SN2'), 'Structure'] = 'C'*9
df.loc[('TG', 'DC', 'SN1'), 'Structure'] = 'C'*10
df.loc[('TG', 'DC', 'SN2'), 'Structure'] = 'C'*10
df.loc[('TG', 'DL', 'SN1'), 'Structure'] = 'C'*11
df.loc[('TG', 'DL', 'SN2'), 'Structure'] = 'C'*11
df.loc[('TG', 'DM', 'SN1'), 'Structure'] = 'C'*13
df.loc[('TG', 'DM', 'SN2'), 'Structure'] = 'C'*13
df.loc[('TG', 'DP', 'SN1'), 'Structure'] = 'C'*15
df.loc[('TG', 'DP', 'SN2'), 'Structure'] = 'C'*15
df.loc[('TG', 'DS', 'SN1'), 'Structure'] = 'C'*17
df.loc[('TG', 'DS', 'SN2'), 'Structure'] = 'C'*17
df.loc[('TG', 'PS', 'SN1'), 'Structure'] = 'C'*17
df.loc[('TG', 'PS', 'SN2'), 'Structure'] = 'C'*15
df.loc[('TG', 'PO', 'SN1'), 'Structure'] = 'C'*15
df.loc[('TG', 'PO', 'SN2'), 'Structure'] = 'CCCCCCC\C=C/CCCCCCCC' #Notice that when adding double bonds, pay attention to how {SN1} and {SN2} lipids are written
df.loc[('TG', 'PL', 'SN1'), 'Structure'] = 'C'* 15
df.loc[('TG', 'PL', 'SN2'), 'Structure'] = 'CCCCCCC\C=C/C\C=C/CCCCC'
df.loc[('TG', 'PE', 'SN1'), 'Structure'] = 'C' *15
df.loc[('TG', 'PE', 'SN2'), 'Structure'] = 'CCCCCCCCCCC\C=C/CCCCCCCC'
df.loc[('TG', 'SO', 'SN1'), 'Structure'] = 'C' * 17
df.loc[('TG', 'SO', 'SN2'), 'Structure'] = 'CCCCCCC\C=/CCCCCCCC'
df.loc[('TG', 'SL', 'SN1'), 'Structure'] = 'C' * 17
df.loc[('TG', 'SL', 'SN2'), 'Structure'] = 'CCCCCCC\C=C/C\C=C/CCCCC'
df.loc[('TG', 'DR', 'SN1'), 'Structure'] = 'CCCC\C=C/CCCCCCC'
df.loc[('TG', 'DR', 'SN2'), 'Structure'] = 'CCCCCCC\C=C/CCCC'
df.loc[('TG', 'DO', 'SN1'), 'Structure'] = 'CCCCCCCC\C=C/CCCCCCC'
df.loc[('TG', 'DO', 'SN2'), 'Structure'] = 'CCCCCCC\C=C/CCCCCCCC'

os.makedirs('Dictionary', exist_ok=True)
df.to_csv("Dictionary/LipidLibrary.csv")
