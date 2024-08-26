import pandas as pd
import os

# Suppress PerformanceWarning for chained assignments
pd.options.mode.chained_assignment = None

# Build Dictrionary
Lipid_split = ['HG', 'TG', 'Sterol']

head_acronyms = ['PA', 'PC', 'PE', 'PG', 'PS', 'PP']
tail_acronyms = ['DD', 'DC', 'DL', 'DM', 'DP', 'DS', 'PS', 'PO', 'PL', 'PE', 'SO', 'SL', 'DR', 'DO']
sterol_acronyms = ['cholesterol','ergosterol', 'beta-Sitosterol','stigmasterol', 'cortisol', 'corticosterone', 'aldosterone','progesterone','beta-estradiol','estradiol']

head_ids = ['hg']
tail_ids = ['sn1', 'sn2']
sterol_ids = ['st']


multi_index_hg = pd.MultiIndex.from_product([head_acronyms, head_ids], names=["acronym", "head_id"])
multi_index_tg = pd.MultiIndex.from_product([tail_acronyms, tail_ids], names=["acronym", "tail_id"])
multi_index_sterol = pd.MultiIndex.from_product([sterol_acronyms, sterol_ids], names=["acronym", 'sterol_id'])

index_hg = pd.MultiIndex.from_product([['HG'], head_acronyms, head_ids], names=["layer", "acronym", "id"])
index_tg = pd.MultiIndex.from_product([['TG'], tail_acronyms, tail_ids], names=["layer", "acronym", "id"])
index_sterol = pd.MultiIndex.from_product([['Sterol'],sterol_acronyms, sterol_ids], names=["layer", "acronym", "id"])
full_index = index_hg.append(index_tg).append(index_sterol)
# full_index = comb_index.append(index_sterol)
df = pd.DataFrame(index=full_index)#, columns=['Smiles'])

# Ensure the index is sorted before setting values
df = df.sort_index()

# add smiles strings
df.at[('HG', 'PA', 'hg'), 'Structure'] = 'O([P@@]([O-])(O)(=O))'       # PA
df.at[('HG', 'PA', 'hg'), 'pull'] = 'P1x'       # PA
df.at[('HG', 'PC', 'hg'), 'Structure'] = 'O([P@@]([O-])(OCC[N+](C)(C)C)(=O))' # PC
df.at[('HG', 'PC', 'hg'), 'pull'] = 'N1x' # PC
df.at[('HG', 'PE', 'hg'), 'Structure'] = 'O([P@@]([O-])(OCC[NH3+])(=O))'  # PE
df.at[('HG', 'PE', 'hg'), 'pull'] = 'N1x'  # PE
df.at[('HG', 'PG', 'hg'), 'Structure'] = 'O([P@@]([O-])(OC[C@@H](O)CO)(=O))' # PG
df.at[('HG', 'PG', 'hg'), 'pull'] = 'O5x' # PG
df.at[('HG', 'PS', 'hg'), 'Structure'] = 'O([P@@]([O-])(OC[C@@H]([NH3+])(C(=O)([O-])))(=O))' # PS
df.at[('HG', 'PS', 'hg'), 'pull'] = 'N1x' # PS
df.at[('HG', 'PP', 'hg'), 'Structure'] = 'O([P@@]([O-])(O[P@@]([O-])(O)(=O))(=O))' # PP
df.at[('HG', 'PP', 'hg'), 'pull'] = 'P2x' # PP
# Figure out smiles for PI 
df.at[('TG', 'DD', 'sn1'), 'Structure'] = 'C'*9
df.at[('TG', 'DD', 'sn2'), 'Structure'] = 'C'*9
df.at[('TG', 'DC', 'sn1'), 'Structure'] = 'C'*10
df.at[('TG', 'DC', 'sn2'), 'Structure'] = 'C'*10
df.at[('TG', 'DL', 'sn1'), 'Structure'] = 'C'*11
df.at[('TG', 'DL', 'sn2'), 'Structure'] = 'C'*11
df.at[('TG', 'DM', 'sn1'), 'Structure'] = 'C'*13
df.at[('TG', 'DM', 'sn2'), 'Structure'] = 'C'*13
df.at[('TG', 'DP', 'sn1'), 'Structure'] = 'C'*15
df.at[('TG', 'DP', 'sn2'), 'Structure'] = 'C'*15
df.at[('TG', 'DS', 'sn1'), 'Structure'] = 'C'*17
df.at[('TG', 'DS', 'sn2'), 'Structure'] = 'C'*17
df.at[('TG', 'PS', 'sn1'), 'Structure'] = 'C'*17
df.at[('TG', 'PS', 'sn2'), 'Structure'] = 'C'*15
df.at[('TG', 'PO', 'sn1'), 'Structure'] = 'C'*15
df.at[('TG', 'PO', 'sn2'), 'Structure'] = 'CCCCCCC\C=C/CCCCCCCC' #Notice that when adding double bonds, pay attention to how {sn1} and {sn2} lipids are written
df.at[('TG', 'PL', 'sn1'), 'Structure'] = 'C'* 15
df.at[('TG', 'PL', 'sn2'), 'Structure'] = 'CCCCCCC\C=C/C\C=C/CCCCC'
df.at[('TG', 'PE', 'sn1'), 'Structure'] = 'C' *15
df.at[('TG', 'PE', 'sn2'), 'Structure'] = 'CCCCCCCCCCC\C=C/CCCCCCCC'
df.at[('TG', 'SO', 'sn1'), 'Structure'] = 'C' * 17
df.at[('TG', 'SO', 'sn2'), 'Structure'] = 'CCCCCCC\C=/CCCCCCCC'
df.at[('TG', 'SL', 'sn1'), 'Structure'] = 'C' * 17
df.at[('TG', 'SL', 'sn2'), 'Structure'] = 'CCCCCCC\C=C/C\C=C/CCCCC'
df.at[('TG', 'DR', 'sn1'), 'Structure'] = 'CCCCCCC\C=C/CCCC' #if written like {SN1}C(=O)O...CCCC\C=C/CCCCCCC
df.at[('TG', 'DR', 'sn2'), 'Structure'] = 'CCCCCCC\C=C/CCCC'
df.at[('TG', 'DO', 'sn1'), 'Structure'] = 'CCCCCCC\C=C/CCCCCCCC' #if written like {SN1}C(=O)O...CCCCCCCC\C=C/CCCCCCC
df.at[('TG', 'DO', 'sn2'), 'Structure'] = 'CCCCCCC\C=C/CCCCCCCC'
#add sterols
df.at[('Sterol', 'cholesterol', 'st'), 'Structure'] = 'C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C'
df.at[('Sterol', 'ergosterol', 'st'), 'Structure'] = 'C[C@H](/C=C/[C@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3C2=CC=C4[C@@]3(CC[C@@H](C4)O)C)C'
df.at[('Sterol', 'beta-Sitosterol', 'st'), 'Structure'] = 'CC[C@H](CC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C)C(C)C'
df.at[('Sterol', 'stigmasterol', 'st'), 'Structure'] = 'CC[C@H](/C=C/[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C)C(C)C'
df.at[('Sterol', 'cortisol', 'st'), 'Structure'] = '[3H]C1[C@H]2[C@@H]3CC[C@@]([C@]3(C[C@@H]([C@@H]2[C@]4(C(C(C(=O)C=C4C1[3H])[3H])[3H])C)O)C)(C(=O)CO)O'
df.at[('Sterol', 'corticosterone', 'st'), 'Structure'] = 'C[C@]12CCC(=O)C=C1CC[C@@H]3[C@@H]2[C@H](C[C@]4([C@H]3CC[C@@H]4C(=O)CO)C)O'
df.at[('Sterol', 'aldosterone', 'st'), 'Structure'] = 'C[C@]12CCC(=O)C=C1CC[C@@H]3[C@@H]2[C@H](C[C@]4([C@H]3CC[C@@H]4C(=O)CO)C=O)O'
df.at[('Sterol', 'progesterone', 'st'), 'Structure'] = 'CC(=O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4=CC(=O)CC[C@]34C)C'
df.at[('Sterol', 'beta-estradiol', 'st'), 'Structure'] = 'C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O'
df.at[('Sterol', 'testosterone', 'st'), 'Structure'] = 'C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C'
df.at[('Sterol', 'estradiol', 'st'), 'Structure'] = 'C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O'

os.makedirs('Dictionary', exist_ok=True)
df.to_csv("Dictionary/LipidLibrary.csv")