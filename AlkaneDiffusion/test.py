import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd

molecules = ['hexane', 'heptane', 'octane', 'decane', 'pentadecane']
sizes = [512, 1024, 2048, 4096]

# Create a multiindex df from product of molecules and sizes
index = pd.MultiIndex.from_product([molecules, sizes], names=['Molecule', 'Size'])

# Create an empty df with the specified index
data = {'Ds': 0, 'Stderr': 0}  
df = pd.DataFrame(data, index=index)

print(df)