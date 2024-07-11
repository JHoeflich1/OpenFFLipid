import os
import numpy as np
import pdb

num_molecules = 900 #number of molecules in sim 
tlen = 2001 #trajectory length
n = 3 #number of atoms in your molecule

msds = np.zeros([num_molecules,tlen])
for i in range(num_molecules):

  # create a new index file
  f = open(f"ndxs/p{i}.ndx","w")
  f.write(f"[ p{i} ]\n")
  nO = n*i+1
  f.write(f"{nO}")
  f.close()

  #run gmx msd with this index file
  command = f"{gmx} msd -f justwater.xtc -s justwater.tpr -o msds/msd{i}.xvg -n ndxs/p{i}.ndx"
  os.system(command) 

  # now read the index file
  f = open(f"msds/msd{i}.xvg")
  lines = f.readlines()
  f.close()
   
  itv = 0
  for l in lines:
    if l[0] !='#' and l[0] != '@':
      vals = l.split()
      msds[i,itv] = float(vals[1])
      itv+=1
np.save("manytraj.npy",msds)