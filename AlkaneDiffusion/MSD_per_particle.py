import os
import numpy as np
import pdb
import subprocess

molecules = {6:'hexane', 7:'heptane'}
num_molecules = [512, 1024]
tlen = 2001 #trajectory length

msd_dict = {}

for key, mol in molecules.items():
    msd_dict[mol] = {}
    for num in num_molecules:
        msd_dict[mol][num] = np.zeros([num, tlen])

for num in num_molecules:
  for key, mol in molecules.items():
    for i in range(num):
      msds = np.zeros([num,tlen])

      f = open(f"ndxs_{mol}_{num}_{i}.ndx","w")
      f.write(f"[ p{i} ]\n")
      n0 = key * i + 1
      f.write(f"{n0}")
      f.close()


      #run gmx msd with each index file
      command = f"gmx msd -f nvt_{mol}_{num}.xtc -s nvt_{mol}_{num}.tpr -o msds/msd_{mol}_{num}_{i}.xvg -n ndxs_{mol}_{num}_{i}.ndx")
      os.system(command)

      # now read the index file
      f = open(f"msds/msd_{mol}_{num}_{i}.xvg")
      lines = f.readlines()
      f.close()
      
      itv = 0
      for l in lines:
        if l[0] !='#' and l[0] != '@':
          vals = l.split()
          msd_dict[mol][num][i,itv] = float(vals[1])
          itv+=1
np.savez("msd_data.npz", **msd_dict)