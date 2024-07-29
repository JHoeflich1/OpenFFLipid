import os
import numpy as np
import pdb

nwaters = 1024
tlen = 7001

msds = np.zeros([nwaters,tlen])
for i in range(nwaters):

  # create a new index file
  f = open(f"ndxs_shirts/p{i}.ndx","w")
  f.write(f"[ p{i} ]\n")
  nO = 3*i+1
  f.write(f"{nO}")
  f.close()

  #run gmx msd with this index file
  command = f"echo 0 | gmx msd -f nvt2_water_1024.xtc -s nvt2_water_1024.tpr -o msds_shirts/msd{i}.xvg -n ndxs_shirts/p{i}.ndx"
  os.system(command) 

  # now read the index file
  f = open(f"msds_shirts/msd{i}.xvg")
  lines = f.readlines()
  f.close()
   
  itv = 0
  for l in lines:
    if l[0] !='#' and l[0] != '@':
      vals = l.split()
      msds[i,itv] = float(vals[1])
      itv+=1
np.save("manytraj.npy",msds)