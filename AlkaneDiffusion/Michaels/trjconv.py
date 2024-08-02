import os
import numpy as np
import pdb

nwaters = 512
tlen = 1001

command = f"gmx trjconv -f nvt2_short_water_512.xtc -s water_512.gro -pbc nojump -o nvt2_short_water_512_nojump.xtc"
os.system(command)

if not os.path.exists('msds_test'):
    os.makedirs('msds_test')

if not os.path.exists('ndxs_test'):
    os.makedirs('ndxs_test') 

msds = np.zeros([nwaters,tlen])
for i in range(nwaters):

  # create a new index file
  f = open(f"ndxs_test/p{i}.ndx","w")
  f.write(f"[ p{i} ]\n")
  nO = 3*i+1
  f.write(f"{nO}")
  f.close()

  #run gmx msd with this index file
  command = f"echo 0 | gmx msd -f nvt2_short_water_512.xtc -s nvt2_water_512.tpr -o msds_test/msd{i}.xvg -n ndxs_test/p{i}.ndx"
  os.system(command) 

  # now read the index file
  f = open(f"msds_test/msd{i}.xvg")
  lines = f.readlines()
  f.close()
   
  itv = 0
  for l in lines:
    if l[0] !='#' and l[0] != '@':
      vals = l.split()
      msds[i,itv] = float(vals[1])
      itv+=1
np.save("manytraj_trjconv.npy",msds)