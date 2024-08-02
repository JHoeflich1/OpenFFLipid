import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.stats as stats

msd_data = np.load("manytraj_trjconv.npy")
nparticles, length = np.shape(msd_data)
for i in range(nparticles):
    plt.plot(msd_data[i,:],'b',alpha=0.02)

avemsd = np.mean(msd_data,axis=0)

plt.plot(avemsd)
plt.xlabel("frames")
plt.ylabel(r"$\langle MSD \rangle$ using trjconv -pbc nojump")
plt.savefig("allMSD_trjconv.png")
plt.clf()  

def func(a, mymsd): #5000 ps time, 5 ps timestep, 1001 frames 
    return a * np.linspace(0, 5000, 1001) - mymsd

slope = scipy.optimize.leastsq(func,0.1, args=avemsd)[0][0]
D = slope/6
print(D) 

finalD = (D / (10**7 * 10**7))*10**12
print(finalD)