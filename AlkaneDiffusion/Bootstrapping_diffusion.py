# Bootstrapping for computing self-diffusion

# in 3D, you can esimate the self diffusion coefficient with the following formula:
# <(x(t) - x(t + tau)**2)> = 6*D*tau   where 6 is 2*n where n is the number of dimensions (3 in this case)

# 1. calculate the mean square displacement of a particle as a function of time 
# 2. average all the square displacements 
# 3. Fit the result to a line with intercept = 0 and find the slope
# 4. Divide the slope by 6 


import numpy as np
import matplotlib.pyplot as plt
import scipy

molecules = {6: 'hexane', 7: 'heptane'}
num_molecules = [512, 1024]


msd_data = np.load("msd_data.npz")

# Reconstruct the dictionary
msd_dict_loaded = {key: msd_data[key] for key in msd_data}

# Accessing a specific array
print(msd_dict_loaded['hexane_512'].shape)

######################################################

nparticles, length = np.shape(msd_data)

for i in range(nparticles):
    plt.plot(msd_data[i,:], 'b', alpha = 0.02)

avemsd = np.mean(msd_data, axis = 0) #average overa all particles

plt.plot(avemsd)
plt.xlabel("frames")
plt.ylabel(r"$\langle MSD \rangle$")
plt.show()

# we want to minimize sum over all i of (a*tau_i - y(tau_i))**2. A linear fit, but with no intercept

#frames are from 0 to 1000, 0.5 ps each 
#this is the function we want to minimize to:
def func(a,mymsd):
    return a*np.linspace(0,1000,2001)-mymsd

slope = scipy.optimize.leastsq(func,0.1,args=avemsd)[0][0]

D = slope/6  # diffusion coefficient is the slope of the msd
print(D) #this is in nm**2/ps

#to change to cm**2 /s 
finalD = (D / (10**7 * 10**7)) * 10**12


## Now we want to perform a bootstrap error analysis by constructing bootstrap samples over particles, and then use to compute the average MSD
# 1. Plot 500 bootstrapped MSDs
# 2. Show the 95% confidence intercal region of MSD trajectories
# 3. Plot the distributions of MSDs. Is it gaussian, if so, what is the standard error in the calculation 

nbootstraps = 500
newmsd_data = np.zeros(np.shape(msd_data)) #create an empty matrix of the same size 
Ds = np.zeros(nbootstraps)
for n in range(nbootstraps):
    for i in range(nparticles):
        #generate bootstrapped data:
        newi = np.rangom.randint(0,nparticles)
        newmsd_data[i,:] = msd_data[newi,:]
    new_avemsd = np.mean(newmsd_data,axis = 0)
    plt.plot(avemsd, 'b', alpha = 0.1)
    slope = scipy.optimize.leastsq(func, 0.1, new_avemsd)[0][0]
    Ds[n] = slope/6 * 0.01 #10^7 * 10^7 /10^14
plt.show()
plt.hist(Ds,bins=50)
#does this look gaussian
stderr_boots = np.std(Ds)


#Now we want to git it to a curve using linear regression and determine D inf