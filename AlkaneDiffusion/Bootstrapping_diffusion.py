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
import pandas as pd

#load in dataframe from MSD_per_particle.py
msd_df = pd.read_pickle("msd_data.pkl")
# print(msd_df.head())

##### make an empty dataframe to store calcualted Ds and Sterr values 
molecules = ['hexane', 'heptane', 'octane', 'decane', 'pentadecane']
sizes = [512, 1024, 2048, 4096]

index = pd.MultiIndex.from_product([molecules, sizes], names=['Molecule', 'Size'])
data = {'Ds': 0, 'Stderr': 0}  
DS_final = pd.DataFrame(index=index, columns=['Ds','Stderr'])
DS_final[:] = np.nan

for molecule in molecules:
    for size in sizes:
        matrix = msd_df.loc[(molecule, size)]
        # print(f'Molecule: {molecule}, Size: {size}')
        # print(matrix.head())


        nparticles, length = np.shape(matrix)

        for i in range(nparticles): 
            plt.plot(matrix[i,:], 'b', alpha = 0.02)
        plt.xlabel('frames')
        plt.ylabel('MSD')
        plt.title(f"All MSD for {molecule} {size}")
        plt.savefig(f"allMSD_{molecule}_{size}")

        avemsd = np.mean(matrix, axis = 0) # mean computed across columns (aka mean across all particles per time point)

        plt.plot(avemsd)
        plt.xlabel("frames")
        plt.ylabel(r"$\langle MSD \rangle$")
        plt.title(f"average MSD for {molecule}_{size}")
        plt.savefig(f"AvgMSD_{molecule}_{size}")
        # plt.show()

        # we want to minimize sum over all i of (a*tau_i - y(tau_i))**2. A linear fit, but with no intercept

        #frames are from 0 to 1000, 0.5 ps each 
        #this is the function we want to minimize to:
        def func(a,mymsd):
            return a*np.linspace(0,1000,2001)-mymsd

        slope = scipy.optimize.leastsq(func,0.1,args=avemsd)[0][0]

        D = slope/6  # diffusion coefficient is the slope of the msd
        # print(D) #this is in nm**2/ps

        #to change to cm**2 /s 
        finalD = (D / (10**7 * 10**7)) * 10**12
        print(f"final diffusion for {molecule} {size} is {finalD} ",r"\frac{cm^{2}}{s}")


        ## Now we want to perform a bootstrap error analysis by constructing bootstrap samples over particles, and then use to compute the average MSD
        # 1. Plot 500 bootstrapped MSDs
        # 2. Show the 95% confidence intercal region of MSD trajectories
        # 3. Plot the distributions of MSDs. Is it gaussian, if so, what is the standard error in the calculation 

        nbootstraps = 500
        newmsd_data = np.zeros(np.shape(matrix)) #create an empty matrix of the same size 
        Ds = np.zeros(nbootstraps)
        for n in range(nbootstraps):
            for i in range(nparticles):
                #generate bootstrapped data:
                newi = np.random.randint(0,nparticles)
                newmsd_data[i,:] = matrix[newi,:]
            new_avemsd = np.mean(newmsd_data,axis = 0)
            plt.plot(avemsd, 'b', alpha = 0.1)
            slope = scipy.optimize.leastsq(func, 0.1, new_avemsd)[0][0]
            Ds[n] = slope/6 * 0.01 #10^7 * 10^7 /10^14 cm^2/s
        plt.show()
        plt.savefig(f"bootstrapped data for {molecule} {size}")

        plt.hist(Ds,bins=50)
        plt.save(f"Ds histogram for {molecule} {size}")
        #does this look gaussian
        stderr_boots = np.std(Ds)

        print(f"Diffusion for {molecule} {size} is {Ds} +/- {stderr_boots} cm^2/s")

        DS_final[(molecule, size), 'Ds'] = Ds
        DS_final[(molecule, size), 'Stderr'] = stderr_boots

    print(DS_final)
    DS_final.to_pickle("DS_final.pkl")
        #Now we want to git it to a curve using linear regression and determine D inf