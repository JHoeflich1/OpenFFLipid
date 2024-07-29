
import numpy as np
import matplotlib.pyplot as plt
import scipy
import pandas as pd
import itertools
import multiprocessing


def BootstrappingMSD(args):
    '''Accepts a moleule and a size. It looks for the .csv file that contains msd data for each particle.
    Diffusion coefficient is calculated, and error bars are calcualted by bootstrap analysis
    returns a dataframe'''  
    molecule, size = args

    #extract msd data in a matrix [particles, tlen]
    data = pd.read_csv(f'msd_data/msd_data_{molecule}_{size}.csv') #add folder in her msd_data/msd_data_...
    msd_columns = [str(i) for i in range(701)]  # Column names from 0 to 700
    msd_data = data[msd_columns]
    # create np array
    msd_matrix = msd_data.values
    nparticles, length = np.shape(msd_matrix)


    for i in range(nparticles): 
        plt.plot(msd_matrix[i,:], 'b', alpha = 0.02)
    plt.xlabel('frames')
    plt.ylabel('MSD')
    plt.title(f"All MSD for {molecule} {size}")
    plt.savefig(f"allMSD_{molecule}_{size}")
    plt.clf() #clear plot

    avemsd = np.mean(msd_matrix, axis = 0) # mean computed across columns (aka mean across all particles per time point)
    # print(avemsd, 'avemsd')

    plt.plot(avemsd)
    plt.xlabel("frames")
    plt.ylabel(r"$\langle MSD \rangle$")
    plt.title(f"average MSD for {molecule}_{size}")
    # plt.savefig(f"AvgMSD_{molecule}_{size}")
    plt.clf() #clear plot

    # we want to minimize sum over all i of (a*tau_i - y(tau_i))**2. A linear fit, but with no intercept

    #frames are from 0 to 7000, 10 ps each 
    #this is the function we want to minimize to:
    def func(a,mymsd):
        return a * np.linspace(0, 7000, 701) - mymsd

    slope = scipy.optimize.leastsq(func,0.1,args=avemsd)[0][0]

    D = slope/6  # diffusion coefficient is the slope of the msd
    # print('D in nm**2/ps',D) #this is in nm**2/ps

    #to change to cm**2 /s 
    finalD = (D / (10**7 * 10**7)) * 10**12
    # print(f"final diffusion for {molecule} {size} is {finalD} ",r"\frac{cm^{2}}{s}")


    ## Now we want to perform a bootstrap error analysis by constructing bootstrap samples over particles, and then use to compute the average MSD
    # 1. Plot 500 bootstrapped MSDs
    # 2. Show the 95% confidence intercal region of MSD trajectories
    # 3. Plot the distributions of MSDs. Is it gaussian, if so, what is the standard error in the calculation 

    nbootstraps = 5000
    newmsd_data = np.zeros(np.shape(msd_matrix)) #create an empty matrix of the same size 
    Ds = np.zeros(nbootstraps)

    for n in range(nbootstraps):
        for i in range(nparticles):
            #generate bootstrapped data:
            newi = np.random.randint(0,nparticles)
            newmsd_data[i, :] = msd_matrix[newi, :]
        new_avemsd = np.mean(newmsd_data,axis = 0)
        plt.plot(avemsd, 'b', alpha = 0.1)
        slope = scipy.optimize.leastsq(func, 0.1, new_avemsd)[0][0]
        Ds[n] = slope/6 * 0.01 #10^7 * 10^7 /10^14 cm^2/s
    plt.show()
    # plt.savefig(f"bootstrapped data for {molecule} {size}")
    plt.clf()

    plt.hist(Ds,bins=50)
    plt.savefig(f"Ds histogram for {molecule} {size}")
    #does this look gaussian
    stderr_boots = np.std(Ds)

    print(f"Diffusion for {molecule} {size} is {finalD:.4g} +/- {stderr_boots:.4g} cm^2/s")

    return finalD, stderr_boots


if __name__ == '__main__':
    molecules = ['water']#,'pentane','hexane','heptane','octane','decane','pentadecane']
    sizes = [512, 1024]#, 2048]
    combinations = list(itertools.product(molecules, sizes))

    pool = multiprocessing.Pool()
    data = pool.map(BootstrappingMSD,combinations)
    summaryDF = pd.DataFrame(data, columns=['Diffusion', 'StdErr'], index=pd.MultiIndex.from_product([molecules, sizes], names=['Molecule', 'Size']))
    print(summaryDF)
    summaryDF.to_csv('summary_diffusion_data.csv')
    pool.close()
    pool.join()




