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
import multiprocessing 

def BootstrappingMSD(args):
    '''Input: molecule and size
    returns a tuple with final diffusion coefficient and standard error of the diffusion coefficient'''  
    molecule, size = args

    # Extract MSD data in a matrix [particles, tlen]
    msd_df = pd.read_pickle("msd_data_new.pkl")

    matrix = msd_df.loc[(molecule, size)]
    matrix = matrix.to_numpy()  # or matrix.values
    nparticles, length = np.shape(matrix)
    # print(np.shape(matrix))
    # Plot all MSDs
    for i in range(nparticles): 
        plt.plot(matrix[i, :], 'b', alpha=0.02)
    plt.xlabel('frames')
    plt.ylabel('MSD')
    plt.title(f"All MSD for {molecule} {size}")
    plt.savefig(f"allMSD_{molecule}_{size}.png")
    plt.clf()  # Clear plot

    avemsd = np.mean(matrix, axis=0)  # Mean computed across columns

    # Plot average MSD
    plt.plot(avemsd)
    plt.xlabel("frames")
    plt.ylabel(r"$\langle MSD \rangle$")
    plt.title(f"Average MSD for {molecule}_{size}")
    # plt.savefig(f"AvgMSD_{molecule}_{size}.png")
    plt.clf()  # Clear plot

    def func(a, mymsd): #2000 ps time, 1 ps timestep, 2001 frames 
        return a * np.linspace(0, 2000, 2001) - mymsd

    slope = scipy.optimize.leastsq(func, 0.1, args=avemsd)[0][0]

    D = slope / 6  # Diffusion coefficient is the slope of the MSD
    finalD = (D / (10**7 * 10**7)) * 10**12  # Convert to cm^2/s

    # Perform bootstrap error analysis
    nbootstraps = 5000
    newmsd_data = np.zeros(np.shape(matrix))  # Create an empty matrix of the same size
    Ds = np.zeros(nbootstraps)

    for n in range(nbootstraps):
        for i in range(nparticles):
            newi = np.random.randint(0, nparticles)
            newmsd_data[i, :] = matrix[newi, :]
        new_avemsd = np.mean(newmsd_data, axis=0)
        slope = scipy.optimize.leastsq(func, 0.1, args=new_avemsd)[0][0]
        Ds[n] = slope / 6 * 0.01  # Convert to cm^2/s
    
    # Plot bootstrapped MSDs
    plt.plot(avemsd, 'b', alpha=0.1)
    plt.title(f"Bootstrapped MSDs for {molecule} {size}")
    # plt.savefig(f"bootstrapped_data_{molecule}_{size}.png")
    plt.clf()

    # Plot histogram of diffusion coefficients
    plt.hist(Ds, bins=50)
    plt.title(f"Histogram of Diffusion Coefficients for {molecule} {size}")
    # plt.savefig(f"Ds_histogram_{molecule}_{size}.png")
    plt.clf()

    stderr_boots = np.std(Ds)

    print(f"Diffusion for {molecule} {size} is {finalD:.4g} +/- {stderr_boots:.4g} cm^2/s")

    return (molecule, size), finalD, stderr_boots


if __name__ == '__main__':
    molecules = ['water']
    sizes = [512, 1024, 2048, 4096]

    index = pd.MultiIndex.from_product([molecules, sizes], names=['Molecule', 'Size'])
    DS_final = pd.DataFrame(index=index, columns=['Ds', 'Stderr'])
    DS_final[:] = np.nan

    # Create a pool of worker processes
    pool = multiprocessing.Pool(multiprocessing.cpu_count())

    results = []

    for molecule in molecules:
        for size in sizes:
            results.append(pool.apply_async(BootstrappingMSD, args=((molecule, size),)))

    pool.close()
    pool.join()

    # Collect results
    for result in results:
        (molecule, size), finalD, stderr_boots = result.get()
        DS_final.loc[(molecule, size), 'Ds'] = finalD
        DS_final.loc[(molecule, size), 'Stderr'] = stderr_boots

    # Save the final DataFrame to a CSV file
    DS_final.to_csv('diffusion_coefficients.csv')
    DS_final.to_pickle('diffusion_coefficients.pkl')