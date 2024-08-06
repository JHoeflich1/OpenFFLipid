import numpy as np
import pandas as pd
import os
import multiprocessing
import matplotlib.pyplot as plt


def calculate_msd(params):
    '''
    Calculates the mean square displacement for a molecule
    inputs:
        params : a tuple containing (alkane, alkane_atoms, alkane_size)
    '''
    alkane, alkane_atoms, alkane_size = params


    ## check pressure profile
    # command_P = f"echo Pressure | gmx energy -f nvt2_{alkane}_{alkane_size}.edr -o {alkane}_{alkane_size}_P.xvg"
    # os.system(command_P)

    # #plot volume 
    # P_len = np.loadtxt(f'{alkane}_{alkane_size}_P.xvg', comments=['#','@'])
    # P_ps = P_len[:,0]
    # P_nm3 = P_len[:,1]

    # plt.figure()
    # plt.plot(P_ps, P_nm3, label=f"{alkane} {alkane_size}")
    # plt.xlabel('Time (ps)')
    # plt.ylabel('Pressure')
    # plt.title(f'Pressure vs Time for {alkane} {alkane_size}')
    # plt.legend()
    # plt.savefig(f"{alkane}_{alkane_size}_pressure.png")
    # plt.close()

    for i in range(alkane_size):
        # Skip calculations that have already done
        msd_file = f"msds/msd_{alkane}_{alkane_size}_{i}.xvg"
        if os.path.exists(msd_file):
            print(f"{msd_file} already exists, skipping this calculation.")
            continue
        with open(f"ndxs/ndxs_{alkane}_{alkane_size}_{i}.ndx", "w") as f:
            f.write(f"[ p{i} ]\n")
            n0 = alkane_atoms * i + 1
            f.write(f"{n0}\n\n")
            #print(f'ndxs_test/ndxs_{alkane}_{alkane_size}_{i}.ndx is written')

        # Run gmx msd with each index file
        #print('moving on to gmx')
        command = f"echo 0 | gmx msd -f nvt2_cut_{alkane}_{alkane_size}.xtc -s nvt2_{alkane}_{alkane_size}.tpr -o msds/msd_{alkane}_{alkane_size}_{i}.xvg -n ndxs/ndxs_{alkane}_{alkane_size}_{i}.ndx"# -rmpbc -pbc" #-selrpos whole_mol_com
        os.system(command)
        # -selrpos whole_mol_com calculates the whole moelcule COM even if only a part of it is selected
        # -rmpbc means that molecules are made whole for each frame
        # -pbc means to use periodic boundary conditions for distance calculation


if __name__ == '__main__':
    molecules =  {17: 'pentane', 20: 'hexane', 23: 'heptane', 26: 'octane', 32: 'decane', 47: 'pentadecane'}
    sizes = [512, 1024, 2048]

    if not os.path.exists('msds'):
        os.makedirs('msds')

    if not os.path.exists('ndxs'):
        os.makedirs('ndxs')

    nprocs = multiprocessing.cpu_count()
    print(f'Running on {nprocs} CPUs:')

    # Create a pool of workers
    pool = multiprocessing.Pool(processes=nprocs)

    # Generate the list of tasks, each task is a tuple with all required parameters
    tasks = [(name, key, size) for key, name in molecules.items() for size in sizes]

    # Map the tasks to the pool
    pool.map(calculate_msd, tasks)
    
    pool.close()
    pool.join()
