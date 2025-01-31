# script is located in root (Heat_of_vap) and assumes you run in the Output_sims folder

import subprocess 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# ### check time series total energy is not drifting in gas phase
# molecules = ['pentane', 'hexane', 'heptane', 'octane', 'decane','pentadecane']
# for molecule in molecules:
#     command_V = f"echo Total-Energy | gmx energy -f ./{molecule}/nvt_{molecule}_gas.edr -o ./{molecule}/{molecule}_gas_U.xvg"
#     subprocess.run(command_V, shell=True, check=True)

#     #plot volume 
#     V_len = np.loadtxt(f'./{molecule}/{molecule}_gas_U.xvg', comments=['#','@'])
#     V_ps = V_len[:,0]
#     V_nm3 = V_len[:,1]

#     plt.figure()
#     plt.plot(V_ps, V_nm3, label=f"{molecule} ")
#     plt.xlabel('Time (ps)')
#     plt.ylabel('Total Energy (kJ/mol)')
#     plt.title(f'Total Energy vs Time for {molecule} ')
#     plt.legend()
#     plt.savefig(f"./{molecule}/{molecule}_total_E.png")
#     plt.close()

# ## check time series total energy is not drifting in liquid phase
# molecules = ['pentane', 'hexane', 'heptane', 'octane', 'decane','pentadecane']
# for molecule in molecules:
#     command_V = f"echo Total-Energy | gmx energy -f ./{molecule}/nvt_{molecule}_liquid.edr -o ./{molecule}/{molecule}_liquid_U.xvg"
#     subprocess.run(command_V, shell=True, check=True)

#     #plot volume 
#     V_len = np.loadtxt(f'./{molecule}/{molecule}_liquid_U.xvg', comments=['#','@'])
#     V_ps = V_len[:,0]
#     V_nm3 = V_len[:,1]

#     plt.figure()
#     plt.plot(V_ps, V_nm3, label=f"{molecule} ")
#     plt.xlabel('Time (ps)')
#     plt.ylabel('Total Energy (kJ/mol)')
#     plt.title(f'Total Energy vs Time for {molecule} ')
#     plt.legend()
#     plt.savefig(f"./{molecule}/{molecule}_total_E.png")
#     plt.close()

#### Calculate Heat of vaporization. Right now I am averaging over the last 2 ns of both gas and liquid simulations


molecules = ['pentane', 'hexane', 'heptane', 'octane', 'decane','pentadecane']

data = []
for mol in molecules:
    #Note that the total run time for the gas phase simulations was 30 ns. only use the last 1 ns. Notice how energy default only prints out every 1000 steps.
    command_gas_U = f"echo Total-Energy | gmx energy -f ./Gas_sims/{mol}/nvt_{mol}_gas.edr -o ./Gas_sims/{mol}/{mol}_gas_U_short.xvg -b 29000 -e 30000"
    subprocess.run(command_gas_U, shell=True, check=True)
    U_gas = np.loadtxt(f'./Gas_sims/{mol}/{mol}_gas_U_short.xvg', comments=['#','@'])
    Ug_ps = U_gas[:,0]
    Ug_kjm = U_gas[:,1]
    print(f'{mol} len(Ug_kjm)',len(Ug_kjm))

    #Note that the total run time for the liquid phase simulations was 2 ns. Only use the last 1 ns. Notice how energy default only prints out every 1000 steps.
    command_liquid_U = f"echo Total-Energy | gmx energy -f ./Liquid_sims/{mol}/nvt_{mol}_liquid.edr -o ./Liquid_sims/{mol}/{mol}_liquid_U_short.xvg -b 100 -e 200"
    subprocess.run(command_liquid_U, shell=True, check=True)
    U_liquid = np.loadtxt(f'./Liquid_sims/{mol}/{mol}_liquid_U_short.xvg', comments=['#','@'])
    Ul_ps = U_liquid[:,0]
    Ul_kjm = U_liquid[:,1]
    print(f'{mol} len(UL_kjm)',len(Ul_kjm))
    #Calcualte the average total energy in the gas and liquid sims, and divide by total moelcules in the simulation 
    Avg_U_gas = np.average(Ug_kjm) / 1
    Avg_U_liq = np.average(Ul_kjm) / 1000 
    #calc stdev 
    Std_U_gas = np.std(Ug_kjm, ddof=1) / 1
    Std_U_liq = np.std(Ul_kjm, ddof=1) / 1000

    N_gas = len(Ug_kjm)
    N_liq = len(Ul_kjm)

    SE_U_gas = Std_U_gas / np.sqrt(N_gas)
    SE_U_liq = Std_U_liq / np.sqrt(N_liq)


    # assuming V gas >> Vliq 
    # dHvap = Ugas - Uliq + RT
    R = 8.314/1000 #kJ/mol.K
    T = 300 #K 

    Hvap = Avg_U_gas - Avg_U_liq  + R*T

    #Uncertainity 
    Hvap_std = np.sqrt(Std_U_gas**2 + Std_U_liq**2)

    #  standard error
    Hvap_SE = np.sqrt(SE_U_gas**2 + SE_U_liq**2)
    
    data.append({
        "molecule": mol,
        "U_gas (kJ/mol)": Avg_U_gas,
        "U_liquid (kJ/mol)": Avg_U_liq,
        "H_vap (kJ/mol)": Hvap,
        "H_vap_std (kJ/mol)": Hvap_std,
        "H_vap_SE (kJ/mol)": Hvap_SE
    })

df = pd.DataFrame(data)

csv = "Hvap_totals.csv"
df.to_csv(csv,index=False)


