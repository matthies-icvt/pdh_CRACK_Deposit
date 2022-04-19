import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

def plot_func(result_dict_gas, result_dict_surf, inactiv, Depo_plots, n_steps, z):

    list_cycle = np.arange(inactiv)
    fig9, (ax20) = plt.subplots(figsize=(8, 5))
    # fig9.suptitle('Moll fraction of compounds at the outlet of reactor')
    ax20.plot(list_cycle, result_dict_gas["end_state"]["C3H8"], label='C3H8')
    ax20.plot(list_cycle, result_dict_gas["end_state"]["C3H6"], label='C3H6')
    ax20.plot(list_cycle, result_dict_gas["end_state"]["CH4"], label='CH4')
    ax20.plot(list_cycle, result_dict_gas["end_state"]["H2"], label='H2')
    ax20.set_xlabel('Cycles')
    ax20.set_ylabel('Moll fraction [mol/mol]')
    ax20.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax20.set_ylim(top=0.6)
    fig9.legend()
    fig9.savefig('Plots/mol_frac.png', dpi=400)

    fig10, ax23 = plt.subplots(figsize=(8, 5))
    # fig10.suptitle('Deposition along the reactor in cycles')
    # Depo_plots = 5  # How many plots should be in the diagramm
    for i in range(Depo_plots):
        pltNR = i * (inactiv / Depo_plots)
        deposition = np.multiply(result_dict_surf["mass_depo"][pltNR], 10e6)
        deposition = deposition.tolist()
        ax23.plot(z[1:], deposition, label='Cycle ' + str(int(pltNR)))
    ax23.set_ylabel('Deposition [mg/s]')
    ax23.set_xlabel('Reactor length [m]')
    ax23.autoscale()
    fig10.legend()
    fig10.savefig('Plots/depo.png', dpi=400)

    fig11, ax24 = plt.subplots(figsize=(8, 5))
    # fig11.suptitle('Surface area per Volume')
    result_dict_surf['cat_area_per_vol'][0].pop(0)
    for i in range(Depo_plots):
        pltNR = pltNR = i * (inactiv / Depo_plots)
        ax24.plot(z[1:], result_dict_surf['cat_area_per_vol'][pltNR], label='Cycle ' + str(int(pltNR)))
    ax24.set_ylabel('Specific Surface [1/m]')
    ax24.set_xlabel('Reactor length [m]')
    fig11.legend()
    fig11.savefig('Plots/surf_vol.png', dpi=400)

    fig12, ax25 = plt.subplots(figsize=(8, 5))
    # fig12.suptitle('Selective surface area per Volume')
    result_dict_surf['cat_area_per_vol2'][0].pop(0)
    z1 = z[1:]
    slice = int(n_steps / 4 - 1)
    for i in range(Depo_plots):
        pltNR = i * (inactiv / Depo_plots)
        ax25.plot(z1[slice:], result_dict_surf['cat_area_per_vol2'][2 * i][slice:], label='Cycle ' + str(2 * i))
    ax25.set_ylabel('Specific Surface [1/m]')
    ax25.set_xlabel('Reactor length [m]')
    fig12.legend()
    fig12.savefig('Plots/surf2_vol.png', dpi=400)

    fig13, axs = plt.subplots(2, 2, sharex=True, figsize=(8, 5))
    # fig13.suptitle('Moll fraction of compounds along the reactor')
    for i in range(Depo_plots):
        pltNR = pltNR = i * (inactiv / Depo_plots)
        axs[0, 0].plot(z[1:], result_dict_gas['conc'][pltNR]['C3H8'], label='Cycle ' + str(int(pltNR)))
        axs[0, 0].set_title('C3H8')
        axs[0, 1].plot(z[1:], result_dict_gas['conc'][pltNR]['C3H6'])
        axs[0, 1].set_title('C3H6')
        axs[1, 0].plot(z[1:], result_dict_gas['conc'][pltNR]['CH4'])
        axs[1, 0].set_title('CH4')
        axs[1, 1].plot(z[1:], result_dict_gas['conc'][pltNR]['H2'])
        axs[1, 1].set_title('H2')
    fig13.text(0.5, 0.04, "Reactor length [m]", ha='center', va='center')
    fig13.text(0.04, 0.5, "Concentration [mol/m^3]", ha='center', va='center', rotation='vertical')
    fig13.legend()
    fig13.savefig('Plots/conc.png', dpi=400)