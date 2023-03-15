# heterogen_reaction_1D.py>

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import time
import Reduction_mech
from itertools import islice
from function_surface_gas_reaction import *

#-------------------- INPUT DATA FROM TERMINAL --------------------#

# n_steps = int(input('Number of reactors in the pfr system. The number must be divisible by 4!\n'))
# while n_steps % 4 != 0:
#     print("The number must be divisible by 4!\n")
#     n_steps = int(input("Number of reactors:\n"))
# iters = int(input('Number of iterations for the pfr-reactor cylcle?\n'))
# mechanism = int(input('Press 0 for automatic gri_30, 1 to choose reduced gri30 or Press 2 to choose mech_13\n'))
# inactiv = int(input('Should surface be decreased due to deposition? 1- yes'))
mechanism = 2
remove = 1
inactiv = 10 # inactiv = 1 no deactivation, inactiv > 1 how many cycles should reactor work with decreasing surface every cycle. %10==0!
Depo_plots = 5
# starts counting time in which program runs
start_time = time.time()

#-------------------- PROCESS PARAMETERS --------------------#
# R
R = 8.314 #J/mol/K
# composition
composition_0 = 'C3H8:10, H2:1'
# relative temperature [K]
T_0_n = 293.0
# inlet temperature [K]
T_0 = 550 + 273.15
# wall temperature [K]
T_wall = 610 + 273.15
# constant pressure [Pa]
pressure = ct.one_atm
# flow rate [m3/s] volumen flow standarized for 273K
vol_flow_rate_N = 55*1e-3*1e-3 #9.166667e-6
# Mol Weight
gas = ct.Solution('reaction_kinetics/Propan_surface.yaml')
initial_state = T_0, pressure, composition_0
gas.TPX = initial_state
MW_0 = gas.mean_molecular_weight/1000 #kg/mol
print(MW_0)
# flow mass [kg/s]
m_dot_0 = vol_flow_rate_N * 10**5 / 8.314 / 273.15 * MW_0
# composition_0 = 'CH4:1, O2:1.5, AR:0.1'
# definition of initial and reactive state
initial_state = T_0, pressure, composition_0
reactive_state = T_wall, pressure, composition_0

#-------------------- REACTOR GEOMETRY --------------------#
length = 40*10**-3 + 60*10**-3 #0.0815  # Kat length + rest of reactor *approximate* PFR length [m]
height = 0.006  # [m]
depth = 0.038  # [m]
area = height*depth  # cross-sectional area of reactor [m**2]
porosity = 0.6 # [-] It's not a porosity per se, it defines what fraction of reactor will be surface reactor.
# Zeolite has a surface of 400-425 m2 per g. Density is ca. 720 g/litre (Datenblatt zsm5 pellets https://www.acsmaterial.com/zsm-5-catalyst.html) -> surface per vol: 405*720*10^3*m2/m3
zsm5_area_per_vol = 291.6*10**6   # Catalyst particle surface area per unit volume [1/m] 
surf_pt = 0.559 #m2 est. surface of inner pt surface
reactor_vol = length * height * depth
cat_area_per_vol = surf_pt / reactor_vol
area_react = area * (1-porosity)
area_cat = area * porosity
K = 10 # Heat transfer coeff. between wall and a reactor

#-------------------- SIM PARAMETERS --------------------#
dt_0 = 0.000001 #s
n_steps = 100
n_reactor_vol = reactor_vol / n_steps
n_area_cat = area_cat / n_steps
reactors = []
for i in range(n_steps):
    reactors.append({'m': 0, 'T': T_0, 'c': composition_0, 'p': 10**5, 'm_dot_in': m_dot_0})

#-------------------- DEPOSITION MODEL PARAMETERS --------------------#
# Surface deactivation coefficient - is used to artificially enlarge the influence of deposition on a surface
alpha1 = 1e5
alpha2 = 2e5

# reaction mechanism for surface reaction
reaction_mech_surf = 'reaction_kinetics/Propan_surface.yaml'
reaction_mech_surf2 = 'reaction_kinetics/Propan_surface2.yaml'

M_depo = 80.0 # Molar mass from which the deposition starts



composition_gas_0 = {'H2': 0.05 , 'C3H8':  0.95}
cov_pt_0 = {'PT(S)': 1., 'H(S)': 1.}

composition_gas = composition_gas_0
mol_gas = {}
cov_pt = cov_pt_0
mol_surf = {}

results = []

for n, reaction_volume in enumerate(reactors):

    #calculation balance
    dt = dt_0

    net_prod = surf_reactions(composition_gas, cov_pt, dt)
    
    #print(net_prod)

    gas_mol_net_production = net_prod['gas'] # in mol/m3/s
    
    cov_mol_net_production = net_prod['cov'] # in mol/m2/s

    sum_mol_gas_n_0 = 5*10**5 * n_reactor_vol / T_0 / R
    sum_mol_gas_n = 0

    sum_cov_n_0 = 5.75e19 / n_steps # mol/m2 // side density 0.01 mol/cm^2 
    sum_mol_surf_n = 0 

    # calculation gas
    for species in gas_mol_net_production:
        species_mol_change = n_reactor_vol * dt * gas_mol_net_production[species] # [mol]

        if species in composition_gas:
            mol_gas[species] = composition_gas[species] * sum_mol_gas_n_0 + species_mol_change
        else:
            mol_gas[species] = species_mol_change

        if mol_gas[species] < 0:
            mol_gas[species] = 0

        sum_mol_gas_n += mol_gas[species]

    for species in mol_gas:
        composition_gas[species] = mol_gas[species] / sum_mol_gas_n

    # calculation coverages
    for species in cov_mol_net_production:
        species_mol_change = n_area_cat * dt * cov_mol_net_production[species] # [mol]
        if species in cov_pt:
            mol_surf[species] = cov_pt[species] * sum_cov_n_0 + species_mol_change
        else:
            mol_surf[species] = species_mol_change

        if mol_surf[species] < 0:
            mol_surf[species] = 0

        sum_mol_surf_n += mol_surf[species]

    for species in cov_pt:
        cov_pt[species] = mol_surf[species] / sum_mol_surf_n

    #print('\n composition next step: \n')

    #storage of results in lists of each reactor
    print("loop nr = {}".format(n))
    print(composition_gas)
    print(cov_pt)
    results.append(composition_gas)

print('####')


