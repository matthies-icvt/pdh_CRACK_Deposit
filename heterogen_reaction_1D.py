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
n_steps = 200
mechanism = 2
remove = 1
inactiv = 10 # inactiv = 1 no deactivation, inactiv > 1 how many cycles should reactor work with decreasing surface every cycle. %10==0!
Depo_plots = 5
# starts counting time in which program runs
start_time = time.time()

#-------------------- PROCESS PARAMETERS --------------------#
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
# composition
composition_0 = 'C3H8:10, H2:1'
# composition_0 = 'CH4:1, O2:1.5, AR:0.1'
# definition of initial and reactive state
initial_state = T_0, pressure, composition_0
reactive_state = T_wall, pressure, composition_0

#-------------------- DEPOSITION MODEL PARAMETERS --------------------#
# Surface deactivation coefficient - is used to artificially enlarge the influence of deposition on a surface
alpha1 = 1e5
alpha2 = 2e5

# reaction mechanism for surface reaction
reaction_mech_surf = 'reaction_kinetics/Propan_surface.yaml'
reaction_mech_surf2 = 'reaction_kinetics/Propan_surface2.yaml'

M_depo = 80.0 # Molar mass from which the deposition starts

#-------------------- REACTOR GEOMETRY --------------------#
length = 40*10**-3 + 60*10**-3 #0.0815  # Kat length + rest of reactor *approximate* PFR length [m]
height = 0.006  # [m]
depth = 0.038  # [m]
area = height*depth  # cross-sectional area of reactor [m**2]
porosity = 0.6 # [-] It's not a porosity per se, it defines what fraction of reactor will be surface reactor.
# Zeolite has a surface of 400-425 m2 per g. Density is ca. 720 g/litre (Datenblatt zsm5 pellets https://www.acsmaterial.com/zsm-5-catalyst.html) -> surface per vol: 405*720*10^3*m2/m3
zsm5_area_per_vol = 291.6*10**6   # Catalyst particle surface area per unit volume [1/m] 
surf_pt = 0.559 #m2 est. surface of inner pt surface
cat_area_per_vol = surf_pt / (length * height * depth)
area_react = area * (1-porosity)
area_cat = area * porosity
beta = 1 # How much from the main flow (gas reactor) will flow to the surface reactor, 1 = whole mass stream
K = 10 # Heat transfer coeff. between wall and a reactor

surf_reactions()

