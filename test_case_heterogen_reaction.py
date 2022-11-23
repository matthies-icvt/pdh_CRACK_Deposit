import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import time
import Reduction_mech
from itertools import islice


# anode = ct.Solution(input_file, "anode")
# cathode = ct.Solution(input_file, "cathode")
# metal = ct.Solution(input_file, "electron")
# electrolyte = ct.Solution(input_file, "electrolyte")
# anode_int = ct.Interface(
#     input_file, "edge_anode_electrolyte", adjacent=[anode, metal, electrolyte])
# cathode_int = ct.Interface(
#     input_file, "edge_cathode_electrolyte", adjacent=[cathode, metal, electrolyte])

reaction_mech = 'gri30.yaml'
reaction_mech_surf = 'reaction_kinetics/Propan_surface_and_gas.yaml'

area_cat = 1 # [m^2] 
print('1!')
gas = ct.Solution(reaction_mech_surf, 'gas')
#pt_cat = ct.Interface(reaction_mech_surf, "pt_surf")
print('2!')
# Import gas model for surface reaction
#gas_surf = ct.Solution(reaction_mech_surf, 'pt_surf')
print('3!')
pt_cat = ct.Interface(reaction_mech_surf, 'pt_surf', adjacent=[gas])
#pt_cat = ct.Interface(reaction_mech_surf, "edge_pt_surf", adjacent=[gas, gas_surf])
print('4!')
#########
print('read data!')
# relative temperature [K]
T_0_n = 293.0
# inlet temperature [K]
T_0 = 800 + 273.15
# wall temperature [K]
T_wall = 610 + 273.15
# constant pressure [Pa]
pressure = ct.one_atm
# composition
composition_0 = 'C3H8:10, H2:1'
# composition_0 = 'CH4:1, O2:1.5, AR:0.1'
# definition of initial and reactive state
initial_state = T_0, pressure, composition_0

gas.TPX = initial_state

pt_cat.TP = T_0, pressure
cov = pt_cat.coverages
print('5!')
gas_r = ct.IdealGasReactor(contents=gas, energy = 'on')
print('6!')

rsurf = ct.ReactorSurface(pt_cat, gas_r, A = area_cat)
print('init reactor!')

sim = ct.ReactorNet([gas_r])

sim.rtol = 1.0e-9
sim.atol = 1.0e-21

print('starting sim!')
sim.advance_to_steady_state()

print('finished sim!')
print('{0:10f}  {1:10f}  {2:10f} {3:10f}'.format(*gas['C3H8', 'H2', 'C3H6', 'CH4'].X))