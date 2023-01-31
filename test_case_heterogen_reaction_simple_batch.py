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

reaction_mech = 'reaction_kinetics/mech_13.yaml'
reaction_mech_surf = 'reaction_kinetics/Propan_surface_and_gas.yaml'

area_cat = 7.85 * 10**(-19) # [m^2]
vol_reactor = 0.05 * 0.006 * 0.038 # katschüttung leervolumen[m3]

gas = ct.Solution(reaction_mech,'gas')
#pt_cat = ct.Interface(reaction_mech_surf, "pt_surf")

# Import gas model for surface reaction
#gas_surf = ct.Solution(reaction_mech_surf, 'pt_surf')

gas_surf = ct.Solution(reaction_mech_surf,'gas')

pt_cat = ct.Interface(reaction_mech_surf, 'pt_surf', adjacent=[gas_surf])

#pt_cat = ct.Interface(reaction_mech_surf, "edge_pt_surf", adjacent=[gas, gas_surf])

#########
print('read data!')
# relative temperature [K]
T_0_n = 293.0
# inlet temperature [K]
T_0 = 600 + 273.15
# wall temperature [K]
T_wall = 610 + 273.15
# constant pressure [Pa]
pressure = ct.one_atm
# composition
composition_0 = 'C3H8:10, H2:0'
# composition_0 = 'CH4:1, O2:1.5, AR:0.1'
# definition of initial and reactive state
initial_state = T_0, pressure, composition_0

gas.TPX = initial_state
conc_0 = gas.concentrations * 1000 #mol/m3

pt_cat.TP = T_0, pressure
cov = pt_cat.coverages
print(cov)

gas_r = ct.IdealGasReactor(contents=gas, energy = 'on')
gas_r.volume = vol_reactor
print('C3H8' +'\t'+ 'H2' +'\t'+ 'C3H6' +'\t'+ 'CH4')
print('{0:10f}  {1:10f}  {2:10f} {3:10f}'.format(*gas['C3H8', 'H2', 'C3H6', 'CH4'].X))
#print(gas.X)

sim_gas = ct.ReactorNet([gas_r])
dt = 0.01 #s
#sim_gas.advance(dt)
print('{0:10f}  {1:10f}  {2:10f} {3:10f}'.format(*gas['C3H8', 'H2', 'C3H6', 'CH4'].X))
#print(gas.X)

##

gas_surf.TPX = initial_state

r_gas_surf = ct.IdealGasReactor(contents=gas_surf, energy = 'on')
r_gas_surf.volume = vol_reactor
rsurf = ct.ReactorSurface(pt_cat, r_gas_surf, A = area_cat)

sim_surf = ct.ReactorNet([r_gas_surf])
print(gas_surf.kinetics_species_names)
print(gas_surf.X)
print(pt_cat.coverages)
sim_surf.advance_to_steady_state()
cov = pt_cat.coverages

gas_surf.TPX = initial_state
pt_cat.coverages = cov
#sim_gas.advance(dt)
print(gas_surf.kinetics_species_names)
print(gas_surf.X)
print(pt_cat.kinetics_species_names)
print(pt_cat.coverages)

#sim_surf.advance(dt)


print('{0:10f}  {1:10f}  {2:10f} {3:10f}'.format(*gas_surf['C3H8', 'H2', 'C3H6', 'CH4'].X))
print(gas_surf.kinetics_species_names)
print(gas_surf.X)
print(pt_cat.coverages)

print('finished sim!')

#print(gas.species())
#print(gas.kinetics_species_names)
#print(gas.net_production_rates)
gas_net_production_rates = gas.net_production_rates # reactionrates of gas in common with mol/m3
#print(gas_surf.species())
#print(gas_surf.net_production_rates)
#print(pt_cat.species())
#print(pt_cat.kinetics_species_names)
#print(pt_cat.net_production_rates)

surf_reaction_spec_names = pt_cat.kinetics_species_names
surf_reaction_spec_idx = []
for species in surf_reaction_spec_names:
    try: 
        surf_reaction_spec_idx.append(gas.kinetics_species_names.index(species))
    except:
        pass

gas_net_production = dict(zip(gas.kinetics_species_names,gas.net_production_rates)) #dict with species kmol/m3/s
surf_net_production = dict(zip(pt_cat.kinetics_species_names,pt_cat.net_production_rates)) #dict with species kmol/m2/s

#print(gas_net_production)
#print(surf_net_production)

#combi_reaction_rates = gas_net_production_rates[surf_reaction_spec_idx] #reactionrates of gas in common with mol/m2

#print(combi_reaction_rates)
#print(surf_reaction_spec_idx)
#production_gas = sim_gas.net_production_rates

dA_surf_cat = rsurf.area

# cV_r_gas = r_gas_surf.V
r_gas_state = r_gas_surf.get_state()
dV_r_gas = r_gas_surf.volume #The volume [m^3] of the reactor.
print('Reactorvolume '+ str(dV_r_gas) + ' m3')
print('Reactorsurface '+ str(dA_surf_cat) + ' m2')
#dV_r_gas = r_gas_state[1]

mol_net_production = dict()
for entry in gas_net_production:
    if entry in surf_net_production:
        mol_prod_entry = 1000*gas_net_production[entry] + 1000*surf_net_production[entry] * dA_surf_cat / dV_r_gas #mol/m3/s
        print(entry)
        #print(mol_prod_entry)
        print(str(gas_net_production[entry]) + ' + ' + str(surf_net_production[entry]))
    else:
         mol_prod_entry = gas_net_production[entry]*1000 #mol//m3/s
    mol_net_production.update({entry:mol_prod_entry})

#print(mol_net_production.items())
# calculation of final concentration
conc_dt = dict()
conc = dict()
for idx, element_mol_prod in enumerate(mol_net_production):
    conc[element_mol_prod] = conc_0[idx]+ mol_net_production[element_mol_prod] * dt #mol/m3
    #print(str(conc) + 'mol/m3')

for spec in ['H2', 'C3H8', 'C3H6', 'CH4']:
    print(str(spec) + ':\t\t' + str(conc[spec]) + ' mol/m3')
    print(str(spec) + '/dt:\t' + str(mol_net_production[spec]) + ' mol/m3/s')

dN_reaction = sum(mol_net_production.values()) * dV_r_gas * dt #mol
print('dN of reactions = ' + str(dN_reaction) + ' mol')

total_N_gas_direct = sum(gas.concentrations) * 1000 * dV_r_gas
print('N of conc direct = ' + str(total_N_gas_direct) + ' mol')
print('N of conc direct + dN = ' + str(total_N_gas_direct+dN_reaction) + ' mol')

total_N_gas = sum(conc.values()) * dV_r_gas
print('N of conc = ' + str(total_N_gas) + ' mol')
###

gas_reaction_spec_names = gas.kinetics_species_names
#calculate mass net production
mass_gas_r = gas_r.mass #kg
print('mass in gas = ' + str(mass_gas_r) + ' kg')

MW_total_gas_r = gas.mean_molecular_weight / 1000 #kg/mol
MW_list = gas.molecular_weights
MW_mean = 0
for idx, x in enumerate(gas.X):
    MW_mean += x * MW_list[idx] / 1000 #kg/mol

print('mean molar cantera = ' + str(MW_total_gas_r) + ' kg/mol')
print('mean molar me = ' + str(MW_mean) + ' kg/mol')

print('mean molar weight in gas = ' + str(MW_mean) + ' kg/mol')
mole_content_total_0 = mass_gas_r / MW_mean # kg/(kg/mol)=mol
print('moles in gas (from mass dens) = ' + str(mole_content_total_0) + ' mol')

dens_moles_gas = gas.density_mole * 1000 #[mol/m3]
print('density in gas = ' + str(dens_moles_gas) + ' mol/m3')

mole_content_total_0 = dV_r_gas * dens_moles_gas
print('moles in gas (from mole dens) = ' + str(mole_content_total_0) + ' mol')

abs_molar_contents_0 = []
for idx, X_i in enumerate(gas.X):
    #print(X_i)
    spec_molar_content_gas_0 = X_i * mole_content_total_0
    abs_molar_contents_0.append(spec_molar_content_gas_0)
    #print(gas_reaction_spec_names[idx] + ' ' + str(spec_molar_content_gas_0))
mass_reactor = r_gas_surf.mass
mol_weights = gas.molecular_weights #g/mol
print(gas.kinetics_species_names) #names (same sorting)

abs_molar_contents_0 = dict(zip(gas_reaction_spec_names,abs_molar_contents_0))
print(abs_molar_contents_0)
print('total n_0 from mass = ' + str(sum(abs_molar_contents_0.values())) + ' mol')

abs_molar_contents = []
for spec_name in gas_reaction_spec_names:
    abs_molar_contents.append(abs_molar_contents_0[spec_name]+ mol_net_production[spec_name] * dV_r_gas * dt) #mol
print(mol_weights)
abs_molar_contents = dict(zip(gas_reaction_spec_names,abs_molar_contents))
print(abs_molar_contents)

print('total n from concentrations = ' + str(total_N_gas) + ' mol')
print('total n from mass = ' + str(sum(abs_molar_contents.values())) + ' mol')

