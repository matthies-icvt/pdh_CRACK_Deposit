"""
Constant-pressure, adiabatic kinetics simulation.

Requires: cantera >= 2.5.0, matplotlib >= 2.0
Keywords: combustion, reactor network, plotting
"""

import sys
import matplotlib.pyplot as plt
import cantera as ct
import numpy as np
import pandas as pd

plot = str(input('Plot biggest components? (y/n)\n'))
n_plots = int(input('How many components to plot (int) \n)'))

# Mechanism?
mech13 = '../reaction_kinetics/mech_13.yaml'
gas = ct.Solution(mech13)

# Other reaction for debug
# gas = ct.Solution('h2o2.yaml')
# gas.TPX = 1001.0, ct.one_atm, 'H2:2,O2:1,N2:4'

# Initial State?
gas.TPX = 1000, ct.one_atm, 'C3H8:10, H2:1'
r = ct.IdealGasConstPressureReactor(gas)

# Reaktor
sim = ct.ReactorNet([r])
sim.verbose = True

# Time step and number of steps
dt_max = 1.e-5
steps = 1e4
t_end = steps * dt_max
# states = ct.SolutionArray(gas, extra=['t'])

# Create ndarrays for results
X = np.array([r.thermo.X])
Names =np.array(r.thermo.species_names)
M = r.thermo.molecular_weights
# Sort acc to Moll Mass
sorted_to = M.argsort()
M = M[sorted_to[::-1]]
Names = Names[sorted_to[::-1]]

# Simulation loop
while sim.time < t_end:
    sim.advance(sim.time + dt_max)
    x = r.thermo.X
    X = np.vstack([X,x])

# # Debug prints
# print(Names)
# print(X[0])
# sort X acc to Moll Mass
for i in range(int(steps) + 1):
    X[i] = X[i][sorted_to[::-1]]
# print(X[0])
# Change ndarray type to str
M_str = M.astype('U')
# Create index array for dataframe
columnI = [Names[i] + ' ' + M_str[i] for i in range(len(Names))]
# Names = np.insert(Names,0, 'time')
times = np.arange(0,int(steps) + 2)
df_result = pd.DataFrame(data=X,
             index=times,
             columns=columnI)
df_result.to_csv('BatchTest.csv')

Id_plot=[]

# Plot biggest compounds that are present in system
for i in range(len(Names)):
    for j in range(int(steps)):
        if X[j][i] > 10e-8:
            if i not in Id_plot:
                Id_plot.append(i)

fig, ax = plt.subplots(figsize=(8,5))
for i in (Id_plot[:n_plots]):
    ax.plot(times[:-1], df_result[(Names[i] + ' ' + M_str[i])][:-1], label=(Names[i]+' '+M_str[i]))
fig.legend()
fig.show()








print("To view a plot of these results, run this script with the option --plot")