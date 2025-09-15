# -*- coding: utf-8 -*-
"""
Created on Mon Sep  8 10:40:39 2025

@author: smfadurm
"""

import cantera as ct

# create a gas-phase object to represent the gas in the pores, with a
# dusty gas transport manager
# g = ct.DustyGas('h2o2.yaml')
g = ct.DustyGas('gri30.yaml')


# set the gas state
composition = {"H2": 1, "H2O":1}
g.TPX = 973.0, ct.one_atm, composition
X = g.X

# set its parameters
g.porosity = 0.5
g.tortuosity = 4.0
g.mean_pore_radius = 1.5e-7
g.mean_particle_diameter = 200e-6  # lengths in meters

# print the multicomponent diffusion coefficients
print(g.multi_diff_coeffs[0:2])
D = g.multi_diff_coeffs
for i in range(len(D)):
    for j in range(len(D)):
        if D[i,j] < 1e-20:
            D[i, j] = 0.

# print the thermal conductivity of the gas phase
# print(g.thermal_conductivity)

# compute molar species fluxes
T1, rho1, Y1 = g.TDY

g.TP = g.T, 1.2 * ct.one_atm
T2, rho2, Y2 = g.TDY
delta = 0.001

# print(g.molar_fluxes(T1, T1, rho1, rho1, Y1, Y1, delta))
# print(g.molar_fluxes(T1, T2, rho1, rho2, Y1, Y2, delta))