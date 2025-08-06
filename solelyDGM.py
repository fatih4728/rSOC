# -*- coding: utf-8 -*-
"""
Created on Mon Aug  4 14:36:42 2025

Here I will try to write a code, just for the DGM
This is written largely with ChatGPT. 
I am not familiar with the source of the DGM.

@author: smfadurm
"""

import cantera as ct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import newton
from CoolProp.CoolProp import PropsSI


# %%

# Physical constants
R = 8.314  # J/(mol·K)
T = 298.15  # K
p = 101325  # Pa

# Geometry
L = 0.01  # Length of the domain [m]
N = 1000  # Number of grid points
dx = L / (N - 1)

# Porous media properties
epsilon = 0.4  # porosity
tau = 0.5      # tortuosity
rp = 1e-6      # pore radius [m]

# Species: H2 (0), N2 (1)
M = np.array([2.016e-3, 28.0134e-3])  # kg/mol
D_12 = 7.8e-5  # binary diffusion coefficient [m²/s], H2-N2 at 298K
D_12_eff = epsilon / tau * D_12

# Knudsen diffusion coefficients
def D_K(M_i):
    return (2 / 3) * rp * np.sqrt(8 * R * T / (np.pi * M_i))

D_K_eff = epsilon / tau * np.array([D_K(M[0]), D_K(M[1])])  # H2, N2

# Boundary conditions for mole fractions
x_H2_left = 0.8
x_H2_right = 0.5

# Initialize mole fraction profile
x_H2 = np.linspace(x_H2_left, x_H2_right, N)
x_N2 = 1 - x_H2

# Initialize flux arrays
J_H2 = np.zeros(N-1)
J_N2 = np.zeros(N-1)

# Solve fluxes using finite differences
for i in range(N - 1):
    dx_H2 = (x_H2[i+1] - x_H2[i]) / dx
    dx_N2 = (x_N2[i+1] - x_N2[i]) / dx

    # Use simplified DGM for binary mixture at constant p
    num = D_12_eff * D_K_eff[0]
    denom = D_12_eff + D_K_eff[0]

    J_H2[i] = - num / denom * dx_H2

    # By mass conservation in binary system
    J_N2[i] = -J_H2[i]

# Average flux
J_H2_avg = np.mean(J_H2)
J_N2_avg = np.mean(J_N2)

print(f"Average H2 flux: {J_H2_avg:.3e} mol/m²·s")
print(f"Average N2 flux: {J_N2_avg:.3e} mol/m²·s")

# Plot mole fraction profile
plt.plot(np.linspace(0, L, N), x_H2, label='H2')
plt.plot(np.linspace(0, L, N), x_N2, label='N2')
plt.xlabel('Position [m]')
plt.ylabel('Mole Fraction')
plt.title('Mole Fraction Profile')
plt.legend()
# plt.grid(True)
plt.show()

plt.figure()
plt.plot(np.linspace(0, L, N - 1), J_H2, label='H2')
# plt.plot(np.linspace(0, L, N - 1), J_N2, label='N2')
plt.xlabel('Position [m]')
plt.ylabel('Mole Flux')
plt.ylim([0, 2e-3])
plt.title('Mole Flux Profile')
plt.legend()
# plt.grid(True)
plt.show()























