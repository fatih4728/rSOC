# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 08:59:08 2025

control volume with an entry, an outlet, and an exchange with the electrode

@author: smfadurm
"""

# %%
import cantera as ct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import newton
from CoolProp.CoolProp import PropsSI

# %%
# universal constants

T = 800 + 273.15        # operating temperature of the soc
P = ct.one_atm          # operating pressure of the soc
F = ct.faraday*1e-3          # faraday constant
R = ct.gas_constant*1e-3     # import universal gas constant

# %%
i = 1.e4           # current density in A/m²
A = 1.1e-4
I = i * A
ndotH2 = I/F/2
ndotH2O = ndotH2
J = np.array([-ndotH2, ndotH2O])     # molar flux in mol/m²/s (A = 1 m²)


# %%
VdotFuel = 140 *1e-6 / 60  # volume flow (mL/min) in m³/s
# VdotAir = 550 *1e-6 / 60   # volume flow in m³/s
Tflow = 150 + 273.15       # Temperature of the volume flows in K
Tstandard = 25 + 273.15


rhoH2 = PropsSI('D', 'P', P, 'T', Tflow, 'H2')
rhoH2O = PropsSI('D', 'P', P, 'T', Tflow, 'H2O')
rho = np.array([rhoH2, rhoH2O])
M = np.array([2e-3, 18e-3])
xIn = np.array([0.9, 0.1])

# mass balance total
mDotIn = VdotFuel * xIn * rho
mDotDiff = J * M
mDotOut = mDotDiff + mDotIn
if mDotOut[0] < 0 or mDotOut[1] < 0:
    print("the outlet is not outletting")
    

# mass balance H2
xOut = mDotOut / mDotOut.sum()

# the utilization factor can be calculated
Uf = 100 - xOut[0]/xIn[0] * 100
print(f"The utilization factor is {Uf:.2f} %")

# %%
import matplotlib.patches as patches


# Create figure and axes
fig, ax = plt.subplots(figsize=(5, 4))

# Draw control volume as a rectangle
control_volume = patches.Rectangle((1, 1), 4, 2, linewidth=2, edgecolor='black', facecolor='none')
ax.add_patch(control_volume)

# Add arrows for mass fluxes
ax.annotate('', xy=(1, 2), xytext=(0, 2), arrowprops=dict(arrowstyle='->', linewidth=2))
ax.annotate('', xy=(5, 2), xytext=(6, 2), arrowprops=dict(arrowstyle='<-', linewidth=2))
ax.annotate('', xy=(3, 3.1), xytext=(3, 4.1), arrowprops=dict(arrowstyle='<-', linewidth=2))

# Labels for arrows
ax.text(0.2, 2.1, r'$\dot{m}_{in}$', fontsize=14, verticalalignment='bottom')
ax.text(5.2, 2.1, r'$\dot{m}_{out}$', fontsize=14, verticalalignment='bottom')
ax.text(3.1, 3.2, r'$\dot{m}_{diff}$', fontsize=14, verticalalignment='bottom')


# Formatting
ax.set_xlim(-1, 7)
ax.set_ylim(0, 4)
ax.set_aspect('equal', adjustable='box')
ax.axis('off')

plt.show()


