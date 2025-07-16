# -*- coding: utf-8 -*-
"""
Created on Mon Jul 14 14:07:55 2025

this will be my code for the dusty gas model

I will first start, by defining a container, where gasses can enter and leave


@author: smfadurm
"""

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
## just a small change to the script


# gasses 
xH2 = 0.9
xH2O = 1 - xH2 
xO2 = 0.21 
pH2 = xH2 * P
pH2O = xH2O * P
pO2 = xO2 * P


# Gas flow
VdotFuel = 140 *1e-6 / 60  # volume flow in m³/s
VdotAir = 550 *1e-6 / 60   # volume flow in m³/s
Tflow = 150 + 273.15       # Temperature of the volume flows in K



# I want to define an environment, that is at first filled
# only with N2, the flush gas. Afterwards, it gets filled with
# the entering gasses. I will do the anode side, for a button cell
# with a diameter of 2 cm² and a thickness of 1 mm

Dia = 2e-2    # diameter of the button cell   / m
Cross = Dia * np.pi      # Cross area of the cell / m²
Depth = 1e-3    # depth of the cell / m
Volume = Cross * Depth   # volume of the anode and the cathode

# I will assume, it is filled with H2O at 150°C and Patm
rhoH2 = PropsSI('D', 'P', P, 'T', Tflow, 'H2')
rhoH2O = PropsSI('D', 'P', P, 'T', Tflow, 'H2O')
mH2O = Volume / rhoH2O
mH2 = 0.

print('The cell as a diameter of ', Dia*1e2 , ' cm and a Volume \
      of ', Volume*1e6 ,' milliliter. \n')
print('The starting mass of H2O is ', mH2O*1e3, ' gram.')      


# I want to calculate the ingoing mass flux
VdotH2 = VdotFuel * xH2
VdotH2O = VdotFuel - VdotH2

mdotH2 = VdotH2 * rhoH2
mdotH2O = VdotH2O * rhoH2O

# for sake of simplicity, I will say, that the outgoing mass flux is equal
mdotOutTotal = mdotH2 + mdotH2O
# mdotOutH2 = mdotOutTotal * 0.5
# mdotOutH2O = mdotOutTotal - mdotOutH2


## here, I will start the loop
# I am not sure, on how to exactly do this, maybe by pressure?
# also, the utilization factor would be helpful
# do I do this steady state or dynamic?
mH2i = []
mH2Oi = []
t = np.arange(0, 1000, 1)
for i in t:
    # loop over 1000 seconds in 1 second steps
    mdotOutH2 = mH2 * 0.1
    mH2 = mdotH2 - mdotOutH2
    
    mdotOutH2O = mH2O * 0.2
    mH2O = mdotH2O - mdotOutH2O
    mH2i.append(mH2)
    mH2Oi.append(mH2O)
mH2i = np.array(mH2i)
mH2Oi = np.array(mH2Oi)


plt.figure()
plt.plot(t, mH2i, label='m_H2')
plt.plot(t, mH2Oi, label='m_H2O')
plt.legend()














