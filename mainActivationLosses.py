# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 10:42:31 2025

temporary test main file for the activation losses.

@author: smfadurm
"""
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 11:59:01 2025

calling the class test

@author: smfadurm
"""


import numpy as np
from matrixFormDGM import DustyGasModelZhou, permeabilityFactorBg
from matrixFormDGM import calculateMolarFlux, x2w, w2x
from controlVolume import ControlVolume
from CoolProp.CoolProp import PropsSI

# %% Parameters %%

# universal constants
T = 700 + 273.15        # operating temperature of the soc
P = 101325.0          # operating pressure of the soc
F = 96485.33212331001          # faraday constant
R = 8.314462618153241     # import universal gas constant
i = np.array([1369.7155, 15468.857, 24999.998])
i = i[2]


# geometry
epsilon = 0.54       # porosity of the electrode (affects mass balance; smaller -> worse)
tau = 4.81            # tortuosity of the electrode
dp = 2.7e-6
rp = 1.5e-7          # radius of the pores (affects mass balance; smaller -> worse)
dEle = 200e-6       # thickness of the electrode (affects mass balance a lot; bigger -> worse)
A = 1.1e-4             # area in cm²

# gas variables [H2, H2O]
M = np.array([2e-3, 18e-3])   
cT = P/R/T 
w = np.array([wH2 := 0.042427, 1 - wH2])
x_ch = w2x(w, M)
c_ch = x_ch * cT



# %% Control volume

# concentration at entrance
x_in = np.array([0.5, 0.5])

# flows
J = calculateMolarFlux(i, A)
VdotFuel = 140 *1e-6 / 60  
Tflow = 150 + 273.15 















