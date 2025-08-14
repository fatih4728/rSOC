# -*- coding: utf-8 -*-
"""
Created on Mon Aug 11 16:22:01 2025

Hier möchte ich einfach mal das Ficksche Gesetz anwenden


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


i = 1.e4           # current density is 1 A/cm²
A = 1e-4
I = i * A
ndotH2 = I/F/2
ndotH2O = ndotH2
J = np.array([-ndotH2, ndotH2O])     # molar flux in mol/m²/s (A = 1 m²)

L = 200e-6      # thickness of the electrode
D_AB = 1e-5

# molar concentrations and fractions
x1 = np.array([0.9, 0.1])
cTotal = P/R/T
c1 = x1 *cTotal

# concentration gradient 
dcdz = -J / D_AB

# concentration at electrolyte
c2 = c1 - dcdz * L
x2 = c2 / cTotal

print(fr"x_H2 at the electrolyte is {x2[0]:.2}")





















