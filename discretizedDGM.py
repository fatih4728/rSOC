# -*- coding: utf-8 -*-
"""
Created on Wed Aug  6 09:36:12 2025

Here I want to use the DGM in the discretized form to solve for x2

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


# necessary coefficients
epsilon = 0.5       # porosity of the electrode
tau = 0.5           # tortuosity of the electrode
M = np.array([2e-3, 28e-3, 18e-3]) # molar weights of H2, N2, H2O
rp = 1e-7           # radius of the pores
N = 3               # number of species
dEle = 1e-5         # thickness of the electrode
xT = P/R/T

# %%
xC = np.array([0.8, 0.1, 0.1])      # the molar concentration at the channel
xM = np.zeros((3))

i = 1e4           # current density is 1 A/cm²
ndotH2 = i/F/2
mdotH2 = ndotH2 * 2

# for k in range(N):
#     for l in range(N):
#         sumTerm = 0
#         if l != k:
#             sumTerm = (xC[l]*J[k] - xC[k]*J[l])/(xT*D[k, l])
#     xM[k] = d * sumTerm + d*J[k]/Dkn[k]+xC[k]*Bg/Dkn[k]/mu*dp



























