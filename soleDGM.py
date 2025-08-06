# -*- coding: utf-8 -*-
"""
Created on Mon Aug  4 15:03:27 2025

Code for the DGM based on the Zhu Paper

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
M_H2 = 2e-3         # molar weight of H2
M_N2 = 28e-3        # molar weight of N2
M_H2O = 18e-3       # molar weight of H2O
rp = 1e-7

# %%
# The matrix will have the components 1: H2; 2: N2; 3: H2O
N = int(3)
components = ["H2", "N2", "H2O"]
M = np.array([2e-3, 28e-3, 18e-3])
H_named = [] 
for k in range(N):
    row = []
    for l in range(N):
        pair = f"{components[k]}-{components[l]}"
        row.append(pair)
    H_named.append(row)
        

## Effective Knudsen diffusion coefficient
def knudsenCoeff(epsilon,tau,M):
    Dkn = 4/3*epsilon/tau*rp*np.sqrt(8*R*T/np.pi/M)
    return Dkn

knudsenMatrix = []
for k in range(N):
    Dknudsen = knudsenCoeff(epsilon, tau, M[k])
    knudsenMatrix.append(Dknudsen)
knudsenMatrix = np.array(knudsenMatrix)

## The effective diffusion coefficient stems from an equation I will get
#  later on. For now I will use a fixed value
diffusionMatrix = np.array([[0.0, 7.8e-5, 6-9e-5], [7.8e-5, 0, 2.1e-5], [6.9e-5, 2.1e-5, 0.0]])

# Permability by the Kozeny-Carman relationship
Bg = epsilon**3 * (2*rp)**2 / 72 / tau / (1 - epsilon)**2



# %%
## I want to have 3 gasses: N2, H2 and H2O. So N = 3
xFC = np.array([0.7, 0.2, 0.1])      # the molar fractions at the flow channel
H = []
for k in range(N):
    row = []
    for l in range(N):
        if k == l:
            sumTerm = 0
            for j in range(N):
                if j != k:
                    sumTerm += xFC[j] / diffusionMatrix[l, j]
            hkl = 1/knudsenMatrix[k] + sumTerm
        else:
            hkl = -xFC[k] / diffusionMatrix[k, l]
        row.append(hkl)
    H.append(row)    
H = np.array(H)
D_DGM = np.linalg.inv(H)


# %%
xMem = np.array([0.5, 0.2, 0.3])      # molar fraction at the membrane
delta = 1e-5                # thickness of the electrode
mu = np.array([1.84e-5, 3.78e-5, 3.26e-5])    # dynamic visosities at 600°C
muMix = xFC*mu                 # dynamic for the mix
muMix = muMix.sum()              
dp = 1e3                       # pressure difference assumed 

Jk = []
for k in range(N):
    firstPart = 0.
    secPart = 0.
    for l in range(N):
        firstPart += D_DGM[k, l]*(xFC[k]-xMem[k]/delta)
        secPart +=D_DGM[k, l]*xFC[l]/knudsenMatrix[l]
    Jkx = - firstPart - secPart*Bg/muMix*dp
    Jk.append(Jkx)
Jk = np.array(Jk)




















