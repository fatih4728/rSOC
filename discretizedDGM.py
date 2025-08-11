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
epsilon = 0.3       # porosity of the electrode
tau = 5.           # tortuosity of the electrode
M = np.array([2e-3, 28e-3, 18e-3]) # molar weights of H2, N2, H2O
rp = 1e-7           # radius of the pores
N = 3               # number of species

dEle = 20e-5         # thickness of the electrode
XT = P/R/T

# %%
xC = np.array([0.8, 0.0, 0.2])      # the molar fraction at the channel
XC = xC*XT                          # the molar concentration



i = 1.e4           # current density is 1 A/cm²
ndotH2 = i/F/2
ndotO2 = ndotH2/2
ndotH2O = ndotH2
J = np.array([-ndotH2, -ndotO2 * 0., ndotH2O])     # molar flux in mol/m²/s (A = 1 m²)

# the diffusion matrix, with arbitrary numbers
D = np.array([[0.0, 7.8e-5, 6.9e-5], [7.8e-5, 0, 2.1e-5], [6.9e-5, 2.1e-5, 0.0]])


# knudsen Array
def knudsenCoeff(epsilon,tau,M):
    Dkn = 4/3*epsilon/tau*rp*np.sqrt(8*R*T/np.pi/M)
    return Dkn
knudsenArray = []
for k in range(N):
    Dknudsen = knudsenCoeff(epsilon, tau, M[k])
    knudsenArray.append(Dknudsen)
Dkn = np.array(knudsenArray)

# Permability by the Kozeny-Carman relationship
Bg = epsilon**3 * (2*rp)**2 / 72 / tau / (1 - epsilon)**2

# viscosity of the mix
mu = np.array([1.84e-5, 3.78e-5, 3.26e-5])    # dynamic visosities at 600°C
muMix = xC*mu                 # dynamic for the mix
muMix = muMix.sum()         

# pressure difference
dp = 1000.       # arbitrary pressure differences in Pa
dpdz = dp/dEle

# calulate the dx/dz
dxdz = np.zeros((3))
for k in range(N):
    sumTerm = 0.
    for l in range(N):  
        if l != k:
            sumTerm += (XC[l]*J[k] - XC[k]*J[l])/(XT*D[k, l])
    dxdz[k] = sumTerm + J[k]/Dkn[k] - XC[k]*Bg/Dkn[k]/muMix*dpdz
dxdz *= -1

XM = XC - dxdz * dEle


# calculate the new x
# for k in range(N):
#     for l in range(N):
#         sumTerm = 0.
#         if l != k:
#             sumTerm = (xC[l]*J[k] - xC[k]*J[l])/(xT*D[k, l])
#     xM[k] = d * sumTerm + d*J[k]/Dkn[k]+xC[k]*Bg/Dkn[k]/muMix*dp

# Np = 1000
# XH2 = np.linspace(XC[0], XM[0], Np)
# XO2 = np.linspace(XC[1], XM[1], Np)
# XH2O = np.linspace(XC[2], XM[2], Np)
# xDEle = np.linspace(0, dEle, Np)
# # plt.plot(xDEle, XH2, label=r'$X_{H2}$')
# # plt.plot(xDEle, XO2, label=r'$X_{O2}$')
# plt.plot(xDEle, XH2O, label=r'$X_{H2O}$')


























