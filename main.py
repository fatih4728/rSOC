# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 11:59:01 2025

calling the class test

@author: smfadurm
"""

def x2w(x, M):
    w = x * M
    return w / w.sum()
def w2x(w, M):
    x = w / M
    return x / x.sum()

import numpy as np
from matrixFormDGM import DustyGasModelZhou
from matrixFormDGM import permeabilityFactorBg
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


# Diffusion Coefficients
D_knudsenEff = np.array([2.8e-5, 9.4e-6])
D_binaryEff = np.array([[0., 8.e-4], [8.e-4, 0.]])

# %% Control volume

# concentration at entrance
x_in = np.array([0.5, 0.5])

# flows
J = calculateMolarFlux(i, A)
VdotFuel = 140 *1e-6 / 60  
Tflow = 150 + 273.15 

# density of the mix
rhoH2 = PropsSI('D', 'P', P, 'T', Tflow, 'H2')
rhoH2O = PropsSI('D', 'P', P, 'T', Tflow, 'H2O')
rho = np.array([rhoH2, rhoH2O])
# Get viscosity of the mix
muH2 = PropsSI('V', 'P', P, 'T', Tflow, 'H2')
muH2O = PropsSI('V', 'P', P, 'T', Tflow, 'H2O')
mu = np.array([muH2, muH2O])    


# calculate the values of the control volume
mDotIn = VdotFuel * x_in * rho
mDotDiff = calculateMolarFlux(i, A) * M
controlVolume = ControlVolume(mDotIn, mDotDiff, x2w(x_in, M))
# this has to give same values to Colin
mDotOut, w_cv, Uf = controlVolume.massBalance()     
x_cv = w2x(w_cv, M)
c_cv = x_cv * cT    # the total concentration may differ though


# %% Dusty Gas Model %%

Bg = permeabilityFactorBg(epsilon, tau, rp)
dgm = DustyGasModelZhou(Bg, c_cv, M, mu, i, dEle, T, 
                        D_binaryEff, D_knudsenEff)

x_tpb, P_tpb = dgm.solveDGM()

w_zhou = x2w(x_tpb, M)
print(f'w [H2, H20] @channel is \t{w_cv}')
print(f'w [H2, H20] @tpb is     \t{w_zhou}')
print(f'The total pressure is {P_tpb*1e-5:.4} bar')














