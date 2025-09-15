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
from DustyGasModel import DustyGasModel
from DustyGasModel import calculateMolarFlux
from controlVolume import ControlVolume
from CoolProp.CoolProp import PropsSI

# universal constants
T = 800 + 273.15        # operating temperature of the soc
P = 101325.0          # operating pressure of the soc
F = 96485.33212331001          # faraday constant
R = 8.314462618153241     # import universal gas constant
i = 2.499                  # current density in A/cm²


# geometry
epsilon = 0.54       # porosity of the electrode (affects mass balance; smaller -> worse)
tau = 4.81            # tortuosity of the electrode
dp = 2.7e-6
rp = 1.5e-7          # radius of the pores (affects mass balance; smaller -> worse)
dEle = 200e-6       # thickness of the electrode (affects mass balance a lot; bigger -> worse)
A = 1.1             # area in cm²

# gas variables [H2, H2O]
M = np.array([2e-3, 18e-3]) # molar weights of H2, H2O
# wC = np.array([0.9, 0.1])
# xC = w2x(wC, M)
xC = np.array([0.5, 0.5])      # the molar fraction at the channel
XC = xC*P/R/T                     # the molar concentration
VdotFuel = 140 *1e-6 / 60  # volume flow (mL/min) in m³/s
# VdotAir = 550 *1e-6 / 60   # volume flow in m³/s
Tflow = 150 + 273.15       # Temperature of the volume flows in K

# viscosity of the mix
mu = np.array([1.84e-5, 3.26e-5])    # dynamic visosities at 600°C
rhoH2 = PropsSI('D', 'P', P, 'T', Tflow, 'H2')
rhoH2O = PropsSI('D', 'P', P, 'T', Tflow, 'H2O')
rho = np.array([rhoH2, rhoH2O])

J = calculateMolarFlux(i, A)

# create the dusty gas model object
dgm = DustyGasModel(porosity = epsilon, 
                    tortuosity = tau, 
                    poreRadius = rp, 
                    c1 = XC, 
                    M = M,
                    mu = mu, 
                    J = J,
                    L = dEle,
                    T = T)
# J = calculateMolarFlux(i, A)

# get the values of the DGM diffusion
xM, dp, addInfo = dgm.calculateMoleFraction() 




# calculate the values of the control volume
mDotIn = VdotFuel * xC * rho
mDotDiff = calculateMolarFlux(i, A) * M
controlVolume = ControlVolume(mDotIn, mDotDiff, xC)
mDotOut, xOut, Uf = controlVolume.massBalance()
p_tpb = P + dp


wC = x2w(xC, M)
wM = x2w(xM, M)
wOut = x2w(xOut, M)
wC_avg = wC * 0.5 + wOut * 0.5

if dp < 1.e5:
    print(f"i = {i*1e4:.4f}")
    print(f"p_tpb = {p_tpb*1e-5:.5} bar ")
    print(f"nH2 = {J[0]:.3}")
    print("w   \t channel \t tpb")
    print(f"H2  \t {wC[0]:.5}  \t {wM[0]:.5}")
    print(f"H2O \t {wC[1]:.5}  \t {wM[1]:.5}")
    # print(f"The sum of the molar fractions is {xM.sum():.4}")
    # print(f"The utilization factor is {Uf*100:.2f} %")







