# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 11:59:01 2025

calling the class test

@author: smfadurm
"""

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
i = 1.                   # current density in A/cm²


# geometry
epsilon = 0.3       # porosity of the electrode (affects mass balance; smaller -> worse)
tau = 5.            # tortuosity of the electrode
rp = 5.e-7          # radius of the pores (affects mass balance; smaller -> worse)
dEle = 200e-6       # thickness of the electrode (affects mass balance a lot; bigger -> worse)
A = 10             # area in cm²

# gas variables [H2, H2O]
M = np.array([2e-3, 18e-3]) # molar weights of H2, H2O
xC = np.array([0.9, 0.1])      # the molar fraction at the channel
XC = xC*P/R/T                     # the molar concentration
VdotFuel = 140 *1e-6 / 60  # volume flow (mL/min) in m³/s
# VdotAir = 550 *1e-6 / 60   # volume flow in m³/s
Tflow = 150 + 273.15       # Temperature of the volume flows in K

# viscosity of the mix
mu = np.array([1.84e-5, 3.26e-5])    # dynamic visosities at 600°C
muMix = (xC*mu).sum()                 # dynamic for the mix
rhoH2 = PropsSI('D', 'P', P, 'T', Tflow, 'H2')
rhoH2O = PropsSI('D', 'P', P, 'T', Tflow, 'H2O')
rho = np.array([rhoH2, rhoH2O])


# create the dusty gas model object
dgm = DustyGasModel(porosity = epsilon, 
                    tortuosity = tau, 
                    poreRadius = rp, 
                    c1 = XC, 
                    M = M, 
                    muMix = muMix, 
                    J = calculateMolarFlux(i, A),
                    L = dEle,
                    T = T)
# get the values of the DGM diffusion
xM, dp, addInfo = dgm.calculateMoleFraction()

# calculate the values of the control volume
mDotIn = VdotFuel * xC * rho
mDotDiff = calculateMolarFlux(i, A) * M
controlVolume = ControlVolume(mDotIn, mDotDiff, xC)
mDotOut, xOut, Uf = controlVolume.massBalance()



if dp < 1.e5:
    print(f"\nxH2 at the electrolyte is {xM[0]:.5}")
    print(f"The sum of the molar fractions is {xM.sum():.4}")
    print(f"dp = {dp*1e-2:.3} mbar")
    print(f"The utilization factor is {Uf*100:.2f} %")







