# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 11:59:01 2025

calling the class test

@author: smfadurm
"""

import matplotlib.pyplot as plt
import numpy as np
from DustyGasModel import DustyGasModelZhou, permeabilityFactorBg
from DustyGasModel import calculateMolarFlux, x2w, w2x
from DustyGasModel import D_Fuller, D_Knudsen
from DustyGasModel import ControlVolume
from CoolProp.CoolProp import PropsSI
from electrochemistryAndThermodynamics_byCurrent import ElectrochemicalSystem
from XY_2D import standardized_plot

def getDensities(components, P, T):
    rhoH2 = PropsSI('D', 'P', P, 'T', T, components[0])
    rhoH2O = PropsSI('D', 'P', P, 'T', T, components[1])
    return np.array([rhoH2, rhoH2O])

# %% Parameters %%

# control parameter
xH2_in = 0.5

# universal constants
T = 700 + 273.15        # operating temperature of the soc
P = 101325.0          # operating pressure of the soc
F = 96485.33212331001          # faraday constant
R = 8.314462618153241     # import universal gas constant


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


## Diffusion Coefficients
D_knudsenEff = D_Knudsen(rp, T, M*1e3) * epsilon / tau
D_binaryEff = np.array([[0., D_Fuller(T)], [D_Fuller(T), 0.]]) *epsilon/tau

# flow parameters
VdotFuel = 140 *1e-6 / 60  
Tflow = 150 + 273.15 

# density of the mix
rho = getDensities(['H2', 'H2O'], P, Tflow)
# Get viscosity of the mix
muH2 = PropsSI('V', 'P', P, 'T', Tflow, 'H2')
muH2O = PropsSI('V', 'P', P, 'T', Tflow, 'H2O')
mu = np.array([muH2, muH2O])  
# concentration at entrance
x_in = np.array([xH2_in, 1 - xH2_in])
w_in = x2w(x_in, M)
mDotIn = VdotFuel * x_in * rho

# parameters for the object Echem
l_tpb = 10e4        # I will have to check this value
xO2 = 0.21


# %% Start the loop here
N = 15
currentDensities =np.linspace(1e-4, 5.e4, N)

# creating lists for storage
voltages = []
pressures=[]
x_tpbList=[]

for currentDensity in currentDensities:
   
    
    # %% Control volume Anode
      
    # flows
    J = calculateMolarFlux(currentDensity, A)
      
    # calculate the values of the control volume
    mDotDiff = J * M
    controlVolume = ControlVolume(mDotIn, mDotDiff, x2w(x_in, M))
    # this has to give same values to Colin
    mDotOut, w_cv, Uf = controlVolume.massBalance()     
    x_cv = w2x(w_cv, M)
    c_cv = x_cv * cT    # the total concentration may differ though
    
    
    # %% Dusty Gas Model Anode
    Bg = permeabilityFactorBg(epsilon, tau, rp) # integrate into class
    dgm = DustyGasModelZhou(Bg, c_cv, M, mu, currentDensity, dEle, T, 
                            D_binaryEff, D_knudsenEff)
    x_tpb, P_tpb = dgm.solveDGM()
    if x_tpb[0] < 0:
        print('##############################')
        print('Hydrogen is  depleted. Try increasing hydrogen \
              concentration or decreasing the current density')
        print(f'Current density: {currentDensity*1e-4:.2}')
        # print(f'Overpotential: {etaAct:.2}')
        print('##############################')
        break
    
    pressures.append(P_tpb)
    x_tpbList.append(x_tpb)
    
    # %% Electrochemistry
    
    xH2, xH2O = x_tpb
    
    # create the object Echem
    Echem = ElectrochemicalSystem(['H2', 'O2', 'H2O'], 
                                  T, P, xH2, xH2O, xO2, 
                                  currentDensity, l_tpb,
                                  k_c3=5.e-14)
    
    # extract values
    dG_R = Echem.calculate_gibbs()
    K1to5 = Echem.calculate_reaction_constants()
    etaAct = Echem.calculateOverpotential()

    # check, whether upper boundary has been reached
    initial_guess = 0.1
    tol = 1e-3
    mask_failed = np.isclose(etaAct, initial_guess, atol=tol)
    if mask_failed:
        break
    
    # calculate the OCV
    U_rev = Echem.calcUrev()
    eta_leak = Echem.calcLeakVoltage()
    OCV = U_rev - eta_leak
    
    
    voltage = OCV - etaAct
    voltages.append(voltage)
    # currentDensities.append(Echem.currentDensity())


voltages = np.array(voltages)
currentDensities = currentDensities[0:len(voltages)]
x_tpbList = np.array(x_tpbList)
xH2 = x_tpbList[:, 0]
xH2O = x_tpbList[:, 1]

standardized_plot(currentDensities*1e-4, [voltages], 
                  r'current density / A cm$^{-2}$', 'voltage / V',
                  savepath=None,
                  labels=['H2-H2O'],
                  marker='.',
                  startAtZero=True,
                  xAxisLim=[0, 5])

# standardized_plot(currentDensities*1e-4, [xH2, xH2O],
#                   r'current density / A cm$^{-2}$', 'molar fraction / 1',
#                   labels=['H2', 'H2O'],
#                   xAxisLim=[0, 5])
















