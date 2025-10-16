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
from electrochemistryAndThermodynamics import ElectrochemicalSystem
from XY_2D import standardized_plot

# %% Parameters %%

# universal constants
T = 700 + 273.15        # operating temperature of the soc
P = 101325.0          # operating pressure of the soc
F = 96485.33212331001          # faraday constant
R = 8.314462618153241     # import universal gas constant
# i = np.array([1369.7155, 15468.857, 24999.998])
# i = i[0]


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
D_knudsen = D_Knudsen(rp, T, M*1e3)
D_knudsenEff = D_knudsen * epsilon / tau
D_binaryEff = np.array([[0., D_Fuller(T)], [D_Fuller(T), 0.]]) *epsilon/tau


# %% Start the loop here

etaActList = np.linspace(0., 0.54, 20)
voltages = []
currentDensities =[]
pressures=[]
for etaAct in etaActList:
        
    # %% Electrochemistry
    # parameters for the object Echem
    
    k_c3 = 0.2e-16
    l_tpb = 10e4        # I will have to check this value
    x_tpb = np.array([0.1, 0.9])
    xH2, xH2O = x_tpb
    xO2 = 0.21
    
    # create the object Echem
    Echem = ElectrochemicalSystem(['H2', 'O2', 'H2O'], 
                                  T, P, xH2, xH2O, xO2, 
                                  etaAct, l_tpb)
    
    # extract values
    dG_R = Echem.calculate_gibbs()
    K1to5 = Echem.calculate_reaction_constants()
    i = Echem.currentDensity()*1e4
    
    # calculate the OCV
    U_rev = Echem.calcUrev()
    eta_leak = Echem.calcLeakVoltage()
    OCV = U_rev - eta_leak
    
    
    # %% Control volume
    
    # concentration at entrance
    x_in = np.array([0.5, 1 - 0.5])
    w_in = x2w(x_in, M)
    
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
    mDotDiff = J * M
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
    if x_tpb[0] < 0:
        print('##############################')
        print('Hydrogen is  depleted. Try increasing hydrogen \
              concentration or decreasing the current density')
        print(f'Current density: {i*1e-4:.2}')
        print(f'Overpotential: {etaAct:.2}')
        print('##############################')
        break
    
    
    w_zhou = x2w(x_tpb, M)
    # print(f'w[H2, H20] @entrance is \t{w_in}')
    # print(f'w[H2, H20] @channel is \t{w_cv}')
    # print(f'w[H2, H20] @tpb is     \t{w_zhou}')
    # print(f'The total pressure is {P_tpb*1e-5:.4} bar')
    pressures.append(P_tpb)
    
    # plt.figure()
    # plt.title('concentration plot')
    # plt.plot(['input', 'channel', 'tpb'], [w_in[0], w_cv[0], w_zhou[0]], 'k--o')
    # plt.ylim([0, 0.11])
    
    # %%
    voltage = OCV - etaAct
    voltages.append(voltage)
    currentDensities.append(Echem.currentDensity())


voltages = np.array(voltages)
currentDensities = np.array(currentDensities)

standardized_plot(currentDensities, [voltages], 
                  r'current density / A cm$^{-2}$', 'voltage / V',
                  savepath=r"C:\users\smfadurm\Desktop\UV",
                  labels=['H2-H2O'],
                  marker='.')

# plt.figure()
# plt.title('UV-Curve @ 800°C')
# plt.plot(currentDensities, 
#          voltages, 'ro' , label = '$x_{H2}$ = ' + str(xH2))
# plt.xlim([0, 5.])
# plt.ylim([0, 1.3])
# plt.xlabel('current / A/cm²')
# plt.ylabel('voltage / V')
# plt.legend()
















