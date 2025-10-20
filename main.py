# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 11:59:01 2025

calling the class test

@author: smfadurm
"""

import numpy as np
from DustyGasModel import DustyGasModelZhou, permeabilityFactorBg
from DustyGasModel import calculateMolarFlux, x2w, w2x
from DustyGasModel import D_Fuller, D_Knudsen
from DustyGasModel import ControlVolume, getDensities, getViscosities
from eChemAndThermo import ElectrochemicalSystem
from XY_2D import standardized_plot


# %% Parameters %%

# control parameter
xH2_in = 0.97
N = 15
currentDensities =np.linspace(1e-4, 4.25e4, N)

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
M_fuel = np.array([2e-3, 18e-3])   
M_air = np.array([32e-3, 28e-3])
cT = P/R/T 
w = np.array([wH2 := 0.042427, 1 - wH2])
x_ch = w2x(w, M_fuel)
c_ch = x_ch * cT


## Diffusion Coefficients
D_knudsenEff = D_Knudsen(rp, T, M_fuel*1e3) * epsilon / tau
D_binaryEff = np.array([[0., D_Fuller(T)], [D_Fuller(T), 0.]]) *epsilon/tau

# flow parameters
VdotFuel = 140 * 1e-6 / 60  
VdotAir = 550 * 1e-6 / 60  
Tflow = 150 + 273.15 

# density of the mix
rho_fuel = getDensities(['H2', 'H2O'], P, Tflow)
rho_air = getDensities(['O2', 'N2'], P, Tflow)
# Get viscosity of the mix
mu_fuel = getViscosities(['H2', 'H2O'], P, Tflow)
mu_air = getViscosities(['O2', 'N2'], P, Tflow)
# concentration at entrance
x_fuel = np.array([xH2_in, 1 - xH2_in])
mDotFuelIn = VdotFuel * x_fuel * rho_fuel
xO2 = 0.21
x_air = np.array([xO2, 1 - xO2])
mDotAirIn = VdotAir * x_air * rho_air

# parameters for the object Echem
l_tpb = 10e4        # I will have to check this value


# %% Start the loop here

# creating lists for storage
voltages = []
pressures=[]
x_tpbList=[]

for currentDensity in currentDensities:
   
    
    # %% Control volume 
      
    # flows
    J_fuel = calculateMolarFlux(currentDensity, A)
    J_air = np.array([0.5 * J_fuel[0], 0])
      
    # fuel electrode
    mDotDiff_fuel = J_fuel * M_fuel
    controlVolume_fuel = ControlVolume(mDotFuelIn, mDotDiff_fuel, x2w(x_fuel, M_fuel))
    mDotOut_fuel, w_fuel_channel, Uf = controlVolume_fuel.massBalance()     
    x_fuel_channel = w2x(w_fuel_channel, M_fuel)
    c_fuel_channel = x_fuel_channel * cT    
    
    # air electrode
    mDotDiff_air = J_air * M_air
    controlVolume_air = ControlVolume(mDotAirIn, mDotDiff_air, x2w(x_air, M_air))
    mDotOut_air, w_air_channel, Uf_air = controlVolume_air.massBalance()     
    x_air_channel = w2x(w_air_channel, M_air)
    c_air_channel = x_air_channel * cT

    
    # %% Dusty Gas Model fuel electrode
    Bg = permeabilityFactorBg(epsilon, tau, rp) # integrate into class
    dgm_fuel = DustyGasModelZhou(Bg, c_fuel_channel, M_fuel, mu_fuel, 
                                 currentDensity, -J_fuel / A,
                                 dEle, T, 
                                 D_binaryEff, D_knudsenEff)
    x_tpb, P_tpb = dgm_fuel.solveDGM()
    if x_tpb[0] < 0:
        print('##############################')
        print('Hydrogen is  depleted. Try increasing hydrogen \
              concentration or decreasing the current density')
        print(f'Current density: {currentDensity*1e-4:.3}')
        # print(f'Overpotential: {etaAct:.2}')
        print('##############################')
        break
   
    pressures.append(P_tpb)
    x_tpbList.append(x_tpb)
   
    # %% Dusty Gas Model Air Electrode
    dgm_air = DustyGasModelZhou(Bg, c_air_channel, M_air, mu_air, 
                                currentDensity, -J_air / A, 
                                dEle, T, 
                                D_binaryEff, D_knudsenEff)
    x_tpb_air, P_air = dgm_air.solveDGM() 
   
    # %% Electrochemistry
    xH2, xH2O = x_tpb
    xO2, xN2 = x_tpb_air
    
    # create the object Echem
    Echem = ElectrochemicalSystem(['H2', 'O2', 'H2O'], 
                                  T, P, xH2, xH2O, xO2, 
                                  currentDensity, l_tpb,
                                  k_c3=7.e-14)
    
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
    eta_leak = Echem.calcLeakVoltage(1.)
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
















