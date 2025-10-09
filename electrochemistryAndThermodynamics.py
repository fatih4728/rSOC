# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 14:18:04 2025

@author: smfadurm
"""

import numpy as np
import matplotlib.pyplot as plt
from getHS import get_thermo_properties_dict

class ElectrochemicalSystem:
    F = 96485.33212331001      # Faraday constant
    R = 8.314462618153241      # Gas constant

    def __init__(self, species_list, T, P, xH2, xH2O, xO, etaAct, l_tpb, a=25, k_c3=0.2e-16):
        self.species_list = species_list
        self.T = T
        self.P = P
        self.a = a
        self.xH2 = xH2
        self.xH2O = xH2O
        self.xO = xO
        self.etaAct = etaAct
        self.l_tpb = l_tpb
        self.k_c3 = k_c3

        # Partial pressures
        self.pH2 = xH2 * P
        self.pH2O = xH2O * P
        self.pO = xO * P
        self.pO2 = 0.5 * xO * P

        # Thermodynamics
        self.TD = get_thermo_properties_dict(species_list, T, P)
        self.calculate_gibbs()
        self.set_adsorption_properties()
        self.calculate_reaction_constants()
        self.set_kinetic_coefficients()

        # Activation energy from equilibrium constants
        K1, K2, K3, K4, K5 = self.K1to5
        self.Ea = np.log(K1*K2*K3*K4/K5*xH2/xH2O) * self.R*T / (2*self.F) * -1

    # ---------------- Thermodynamics ----------------
    def calculate_gibbs(self):
        self.dG_R = self.TD['H2O']['G'] - self.TD['H2']['G'] - 0.5*self.TD['O2']['G']
        self.H2 = np.array([self.TD['H2']['H'], self.TD['H2']['S']])
        self.H2O = np.array([self.TD['H2O']['H'], self.TD['H2O']['S']])
        return self.dG_R

    def set_adsorption_properties(self):
        self.Ni = np.array([0.0e3, 0.0])
        self.HNi = np.array([-31.19e3, 41.14])
        self.YSZ = np.array([0.0e3, 0.0])
        self.O2YSZ = np.array([-85.6e3, 139.6])
        self.OHYSZ = np.array([-165.23e3, 175.1])
        self.H2OYSZ = np.array([-254.45e3, 54.34])
        self.OxYSZ = np.array([-85.6e3, 139.6])
        self.VoYSZ = np.array([0.0e3, 0.0])

    def K(self, dHS):
        dH, dS = dHS
        return np.exp(-dH/self.R/self.T + dS/self.R)

    def calculate_reaction_constants(self):
        HS1 = 2*self.HNi - self.H2 - 2*self.Ni
        HS2 = self.Ni + self.OHYSZ - self.HNi - self.O2YSZ
        HS3 = self.Ni + self.H2OYSZ - self.HNi - self.OHYSZ
        HS4 = self.H2O + self.YSZ - self.H2OYSZ
        HS5 = self.O2YSZ + self.VoYSZ - self.OxYSZ - self.YSZ

        self.K1 = self.K(HS1)
        self.K2 = self.K(HS2)
        self.K3 = self.K(HS3)
        self.K4 = self.K(HS4)
        self.K5 = self.K(HS5) / self.a

        self.K1to5 = np.array([self.K1, self.K2, self.K3, self.K4, self.K5])
        return self.K1to5

    def set_kinetic_coefficients(self):
        self.A = 2.4e-7
        self.B = 3.55
        self.C = 0.144
        self.D = 0.437
        self.AtoD = np.array([self.A, self.B, self.C, self.D])

    # ---------------- Electrochemistry ----------------
    def calcUrev(self):
        return (-self.dG_R/2/self.F - self.R*self.T/2/self.F *
                np.log(self.pH2O/self.pH2/(self.pO2**0.5) * self.P**0.5))

    def calcLeakVoltage(self, i_A=0., i_A_lim=1.):
        A, B, C, D = self.AtoD
        eta_leak_0 = -A * np.log(self.pO2)**B + C * (self.pH2/self.P)**D
        return eta_leak_0 * (1 - np.tanh(i_A/i_A_lim))

    def currentDensity(self):
        K1, K2, K3, K4, K5 = self.K1to5
        pH2ad = 1 / K1
        i03star = self.l_tpb * self.F * self.k_c3 * ((K2*K4/K5)**0.25) * K4**0.75
        i03cross = i03star * ((self.pH2/pH2ad)**0.25 * self.pH2O**0.75) / (1 + (self.pH2/pH2ad)**0.5)
        i3 = i03cross * (np.exp(0.5*self.F*self.etaAct/self.R/self.T) -
                         np.exp(-1.5*self.F*self.etaAct/self.R/self.T)) / \
             (np.exp(-self.F*self.etaAct/self.R/self.T) * (1 + 1/K5) + (K2/K3/K4/K5*self.xH2O)**0.5)
        return i3 / 10000

    # ---------------- Surface Coverages ----------------
    def DEN_Ni(self, N):
        xH2 = np.linspace(0., 1., N)
        return 1 + (self.K1to5[0]*xH2)**0.5

    def DEN_YSZ(self, N):
        K1, K2, K3, K4, K5 = self.K1to5
        xH2 = np.linspace(0, 1, N)
        return 1 + 1/K5 + (K2/K3/K4/K5*(1-xH2))**0.5 * np.exp(self.F*self.etaAct[-1]/self.R/self.T) + 1/K4*(1-xH2)

    def coveragesNickel(self, N):
        xH2 = np.linspace(0., 1., N)
        covNi = 1/self.DEN_Ni(N)
        covH = (self.K1to5[0]*xH2)**0.5/self.DEN_Ni(N)
        return covNi, covH

    def coveragesYSZ(self, N):
        K1, K2, K3, K4, K5 = self.K1to5
        covYSZ = 1 / self.DEN_YSZ(N)
        xH2 = np.linspace(0, 1, N)
        xH2O = 1 - xH2
        covH2O = xH2O / K4 / self.DEN_YSZ(N)
        covO2m = 1 / K5 / self.DEN_YSZ(N)
        covOHm = 1 - covYSZ - covO2m - covH2O
        return covYSZ, covO2m, covH2O, covOHm
    
    
    
if __name__=="__main__":
    
    # universal constants
    
    T = 800 + 273.15        # operating temperature of the soc
    P = 101325.0          # operating pressure of the soc
    F = 96485.33212331001          # faraday constant
    R = 8.314462618153241     # import universal gas constant
    
    # gasses 
    xH2 = 0.9
    xH2O = 1 - xH2 
    xO2 = 0.21
    pH2 = xH2 * P
    pH2O = xH2O * P
    pO2 = xO2 * P
       
    
    # %%
    # parameters for the object Echem
    N = 10000
    species_list = ['H2', 'O2', 'H2O']
    etaAct = np.linspace(0., 0.8, N)
    # beta3 = 0.5
    k_c3 = 0.2e-16
    l_tpb = 10e4        # I will have to check this value
    
    # create the object Echem
    Echem = ElectrochemicalSystem(species_list, T, P, xH2, xH2O, 
                                  xO2, etaAct, l_tpb)
    
    # extract values
    dG_R = Echem.calculate_gibbs()
    K1to5 = Echem.calculate_reaction_constants()
    
    # calculate the OCV
    U_rev = Echem.calcUrev()
    eta_leak = Echem.calcLeakVoltage()
    OCV = U_rev - eta_leak
    
    
    # %%
    xH2_array = np.linspace(0., 1., N)
    
    
    # %%
    # # plot Nickel
    
    covNi, covH = Echem.coveragesNickel(N)
    
    plt.figure()
    plt.title('Surface Coverage Ni (p = 1 atm)')
    plt.plot(xH2_array, covNi, label = r'$\theta$Ni')
    plt.plot(xH2_array, covH, label = r'$\theta$H(Ni)')
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.xlabel(r'$x_{H2}$ / mol/mol')
    plt.ylabel('coverage / 1')
    plt.legend()
    
    # plot YSZ
    covYSZ, covO2m, covH2O, covOHm = Echem.coveragesYSZ(N)
    # plt.figure()
    # plt.title('Surface Coverage YSZ (p = 1 atm)')
    # plt.plot(xH2_array, covYSZ, label = r'$\theta$YSZ')
    # plt.plot(xH2_array, covO2m, label = r'$\theta O^{2-}(YSZ)$')
    # plt.plot(xH2_array, covH2O, label = r'$\theta H2O(YSZ)$', )
    # plt.plot(xH2_array, covOHm, label = r'$\theta OH^{-}(YSZ)$')
    # plt.xlim([0, 1])
    # plt.ylim([0, 1])
    # plt.xlabel(r'$x_{H2}$ / mol/mol')
    # plt.ylabel('coverage / 1')
    # plt.legend()
    
    
    # why are there two Ea??
    K1, K2, K3, K4, K5 = K1to5
    Ea = -np.log(K1*K2*K3*K4/K5)*R*T/2/F  # why not partial pressure?
    
    
    # ## calculation of i3
    voltage = OCV - etaAct
    print(voltage)
    
    plt.figure()
    plt.title('UV-Curve @ 800°C')
    plt.plot(Echem.currentDensity(), 
             voltage, label = '$x_{H2}$ = ' + str(xH2))
    plt.xlim([0, 3.])
    plt.ylim([0, 1.3])
    plt.xlabel('current / A/cm²')
    plt.ylabel('voltage / V')
    plt.legend()
    
    
    
    
    
    
    
    
    
    
    
    
