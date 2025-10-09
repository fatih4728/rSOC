# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 09:41:06 2025


In this model I will try to recreate the SOFC model of Colin

@author: smfadurm
"""
import numpy as np
import matplotlib.pyplot as plt
from getHS import get_thermo_properties_dict


# %%
class Echem:
    def __init__(self, T, P, xH2, xH2O, xO, AtoD, etaAct, l_tpb, 
                  K1to5, k_c3 = 0.2e-16):
        self.T = T
        self.P = P
        self.xH2 = xH2
        self.xH2O = xH2O
        self.xO = xO
        self.A_D = AtoD
        self.etaAct = etaAct
        self.l_tpb = l_tpb
        self.k_c3 = k_c3
        self.pH2 = xH2 * P
        self.pH2O = xH2O * P
        self.pO = xO * P
        self.pO2 = xO / 2* P
        self.F = 96485.33212331001          # faraday constant
        self.R = 8.314462618153241
        self.K1to5 = K1to5
        K1, K2, K3, K4, K5 = K1to5
        Ea = np.log(K1*K2*K3*K4/K5*xH2/xH2O)*R*T/2/F*-1
        self.Ea = Ea

    def calcUrev(self, dG_R):
        """
        calculate the reversible voltage (Nernst voltage)
    
        Parameters
        ----------
        dG_R : float
            reaction enthalpy.
    
        Returns
        -------
        float
            reversible voltage Urev.
    
        """
        return (-dG_R/2/self.F - self.R*self.T/2/self.F*
                np.log(self.pH2O/self.pH2/(self.pO2**0.5) * 
                       (self.P**0.5)))

    def calcLeakVoltage(self, i_A= 0., i_A_lim = 1.):
        """
        calculates leak voltage. I do not know, whether i_A and i_A_lim will be changed ever
    
        Parameters
        ----------
        i_A : float
            dunno.
        i_A_lim : float
            dunno.
      
        Returns
        -------
        eta_leak : TYPE
            DESCRIPTION.
    
        """
        A, B, C, D = self.A_D
        eta_leak_0 = (-A * np.log(self.pO2)**B + 
                      C * (self.pH2/self.P)**D)
        eta_leak = eta_leak_0 * (1 - np.tanh(i_A/i_A_lim))
        return eta_leak
    
    def currentDensity(self):
        K1, K2, K3, K4, K5 = self.K1to5
        pH2ad = 1/K1
        i03star = (self.l_tpb * self.F * self.k_c3 * (
                  K2*K4/K5)**0.25*K4**0.75)
        i03cross = i03star * ((self.pH2/pH2ad)**0.25 * pH2O**0.75)/(1+(pH2/pH2ad)**0.5)
        i3 = (i03cross * (np.exp(0.5*self.F*self.etaAct/self.R/self.T)-
                          np.exp(-1.5*self.F*self.etaAct/self.R/self.T)) / 
              (np.exp(-self.F*self.etaAct/self.R/self.T)*
               (1+1/K5)+(K2/K3/K4/K5*self.xH2O)**0.5))
        return i3 / 10000

    def DEN_Ni(self, N):
         xH2 = np.linspace(0., 1., N)
         return 1 + (self.K1to5[0]*xH2)**0.5
 
    def DEN_YSZ(self, N):
         # 1 + 1/K5 + K2/K5*(K1*xH2)**0.5 * np.exp(F*Ea/R/T) + 1/K4*(1-xH2)
         K1, K2, K3, K4, K5 = self.K1to5
         xH2 = np.linspace(0, 1, N)
         return (1 + 1/K5 + (K2/K3/K4/K5*(1-xH2))**0.5 *
                 np.exp(self.F*self.etaAct[-1]/self.R/self.T) + 
                 1/K4*(1-xH2))

    def coveragesNickel(self, N):
        xH2 = np.linspace(0., 1., N)
        covNi = 1/Echem.DEN_Ni(N)
        covH = (K1to5[0]*xH2)**0.5/Echem.DEN_Ni(N)
        return covNi, covH

    def coveragesYSZ(self, N):
        K1, K2, K3, K4, K5 = self.K1to5
        covYSZ = 1 / Echem.DEN_YSZ(N)
        xH2 = np.linspace(0, 1, N)
        xH2O = 1 - xH2 # isnt this necessary?
        covH2O = 1 / K4*(xH2O)/Echem.DEN_YSZ(N)    
        covO2m = 1 / K5 / Echem.DEN_YSZ(N)
        # covOHm = K2/K5*(K1*xH2)**0.5*np.exp(F*Ea1/R/T)/DEN_YSZ(xH2, Ea)
        covOHm = 1 - covYSZ - covO2m - covH2O
        # covSum = covYSZ + covO2m + covH2O + covOHm
        # if abs(covSum-0.01) != 1:
        #     print(covYSZ + covO2m + covH2O + covOHm)
        #     print('coverages not adding up to 1')
        return covYSZ, covO2m, covH2O, covOHm
    # coveragesYSZ = np.vectorize(coveragesYSZ)

class Reaction:
    def __init__(self, P, T):
        self.P = P
        self.T = T  


def GibbsFreeEnergy(HS):
    """
    simple function to calculate the free Gibbs energy

    Parameters
    ----------
    HS : Array
        Enthalpy and entropie.

    Returns
    -------
    float
        Free Gibbs enthalpy.

    """
    H, S = HS
    return H - T * S

def ReactionH2O2_H2O(H2, O2, H2O):
    dG_R = (GibbsFreeEnergy(H2O) - 
            (GibbsFreeEnergy(H2) + 
             GibbsFreeEnergy(O2)*0.5))
    return dG_R

def K(dHS):
    dH, dS = dHS
    return np.exp(-dH/R/T + dS/R)

# %%

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
       

    # %%  Automatically calling the thermodynamic states
    species_list = ['H2', 'O2', 'H2O']  
    TD = get_thermo_properties_dict(species_list, T, P)
    dG_R = TD['H2O']['G'] - TD['H2']['G'] - TD['O2']['G'] * 0.5

    H2 = np.array([TD['H2']['H'], TD['H2']['S']])
    # O2 = np.array([TD['O2']['H'], TD['O2']['S']])
    H2O = np.array([TD['H2O']['H'], TD['H2O']['S']])
    
    # %%
            
    # adsorption enthalpies and entropies
    Ni = np.array([0.0e3, 0.0])
    HNi = np.array([-31.19e3, 41.14])
    YSZ = np.array([0.0e3, 0.0])
    O2YSZ = np.array([-85.6e3, 139.6])
    OHYSZ = np.array([-165.23e3, 175.1])
    H2OYSZ = np.array([-254.45e3, 54.34])
    OxYSZ = np.array([-85.6e3, 139.6])
    VoYSZ = np.array([0.0e3, 0.0])
    
    # reaction 1: H2 + 2 Ni <-> 2H(Ni)
    HS1 = 2*HNi - H2 - 2 *Ni
    K1 = K(HS1)
    # reaction 2: H(ni) + O2-(YSZ) <-> (Ni) + OH-(YSZ) + e-(Ni)
    HS2 = Ni + OHYSZ - HNi - O2YSZ
    K2 = K(HS2)
    # # reaction 3: H(Ni) + OH-(YSZ) <-> (Ni) + H2O(YSZ) + e-(Ni) -- RDS
    HS3 = Ni + H2OYSZ - HNi - OHYSZ
    K3 = K(HS3)
    # # reaction 4: H2O(YSZ) <-> H2O(g) + (YSZ)
    HS4 = H2O + YSZ - H2OYSZ
    K4 = K(HS4)
    # # reaction 5: Ox(YSZ) + (YSZ) <-> O2-(YSZ) + Vö(YSZ)
    HS5 = O2YSZ + VoYSZ - OxYSZ - YSZ
    a = 25
    K5 = K(HS5) / a
    
    # %%
    # Kinetic coefficients
    # i_A = 0.
    # i_A_lim = 1.0
    A = 2.4e-7; B = 3.55; C = 0.144; D = 0.437
    AtoD = np.array([A, B, C, D])
    K1to5 = np.array([K1, K2, K3, K4, K5])
    
    # %%
    # create the object Echem

    N = 10000
    etaAct = np.linspace(0., 0.8, N)
    # beta3 = 0.5
    k_c3 = 0.2e-16
    l_tpb = 10e4        # I will have to check this value
    Echem = Echem(T, P, xH2, xH2O, xO2, AtoD, etaAct, l_tpb, K1to5)
    
    # calculate the OCV
    U_rev = Echem.calcUrev(dG_R)
    
    eta_leak = Echem.calcLeakVoltage()
    
    OCV = U_rev - eta_leak
    
    
    # %%
    
    xH2_array = np.linspace(0., 1., N)
    
    
    # %%
    # # plot Nickel
    
    covNi, covH = Echem.coveragesNickel(N)
    
    # plt.figure()
    # plt.title('Surface Coverage Ni (p = 1 atm)')
    # plt.plot(xH2_array, covNi, label = r'$\theta$Ni')
    # plt.plot(xH2_array, covH, label = r'$\theta$H(Ni)')
    # plt.xlim([0, 1])
    # plt.ylim([0, 1])
    # plt.xlabel(r'$x_{H2}$ / mol/mol')
    # plt.ylabel('coverage / 1')
    # plt.legend()
    
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




















