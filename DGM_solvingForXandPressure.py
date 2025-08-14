# -*- coding: utf-8 -*-
"""
Created on Wed Aug  6 09:36:12 2025

Here I want to use the DGM in the discretized form to solve for x2

I further want to solve a fitting pressure gradient

@author: smfadurm
"""


import cantera as ct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import newton
from CoolProp.CoolProp import PropsSI

# %%

# Functions used here
def permeabilityFactorBg(epsilon, tau, rp):
    """
    This calculates the permeability factor commonly noted as Bg
    Here the Kozeny-Carman relationship is used, according to the
    Zhou paper.
    ___
    epsilon = porosity
    tau = tortuosity
    rp = pore radius
    ___
    
    """
    return epsilon**3 * (2*rp)**2 / 72 / tau / (1 - epsilon)**2

def knudsenCoeff(porosity,tortousity, molarWeight):
    """
    This function calculates the Knudsen Coefficient
    """
    return 4/3*porosity/tortousity*rp*np.sqrt(8*R*T/np.pi/molarWeight)

def calculateMolarFlux(i, A):
    """
    This function calculates the anode molar flux of a SOC
    Only valid for pure hydrogen and steam

    Parameters
    ----------
    i : float
        current density in A/cm².
    A : Area
        Area in cm².

    Returns
    -------
    J : Array
        molar flux of H2 and H2O.

    """
    ndotH2 = i/F/2
    ndotH2O = ndotH2
    j = np.array([-ndotH2, ndotH2O])
    return j * A

def Fdxdz(J, XC, D, Dkn, Bg, muMix, dpdz):
    """
    this calculates the derivative of the molar concentration

    Parameters
    ----------
    J : Float/Array
        Molar Flux.
    XC : Float/Array
        Concentration at channel.
    D : Float/Array
        Diffusion Coefficient.
    Dkn : Float/Array
        Knudsen Coefficient.
    Bg : Float
        Permeability factor.
    muMix : Float
        Kinetic viscosity of the mix.
    dpdz : Float
        pressure gradient.

    Returns
    -------
    dxdz : Float/Array
        concentration gradient.

    """
    XT = XC.sum()
    N = len(J)
    dxdz = np.zeros((N))
    for k in range(N):
        sumTerm = 0.
        for l in range(N):  
            if l != k:
                sumTerm += (XC[l]*J[k] - XC[k]*J[l])/(XT*D[k, l])
        dxdz[k] = sumTerm + J[k]/Dkn[k] - XC[k]*Bg/Dkn[k]/muMix*dpdz
    dxdz *= -1
    return dxdz

def calculateMoleFraction(L, acc, maxSteps):
        """
        calculate the mole fraction at L

        Parameters
        ----------
        L : float
            length of solid phase.
        acc : float
            accuracy.
        maxSteps : float/integer
            max steps before simulation breaks.

                

        Returns
        -------
        xM : float/array
            molar fraction at L.
        dp : float
            pressure difference required.

        """
        dp = 1
        counter = 0
        xM = np.array([0., 0.])    
        while abs(xM.sum()-1) > acc:
            counter += 1
            if xM.sum()-1 > 0:
                dp += 1
            else:
                dp -= 1
            dpdz = dp/L        
            dxdz = Fdxdz(J, XC, D, Dkn, Bg, muMix, dpdz)
            XM = XC - dxdz * L    
            xM = XM / XT
                        
            if  counter > maxSteps:
                break
            if dp >= 1.e5:
                print("pressure difference required is > 1 bar")
                break
        return xM, dp, [dpdz, dxdz, counter]

# %%
# universal constants

if __name__=="__main__":

    # universal constants
    T = 800 + 273.15        # operating temperature of the soc
    P = ct.one_atm          # operating pressure of the soc
    F = ct.faraday*1e-3          # faraday constant
    R = ct.gas_constant*1e-3     # import universal gas constant
    i = 5.                   # current density in A/cm²

    
    # geometry
    epsilon = 0.3       # porosity of the electrode (affects mass balance; smaller -> worse)
    tau = 5.            # tortuosity of the electrode
    rp = 5.e-7          # radius of the pores (affects mass balance; smaller -> worse)
    dEle = 200e-6       # thickness of the electrode (affects mass balance a lot; bigger -> worse)
    A = 100             # area in cm²
    
    # gas variables
    M = np.array([2e-3, 18e-3]) # molar weights of H2, H2O
    XT = P/R/T
    xC = np.array([0.9, 0.1])      # the molar fraction at the channel
    XC = xC*XT                     # the molar concentration
    
    # viscosity of the mix
    mu = np.array([1.84e-5, 3.26e-5])    # dynamic visosities at 600°C
    muMix = (xC*mu).sum()                 # dynamic for the mix
    
    # %%
    
    # molar flux of the species
    J = calculateMolarFlux(i, A)
    
    # diffusion matrix
    D = np.zeros((2, 2)) + 1.e-5
    
    # knudsen Array
    Dkn = knudsenCoeff(epsilon, tau, M)
    
    # Permeability by the Kozeny-Carman relationship
    Bg = permeabilityFactorBg(epsilon, tau, rp)
    
    # start the loop to calculate mole fraction and pressure difference            
    xM, dp, addInfo = calculateMoleFraction(L = dEle, acc = 1.e-5, maxSteps = 1e6)

    
    if dp < 1.e5:
        print(f"\nxH2 at the electrolyte is {xM[0]:.5}")
        print(f"The sum of the molar fractions is {xM.sum():.4}")
        print(f"dp = {dp*1e-2:.3} mbar")























