# -*- coding: utf-8 -*-
"""
Created on Wed Sep  3 14:03:41 2025

Here I want the Ficks Diffusion

@author: smfadurm
"""

import numpy as np


def calculateConcentration(J, L, D, c1):
    return c1 - J*L/D

if __name__=="__main__":

    # universal constants
    T = 800 + 273.15        # operating temperature of the soc
    P = 101325.0          # operating pressure of the soc
    F = 96485.33212331001          # faraday constant
    R = 8.314462618153241     # import universal gas constant
    i = 5.                   # current density in A/cm²

    
    # geometry
    epsilon = 0.3       # porosity of the electrode (affects mass balance; smaller -> worse)
    tau = 5.            # tortuosity of the electrode
    rp = 5.e-7          # radius of the pores (affects mass balance; smaller -> worse)
    dEle = 200e-6       # thickness of the electrode (affects mass balance a lot; bigger -> worse)
    A = 100             # area in cm²
    
    # gas variables
    M = np.array([2e-3, 18e-3]) # molar weights of H2, H2O

    J = np.array([-7.81e-7, 7.81e-7])
    L = 200*1e-6
    D = 1e-9
    Deff = D * epsilon / tau
    
    cT = P/R/T
    x1 = np.array([0.5, 0.5])
    c1 = x1 * cT
    c2 = calculateConcentration(J, L, Deff, c1)
    x2 = c2 / cT
    
    print(x1)
    print(x2)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    