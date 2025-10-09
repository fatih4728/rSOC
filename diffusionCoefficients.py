# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 10:59:55 2025

In this file I want to calculate the diffusion coeffcients

@author: smfadurm
"""


# %%
import numpy as np

import cantera as ct

p = ct.one_atm
R = ct.gas_constant*1e-3
T = 900

cT = p/R/T
M_A = 2e-3
M_B = 18e-3
D_AB_Colin = 8.e-4
D_Kn_Colin = 2.8e-5, 9.4e-6

# %%
def D_Fuller(T, P = 1.):
    """    
    Fuller equation for the binary diffusion coefficient

    Parameters
    ----------
    T : float
        Temperature / K.
    P : float, optional
        Pressure / atm. The default is 1..

    Returns
    -------
    Binary diffusion coefficient.

    """
    M_A = 2.01588
    M_B = 18.01528
    v_A = 7.07
    v_B = 12.7
    return (1e-3*T**1.75*(1/M_A+1/M_B)**0.5/
            (P*(v_A**(1/3)+v_B**(1/3))**2))

# D_Full = D_Fuller(800+273.15)
# print(D_Full*1.e-4)

# %%
def D_Knudsen(rp, T, M):
    """
    Source is Transport Phenomena in materials processing

    Parameters
    ----------
    rp : float
        pore radius / cm.
    T : float
        temperature / K.
    M : float
        molar mass / g/mol.

    Returns
    -------
    float
        Knudsen Diffusion Coefficient.

    """
    return 9700*rp*(T/M)**0.5

rp = 1.5e-6*1e0
D_knudsen = np.array([D_Knudsen(rp, 1000, 2), D_Knudsen(rp, 1000, 18)])
print(D_knudsen*1e-4/2)


# %%
def D_ChapmanEnskog(cT, M_A, M_B, T):

    sigmaH2 = 2.827
    sigmaH2O = 2.641
    sigmaAB = (sigmaH2 + sigmaH2O) / 2
    kB = ct.boltzmann
    bruch = 3 / (32 * cT*1e-8 * sigmaAB**2) 
    klammer = (8 * kB * T / np.pi* (1/M_A + 1/M_B))**0.5
    return bruch*klammer


# D_AB = D_ChapmanEnskog(cT, M_A, M_B, T)
# print(D_AB*1.e-4)
