# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 10:59:55 2025

In this file I want to calculate the diffusion coeffcients

@author: smfadurm
"""


import numpy as np
import cantera as ct

p = ct.one_atm
R = 8.3145
T = 900

cT = p/R/T
M_A = 2e-3
M_B = 18e-3
D_Colin = 3.23e-8

def D_ChapmanEnskog(cT, M_A, M_B, T):
    sigmaH2 = 2.827
    sigmaH2O = 2.641
    sigmaAB = (sigmaH2 + sigmaH2O) / 2
    kB = ct.boltzmann
    bruch = 3 / (32 * cT*1e-8 * sigmaAB**2) 
    klammer = (8 * kB * T / np.pi* (1/M_A + 1/M_B))**0.5
    return bruch*klammer


D_AB = D_ChapmanEnskog(cT, M_A, M_B, T)
print(D_AB*1.e-4)

# def D_Fuller(T, P = 1.):
#     M_A = 2.01588
#     M_B = 18.01528
#     v_A = 7.07
#     v_B = 12.7
#     return (1e-3*T**1.75*(1/M_A+1/M_B)**0.5/
#             (P*(v_A**(1/3)+v_B**(1/3))**2))

# D_Full = D_Fuller(800+273.15)
# print(D_Full*1.e-4)


